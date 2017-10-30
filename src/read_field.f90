



module read_field

public get_field


contains

!------------------------------------------------------------------------------
subroutine get_field(iter)
     use globals, only: input_dir, type, type_fmt
     implicit none

     integer        :: iter
     integer        :: Ntemp(3), hcheck, dcheck
     character(100) :: input_name
     character(8)   :: it_number


     ! setup field for reading
     write(*,*) " Step number ", iter
     write(it_number,"(I8.8)") iter
     input_name = trim(adjustl(input_dir))
     !print *, " Reading field from", trim(adjustl(input_name))

     ! select filetype
     if(index(type,"h5").ne.0 .OR. index(type,"H5").ne.0) then
        if(index(type_fmt,"complex").ne.0 .OR. index(type_fmt,"COMPLEX").ne.0) then
           call read_esio_cfield(input_name)
        elseif(index(type_fmt,"real").ne.0 .OR. index(type_fmt,"REAL").ne.0) then
           call read_esio_rfield(input_name)
        elseif(index(type_fmt,"jhtdb").ne.0 .OR. index(type_fmt,"JHTDB").ne.0) then
        !   call jhtdb_slabs() 
           call jhtdb_slices()
        else
           write(*,*) "UNRECOGNIZED INPUT DATA TYPE, USE complex OR real.  EXITING..."
           stop
        endif

     elseif(index(type,"dat").ne.0 .OR. index(type,"DAT").ne.0) then
        if(index(type_fmt,"complex").ne.0 .OR. index(type_fmt,"COMPLEX").ne.0) then
           write(*,*) "DAT STYLE FILE ONLY AVAILABLE FOR real.  EXITING..."
           stop
        elseif(index(type_fmt,"real").ne.0 .OR. index(type_fmt,"REAL").ne.0) then
           call read_basic_field(input_name)
        else
           write(*,*) "DAT STYLE FILE ONLY AVAILABLE FOR real.  EXITING..."
           stop
        endif

     else
        write(*,*) "UNRECOGNIZED INPUT FILE TYPE, EXITING..."
        stop
     endif


end subroutine get_field


!------------------------------------------------------------------------------
subroutine read_basic_field(filename)
     use globals
     use, intrinsic :: iso_c_binding
     implicit none
     include 'fftw3.f03'


        integer            :: N_file, i, j, k, oo
        integer            :: Ntemp(3)
        integer, parameter :: field_input = 90
        integer, parameter :: input_unit = 10
        real(8)            :: kol
        character(*)       :: filename
        character(8)       :: it_number


        print *, " Reading field from", trim(adjustl(filename))

        ! Open field and check size
        open(unit=field_input,file=trim(adjustl(filename)),form="UNFORMATTED",status="OLD")
        read(field_input) Ntemp(1), Ntemp(2), Ntemp(3), t
        print *, " Time from field file", t
        if(Ntemp(1).ne.N(1) .OR. Ntemp(2).ne.N(2) .OR. Ntemp(3).ne.N(3)) then
           print*, " WARNING: Field size does not match input (this probably won't work)"
           !stop
        endif


        ! read field (new format from fvles)
        do k = 1, N(3)
           do j = 1, N(2)
              do i = 1, N(1)
                 read(field_input) u(i,j,k)
              enddo
           enddo
        enddo

        do k = 1, N(3)
           do j = 1, N(2)
              do i = 1, N(1)
                 read(field_input) v(i,j,k)
              enddo
           enddo
        enddo

        do k = 1, N(3)
           do j = 1, N(2)
              do i = 1, N(1)
                 read(field_input) w(i,j,k)
              enddo
           enddo
        enddo

        close(unit=field_input)


        ! transform fields
        call dfftw_execute_dft_r2c(plan_r2c,u, s_u)
        call dfftw_execute_dft_r2c(plan_r2c,v, s_v)
        call dfftw_execute_dft_r2c(plan_r2c,w, s_w)
        s_u = s_u*fieldnorm
        s_v = s_v*fieldnorm
        s_w = s_w*fieldnorm

        ! nu and diss
        if (diss_flag == 0) then
            ! diss_flag = 0 means dissipation is read in from file
            diss = 62.3468D0
            print*, " WARNING: no nu and dissipation specified with old format"
            print*, "          assuming diss =",diss
        else if (diss_flag == 1) then
            ! diss_flag = 1 means dissipation is calculated from field
            print*, " ERROR: Cannot calculate dissipation without nu, which is not specified for the old format"
            stop
        else if (diss_flag == 2) then
            ! diss_flag = 2 means dissipation is specified in input file
            print*, " WARNING: no nu specified with old format"
            print*, "          calculating nu using specified diss =", diss
        endif
        kol = 3.0D0 * maxval(delta*dble(f3))/(2.0D0*pi)
        nu = (kol**4 * diss)**(1D0/3D0)
        print*, "          calculated nu =",nu

    return

    end subroutine read_basic_field


!------------------------------------------------------------------------------
subroutine read_esio_rfield(filename)
use mpi
use esio
use globals
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

  integer :: i, j, k
  integer :: ierr, idummy(1), Nin(3), rflag
  real(8) :: rdummy(1), imax, rmax, Lin(3)
  character(*) :: filename
  real(C_DOUBLE), dimension(1:N(1),1:N(2),1:N(3)) :: ubuffer
  type(esio_handle) :: ehandle


   print*, " Reading data from HDF5 file: ",trim(adjustl(filename))//".h5"
   print*, " "

   ! initialize esio file/fields
   call esio_handle_initialize(ehandle,MPI_COMM_WORLD)
   call esio_file_open(ehandle,trim(adjustl(filename))//".h5",.true.)
   call esio_line_establish(ehandle,1,1,1,ierr)

   call esio_line_read_double(ehandle,"time",rdummy,0)
   t = rdummy(1)

   ! consistency checks
   call esio_line_read_double(ehandle,"nu",rdummy,0)
   nu = rdummy(1)

   if (diss_flag == 0) then
       call esio_line_read_double(ehandle,"targetDiss",rdummy,0)
       diss = rdummy(1)
   endif

   call esio_line_read_integer(ehandle,"Nx",idummy,0)
   Nin(1) = idummy(1)

   call esio_line_read_integer(ehandle,"Ny",idummy,0)
   Nin(2) = idummy(1)

   call esio_line_read_integer(ehandle,"Nz",idummy,0)
   Nin(3) = idummy(1)

   call esio_line_read_double(ehandle,"Lx",rdummy,0)
   Lin(1) = rdummy(1)

   call esio_line_read_double(ehandle,"Ly",rdummy,0)
   Lin(2) = rdummy(1)

   call esio_line_read_double(ehandle,"Lz",rdummy,0)
   Lin(3) = rdummy(1)

   call check_box(Lin,L,Nin,N)

   print*, " viscosity from file: ", nu
   if (diss_flag == 0) print*, " dissipation from file: ", diss

   ! velocity fields
   call esio_field_establish (ehandle, Nin(3), 1, Nin(3), Nin(2), 1, Nin(2), Nin(1), 1, Nin(1), ierr)
print*, " esio field established"
   call esio_field_readv_double (ehandle, "u", ubuffer(:,:,:), 1)
print*, " esio u read"
   u(1:N(1),1:N(2),1:N(3)) = ubuffer(1:N(1),1:N(2),1:N(3))
   call esio_field_readv_double (ehandle, "v", ubuffer(:,:,:), 1)
   v(1:N(1),1:N(2),1:N(3)) = ubuffer(1:N(1),1:N(2),1:N(3))
   call esio_field_readv_double (ehandle, "w", ubuffer(:,:,:), 1)
   w(1:N(1),1:N(2),1:N(3)) = ubuffer(1:N(1),1:N(2),1:N(3))
print*, " esio all fields read"

   ! transform fields
   call dfftw_execute_dft_r2c(plan_r2c,u, s_u)
   call dfftw_execute_dft_r2c(plan_r2c,v, s_v)
   call dfftw_execute_dft_r2c(plan_r2c,w, s_w)
   s_u = s_u*fieldnorm
   s_v = s_v*fieldnorm
   s_w = s_w*fieldnorm

   ! clean-up
   call esio_file_close(ehandle)
   call esio_handle_finalize(ehandle)


return

end subroutine read_esio_rfield


!------------------------------------------------------------------------------
subroutine read_esio_cfield(filename)
use mpi
use esio
use globals
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

  integer :: i, j, k
  integer :: ierr, idummy(1), Nin(3), rflag
  real(8) :: rdummy(1), imax, rmax, Lin(3)
  real(8) :: real_part, imag_part, kfmax, kfmin
  character(*) :: filename
  real(C_DOUBLE), dimension(2,1:N(2),1:N(1)/2,1:N(3)) :: ubuffer ! y,x,z
!  real(C_DOUBLE), dimension(2,1:N(3),1:N(1)/2,1:N(2)) :: ubuffer
  type(esio_handle) :: ehandle


   print*, " Reading data from HDF5 file: ",trim(adjustl(filename))//".h5"
   print*, " (note: set for POONGBACK-style files...)"
   print*, " "

   ! initialize esio file/fields
   call esio_handle_initialize(ehandle,MPI_COMM_WORLD)
   call esio_file_open(ehandle,trim(adjustl(filename))//".h5",.true.)
   call esio_line_establish(ehandle,1,1,1,ierr)

   call esio_line_read_double(ehandle,"time",rdummy,0)
   t = rdummy(1)

   call esio_line_read_integer(ehandle,"timestep",idummy,0)
   ts = idummy(1)

   ! consitency checks
   call esio_line_read_double(ehandle,"nu",rdummy,0)
   nu = rdummy(1)

   if (diss_flag.eq.0) then
       call esio_line_read_double(ehandle,"targetDiss",rdummy,0)
       diss = rdummy(1)

       call esio_line_read_double(ehandle,"kfmax",rdummy,0)
       kfmax = rdummy(1)

       call esio_line_read_double(ehandle,"kfmin",rdummy,0)
       kfmin = rdummy(1)
   endif

   if ((N(1).eq.N(2)).and.(N(2).eq.N(3))) then
       call esio_line_read_integer(ehandle,"nx",idummy,0)
       Nin(1) = idummy(1)
       Nin(2) = idummy(1)
       Nin(3) = idummy(1)
       Lin = L
   else
       call esio_line_read_integer(ehandle,"nx",idummy,0)
       Nin(1) = idummy(1)
       call esio_line_read_integer(ehandle,"ny",idummy,0)
       Nin(2) = idummy(1)
       call esio_line_read_integer(ehandle,"nz",idummy,0)
       Nin(3) = idummy(1)
       
       call esio_line_read_double(ehandle,"LMFx",rdummy,0)
       Lin(1) = 2*rdummy(1)
       call esio_line_read_double(ehandle,"LMFy",rdummy,0)
       Lin(2) = 2*rdummy(1)
       call esio_line_read_double(ehandle,"LMFz",rdummy,0)
       Lin(3) = 2*rdummy(1)
   endif

   print*, " File info... "
   print*, " + viscosity: ", nu
   print*, " + N: ", Nin
   print*, " + LM: ", Lin
   if (diss_flag.eq.0) then
        print*, " + dissipation: ", diss
        print*, " + min forcing k: ", kfmin
        print*, " + max forcing k: ", kfmax
   endif

   Lin = L
   call check_box(Lin,L,Nin,N)

   ! initialize field read :: z, x, y
!   call esio_field_establish (ehandle, N(3), 1, N(3), N(1)/2, 1, N(1)/2, N(2), 1, N(2), ierr)
   call esio_field_establish (ehandle, N(2), 1, N(2), N(1)/2, 1, N(1)/2, N(3), 1, N(3), ierr)

   ! velocity fields :: z, x, y
   call esio_field_readv_double (ehandle, "u", ubuffer(:,:,:,:),2)
   s_u = 0.0D0
   do i = 1,N(1)/2
      do j = 1,N(2)
         do k = 1,N(3)
!            real_part = ubuffer(1,k,i,j)
!            imag_part = ubuffer(2,k,i,j)
            real_part = ubuffer(1,j,i,k)
            imag_part = ubuffer(2,j,i,k)
            s_u(i-1,j-1,k-1) = cmplx(real_part, imag_part)
         enddo
      enddo
   enddo

   call esio_field_readv_double (ehandle, "v", ubuffer(:,:,:,:),2)
   s_v = 0.0D0
   do i = 1,N(1)/2
      do j = 1,N(2)
         do k = 1,N(3)
!            real_part = ubuffer(1,k,i,j)
!            imag_part = ubuffer(2,k,i,j)
            real_part = ubuffer(1,j,i,k)
            imag_part = ubuffer(2,j,i,k)
            s_v(i-1,j-1,k-1) = cmplx(real_part, imag_part)
         enddo
      enddo
   enddo

   call esio_field_readv_double (ehandle, "w", ubuffer(:,:,:,:),2)
   s_w = 0.0D0
   do i = 1,N(1)/2
      do j = 1,N(2)
         do k = 1,N(3)
!            real_part = ubuffer(1,k,i,j)
!            imag_part = ubuffer(2,k,i,j)
            real_part = ubuffer(1,j,i,k)
            imag_part = ubuffer(2,j,i,k)
            s_w(i-1,j-1,k-1) = cmplx(real_part, imag_part)
         enddo
      enddo
   enddo


   ! clean-up
   call esio_file_close(ehandle)
   call esio_handle_finalize(ehandle)


return

end subroutine read_esio_cfield


!------------------------------------------------------------------------------
subroutine jhtdb_slices()
     use globals
     use, intrinsic :: iso_c_binding
     implicit none
     include 'fftw3.f03'


        integer            :: N_file, i, j, k, oo
        integer            :: Ntemp(3), inc
        integer, parameter :: ufield_input = 90
        integer, parameter :: vfield_input = 91
        integer, parameter :: wfield_input = 92
        integer, parameter :: input_unit = 10
        real(8)            :: kol
        real(4)            :: ubuffer(N(1),N(2)*N(3))
        character(50)      :: ufilename, vfilename, wfilename
        character(8)       :: it_number


        !print *, " Reading field from", trim(adjustl(filename))
        ufilename = './data/JHTDB_DNS_1024/jhtdb_dns1024_u.dat'
        vfilename = './data/JHTDB_DNS_1024/jhtdb_dns1024_v.dat'
        wfilename = './data/JHTDB_DNS_1024/jhtdb_dns1024_w.dat'

        ! Open field and check size
        open(unit=ufield_input,file=trim(adjustl(ufilename)),form="UNFORMATTED",status="OLD")
        open(unit=vfield_input,file=trim(adjustl(vfilename)),form="UNFORMATTED",status="OLD")
        open(unit=wfield_input,file=trim(adjustl(wfilename)),form="UNFORMATTED",status="OLD")

        ! read field
        do i = 1, N(1)
           read(ufield_input) ubuffer(i,1:N(2)*N(3))
        enddo

        do i = 1,N(1)
           inc = 0
           do j = 1,N(2)
              do k = 1,N(3)
                 inc = inc + 1
!                 u(i,j,k) = dble(ubuffer(i,inc))
                 u(i,j,k) = 1.0d0*ubuffer(i,inc)
              enddo
           enddo
        enddo

        do i = 1, N(1)
           read(vfield_input) ubuffer(i,1:N(2)*N(3))
        enddo

        do i = 1,N(1)
           inc = 0
           do j = 1,N(2)
              do k = 1,N(3)
                 inc = inc + 1
!                 v(i,j,k) = dble(ubuffer(i,inc))
                 v(i,j,k) = 1.0d0*ubuffer(i,inc)
              enddo
           enddo
        enddo


        do i = 1, N(1)
           read(wfield_input) ubuffer(i,1:N(2)*N(3))
        enddo

        do i = 1,N(1)
           inc = 0
           do j = 1,N(2)
              do k = 1,N(3)
                 inc = inc + 1
!                 w(i,j,k) = dble(ubuffer(i,inc))
                 w(i,j,k) = 1.0d0*ubuffer(i,inc)
              enddo
           enddo
        enddo

        close(unit=ufield_input)
        close(unit=vfield_input)
        close(unit=wfield_input)

        ! transform fields
        call dfftw_execute_dft_r2c(plan_r2c,u, s_u)
        call dfftw_execute_dft_r2c(plan_r2c,v, s_v)
        call dfftw_execute_dft_r2c(plan_r2c,w, s_w)
        s_u = s_u*fieldnorm
        s_v = s_v*fieldnorm
        s_w = s_w*fieldnorm

        ! nu and diss
        ! We assume nu is known, but it could also be calculated.
        ! See read_basic_field() for implementation of nu calculation
        if (diss_flag == 0) then
            ! diss_flag = 0 means dissipation is read in from file
            ! diss = 0.0928d0
            diss = 0.103d0
            nu = 0.000185d0
            print*, " WARNING: no nu and dissipation specified with old format"
            print*, "          assuming diss =",diss
            print*, "          assuming nu   =",nu
        else if (diss_flag == 1) then
            ! diss_flag = 1 means dissipation is calculated from field
            print*, " ERROR: Cannot calculate dissipation without nu, which is not specified for the old format"
            stop
        else if (diss_flag == 2) then
            ! diss_flag = 2 means dissipation is specified in input file
            nu = 0.000185d0
            print*, " WARNING: Using the dissipation specified in the input file."
            print*, " WARNING: no nu specified with old format"
            print*, "          assuming nu   =",nu
        endif

    return

    end subroutine jhtdb_slices


!------------------------------------------------------------------------------
subroutine jhtdb_slabs()
     use globals
     implicit none

     integer        :: iter, nslab, i, j, k, nparts
     integer        :: Ntemp(3), hcheck, dcheck
     character(100) :: data_dir, fbase
     character(100) :: input_name
     character(8)   :: it_number
     character(1)   :: c_nslab1
     character(3)   :: c_nslab3

     ! nu and diss
     ! We assume nu is known, but it could also be calculated.
     ! See read_basic_field() for implementation of nu calculation
     if (diss_flag == 0) then
         ! diss_flag = 0 means dissipation is read in from file
         diss = 0.000185d0
         nu = 0.0928d0
         print*, " WARNING: no nu and dissipation specified with old format"
         print*, "          assuming diss =",diss
         print*, "          assuming nu   =",nu
     else if (diss_flag == 1) then
         ! diss_flag = 1 means dissipation is calculated from field
         print*, " ERROR: Cannot calculate dissipation without nu, which is not specified for the old format"
         stop
     else if (diss_flag == 2) then
         ! diss_flag = 2 means dissipation is specified in input file
         nu = 0.0928d0
         print*, " WARNING: Using the dissipation specified in the input file."
         print*, " WARNING: no nu specified with old format"
         print*, "          assuming nu   =",nu
     endif

     !data_dir = "./data/jhtdb/t0/"
     data_dir = "./data/JHTDB_DNS_1024/t0/"
     fbase = "isotropic1024coarse"
     nparts = 8

     do i = 1,nparts

        ! construct file name
        nslab = N(3)/nparts * (i-1)
        if(i.eq.1) write(c_nslab1,'(I1)') nslab
        if(i.gt.1) write(c_nslab3,'(I3)') nslab

        if(i.eq.1) input_name = trim(adjustl(data_dir))//trim(adjustl(fbase))//"_"//c_nslab1//"z.h5"
        if(i.gt.1) input_name = trim(adjustl(data_dir))//trim(adjustl(fbase))//"_"//c_nslab3//"z.h5"

        ! read slab
        call read_h5_field(input_name,nparts,i)

     enddo

     print*, " "     
     print*, " Slabs read and assembled"
     print*, "  - max u:", maxval(u)
     print*, "  - max v:", maxval(v)
     print*, "  - max w:", maxval(w)

     ! test slice for plotting
     open(unit=11,file="slice.txt",status='new')
     i = 512
     do j = 1,N(2)
        do k = 1,N(3)
        write(11,*) i, j, k, u(i,j,k)
        enddo
     enddo
     close(11)

     ! transform fields 
     call dfftw_execute_dft_r2c(plan_r2c,u, s_u)
     call dfftw_execute_dft_r2c(plan_r2c,v, s_v)
     call dfftw_execute_dft_r2c(plan_r2c,w, s_w)
     s_u = s_u*fieldnorm
     s_v = s_v*fieldnorm
     s_w = s_w*fieldnorm

     ! swap storage for compatibility with POONGBACK
     !call swap()

     return

end subroutine jhtdb_slabs
!------------------------------------------------------------------------------


subroutine read_h5_field(filename,nparts,part)
use mpi
use hdf5
use esio
use globals
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

  integer :: i, j, k, part, nslab, nstart, nend, nparts
  integer :: ierr, idummy(1), Nin(3), rflag
  real(8) :: rdummy(1), imax, rmax, Lin(3)
  character(*) :: filename
!  real(C_DOUBLE), dimension(1:N(3)/nparts,1:N(2),1:N(1),3) :: ubuffer
  real(C_FLOAT), dimension(1:N(3)/nparts,1:N(2),1:N(1),3) :: ubuffer
  type(esio_handle) :: ehandle

  ! hdf5 interface 
  INTEGER(HID_T) :: file_id
  INTEGER(HID_T) :: dset_id
  INTEGER(HSIZE_T), DIMENSION(4) :: data_dims


   print*, " Reading data from HDF5 file: ",trim(adjustl(filename))
   print*, " "

   nslab = N(3)/nparts
   nstart = nslab * (part-1) + 1
   nend = nslab * part
   data_dims(1) = nslab
   data_dims(2) = N(2)
   data_dims(3) = N(1)
   data_dims(4) = 3

   ! Initialize FORTRAN interface
   CALL h5open_f(ierr)

   ! Open an existing file
   CALL h5fopen_f (trim(adjustl(filename)), H5F_ACC_RDWR_F, file_id, ierr)

   ! Open an existing dataset
   CALL h5dopen_f(file_id, "u00000", dset_id, ierr)

   ! Read the dataset
   !CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, ubuffer, data_dims, ierr)
   CALL h5dread_f(dset_id, H5T_NATIVE_REAL, ubuffer, data_dims, ierr)

   ! Close the dataset
   CALL h5dclose_f(dset_id, ierr)

   ! Close the file
   CALL h5fclose_f(file_id, ierr)

   ! Close FORTRAN interface
   CALL h5close_f(ierr)

   u(1:N(1),1:N(2),nstart:nend) = dble(ubuffer(1:nslab,1:N(2),1:N(1),1))
   v(1:N(1),1:N(2),nstart:nend) = dble(ubuffer(1:nslab,1:N(2),1:N(1),2))
   w(1:N(1),1:N(2),nstart:nend) = dble(ubuffer(1:nslab,1:N(2),1:N(1),3))


return

end subroutine read_h5_field
!------------------------------------------------------------------------------


!!$subroutine swap()
!!$use mpi
!!$use esio
!!$use globals
!!$use, intrinsic :: iso_c_binding
!!$implicit none
!!$include 'fftw3.f03'
!!$
!!$  integer :: i, j, k
!!$
!!$  do i = 0,N(1)/2
!!$     do j = 0,N(2)-1
!!$        do k = 0,N(3)-1
!!$           s_u(j,i,k) = s_u(i,j,k)
!!$           s_v(j,i,k) = s_v(i,j,k)
!!$           s_w(j,i,k) = s_w(i,j,k)
!!$        enddo
!!$     enddo
!!$  enddo
!!$
!!$return
!!$
!!$end subroutine swap
!------------------------------------------------------------------------------


subroutine check_box(Lin,Lsim,Nin,Nsim)
use mpi
use, intrinsic :: iso_c_binding
implicit none

   integer, INTENT(in)  :: Nin(3), Nsim(3)
   real(C_DOUBLE), INTENT(in)  :: Lin(3), Lsim(3)

  if(Nin(1).ne.Nsim(1) .OR. Nin(2).ne.Nsim(2) .OR. Nin(3).ne.Nsim(3)) then
     write(*,*) " WARNING:: input N differs from simulation N"
     write(*,*) "           input N:", Nin
     write(*,*) "           simulation N:", Nsim
     !stop
  elseif(Lin(1).ne.Lsim(1) .OR. Lin(2).ne.Lsim(2) .OR. Lin(3).ne.Lsim(3)) then
     write(*,*) " WARNIGN:: input box size differs from simulation"
     !stop
  endif


  return

end subroutine check_box

end module read_field
!------------------------------------------------------------------------------
