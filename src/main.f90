
!------------------------------------------------------------------------------
!
!  SWH June 2012: rewrite of analyze.c to correct data structure problems,
!                 update to fftw3.3, and remove the need to write out
!                 transformed fields beforehand
!
!                 fieldconv.c is no longer necessary as is the whole directory
!
!  SWH Feb 2016:  - adding H5/ESIO file support
!                 - changing filter to be real-box
!                 - divided analye.f90 into multiple files for 
!                   clarity (hopefully)
!
!-------------------------------------------------------------------------------


program main
  use mpi
  use globals
  use read_field
  use slice
  use statistics
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'

  integer            :: samples, it, Nmax, i, ones(3), fiso(3)
  integer            :: world_size, world_rank, ierr
  integer, parameter :: field_input = 90
  integer, parameter :: input_unit = 10
  real(8)            :: diss_calc
  character(100)     :: input_name, ifile
  character(100)     :: buffer, flag
  real(8)            :: f3sq(3)


  ! initialize MPI
  call MPI_Init (ierr)
  call MPI_Comm_size (MPI_COMM_WORLD, world_size, ierr)
  call MPI_Comm_rank (MPI_COMM_WORLD, world_rank, ierr)

  ! read input parameters
  call get_input()

  ! initialize fields
  call initialize()

  ! loop over specified restart field data
  Nmax = maxval(N)
  !Nmax = maxval(N/f3)+1
  samples = 0
  do it = it_ini, it_fin, it_step
 
     samples = samples + 1

     ! read-in field (output is in wavespace)
     call get_field(it)

     ! write a slice
     !call write_slice(N,Nmax,f3,ftype,L,s_u,s_v,s_w,N(3)/(2*f3(3)))

     if (stretch) then
         f3sq = 1.0
         if (sqSt(1) .gt. 1) f3sq(1) = sqSt(1)
         if (sqSt(2) .gt. 1) f3sq(2) = sqSt(2)
         if (sqSt(3) .gt. 1) f3sq(3) = sqSt(3)
         call stretchSqueeze(N,L,sqSt,f3sq,s_u,s_v,s_w)
         stop
     end if

     ! FV filter
     if(FV) then
        do i = 1,3
           if(f3(i).gt.1) call lopass(N,f3,L,i,s_u,s_v,s_w)
        enddo
     endif

     if(OutField) then
        call print_realField(N,L,f3,s_u,s_v,s_w)
        stop 
     endif

     ! down sample?
     if(DS) then

        if(f3(1).gt.1 .OR. f3(2).gt.1 .OR. f3(3).gt.1) then
           call sample_co(N,f3,L,s_u,s_v,s_w,s_uo,s_vo,s_wo)
           !call sample_stag(N,f3,L,s_u,s_v,s_w,s_uo,s_vo,s_wo)
        else
           s_uo = s_u
           s_vo = s_v
           s_wo = s_w
        endif
     
        if (diss_flag == 1) then
           call calc_dissipation(N/f3,L,f3,nu,s_uo,s_vo,s_wo,diss)
           write(*,*) "Calculated dissipation: ", diss
        else if (diss_flag == 2) then
           write(*,*) "Dissipation manually set to: ", diss
        endif

        ! calculate spectrum for down-sampled field
        ones = 1
        call spectrum_3D(N/f3,Nmax,ones,ftype,L,s_uo,s_vo,s_wo,spec,nu,diss)
        call spectrum_1D(N/f3,Nmax,ones,ftype,L,s_uo,s_vo,s_wo,spec_x,1)
        call spectrum_1D(N/f3,Nmax,ones,ftype,L,s_uo,s_vo,s_wo,spec_y,2)
        call spectrum_1D(N/f3,Nmax,ones,ftype,L,s_uo,s_vo,s_wo,spec_z,3)


     else

        if (diss_flag == 1) then
           call calc_dissipation(N,L,f3,nu,s_u,s_v,s_w,diss)
           write(*,*) "Calculated dissipation: ", diss
        else if (diss_flag == 2) then
           write(*,*) "Dissipation manually set to: ", diss
        endif

        ! if not down-sampling...
        call spectrum_3D(N,Nmax,f3,ftype,L,s_u,s_v,s_w,spec,nu,diss)
        call spectrum_1D(N,Nmax,f3,ftype,L,s_u,s_v,s_w,spec_x,1)
        call spectrum_1D(N,Nmax,f3,ftype,L,s_u,s_v,s_w,spec_y,2)
        call spectrum_1D(N,Nmax,f3,ftype,L,s_u,s_v,s_w,spec_z,3)

        if (it.eq.it_ini) then
            call save_aniso_stats(N/f3,Nmax,L,f3,nu,s_u,s_v,s_w,it,t,.true.)
        else
            call save_aniso_stats(N/f3,Nmax,L,f3,nu,s_u,s_v,s_w,it,t,.false.)
        endif

     endif

     ! calculate structure function for field
     !call structure(N,Nmax,u,v,w,strct)

!!$     if (it.eq.it_ini) then
!!$        call save_aniso_spectral(N/f3,Nmax,s_uo,s_vo,s_wo,t,"./results/"//input_dir,.true.)
!!$     else
!!$        call save_aniso_spectral(N/f3,Nmax,s_uo,s_vo,s_wo,t,"./results/"//input_dir,.false.)
!!$     end if
!!$
!!$     call save_aniso_stats(N/f3,Nmax,L,f3,nu,s_uo,s_vo,s_wo,"./results/"//input_dir,it,t)
!!$ 

  enddo


  if(DS) then
     ! save stats, down-sampled!
     call save_spectrum(N/f3,Nmax,ones,L,nu,diss,spec,spec_x,spec_y,spec_z,samples,output_prefix)

  else
     ! not down-sampled
     call save_spectrum(N,Nmax,f3,L,nu,diss,spec,spec_x,spec_y,spec_z,samples,output_prefix)

  endif

  call lind_spec1D(N,Nmax,L,s_u,s_v,s_w,output_prefix)

  ! strct, update later...
  !call save_structure(N,Nmax,f3,nu,diss,strct,samples,input_dir)

  ! clean up
  call cleanup()
  call MPI_Finalize(ierr)

end program main


!--------------------------------------------------------------------------------
subroutine get_input()
use globals
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

  
  integer        :: fcount
  character(100) :: buffer, flag, ifile
  integer, parameter :: input_unit = 10

  ! Defaults:
  N = 32
  it_ini = 0
  it_fin = 1000
  it_step = 100
  ifile = "./input"
  pad = 0

  ! get input file
  fcount = 1
  do while (fcount .lt. IARGC())
     call getarg(fcount,buffer)
     read (buffer, *) flag
     if (flag(1:1) =='-') then
        if (flag(2:2) == 'i') then
           call getarg(fcount+1,buffer)
           ifile = buffer
           fcount=fcount+2
        else
           print *, "Bad option: use (i) to indicate input file"
           stop
        endif
     else
        print *, "Bad Flag: use (-) to indicate"
        stop
     endif
  enddo

  ! read inputs
  open (unit=input_unit,file=ifile,form="FORMATTED",status="OLD")
  read(input_unit,*)
  read(input_unit,'(A)') input_dir
  read(input_unit,*)
  read(input_unit,'(A)') output_prefix
  read(input_unit,*)
  read(input_unit,*) type
  read(input_unit,*)
  read(input_unit,*) type_fmt
  read(input_unit,*)
  read(input_unit,*) N(1),N(2),N(3)
  read(input_unit,*)
  read(input_unit,*) f3(1),f3(2),f3(3)
  read(input_unit,*)
  read(input_unit,*) L(1),L(2),L(3)
  read(input_unit,*)
  read(input_unit,*) FV
  read(input_unit,*)
  read(input_unit,*) DS
  read(input_unit,*)
  read(input_unit,*) ftype
  read(input_unit,*)
  read(input_unit,*) it_ini
  read(input_unit,*)
  read(input_unit,*) it_fin
  read(input_unit,*)
  read(input_unit,*) it_step
  read(input_unit,*)
  read(input_unit,*) diss_flag
  read(input_unit,*)
  read(input_unit,*) diss
  read(input_unit,*)
  read(input_unit,*) OutField
  read(input_unit,*)
  read(input_unit,*) stretch
  read(input_unit,*)
  read(input_unit,*) sqSt(1),sqSt(2),sqSt(3)

  close (unit=input_unit)

  print*, " "
  write(*,'(a,I4,1X,I4,1X,I4)') "  Grid size: ", N
  print*, " "
  delta = L*pi/dble(N)*dble(f3)

  if(ftype.ne.0 .AND. ftype.ne.1) then
     print*, " WARNING, bad 1D spectra type (setting to box)"
     ftype = 0
  endif

  return
end subroutine get_input


!--------------------------------------------------------------------------------
subroutine initialize()
use globals
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'


  integer :: Nmax


  Nmax = maxval(N)
  allocate( spec(0:Nmax-1,0:1) )
  allocate( spec_x(0:N(1)-1,0:1) )
  allocate( spec_y(0:N(2)-1,0:1) )
  allocate( spec_z(0:N(3)-1,0:1) )
  allocate( strct(0:Nmax/2) )
  spec = 0.0D0
  spec_x = 0.0D0
  spec_y = 0.0D0
  spec_z = 0.0D0
  strct = 0.0D0

  fieldsize = N(1)*N(2)*N(3)
  sfieldsize = ((N(1)+2)/2)*N(2)*N(3)
  fieldnorm = 1.0D0/dble(fieldsize)

  ! Allocate velocities
  allocate( u_temp(1-pad:N(1)+pad, 1-pad:N(2)+pad, 1-pad:N(3)+pad) )
  allocate( u(1:N(1),1:N(2),1:N(3)) )
  allocate( v(1:N(1),1:N(2),1:N(3)) )
  allocate( w(1:N(1),1:N(2),1:N(3)) )
  allocate( uo(1:N(1)/f3(1),1:N(2)/f3(2),1:N(3)/f3(3)) )
  allocate( s_u(0:N(1)/2,0:N(2)-1,0:N(3)-1) )
  allocate( s_v(0:N(1)/2,0:N(2)-1,0:N(3)-1) )
  allocate( s_w(0:N(1)/2,0:N(2)-1,0:N(3)-1) )
  allocate( s_uo(0:N(1)/(2*f3(1)),0:N(2)/f3(2)-1,0:N(3)/f3(3)-1) )
  allocate( s_vo(0:N(1)/(2*f3(1)),0:N(2)/f3(2)-1,0:N(3)/f3(3)-1) )
  allocate( s_wo(0:N(1)/(2*f3(1)),0:N(2)/f3(2)-1,0:N(3)/f3(3)-1) )

  ! create plans
  call dfftw_plan_dft_r2c_3d(plan_r2c,N(1),N(2),N(3), u, s_u, FFTW_ESTIMATE)
  call dfftw_plan_dft_r2c_3d(plan_r2c2,N(1)/f3(1),N(2)/f3(2),N(3)/f3(3), uo, s_uo, FFTW_ESTIMATE)
  call dfftw_plan_dft_c2r_3d(plan_c2r,N(1),N(2),N(3), s_u, u, FFTW_ESTIMATE)
  call dfftw_plan_dft_c2r_3d(plan_c2r2,N(1)/f3(1),N(2)/f3(2),N(3)/f3(3), s_uo, uo, FFTW_ESTIMATE)



end subroutine initialize


!--------------------------------------------------------------------------------
subroutine cleanup()
use globals
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

  deallocate(u,v,w,u_temp,uo)
  deallocate(s_u,s_v,s_w,s_uo)
  deallocate(spec, spec_x, spec_y, spec_z, strct)
  call dfftw_destroy_plan(plan_r2c)
  call dfftw_destroy_plan(plan_r2c2)
  call dfftw_destroy_plan(plan_c2r)
  call dfftw_destroy_plan(plan_c2r2)

end subroutine cleanup
