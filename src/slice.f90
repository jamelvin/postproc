
module slice

public write_slice


contains
  

!------------------------------------------------------------------------------
subroutine write_slice(N,Nmax,f3,ftype,L,u,v,w,z_slice)
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

  integer :: N(3), Nmax, f3(3), ftype
  integer :: count(0:Nmax-1), ii, jj, kk
  integer :: ix, iy, iz, ix2, iy2, iz2, z_slice
  real(8) :: kx, ky, kz, k2, L(3), delta(3), e_k(3), vol
  real(8) :: kx_max, ky_max, kz_max, k_mag, k_max
  real(8) :: e, tot_energy, ave_energy, k, area, max, check
  real(8) :: pi = 3.14159265
  complex(C_DOUBLE_COMPLEX) :: u(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX) :: v(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX) :: w(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX) :: f_u(0:N(1)/(2*f3(1)),0:N(2)/f3(2)-1,0:N(3)/f3(3)-1)
  complex(C_DOUBLE_COMPLEX) :: f_v(0:N(1)/(2*f3(1)),0:N(2)/f3(2)-1,0:N(3)/f3(3)-1)
  complex(C_DOUBLE_COMPLEX) :: f_w(0:N(1)/(2*f3(1)),0:N(2)/f3(2)-1,0:N(3)/f3(3)-1)
  real(C_DOUBLE) :: r_u(0:N(1)/f3(1)-1,0:N(2)/f3(2)-1,0:N(3)/f3(2)-1)
  real(C_DOUBLE) :: r_v(0:N(1)/f3(1)-1,0:N(2)/f3(2)-1,0:N(3)/f3(2)-1)
  real(C_DOUBLE) :: r_w(0:N(1)/f3(1)-1,0:N(2)/f3(2)-1,0:N(3)/f3(2)-1)
  type(C_PTR) :: plan_c2r_small
  

  write(*,*) "  ...extracting data slice"

  ! create plans for filtered field
  call dfftw_plan_dft_c2r_3d(plan_c2r_small,N(1)/f3(1),N(2)/f3(2),N(3)/f3(3), f_u, r_u, FFTW_ESTIMATE)

  ! initialize
  f_u = 0.0d0
  f_v = 0.0d0
  f_w = 0.0d0

  r_u = 0.0d0
  r_v = 0.0d0
  r_w = 0.0d0

  ! filtered grid spacing
  delta = L*pi/dble(N) * dble(f3)

  ! max wavenumber (can toss nyquist here)
  kx_max = 2.0D0/L(1) * dble(N(1)-1)/dble(2*f3(1))
  ky_max = 2.0D0/L(2) * dble(N(2)-1)/dble(2*f3(2))
  kz_max = 2.0D0/L(3) * dble(N(3)-1)/dble(2*f3(3))

  ! 3D energy spectrum
  do ix = 0, N(1)/2
     kx = dble(ix) * 2.0D0/L(1)

     do iy = 0, N(2)-1
        ky = dble(iy) * 2.0D0/L(2)
        if (iy.gt.N(2)/2) ky = dble(iy - N(2)) * 2.0D0/L(2)

        do iz = 0, N(3)-1
           kz = dble(iz) * 2.0D0/L(3)
           if (iz.gt.N(3)/2) kz = dble(iz - N(3)) * 2.0D0/L(3)


           if(ftype.eq.0) then

              ! sharp spectral box
              if( kx.gt.kx_max ) cycle
              if( abs(ky).gt.ky_max ) cycle
              if( abs(kz).gt.kz_max ) cycle

           elseif(ftype.eq.1) then

              ! ellipsoidal
              check = (kx/kx_max)**2 + (ky/ky_max)**2 + (kz/kz_max)**2
              if(check.gt.1.0D0) cycle

           endif

           ! filtered array index
           ix2 = ix
           iy2 = ky * L(2)/2.0d0 + N(2)/f3(2)
           iz2 = kz * L(3)/2.0d0 + N(3)/f3(3)

           ! copy filtered field
           f_u(ix2,iy2,iz2) = u(ix,iy,iz)
           f_v(ix2,iy2,iz2) = v(ix,iy,iz)
           f_w(ix2,iy2,iz2) = w(ix,iy,iz)
 
        enddo
     enddo
  enddo

  ! transform fields to realspace
  call dfftw_execute_dft_c2r(plan_c2r_small,f_u,r_u)
  call dfftw_execute_dft_c2r(plan_c2r_small,f_v,r_v)
  call dfftw_execute_dft_c2r(plan_c2r_small,f_w,r_w)

  ! write slice
  open(unit=35, file="slice.dat", form="FORMATTED",status="UNKNOWN")
  kk = z_slice !N(3)/(2*f3(3))
  do ii = 1,N(1)/f3(1)
     do jj = 1,N(2)/f3(2)
        write(35,'(E16.7,1X,E16.7,1X,E16.7)') r_u(ii,jj,kk), r_v(ii,jj,kk), r_w(ii,jj,kk) 
     enddo
  enddo

  ! clean up
  close(35)
  call dfftw_destroy_plan(plan_c2r_small)


end subroutine

end module slice
!------------------------------------------------------------------------------






