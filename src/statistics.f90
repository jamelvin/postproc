

module statistics

public spectrum_1D, structure, save_spectrum, save_structure, lopass, &
       spectrum_3D, sample_co, sample_stag, lopass_onepass


contains

!------------------------------------------------------------------------------
subroutine lopass(N,f3,L,ind,phi1,phi2,phi3)
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

  integer :: N(3), f3(3), ind
  integer :: ix, iy, iz, i
  real(8) :: kx, ky, kz, L(3), delta, e_k(3)
  real(8) :: k(3), filter, arg, pw
  real(8) :: pi = 3.14159265
  complex(C_DOUBLE_COMPLEX), INTENT(INOUT) :: phi1(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(INOUT) :: phi2(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(INOUT) :: phi3(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  

  write(*,'(a,I1,a)') "   ...applying k",ind," FV filter"

  ! filtered grid size, N/2 * 2pi/N*f = pi*f : if f>1, filter becomes negative...
  delta = L(ind)*pi/dble(N(ind)) * dble(f3(ind))


  ! 1D filter
  do ix = 0, N(1)/2
     k(1) = dble(ix) * 2.0D0/L(1)
     do iy = 0, N(2)-1
        k(2) = dble(iy) * 2.0D0/L(2)
        if (iy.gt.N(2)/2) k(2) = dble(iy - N(2)) * 2.0D0/L(2)
        do iz = 0, N(3)-1
           k(3) = dble(iz) * 2.0D0/L(3)
           if (iz.gt.N(3)/2) k(3) = dble(iz - N(3)) * 2.0D0/L(3)

           ! real box filter, k_max*delta = N/(2f)*2pi/N*f=pi
           filter = 1.0D0
           arg = 0.5D0*k(ind)*delta
           if(abs(arg).gt.0.0D0) filter = sin(arg)/arg
           !this is wrong... if(arg.gt.0.0D0) filter = sin(arg)/arg

           phi1(ix,iy,iz) = filter * phi1(ix,iy,iz)
           phi2(ix,iy,iz) = filter * phi2(ix,iy,iz)
           phi3(ix,iy,iz) = filter * phi3(ix,iy,iz)

        enddo
     enddo
  enddo

return

end subroutine lopass


!------------------------------------------------------------------------------
subroutine lopass_onepass(N,f3,L,phi1,phi2,phi3)
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

  integer :: N(3), f3(3), ind
  integer :: ix, iy, iz, i
  real(8) :: kx, ky, kz, L(3), delta(3), e_k(3)
  real(8) :: k(3), filter, arg, pw
  real(8) :: pi = 3.14159265
  complex(C_DOUBLE_COMPLEX), INTENT(INOUT) :: phi1(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(INOUT) :: phi2(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(INOUT) :: phi3(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  

  write(*,'(a,I1,a)') "   ...applying k",ind," filter"

  ! filtered grid size, N/2 * 2pi/N*f = pi*f : if f>1, filter becomes negative...
  delta = L*pi/dble(N) * dble(f3)


  ! 1D filter
  do ix = 0, N(1)/2
     k(1) = dble(ix) * 2.0D0/L(1)
     do iy = 0, N(2)-1
        k(2) = dble(iy) * 2.0D0/L(2)
        if (iy.gt.N(2)/2) k(2) = dble(iy - N(2)) * 2.0D0/L(2)
        do iz = 0, N(3)-1
           k(3) = dble(iz) * 2.0D0/L(3)
           if (iz.gt.N(3)/2) k(3) = dble(iz - N(3)) * 2.0D0/L(3)

           ! real box filter, k_max*delta = N/(2f)*2pi/N*f=pi
           filter = 1.0D0
!           arg = 0.5D0*2.0D0/3.0D0*dot_product(k,delta)
           arg = 0.5D0*1.0D0/3.0D0*sqrt((k(1)*delta(1))**2+(k(2)*delta(2))**2+(k(3)*delta(3))**2)
           if(abs(arg).gt.0.0D0) filter = sin(arg)/arg
           !this is wrong... if(arg.gt.0.0D0) filter = sin(arg)/arg

           phi1(ix,iy,iz) = filter * phi1(ix,iy,iz)
           phi2(ix,iy,iz) = filter * phi2(ix,iy,iz)
           phi3(ix,iy,iz) = filter * phi3(ix,iy,iz)

        enddo
     enddo
  enddo

return

end subroutine lopass_onepass


!------------------------------------------------------------------------------
subroutine sample_co(N,f3,L,phi1,phi2,phi3,phi1o,phi2o,phi3o)
use, intrinsic :: iso_c_binding
use globals, only : plan_c2r, plan_r2c2
implicit none
include 'fftw3.f03'

  integer :: N(3), f3(3), ind
  integer :: ix, iy, iz, i, j, k
  real(8) :: kx, ky, kz, L(3), delta, e_k(3)
  real(8) :: filter, arg, fnorm
  real(8) :: pi = 3.14159265
  complex(C_DOUBLE_COMPLEX), INTENT(IN) :: phi1(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(IN) :: phi2(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(IN) :: phi3(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(OUT) :: phi1o(0:N(1)/(2*f3(1)),0:N(2)/f3(2)-1,0:N(3)/f3(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(OUT) :: phi2o(0:N(1)/(2*f3(1)),0:N(2)/f3(2)-1,0:N(3)/f3(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(OUT) :: phi3o(0:N(1)/(2*f3(1)),0:N(2)/f3(2)-1,0:N(3)/f3(3)-1)
  real(C_DOUBLE), dimension(1:N(1),1:N(2),1:N(3)) :: u,v,w
  real(C_DOUBLE), dimension(1:N(1)/f3(1),1:N(2)/f3(2),1:N(3)/f3(3)) :: uo,vo,wo
  

  write(*,'(a,I1,a)') "   ...down sampling to a colocated grid"

  ! transform to realspace
  call dfftw_execute_dft_c2r(plan_c2r,phi1, u)
  call dfftw_execute_dft_c2r(plan_c2r,phi2, v)
  call dfftw_execute_dft_c2r(plan_c2r,phi3, w)

  ! sample
  do i = 1, N(1)/f3(1)
     do j = 1, N(2)/f3(2)
        do k = 1, N(3)/f3(3)

           ix = (i-1)*f3(1)+1
           iy = (j-1)*f3(2)+1
           iz = (k-1)*f3(3)+1

           uo(i,j,k) = u(ix,iy,iz)
           vo(i,j,k) = v(ix,iy,iz)
           wo(i,j,k) = w(ix,iy,iz)

        enddo
     enddo
  enddo

  ! transform back to wavespace
  call dfftw_execute_dft_r2c(plan_r2c2,uo, phi1o)
  call dfftw_execute_dft_r2c(plan_r2c2,vo, phi2o)
  call dfftw_execute_dft_r2c(plan_r2c2,wo, phi3o)
  fnorm = 1.0D0/dble(N(1)/f3(1)*N(2)/f3(2)*N(3)/f3(3))
  phi1o = phi1o*fnorm
  phi2o = phi2o*fnorm
  phi3o = phi3o*fnorm

return

end subroutine sample_co
  
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine sample_stag(N,f3,L,phi1,phi2,phi3,phi1o,phi2o,phi3o)
use, intrinsic :: iso_c_binding
use globals, only : plan_c2r, plan_r2c2
implicit none
include 'fftw3.f03'

  integer :: N(3), f3(3), ind
  integer :: ix, iy, iz, i, j, k
  real(8) :: kx, ky, kz, L(3), delta, e_k(3)
  real(8) :: filter, arg, fnorm
  real(8) :: pi = 3.14159265
  complex(C_DOUBLE_COMPLEX), INTENT(IN) :: phi1(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(IN) :: phi2(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(IN) :: phi3(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(OUT) :: phi1o(0:N(1)/(2*f3(1)),0:N(2)/f3(2)-1,0:N(3)/f3(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(OUT) :: phi2o(0:N(1)/(2*f3(1)),0:N(2)/f3(2)-1,0:N(3)/f3(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(OUT) :: phi3o(0:N(1)/(2*f3(1)),0:N(2)/f3(2)-1,0:N(3)/f3(3)-1)
  real(C_DOUBLE), dimension(1:N(1),1:N(2),1:N(3)) :: u,v,w
  real(C_DOUBLE), dimension(1:N(1)/f3(1),1:N(2)/f3(2),1:N(3)/f3(3)) :: uo,vo,wo
  

  write(*,'(a,I1,a)') "   ...down sampling to a staggered grid"

  ! transform to realspace
  call dfftw_execute_dft_c2r(plan_c2r,phi1, u)
  call dfftw_execute_dft_c2r(plan_c2r,phi2, v)
  call dfftw_execute_dft_c2r(plan_c2r,phi3, w)

  ! sample
  do i = 1, N(1)/f3(1)
     do j = 1, N(2)/f3(2)
        do k = 1, N(3)/f3(3)

           ix = (i-1)*f3(1)+1
           iy = (j-1)*f3(2)+1
           iz = (k-1)*f3(3)+1

           uo(i,j,k) = u(ix,modulo(iy+f3(2)-1,N(2))+1,modulo(iz+f3(3)-1,N(3))+1)
           vo(i,j,k) = v(modulo(ix+f3(1)-1,N(1))+1,iy,modulo(iz+f3(3)-1,N(3))+1)
           wo(i,j,k) = w(modulo(ix+f3(1)-1,N(1))+1,modulo(iy+f3(2)-1,N(2))+1,iz)

        enddo
     enddo
  enddo

  ! transform back to wavespace
  call dfftw_execute_dft_r2c(plan_r2c2,uo, phi1o)
  call dfftw_execute_dft_r2c(plan_r2c2,vo, phi2o)
  call dfftw_execute_dft_r2c(plan_r2c2,wo, phi3o)
  fnorm = 1.0D0/dble(N(1)/f3(1)*N(2)/f3(2)*N(3)/f3(3))
  phi1o = phi1o*fnorm
  phi2o = phi2o*fnorm
  phi3o = phi3o*fnorm

return

end subroutine sample_stag
  

!------------------------------------------------------------------------------
subroutine spectrum_3D(N,Nmax,f3,ftype,L,u,v,w,spec,nu,diss)
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

  integer :: N(3), Nmax, f3(3), ftype
  integer :: count(0:Nmax-1)
  integer :: ix, iy, iz, weight, bin, i
  real(8) :: kx, ky, kz, k2, L(3), delta(3), e_k(3), vol
  real(8) :: kx_max, ky_max, kz_max, k_mag, k_max
  real(8) :: e, tot_energy, ave_energy, k, area, max, check
  real(8) :: diss, rms_u, taylor, Re_taylor, L_int, E_int, nu, C_ls
  real(8) :: pi = 3.14159265
  real(8), dimension(0:Nmax-1)                  :: ave_k, energy
  real(8), dimension(0:Nmax-1,0:1), INTENT(OUT) :: spec
  complex(C_DOUBLE_COMPLEX)        :: u(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX)        :: v(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX)        :: w(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  

  write(*,*) "  ...calculating 3d spectrum"

  ! initialize
  ave_k = 0.0D0
  energy = 0.0D0
  count = 0
  tot_energy = 0.0D0
  ave_energy = 0.0D0

  ! filtered grid spacing
  delta = L*pi/dble(N) * dble(f3)

  ! max wavenumber (can toss nyquist here)
  !kx_max = 2.0D0/L(1) * dble(N(1)-1)/dble(2*f3(1))
  !ky_max = 2.0D0/L(2) * dble(N(2)-1)/dble(2*f3(2))
  !kz_max = 2.0D0/L(3) * dble(N(3)-1)/dble(2*f3(3))
  kx_max = 2.0D0/L(1) * (dble(N(1))/dble(2*f3(1)) - 1.0d0)
  ky_max = 2.0D0/L(2) * (dble(N(2))/dble(2*f3(2)) - 1.0d0)
  kz_max = 2.0D0/L(3) * (dble(N(3))/dble(2*f3(3)) - 1.0d0)


  ! 3D energy spectrum
  do ix = 0, N(1)/2
     kx = dble(ix) * 2.0D0/L(1)
     weight = 2
     if (kx.eq.0 .or. kx.eq.dble(N(1))/2.0D0) weight = 1
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

           e = real(u(ix,iy,iz))*real(u(ix,iy,iz)) + imag(u(ix,iy,iz))*imag(u(ix,iy,iz)) + &
               real(v(ix,iy,iz))*real(v(ix,iy,iz)) + imag(v(ix,iy,iz))*imag(v(ix,iy,iz)) + &
               real(w(ix,iy,iz))*real(w(ix,iy,iz)) + imag(w(ix,iy,iz))*imag(w(ix,iy,iz))
           e = 0.5D0*e

           k_mag = sqrt(kx*kx+ky*ky+kz*kz)
           bin = nint(k_mag)
           ave_k(bin) = ave_k(bin) + k_mag*dble(weight)
           area = 4.0D0*pi*k_mag**2 ! area/count gives average dA later
           energy(bin) = energy(bin) + area * e * dble(weight)
           count(bin) = count(bin) + weight 
           tot_energy = tot_energy + e*dble(weight);

        enddo
     enddo
  enddo

!  vol = 4.0d0/3.0d0 * pi * kx_max*ky_max*kz_max

  ! really dV
  vol = 2.0d0/L(1) * 2.0d0/L(2) * 2.0d0/L(3)
  !ave_energy = ave_energy + tot_energy

  do i = 0, Nmax-1
     if (count(i).gt.0) then

        spec(i,0) = spec(i,0) + ave_k(i)/dble(count(i))

        ! energy per bin
        spec(i,1) = spec(i,1) + energy(i)/(dble(count(i))*vol)

     else
        spec(i,0) = 0.0D0
        spec(i,1) = 0.0D0
     endif
  enddo

  max = 0.0D0;
  do i = 0, Nmax-1
     if (count(i).gt.0) then
        if (ave_k(i)/count(i).gt.max) max = ave_k(i)/dble(count(i))
     endif
  enddo

  ! report
  E_int = 0.0d0
  do i = 0, Nmax-2
     if (count(i).gt.0) then
        E_int =  E_int + (spec(i,1)+spec(i+1,1))/(spec(i,0)+spec(i+1,0)) * (spec(i+1,0)-spec(i,0))
     endif
  enddo
  !E_int = E_int + spec(0,1)/spec(0,0) * spec(1,0)*0.5d0
  
  rms_u = sqrt(2.0d0/3.0d0*tot_energy)
  taylor = sqrt(15.0d0*nu*rms_u**2/diss)
  Re_taylor = rms_u*taylor/nu
  L_int = pi/(2.0d0*rms_u**2)*E_int
  C_ls = diss * L_int/(tot_energy**1.5d0)
     
  write(*,*) " "
  print *, "   - Total energy = ", tot_energy
  print *, "   - u rms = ", rms_u
  print *, "   - Taylor lengthscale = ", taylor
  print *, "   - Re_taylor = ", Re_taylor
  print *, "   - Integral lengthscale = ", L_int
  print *, "   - 15*C_ls = ", 15.0d0*C_ls
  write(*,*) " "

  open(unit=51, file="stats_basic.txt", form="FORMATTED",access="APPEND",status="UNKNOWN")
  write(51,'(F15.8,1X,F15.8,1X,F15.8,1X,F15.8,1X,F15.8)') tot_energy, rms_u, taylor, Re_taylor, L_int
  close(51)

return

end subroutine spectrum_3D


!------------------------------------------------------------------------------
subroutine spectrum_1D(N,Nmax,f3,ftype,L,u,v,w,spec_dir,dir)
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

  integer :: N(3), Nmax, f3(3), dir, ftype
  integer :: count(0:Nmax-1)
  integer :: ix, iy, iz, weight, bin, i, indp, indm
  real(8) :: k(3), k2, L(3), delta(3), e_k(3)
  real(8) :: k3_max(3), k_mag, k_max, check
  real(8) :: e, tot_energy, ave_energy, area, vol, factor
  real(8) :: pi = 3.14159265
  real(8), dimension(0:Nmax-1)                    :: ave_k, energy
  real(8), dimension(0:N(dir)-1,0:1), INTENT(OUT) :: spec_dir
  complex(C_DOUBLE_COMPLEX)        :: u(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX)        :: v(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX)        :: w(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  


  write(*,'(a,I1)') "  ...calculating 1d spectrum: ", dir

  ! intialize
  ave_k = 0.0D0
  energy = 0.0D0
  count = 0
  tot_energy = 0.0D0
  ave_energy = 0.0D0

  ! filtered grid spacing
  delta = L*pi/dble(N) * dble(f3)

  ! max wavenumber (can toss nyquist here)
!  k3_max(1) = 2.0D0/L(1) * dble(N(1)-1)/dble(2*f3(1))
!  k3_max(2) = 2.0D0/L(2) * dble(N(2)-1)/dble(2*f3(2))
!  k3_max(3) = 2.0D0/L(3) * dble(N(3)-1)/dble(2*f3(3))
  k3_max(1) = 2.0D0/L(1) * (dble(N(1))/dble(2*f3(1)) - 1.0d0)
  k3_max(2) = 2.0D0/L(2) * (dble(N(2))/dble(2*f3(2)) - 1.0d0)
  k3_max(3) = 2.0D0/L(3) * (dble(N(3))/dble(2*f3(3)) - 1.0d0)


  ! other indices
  indp = modulo(dir,3)+1
  indm = modulo(dir-2,3)+1

  ! 1D energy spectrum 
  do ix = 0, N(1)/2
     k(1) = dble(ix) * 2.0D0/L(1)
     weight = 2
     if (k(1).eq.0 .OR. k(1).eq.dble(N(1)/2.0D0)) weight = 1

     do iy = 0, N(2)-1
        k(2) = dble(iy) * 2.0D0/L(2)
        if (iy.gt.N(2)/2) k(2) = dble(iy - N(2)) * 2.0D0/L(2)

        do iz = 0, N(3)-1
           k(3) = dble(iz) * 2.0D0/L(3)
           if (iz.gt.N(3)/2) k(3) = dble(iz - N(3)) * 2.0D0/L(3)

           k_mag = sqrt(dot_product(k,k))


           ! area/count gives average dA later
           if(ftype.eq.0) then

              ! sharp spectral box (if down-sampled, this does nothing)
              if( k(1).gt.k3_max(1) ) cycle
              if( abs(k(2)).gt.k3_max(2) ) cycle
              if( abs(k(3)).gt.k3_max(3) ) cycle

              ! rectangle
              area = k3_max(indp) * k3_max(indm)
              area = 4.0d0 * area

           elseif(ftype.eq.1) then

              ! spherical/ellipsoidal filter
              check = (k(1)/k3_max(1))**2 + (k(2)/k3_max(2))**2 + (k(3)/k3_max(3))**2
              if(check.gt.1.0D0) cycle

              ! ellipse
              area = pi * k3_max(indp) * k3_max(indm) * sqrt(1.0D0-(k(dir)/k3_max(dir))**2)

           endif

           e = real(u(ix,iy,iz))*real(u(ix,iy,iz)) + imag(u(ix,iy,iz))*imag(u(ix,iy,iz)) + &
               real(v(ix,iy,iz))*real(v(ix,iy,iz)) + imag(v(ix,iy,iz))*imag(v(ix,iy,iz)) + &
               real(w(ix,iy,iz))*real(w(ix,iy,iz)) + imag(w(ix,iy,iz))*imag(w(ix,iy,iz))
           e = 0.50D0*e

           bin = nint(abs(k(dir))*L(dir)/2.0d0)
           ave_k(bin) = ave_k(bin) + abs(k(dir))*dble(weight)
           energy(bin) = energy(bin) + area * e * dble(weight)

           ! Eii(ki) (alternate, see Pope)
           !factor = 1.0D0/(pi*k_mag**2) * (1.0D0-(k(dir)/k_mag)**2) * 2.0D0/L(indp) * 2.0D0/L(indm) !* area !k3_max(indp) * k3_max(indm)
           !energy(bin) = energy(bin) + factor * e

           count(bin) = count(bin) + weight 

        enddo
     enddo
  enddo

  !if(ftype.eq.0) then
  !   vol = k3_max(1)*k3_max(2)*k3_max(3)
  !elseif(ftype.eq.1) then
  !   vol = 4.0d0/3.0d0 * pi * k3_max(1)*k3_max(2)*k3_max(3)
  !endif

  ! really dV
  vol = 2.0d0/L(1) * 2.0d0/L(2) * 2.0d0/L(3)

  !print*,"  Max wavenumbers:", k3_max
  !print*,"  Wavespace volume:", vol

  do i = 0, N(dir)-1
     if (count(i).gt.0) then

        spec_dir(i,0) = spec_dir(i,0) + ave_k(i) / dble(count(i))

        ! energy per bin
        spec_dir(i,1) = spec_dir(i,1) + energy(i) / (dble(count(i))*vol)

        !print*,spec_dir(i,0),count(i)

     else
        spec_dir(i,0) = 0.0D0
        spec_dir(i,1) = 0.0D0
     endif


  enddo



return

end subroutine spectrum_1D


!------------------------------------------------------------------------------
subroutine lind_spec1D(N,Nmax,L,u,v,w,prefix)
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

  integer :: N(3), Nmax, ix, iy, iz, weight, i, binx, biny, binz
  real(8) :: k(3), L(3), e
  real(8), dimension(0:(Nmax/2)-1) :: energyX, energyY, energyZ, energyH, energyV
  integer, parameter :: spec_output_h  = 193
  integer, parameter :: spec_output_v  = 194
  character(100)     :: spec_name_h, spec_name_v
  character(len=*)   :: prefix
  complex(C_DOUBLE_COMPLEX)        :: u(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX)        :: v(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX)        :: w(0:N(1)/2,0:N(2)-1,0:N(3)-1)

  write(*,'(a)') "  ...calculating 1d spectrums via Lindborg"

  ! intialize
  energyX = 0.0D0
  energyY = 0.0D0
  energyZ = 0.0D0

  ! 1D energy spectrum 
  do ix = 0, N(1)/2
     weight = 2
     k(1) = dble(ix) * 2.0D0/L(1)
     if ((ix .eq. 0) .OR. (ix .eq. N(1)/2)) weight = 1

     do iy = 0, N(2)-1
        k(2) = dble(iy) * 2.0D0/L(2)
        if (iy.gt.N(2)/2) k(2) = dble(iy - N(2)) * 2.0D0/L(2)

        do iz = 0, N(3)-1
           k(3) = dble(iz) * 2.0D0/L(3)
           if (iz.gt.N(3)/2) k(3) = dble(iz - N(3)) * 2.0D0/L(3)
           
           e = real(u(ix,iy,iz))*real(u(ix,iy,iz)) + imag(u(ix,iy,iz))*imag(u(ix,iy,iz)) + &
               real(v(ix,iy,iz))*real(v(ix,iy,iz)) + imag(v(ix,iy,iz))*imag(v(ix,iy,iz)) + &
               real(w(ix,iy,iz))*real(w(ix,iy,iz)) + imag(w(ix,iy,iz))*imag(w(ix,iy,iz))
           e = 0.5D0*e

           binx = nint(abs(k(1))*L(1)/2.0D0)
           biny = nint(abs(k(2))*L(2)/2.0D0)
           binz = nint(abs(k(3))*L(3)/2.0D0)

           energyX(binx) = energyX(binx) + e * dble(weight)
           energyY(biny) = energyY(biny) + e * dble(weight)
           energyZ(binz) = energyZ(binz) + e * dble(weight)
        enddo
     enddo
  enddo

  spec_name_h = trim(adjustl(prefix))//"_h.spec"
  spec_name_v = trim(adjustl(prefix))//"_v.spec"
  open(unit=spec_output_h, file=trim(adjustl(spec_name_h)), form="FORMATTED",status="UNKNOWN")
  open(unit=spec_output_v, file=trim(adjustl(spec_name_v)), form="FORMATTED",status="UNKNOWN")
  print *, " Saving energy spectrum to: ", spec_name_h

  ! This assumes that the vertical (stratified) direction is y and also that the resolution and box
  ! size is the same in the horizontal directions (all of this should be the case)
  do i = 0, (Nmax/2) - 1
     energyH = L(1)/2.0D0 * energyX(i) + L(3)/2.0D0 * energyZ(i)
     energyV = L(2)/2.0D0 * energyY(i)
     write(spec_output_h,"(2E13.4)") dble(i) * 2.0D0/L(2), energyH
     write(spec_output_v,"(2E13.4)") dble(i) * 2.0D0/L(2), energyV
  end do

  close(spec_output_h)
  close(spec_output_v)

return

end subroutine lind_spec1D

!------------------------------------------------------------------------------
subroutine structure(N,Nmax,f3,L,u,v,w,strct)
implicit none

  integer :: N(3), Nmax, f3(3)
  integer :: ix, iy, iz, sep, factor, i
  real(8) :: max, L(3), delta(3)
  real(8) :: pi = 3.14159265
  real(8), dimension(0:Nmax/2) :: strct, sum3, num
  real(8)    :: u(0:N(1)-1,0:N(2)-1,0:N(3)-1)
  real(8)    :: v(0:N(1)-1,0:N(2)-1,0:N(3)-1)
  real(8)    :: w(0:N(1)-1,0:N(2)-1,0:N(3)-1)


  factor = 1
  num = 0.0D0
  sum3 = 0.0D0
  delta = L*pi*dble(f3)/dble(N)

  ! 3rd-order structure function
  do ix = 0, N(1)-1
     do iy = 0, N(2)-1
        do iz = 0, N(3)-1
           do sep = -N(1)/2, N(1)/2 ! should be min N for aniso grids
              i = abs(sep)
              factor = 1
              if (sep.lt.0) factor = -1
              sum3(i) = sum3(i) + dble(factor) * ( (u(ix-sep,iy,iz)-u(ix,iy,iz))**3 + &
                                                   (v(ix,iy-sep,iz)-v(ix,iy,iz))**3 + &
                                                   (w(ix,iy,iz-sep)-w(ix,iy,iz))**3 )
              num(i) = num(i) + 3.0D0
           enddo
        enddo
     enddo
  enddo

  do i = 0, Nmax/2 !N(1)/2
     strct(i) = strct(i) + sum3(i)/num(i)
  enddo

return

end subroutine structure


!------------------------------------------------------------------------------
    subroutine save_aniso_spectral(N,Nmax,s_u,s_v,s_w,t,prefix,flag)
        use, intrinsic :: iso_c_binding
        implicit none
        !Inputs/Outputs
        integer, intent(in) :: N(3), Nmax
        real(8), intent(in) :: t
        complex(C_DOUBLE_COMPLEX)   :: s_u(0:N(1)/2,0:N(2)-1,0:N(3)-1)
        complex(C_DOUBLE_COMPLEX)   :: s_v(0:N(1)/2,0:N(2)-1,0:N(3)-1)
        complex(C_DOUBLE_COMPLEX)   :: s_w(0:N(1)/2,0:N(2)-1,0:N(3)-1)
        logical, intent(in) :: flag
        character(50), intent(in)   :: prefix
        !Dummy variables
        integer :: ix, iy, iz, i, j
        integer :: weight
        real(8) :: delta(3,3), trace_tau
        complex(C_DOUBLE_COMPLEX)   :: s_vel(0:N(1)/2,0:N(2)-1,0:N(3)-1,3)

        !Fourier coefficients of two-point velocity correlation:
        complex(C_DOUBLE_COMPLEX)   :: R_ij(0:N(1)/2,0:N(2)-1,0:N(3)-1,3,3)
        !Reynolds stress tensor:
        complex(C_DOUBLE_COMPLEX)   :: tau(3,3)
        real(8), dimension(3,3)   :: a, b, b2, b3
        real(8)    :: II_b, III_b, eta, ksi, variance
        !Saving variables:
        character(100)      :: aniso_name
        integer, parameter  :: aniso_output = 55

        s_vel(:,:,:,1) = s_u
        s_vel(:,:,:,2) = s_v
        s_vel(:,:,:,3) = s_w

        tau = dcmplx(0.0D0, 0.0D0)
        R_ij = dcmplx(0.0D0, 0.0D0)
        ! 3D energy spectrum
        do ix = 0, N(1)/2
            weight = 2
            if (ix.eq.0) weight = 1  ! old version had kx=N/2 weight also at 1 but, this doesn't make sense
            do iy = 0, N(2)-1
                do iz = 0, N(3)-1

                    do i = 1, 3
                        do j = 1,3
                            R_ij(ix,iy,iz,i,j) = dconjg(s_vel(ix,iy,iz,i))*s_vel(ix,iy,iz,j)
                            tau(i,j) = tau(i,j) + weight*R_ij(ix,iy,iz,i,j)
                        enddo
                    enddo

                enddo
            enddo
        enddo

        !tau = real(tau,8)
        trace_tau = (tau(1,1) + tau(2,2) + tau(3,3))

        delta = 0.0D0
        delta(1,1) = 1.0D0
        delta(2,2) = 1.0D0
        delta(3,3) = 1.0D0

        !Reynolds stress anisotropy:
        a = tau - trace_tau*delta/(3.0D0)

        !Normalized reynolds stress anisotropy and its powers
        b = a/trace_tau
        b2 = matmul(b,b)
        b3 = matmul(b,b2)

        !Invariants of normalized reynolds stress anisotropy
        !    write(*,*) " Trace of tensor: ", (a(1,1) + a(2,2) + a(3,3))
        II_b = -0.5D0*(b2(1,1) + b2(2,2) + b2(3,3))
        III_b = 1.0d0/3.0d0*(b3(1,1)+b3(2,2)+b3(3,3))

        !Normalized invariants
        eta = sqrt(-II_b/3.0D0)
        ksi = (0.5D0*III_b)**(1.0D0/3.0D0)

        !Write the output to a file each time step:
        aniso_name = trim(adjustl(prefix))//".spectral.anis"
        if (flag) then
            open(unit=aniso_output, file=trim(adjustl(aniso_name)),form="FORMATTED",status="UNKNOWN")
            write(aniso_output,"(A,TR25,A,TR17,A,TR22,A,TR22,A,TR22,A)") "#            t","II_b","III_B","eta","ksi","kin_en"
        else
            open(unit=aniso_output, file=trim(adjustl(aniso_name)), access="APPEND",form="FORMATTED",status="OLD")
        end if
        !print *, "Saving anisotropy data to: ", aniso_name
        write(aniso_output,"(6E13.4)") t,II_b, III_b, eta, ksi, 0.5D0*trace_tau
        close(aniso_output)

    end subroutine

    !------------------------------------------------------------------------------
    subroutine save_aniso_stats(N,Nmax,L,f3,nu,s_u,s_v,s_w,timestep,time,create_new_file)
        use, intrinsic :: iso_c_binding
        use globals, only : pi, ftype, diss
        implicit none
        include 'mpif.h'

        !Inputs/Outputs
        integer, intent(in) :: N(3), Nmax, f3(3), timestep
        real(8), intent(in) :: nu, L(3), time
        complex(C_DOUBLE_COMPLEX)   :: s_u(0:N(1)/2,0:N(2)-1,0:N(3)-1)
        complex(C_DOUBLE_COMPLEX)   :: s_v(0:N(1)/2,0:N(2)-1,0:N(3)-1)
        complex(C_DOUBLE_COMPLEX)   :: s_w(0:N(1)/2,0:N(2)-1,0:N(3)-1)
        logical, intent(in) :: create_new_file

        !Temporary variables
        real(c_double)    :: epslonp(3),epslonr(3),var(3), varp(3)
        real(C_DOUBLE)    :: intsc1d(3),Eii_at_0(3),urms,ksq,Eii_at_0p(3)
        real(c_double)    :: ekkr,ekkp,uvar
        real(c_double)    :: tke, epslon
        real(c_double)    :: Lavg
        complex(C_DOUBLE_COMPLEX) :: Rij_p(3,3), Rij(3,3), Oij_p(3,3), Oij(3,3)
        real(C_DOUBLE), dimension(3,3)   :: a, b, b2, b3, dij
        real(C_DOUBLE)    :: II_b, III_b, II_v, III_v, eta, ksi, trace_tau
        real(C_DOUBLE)    :: delta(3), area
        real(C_DOUBLE)    :: k3_max(3), k(3)
        integer :: n_total,xsz_dummy,ysz_dummy
        integer, dimension(3) :: counter = 0
        complex(C_DOUBLE_COMPLEX) :: un_data(N(2),N(1)/2+1,N(3),3)
        integer :: x,y,z, i, j, m
        integer :: indp, indm
        integer :: ix, iy, iz
        integer :: ierr
        real(C_DOUBLE) :: xfactor
        real(C_DOUBLE) :: check
        character(len=3) :: status_flag

        n_total = (N(1)/2)*N(2)*N(3)
        xsz_dummy = N(1)/2
        ysz_dummy = N(2)

        delta = L*pi/dble(N) * dble(f3)
        ! toss nyquist
        k3_max(1) = 2.0D0/L(1) * dble(N(1)/2-1)/dble(f3(1))
        k3_max(2) = 2.0D0/L(2) * dble(N(2)/2-1)/dble(f3(2))
        k3_max(3) = 2.0D0/L(3) * dble(N(3)/2-1)/dble(f3(3))
        print*, "k max: ", k3_max

        ! Transform the data back into a POONGBACK-style field:
        un_data = 0d0
        do ix = 1,N(1)/2
            do iy = 1,N(2)
                do iz = 1,N(3)
                    un_data(iy,ix,iz,1) = s_u(ix-1,iy-1,iz-1)
                    un_data(iy,ix,iz,2) = s_v(ix-1,iy-1,iz-1)
                    un_data(iy,ix,iz,3) = s_w(ix-1,iy-1,iz-1)
                enddo
            enddo
        enddo

        !
        ! Root Mean Square Velocity from two pts. corr
        !
        ! first, we calculate R_ij(s) in wavespace
        !
        ! R_ij(s) = 2 * Sum_n <|c_n|^2> cos(w_n s)
        ! w/ w_n = 2*pi*n/N
        !
        ! Thus,
        ! R_ij(0) = 2 * Sum_n <|c_n|^2>
        !         = u'^2 (u rms squared)
        !
        ! dissipation = viscosity * k^2 * |c_n|^2
        !
        ekkp = 0d0
        Eii_at_0p = 0d0
        Rij = 0d0
        Rij_p = 0d0
        counter = 0

        do i = 1,3
            indp = modulo(i,3)+1
            indm = modulo(i-2,3)+1
            varp(i)=0d0
            epslonp(i)=0d0

            do j=1,n_total
                ! Indices start at 1 (not 0)
                ! This is different than other functions
                y = mod(int(j-1),ysz_dummy) + 1
                x = mod(int(j-1)/ysz_dummy ,xsz_dummy) + 1
                z = (j-1)/(ysz_dummy*xsz_dummy) + 1

                !Get the wavenumbers:
                k(1) = dble(x-1) * 2.0D0/L(1)
                k(2) = dble(y-1) * 2.0D0/L(2)
                if ((y-1)>N(2)/2) k(2) = dble(y-1 - N(2)) * 2.0D0/L(2)
                k(3) = dble(z-1) * 2.0D0/L(3)
                if ((z-1)>N(3)/2) k(3) = dble(z-1 - N(3)) * 2.0D0/L(3)

                ! Area adjustment used for non-2pi planes of integration
                ! area/count gives average dA later
                if(ftype==0) then
                    ! sharp spectral box (if down-sampled, this does nothing)
                    if( k(1)>k3_max(1) ) cycle
                    if( abs(k(2))>k3_max(2) ) cycle
                    if( abs(k(3))>k3_max(3) ) cycle

                    !the rectangular area
                    area = 4.0d0* k3_max(indp) * k3_max(indm)

                elseif(ftype==1) then
                    ! spherical/ellipsoidal filter
                    check = (k(1)/k3_max(1))**2 + (k(2)/k3_max(2))**2 + (k(3)/k3_max(3))**2
                    if(check>1.0D0) cycle

                    !1/4 the elliptical area
                    area = pi* k3_max(indp) * k3_max(indm) * sqrt(1.0D0-(k(i)/k3_max(i))**2)
                endif


                ! Weight nonzero numbers to account for unrepresented wavenumbers
                if (k(1)==0) then
                    xfactor = 1.0d0
                else
                    xfactor = 2.0d0
                endif

                ksq = k(1)**2 + k(2)**2 + k(3)**2

                ! Fourier Transform of the covariance of velocity (Pope 6.152)
                do m = 1, 3
                    Rij_p(i,m) = Rij_p(i,m) + xfactor*un_data(y,x,z,i)*dconjg(un_data(y,x,z,m))
                enddo

                ! Everything is the same for kx=0 except for a factor of 2, or not
                epslonp(i) = epslonp(i) + xfactor*nu*ksq*un_data(y,x,z,i)*dconjg(un_data(y,x,z,i))
                varp(i) = varp(i) + xfactor*un_data(y,x,z,i)*dconjg(un_data(y,x,z,i))
                ! To compute integral scale as integral of E(k)/k
                if(ksq/=0) ekkp = ekkp + xfactor*un_data(y,x,z,i)*dconjg(un_data(y,x,z,i))/sqrt(ksq)
                ! To compute integral scale as the k=0 mode of 1-D energy spectrum
                if(k(i)==0) then
                    Eii_at_0p(i) = Eii_at_0p(i) + &
                        area*xfactor*un_data(y,x,z,i)*dconjg(un_data(y,x,z,i))
                    counter(i) = counter(i) + int(xfactor)
                endif
            enddo
            var(i) = 0d0
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (varp(i),var(i),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

            epslonr(i) = 0d0
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (epslonp(i),epslonr(i),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
            do m = 1,3
                call MPI_BARRIER(MPI_COMM_WORLD,ierr)
                call MPI_ALLREDUCE (Rij_p(i,m), Rij(i,m),1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
            enddo
        enddo

        ekkr = 0d0
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE (ekkp,ekkr,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

        ! Overall RMS value (average of the 1D rms values)
        uvar = (var(1)+var(2)+var(3))/3.0
        urms = dsqrt(uvar)

        ! Total Dissipation epslon = epslon_ii
        epslon = epslonr(1)+epslonr(2)+epslonr(3)

        ! Turbulent Kinetic Energy
        ! K = .5 * q^2
        ! q = sqrt(3) * u' (isotropic)
        ! Thus, K = .5 * 3 * (u')^2
        tke = .5d0 * 3 * urms**2

        Eii_at_0 = 0d0
        do i = 1,3
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (Eii_at_0p(i),Eii_at_0(i),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        enddo

        ! Integral Length Scales in the 11, 22, 33 directions
        ! See "Distortion of Homogeneous Turbulence by axisymmetric strain
        !   and dilatation," Lee, 1989, Formula 34
        ! Or Pope 6.206, 6.213
        intsc1d = 0d0
        do i = 1,3
            j = modulo(i,3)+1
            m = modulo(i-2,3)+1

            ! Note the factors of 2.0 and 0.5 are used in Pope 6.206 & 6.213
            ! but not in Lee 1989
            Eii_at_0(i) = 2.0d0*Eii_at_0(i)/counter(i)
            ! From Pope 6.213:
            intsc1d(i) = 0.5*pi*Eii_at_0(i)/(var(i))
        enddo

        Lavg = sum(intsc1d)/3.0
        uvar = (var(1)+var(2)+var(3))/3.0

        !-----------------------------------------------------------------------
        ! Anisotropy Stats
        !-----------------------------------------------------------------------
        trace_tau = (Rij(1,1) + Rij(2,2) + Rij(3,3))

        dij = 0.0D0
        dij(1,1) = 1.0D0
        dij(2,2) = 1.0D0
        dij(3,3) = 1.0D0

        ! Reynolds stress anisotropy:
        a = Rij - trace_tau*dij/(3.0D0)

        ! Normalized reynolds stress anisotropy and its powers
        b = a/trace_tau
        b2 = matmul(b,b)
        b3 = matmul(b,b2)

        ! Invariants of normalized reynolds stress anisotropy
        II_b = -0.5D0*(b2(1,1) + b2(2,2) + b2(3,3))
        III_b = 1.0d0/3.0d0*(b3(1,1)+b3(2,2)+b3(3,3))

        ! Normalized invariants
        eta = sqrt(-II_b/3.0D0)
        ksi = (0.5D0*III_b)**(1.0D0/3.0D0)
         !----------------------------------------------------------------------

        ! Output the stats to text files
        call write_stats(create_new_file, "var", timestep, time, &
            ["<u^2>","<v^2>","<w^2>"], var)
        call write_stats(create_new_file, "avg_var", timestep, time, &
            ["<u^2>/<uii^2>","<v^2>/<uii^2>","<w^2>/<uii^2>"], var/uvar)
        call write_stats(create_new_file, "instsc_ii", timestep, time, &
            ["L_11","L_22","L_33"], intsc1d)
        call write_stats(create_new_file, "avg_Lii", timestep, time, &
            ["L_11/L_ii","L_22/L_ii","L_33/L_ii"], intsc1d/Lavg)
        call write_stats(create_new_file, "E_ii_at_0", timestep, time, &
            ["E_11(k_1=0)","E_22(k_2=0)","E_33(k_3=0)"], Eii_at_0)
        call write_stats(create_new_file, "b_tensor", timestep, time, &
            ["b_11","b_22","b_33"], [b(1,1), b(2,2), b(3,3)])
        call write_stats(create_new_file, "aniso", timestep, time, &
            ["II_b","III_b","k"], [II_b,III_b,(0.5d0*trace_tau)])
        call write_stats(create_new_file, "urms", timestep, time, &
            ["urms"], [urms])
        call write_stats(create_new_file, "tke", timestep, time, &
            ["tke"], [tke])
        call write_stats(create_new_file, "diss", timestep, time, &
            ["dissipation rate"], [diss])

    end subroutine

    !--------------------------------------------------------------------------
    ! Outputs the statistics with a header to a text file
    !
    ! This function is meant to cut down on repitition when saving data,
    ! while still providing useful features like headers.
    ! This function will overwrite any existing files if "create_new_file"
    ! is set to true.
    ! By setting "create_new_file" initially to "true", and then subsequently
    ! to "false", you can call this function for each timestep to store all
    ! of the timesteps in a single set of files.
    !
    ! INPUTS:
    ! \param[in] create_new_file - True if the file should be created
    ! \param[in] title - A short postfix contained the title of the data
    ! \param[in] timestep - The timestep for the data being saved
    ! \param[in] time - The simulation time for the data being saved
    ! \param[in] headers - An array of character arrays containing table headers
    ! \param[in] data - An array of arbitrary size containing the data to be saved.
    !--------------------------------------------------------------------------
    subroutine write_stats(create_new_file, title, timestep, time, headers, data)
        use, intrinsic :: iso_c_binding
        use globals, only: output_prefix, N, f3, L, diss, FV, DS
        implicit none

        ! Input/output
        logical, intent(in)        :: create_new_file
        character(*), intent(in)   :: title
        integer, intent(in)        :: timestep
        real(C_DOUBLE), intent(in) :: time
        character(*) :: headers(*) ! An array of character arrays
        real(C_DOUBLE) :: data(:)  ! Array size not specified
        ! Local variables
        integer :: unit
        integer :: i, num_columns
        character(255) :: filename
        character(255) :: data_format, header_format
        character(8)  :: date
        character(10) :: savetime
        integer,dimension(8) :: values

        filename = trim(adjustl(output_prefix))//"_"//trim(adjustl(title))//".txt"

        data_format   = "(1I7, 1e24.12, 1x, *(e24.15))"

        if (create_new_file) then
            open(newunit=unit, form="formatted", status="UNKNOWN", &
                file=filename)
            call date_and_time(DATE=date)
            call date_and_time(TIME=savetime)
            call date_and_time(VALUES=values)
            write(unit,"(A)") "#------------------------------------------------------------------------------"
            write(unit,"(A)") "# POSTPROCESSED DATA"
            write(unit,"(A)") "#"
            write(unit,"(A)") "# Postprocessing performed at:"
            write(unit,"(A,A,2x,A)") "# ", &
                date(5:6)//"/"//date(7:8)//"/"//date(1:4), &
                savetime(1:2)//":"//savetime(3:4)//":"//savetime(5:6)
            write(unit,"(A)") "#"
            write(unit,"(A)") "# Field description:"
            write(unit,"(A,3I8)")   "# N      = ", N
            write(unit,"(A,3I8)")   "# filter = ", f3
            write(unit,"(A,3f8.3)") "# length = ", L
            write(unit,"(A,f8.3)")  "# diss   = ", diss
            write(unit,"(A,L)")     "# FV filter?   ", FV
            write(unit,"(A,L)")     "# Downsampled? ", DS
            write(unit,"(A)") "#------------------------------------------------------------------------------"
            ! Write the column headings
            header_format = "(A7, 1A24,    1x)"
            write(unit, header_format, advance='no') "# step ", "time"
            ! Separating this out into a loop keeps the compiler happy
            num_columns = size(data,1)
            do i = 1, num_columns-1
                write(unit, "(A24)", advance='no') headers(i)
            enddo
            write(unit, "(A24)") headers(num_columns)
        else
            open(newunit=unit, form="formatted", position="APPEND", &
                status="OLD", file=filename)
        endif

        data_format   = "(1I7, 1e24.12, 1x, *(e24.15))"
        write(unit,data_format) timestep, time, data
        close(unit)


end subroutine


!---------------------------------------------------------------------
subroutine calc_dissipation(N,L,f3,nu,s_u,s_v,s_w,eps)
    use, intrinsic :: iso_c_binding
    implicit none
    !Inputs/Outputs
    integer, intent(in) :: N(3), f3(3)
    real(C_DOUBLE), intent(in) :: L(3)
    real(C_DOUBLE), intent(in) :: nu
    real(C_DOUBLE), intent(out) :: eps
    complex(C_DOUBLE_COMPLEX)   :: s_u(0:N(1)/2,0:N(2)-1,0:N(3)-1)
    complex(C_DOUBLE_COMPLEX)   :: s_v(0:N(1)/2,0:N(2)-1,0:N(3)-1)
    complex(C_DOUBLE_COMPLEX)   :: s_w(0:N(1)/2,0:N(2)-1,0:N(3)-1)
        !Dummy variables
        integer :: ix, iy, iz
        integer :: weight
        real(C_DOUBLE) :: e, k2
        real(C_DOUBLE) :: kx,ky,kz,kx_max, ky_max, kz_max


        eps = 0.0d0
        ! toss nyquist
        ! max wavenumber (can toss nyquist here)
        kx_max = 2.0D0/L(1) * (dble(N(1)/2)/dble(f3(1))-1D0)
        ky_max = 2.0D0/L(2) * (dble(N(2)/2)/dble(f3(2))-1D0)
        kz_max = 2.0D0/L(3) * (dble(N(3)/2)/dble(f3(3))-1D0)

        ! 3D energy spectrum
        do ix = 0, N(1)/2
            kx = dble(ix) * 2.0D0/L(1)
            weight = 2
            if (kx.eq.0 .or. kx.eq.dble(N(1))/2.0D0) weight = 1
            do iy = 0, N(2)-1
                ky = dble(iy) * 2.0D0/L(2)
                if (ky.gt.dble(N(2))/2.0D0) ky = dble(iy - N(2)) * 2.0D0/L(2)
                do iz = 0, N(3)-1
                    kz = dble(iz) * 2.0D0/L(3)
                    if (kz.gt.dble(N(3))/2.0D0) kz = dble(iz - N(3)) * 2.0D0/L(3)

                    ! sharp spectral (discard any wavenumbers unsupported by the filtered grid)
                    if( kx.gt.kx_max ) cycle
                    if( abs(ky).gt.ky_max ) cycle
                    if( abs(kz).gt.kz_max ) cycle

                    k2 = kx**2+ky**2+kz**2

                    e = 0.5d0*dble(weight) * &
                        (real(s_u(ix,iy,iz))*real(s_u(ix,iy,iz)) + imag(s_u(ix,iy,iz))*imag(s_u(ix,iy,iz)) + &
                         real(s_v(ix,iy,iz))*real(s_v(ix,iy,iz)) + imag(s_v(ix,iy,iz))*imag(s_v(ix,iy,iz)) + &
                         real(s_w(ix,iy,iz))*real(s_w(ix,iy,iz)) + imag(s_w(ix,iy,iz))*imag(s_w(ix,iy,iz)))

                    !See Pope  6.160, 6.191:
                    eps = eps + (2.0d0*nu*k2*e)

                enddo
            enddo
        enddo

    return

    end subroutine


!---------------------------------------------------------------------
subroutine save_spectrum(N,Nmax,f3,L,nu,diss,spec,spec_x,spec_y,spec_z,samples,prefix)
implicit none

  integer                :: N(3), Nmax, f3(3)
  integer                :: i, samples
  integer    , parameter :: spec_output  = 91
  integer    , parameter :: spec_output_x  = 93
  integer    , parameter :: spec_output_y  = 94
  integer    , parameter :: spec_output_z  = 95
  real(8)                :: delta(3), eta, max, kmin, diss, L(3), nu, kmax(3)
  real(8)                :: pi = 3.14159265
  character(len=*)       :: prefix
  character(100)         :: spec_name
  character(100)         :: spec_name_x, spec_name_y, spec_name_z
  character(3)           :: cN
  real(8), dimension(0:Nmax-1,0:1) :: spec
  real(8), dimension(0:N(1)-1,0:1) :: spec_x
  real(8), dimension(0:N(2)-1,0:1) :: spec_y
  real(8), dimension(0:N(3)-1,0:1) :: spec_z



  if (N(1) .ge. 100) then
     write(cN,"(I3.3)") N(1) ! for now, use N(1) only
  else
     write(cN,"(I2.2)") N(1)
  endif

  spec_name = trim(adjustl(prefix))//".spec"
  spec_name_x = trim(adjustl(prefix))//"_x.spec"
  spec_name_y = trim(adjustl(prefix))//"_y.spec"
  spec_name_z = trim(adjustl(prefix))//"_z.spec"
  open(unit=spec_output, file=trim(adjustl(spec_name)), form="FORMATTED",status="UNKNOWN")
  open(unit=spec_output_x, file=trim(adjustl(spec_name_x)), form="FORMATTED",status="UNKNOWN")
  open(unit=spec_output_y, file=trim(adjustl(spec_name_y)), form="FORMATTED",status="UNKNOWN")
  open(unit=spec_output_z, file=trim(adjustl(spec_name_z)), form="FORMATTED",status="UNKNOWN")

  print *, " Saving energy spectrum to: ", spec_name

  delta = L*pi/dble(N)*dble(f3)
  kmin = 2.0D0/maxval(L) * sqrt(3.0D0)
  print *, " Minimum average wavenumber:", kmin, "samples no:",samples
  max = 0.0D0

  ! max wavenumber (can toss nyquist here)
  kmax(1) = 2.0D0/L(1) * (dble(N(1))/dble(2*f3(1)) - 1.0d0)
  kmax(2) = 2.0D0/L(2) * (dble(N(2))/dble(2*f3(2)) - 1.0d0)
  kmax(3) = 2.0D0/L(3) * (dble(N(3))/dble(2*f3(3)) - 1.0d0)


  do i = 0, Nmax/2

         ! unnormalized
         !write(spec_output,*) spec(i,0)/dble(samples), spec(i,1)/dble(samples)

         ! normalized
        if(i.le.maxval(N/2)) then
           if(spec(i,1).gt.0.0D0 .OR. i.eq.0) write(spec_output,"(2E13.4)") spec(i,0)/dble(samples), &
              spec(i,1)/(dble(samples))
        endif

        ! directional
        if(i.le.N(1)/(2*f3(1))-1) then
           if(spec_x(i,1).gt.0.0D0) write(spec_output_x,"(2E13.4)") spec_x(i,0)/dble(samples), &
              spec_x(i,1)/(dble(samples))*(diss**(-2.0D0/3.0D0)*(2.0D0/L(1))**(5.0D0/3.0D0))
        endif
        if(i.le.N(2)/(2*f3(2))-1) then
           if(spec_y(i,1).gt.0.0D0) write(spec_output_y,"(2E13.4)") spec_y(i,0)/dble(samples), &
              spec_y(i,1)/(dble(samples))*(diss**(-2.0D0/3.0D0)*(2.0D0/L(2))**(5.0D0/3.0D0))
        endif
        if(i.le.N(3)/(2*f3(3))-1) then
           if(spec_z(i,1).gt.0.0D0) write(spec_output_z,"(2E13.4)") spec_z(i,0)/dble(samples), &
              spec_z(i,1)/(dble(samples))*(diss**(-2.0D0/3.0D0)*(2.0D0/L(3))**(5.0D0/3.0D0))
        endif

  enddo

  close(spec_output)
  close(spec_output_x)
  close(spec_output_y)
  close(spec_output_z)

end subroutine save_spectrum


!------------------------------------------------------------------------------
subroutine save_structure(N,Nmax,f3,L,nu,diss,strct,samples,prefix)
implicit none

  integer                :: N(3), Nmax, f3(3)
  integer                :: i, samples
  integer    , parameter :: strct_output = 92
  real(8)                :: eta, max, kmin, nu, diss, L(3), delta(3)
  real(8)                :: pi = 3.14159265
  character(100)         :: strct_name, prefix
  character(3)           :: cN
  real(8), dimension(0:Nmax/2)     :: strct



  if (N(1) .ge. 100) then
     write(cN,"(I3.3)") N(1) ! for now, use N(1) only
  else
     write(cN,"(I2.2)") N(1)
  endif

  strct_name = trim(adjustl(prefix))//".strct"
  open(unit=strct_output,file=trim(adjustl(strct_name)),form="FORMATTED",status="UNKNOWN")
  print *, " Saving structure function to: ", strct_name

  delta = L*pi/dble(N)*dble(f3)

  do i = 1, Nmax/2 !N(1)/2 ! Nmin ! in this case, i is r/delta
     ! unnormalized
     !write(strct_output,*) i, strct(i)/dble(samples)
     ! old scaling, doesn't make sense
     !write(strct_output,*) i * delta/eta, strct(i)/dble(samples)

     ! normalized
     write(strct_output,*) i, strct(i)/(dble(samples))*1.0D0/(diss*maxval(delta))

  enddo


  close(strct_output)

end subroutine save_structure


!------------------------------------------------------------------------------
subroutine print_realField(N,L,f3,phi1,phi2,phi3)
use, intrinsic :: iso_c_binding
use globals, only : plan_c2r, plan_r2c2
implicit none
include 'fftw3.f03'

  integer :: N(3), f3(3)
  integer :: ix, iy, iz, i, j, k
  real(8) :: x, y, z
  real(8) :: L(3)
  real(8) :: pi = 3.14159265
  complex(C_DOUBLE_COMPLEX), INTENT(IN) :: phi1(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(IN) :: phi2(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(IN) :: phi3(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX) :: phi_f1(0:N(1)/(2*f3(1)),0:(N(2)/f3(2))-1,0:(N(3)/f3(3))-1)
  complex(C_DOUBLE_COMPLEX) :: phi_f2(0:N(1)/(2*f3(1)),0:(N(2)/f3(2))-1,0:(N(3)/f3(3))-1)
  complex(C_DOUBLE_COMPLEX) :: phi_f3(0:N(1)/(2*f3(1)),0:(N(2)/f3(2))-1,0:(N(3)/f3(3))-1)
  real(C_DOUBLE), dimension(1:N(1)/f3(1),1:N(2)/f3(2),1:N(3)/f3(3)) :: u,v,w
  type(C_PTR) :: planf_c2r

  write(*,'(a,I1,a)') "   ...exporting field in real space"

  if (f3(1) + f3(2) + f3(3) > 1.0) then
     do i = 0, N(1)/(2.0*f3(1))
        ix = i
        do j = 0, N(2)/f3(2)-1
           iy = j
           if (j .gt. N(2)/(2.0*f3(2))) then
              iy = N(2) - (N(2)/f3(2) - j)
           end if 
           do k = 0, N(3)/f3(3)-1
              iz = k
              if (k .gt. N(3)/(2.0*f3(3))) then
                 iz = N(3) - (N(3)/f3(3) - k)
              end if
   
              phi_f1(i,j,k) = phi1(ix,iy,iz)
              phi_f2(i,j,k) = phi2(ix,iy,iz)
              phi_f3(i,j,k) = phi3(ix,iy,iz)
   
           end do
        end do
     end do
  end if

  call dfftw_plan_dft_c2r_3d(planf_c2r,N(1)/f3(1),N(2)/f3(2),N(3)/f3(3), phi_f1, u, FFTW_ESTIMATE)

  ! transform to realspace
  call dfftw_execute_dft_c2r(planf_c2r,phi_f1, u)
  call dfftw_execute_dft_c2r(planf_c2r,phi_f2, v)
  call dfftw_execute_dft_c2r(planf_c2r,phi_f3, w)

  ! sample
  do i = 1, N(1)/f3(1)
     do j = 1, N(2)/f3(2)
        do k = 1, N(3)/f3(3)

           ix = (i-1)+1
           iy = (j-1)+1
           iz = (k-1)+1

           x = (f3(1)*L(1)/N(1))*(i-1)*pi
           y = (f3(2)*L(2)/N(2))*(j-1)*pi
           z = (f3(3)*L(3)/N(3))*(k-1)*pi

           write (*,"(6e24.16)") x, y, z, u(ix,iy,iz), v(ix,iy,iz), w(ix,iy,iz)

        enddo
     enddo
  enddo
  
  call dfftw_destroy_plan(planf_c2r)
end subroutine print_realField



!------------------------------------------------------------------------------
subroutine stretchSqueeze(N,L,f3,f3sq,phi1,phi2,phi3)
use, intrinsic :: iso_c_binding
use globals, only : plan_c2r, plan_r2c2, plan_c2r2
implicit none
include 'fftw3.f03'

  integer :: N(3)
  integer :: ix, iy, iz, i, j, k
  real(8) :: x, y, z
  real(8) :: L(3), f3(3), f3sq(3), scale
  real(8) :: pi = 3.14159265
  complex(C_DOUBLE_COMPLEX), INTENT(IN) :: phi1(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(IN) :: phi2(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(IN) :: phi3(0:N(1)/2,0:N(2)-1,0:N(3)-1)

  complex(C_DOUBLE_COMPLEX) :: phisq1(0:N(1)/(2*f3sq(1)),0:N(2)/f3sq(2)-1,0:N(3)/f3sq(3)-1)
  complex(C_DOUBLE_COMPLEX) :: phisq2(0:N(1)/(2*f3sq(1)),0:N(2)/f3sq(2)-1,0:N(3)/f3sq(3)-1)
  complex(C_DOUBLE_COMPLEX) :: phisq3(0:N(1)/(2*f3sq(1)),0:N(2)/f3sq(2)-1,0:N(3)/f3sq(3)-1)

  real(C_DOUBLE) :: phif1(0:1,0:N(2)/f3(2)-1,0:N(1)/(2*f3(1))-1,0:N(3)/f3(3)-1)
  real(C_DOUBLE) :: phif2(0:1,0:N(2)/f3(2)-1,0:N(1)/(2*f3(1))-1,0:N(3)/f3(3)-1)
  real(C_DOUBLE) :: phif3(0:1,0:N(2)/f3(2)-1,0:N(1)/(2*f3(1))-1,0:N(3)/f3(3)-1)

  real(C_DOUBLE), dimension(1:N(1),1:N(2),1:N(3)) :: u,v,w
  real(C_DOUBLE), dimension(1:N(1)/f3sq(1),1:N(2)/f3sq(2),1:N(3)/f3sq(3)) :: usq,vsq,wsq

  real :: debug_help, debug_help1

  !type(C_PTR) :: plansq_r2c

  !  call dfftw_plan_dft_r2c_3d(plansq_r2c,N(1)/f3sq(1),N(2)/f3sq(2),N(3)/f3sq(3),usq, phisq1, FFTW_ESTIMATE)

  ! First we need to deal with the squeezed directions.  To do that, we
  ! transform to real space, downsample the domain and then transform back to
  ! fourier space

  debug_help = u(5,3,7)
  debug_help1 = u(7,4,11)

  ! transform to realspace
  call dfftw_execute_dft_c2r(plan_c2r,phi1,u)
  call dfftw_execute_dft_c2r(plan_c2r,phi2,v)
  call dfftw_execute_dft_c2r(plan_c2r,phi3,w)

  debug_help = real(phi1(6,5,9))
  debug_help1 = imag(phi1(5,6,7))

  ! downsample the domain in the squeezed directions (i.e. take every nth entry)
  do i = 1, N(1)/f3sq(1)
     do j = 1, N(2)/f3sq(2)
        do k = 1, N(3)/f3sq(3)
           ix=i*f3sq(1)
           iy=j*f3sq(2)
           iz=k*f3sq(3)
           usq(i,j,k) = u(ix,iy,iz)
           vsq(i,j,k) = v(ix,iy,iz)
           wsq(i,j,k) = w(ix,iy,iz)
        enddo
     enddo
  enddo

  debug_help = usq(5,6,7)
  debug_help1 = usq(3,4,7)

  ! transform back to fourier space
  call dfftw_execute_dft_r2c(plan_r2c2,usq,phisq1)
  call dfftw_execute_dft_r2c(plan_r2c2,vsq,phisq2)
  call dfftw_execute_dft_r2c(plan_r2c2,wsq,phisq3)

  debug_help = real(phisq1(5,6,7))
  debug_help1 = imag(phisq1(3,4,7))

  ! Next, we handle the stretched directions.  To do that, we simply add 0's
  ! into the new additional wave numbers, so we copy over the first half of the
  ! vector, and copy the last half to the end of the new vector.

  phif1 = 0.0
  phif2 = 0.0
  phif3 = 0.0
  scale = 1.0/(N(1)/f3sq(1) * N(2)/f3sq(2) * N(3)/f3sq(3))
  do i = 0, N(1)/(2*f3sq(1)) - 1
     ix = i
     do j = 0, N(2)/f3sq(2) - 1
        if(f3sq(2) == f3(2)) then
           iy=j
        else
           iy=j
           if(j > N(2)/(2*f3sq(2))) then
              iy=N(2)/f3(2) - N(2)/f3sq(2) + j
           endif
        endif
        do k = 0, N(3)/f3sq(3) - 1
           if(f3sq(3) == f3(3)) then 
              iz=k
           else
              iz=k
              if(k > N(3)/(2*f3sq(3))) then
                 iz=N(3)/f3(3) - N(3)/f3sq(3) + k
              endif
           endif
           phif1(0,iy,ix,iz) = real(phisq1(i,j,k))*scale 
           phif2(0,iy,ix,iz) = real(phisq2(i,j,k))*scale
           phif3(0,iy,ix,iz) = real(phisq3(i,j,k))*scale
           phif1(1,iy,ix,iz) = imag(phisq1(i,j,k))*scale
           phif2(1,iy,ix,iz) = imag(phisq2(i,j,k))*scale
           phif3(1,iy,ix,iz) = imag(phisq3(i,j,k))*scale
        enddo
     enddo
  enddo

  debug_help = phif1(1,6,2,7)
  debug_help1 = phif1(0,5,8,7)

  call exportField(phif1,phif2,phif3)

  !call dfftw_destroy_plan(plansq_r2c)
end subroutine stretchSqueeze

!------------------------------------------------------------------------------
subroutine stretchSqueezeAlt(N,L,f3,f3sq,phi1,phi2,phi3)
use, intrinsic :: iso_c_binding
use globals, only : plan_c2r, plan_r2c2, plan_c2r2
implicit none
include 'fftw3.f03'

  integer :: N(3)
  integer :: ix, iy, iz, i, j, k
  real(8) :: x, y, z
  real(8) :: L(3), f3(3), f3sq(3), scale
  real(8) :: pi = 3.14159265
  complex(C_DOUBLE_COMPLEX), INTENT(IN) :: phi1(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(IN) :: phi2(0:N(1)/2,0:N(2)-1,0:N(3)-1)
  complex(C_DOUBLE_COMPLEX), INTENT(IN) :: phi3(0:N(1)/2,0:N(2)-1,0:N(3)-1)

  complex(C_DOUBLE_COMPLEX) ::
phisq1(0:N(1)/(2*f3sq(1)),0:N(2)/f3sq(2)-1,0:N(3)/f3sq(3)-1)
  complex(C_DOUBLE_COMPLEX) ::
phisq2(0:N(1)/(2*f3sq(1)),0:N(2)/f3sq(2)-1,0:N(3)/f3sq(3)-1)
  complex(C_DOUBLE_COMPLEX) ::
phisq3(0:N(1)/(2*f3sq(1)),0:N(2)/f3sq(2)-1,0:N(3)/f3sq(3)-1)

  real(C_DOUBLE) :: phif1(0:1,0:N(2)/f3(2)-1,0:N(1)/(2*f3(1))-1,0:N(3)/f3(3)-1)
  real(C_DOUBLE) :: phif2(0:1,0:N(2)/f3(2)-1,0:N(1)/(2*f3(1))-1,0:N(3)/f3(3)-1)
  real(C_DOUBLE) :: phif3(0:1,0:N(2)/f3(2)-1,0:N(1)/(2*f3(1))-1,0:N(3)/f3(3)-1)

  real(C_DOUBLE), dimension(1:N(1),1:N(2),1:N(3)) :: u,v,w
  real(C_DOUBLE), dimension(1:N(1)/f3sq(1),1:N(2)/f3sq(2),1:N(3)/f3sq(3)) ::
usq,vsq,wsq

  real :: debug_help, debug_help1

  !type(C_PTR) :: plansq_r2c

  !  call
  !  dfftw_plan_dft_r2c_3d(plansq_r2c,N(1)/f3sq(1),N(2)/f3sq(2),N(3)/f3sq(3),usq,
  !  phisq1, FFTW_ESTIMATE)

  ! First we need to deal with the squeezed directions.  To do that, we
  ! transform to real space, downsample the domain and then transform back to
  ! fourier space

  debug_help = u(5,3,7)
  debug_help1 = u(7,4,11)

  ! transform to realspace
  call dfftw_execute_dft_c2r(plan_c2r,phi1,u)
  call dfftw_execute_dft_c2r(plan_c2r,phi2,v)
  call dfftw_execute_dft_c2r(plan_c2r,phi3,w)

  debug_help = real(phi1(6,5,9))
  debug_help1 = imag(phi1(5,6,7))

  ! downsample the domain in the squeezed directions (i.e. take every nth entry)
  do i = 1, N(1)/f3sq(1)
     do j = 1, N(2)/f3sq(2)
        do k = 1, N(3)/f3sq(3)
           ix=i*f3sq(1)
           iy=j*f3sq(2)
           iz=k*f3sq(3)
           usq(i,j,k) = u(ix,iy,iz)
           vsq(i,j,k) = v(ix,iy,iz)
           wsq(i,j,k) = w(ix,iy,iz)
        enddo
     enddo
  enddo

  debug_help = usq(5,6,7)
  debug_help1 = usq(3,4,7)

  ! transform back to fourier space
  call dfftw_execute_dft_r2c(plan_r2c2,usq,phisq1)
  call dfftw_execute_dft_r2c(plan_r2c2,vsq,phisq2)
  call dfftw_execute_dft_r2c(plan_r2c2,wsq,phisq3)

  debug_help = real(phisq1(5,6,7))
  debug_help1 = imag(phisq1(3,4,7))

  ! Next, we handle the stretched directions.  To do that, we simply add 0's
  ! into the new additional wave numbers, so we copy over the first half of the
  ! vector, and copy the last half to the end of the new vector.

  phif1 = 0.0
  phif2 = 0.0
  phif3 = 0.0
  scale = 1.0/(N(1)/f3sq(1) * N(2)/f3sq(2) * N(3)/f3sq(3))
  do i = 0, N(1)/(2*f3sq(1)) - 1
     ix = i
     do j = 0, N(2)/f3sq(2) - 1
        if(f3sq(2) == f3(2)) then
           iy=j
        else
           iy=j
           if(j > N(2)/(2*f3sq(2))) then
              iy=N(2)/f3(2) - N(2)/f3sq(2) + j
           endif
        endif
        do k = 0, N(3)/f3sq(3) - 1
           if(f3sq(3) == f3(3)) then
              iz=k
           else
              iz=k
              if(k > N(3)/(2*f3sq(3))) then
                 iz=N(3)/f3(3) - N(3)/f3sq(3) + k
              endif
           endif
           phif1(0,iy,ix,iz) = real(phisq1(i,j,k))*scale
           phif2(0,iy,ix,iz) = real(phisq2(i,j,k))*scale
           phif3(0,iy,ix,iz) = real(phisq3(i,j,k))*scale
           phif1(1,iy,ix,iz) = imag(phisq1(i,j,k))*scale
           phif2(1,iy,ix,iz) = imag(phisq2(i,j,k))*scale
           phif3(1,iy,ix,iz) = imag(phisq3(i,j,k))*scale
        enddo
     enddo
  enddo

  debug_help = phif1(1,6,2,7)
  debug_help1 = phif1(0,5,8,7)

  call exportField(phif1,phif2,phif3)

  !call dfftw_destroy_plan(plansq_r2c)
end subroutine stretchSqueezeAlt

!--------------------------------------------------------------------------------
subroutine exportField(uoutf,voutf,woutf)
  use mpi
  use globals
  use esio
  use, intrinsic :: iso_c_binding
  implicit none

  integer           :: world_size, world_rank, ierr
  integer           :: nvFilterSize, nvIndex, nvCounter
  real(C_DOUBLE)    :: alphaNV,kfmax,kfmin,nvFilter,targetDiss
  real(C_DOUBLE), allocatable :: energyInpStore(:)
  type(esio_handle) :: ohandle
  integer           :: timestep, forcing_on
  real(C_DOUBLE)    :: Re, time
  real(C_DOUBLE), INTENT(IN) :: uoutf(0:1,0:N(2)/sqSt(2)-1,0:N(1)/(2*sqSt(1))-1,0:N(3)/sqSt(3)-1)
  real(C_DOUBLE), INTENT(IN) :: voutf(0:1,0:N(2)/sqSt(2)-1,0:N(1)/(2*sqSt(1))-1,0:N(3)/sqSt(3)-1)
  real(C_DOUBLE), INTENT(IN) :: woutf(0:1,0:N(2)/sqSt(2)-1,0:N(1)/(2*sqSt(1))-1,0:N(3)/sqSt(3)-1)

  ! Write new field...

  ! initialize esio resized file
  call esio_handle_initialize(ohandle,MPI_COMM_WORLD)
  call esio_file_create (ohandle,"filteredField.h5", .true.)
  call esio_line_establish(ohandle,1,1,1,ierr)

  ! Write out field
  call esio_field_establish (ohandle, INT(N(2)/sqSt(2)), 1, INT(N(2)/sqSt(2)), INT(N(1)/(2*sqSt(1))), 1, INT(N(1)/(2*sqSt(1))), INT(N(3)/sqSt(3)), 1, INT(N(3)/sqSt(3)), ierr)
  call esio_field_writev_double (ohandle, "u",  uoutf(:,:,:,:),2)
  call esio_field_writev_double (ohandle, "v",  voutf(:,:,:,:),2)
  call esio_field_writev_double (ohandle, "w",  woutf(:,:,:,:),2)

  ! write scalar values
  call esio_line_establish (ohandle, 1, 1, 1,ierr)
  call esio_line_write_integer(ohandle,"timestep", [INT(ts)],0,"Current step in time")
  call esio_line_write_double(ohandle,"time",[DBLE(t)],0,"Simulation Time")
  call esio_line_write_double(ohandle,"nu",[DBLE(nu)],0,"Viscosity or 1/Re")
  call esio_line_write_double(ohandle,"Re",[DBLE(1./nu)],0,"Reynold Number (not Re_tau) = 1/Re")
  call esio_line_write_double(ohandle,"LMFx",[DBLE(L(1)/(2.0*sqSt(1)))],0,"Length Multiplication Factor: Lx=LMFx*2pi")
  call esio_line_write_double(ohandle,"LMFy",[DBLE(L(2)/(2.0*sqSt(2)))],0,"Length Multiplication Factor: Ly=LMFy*2pi")
  call esio_line_write_double(ohandle,"LMFz",[DBLE(L(3)/(2.0*sqSt(3)))],0,"Length Multiplication Factor: Lz=LMFz*2pi")
  L = L*pi

  call esio_line_write_double(ohandle,"Lx",[DBLE(L(1)/sqSt(1))],0,"Length in streamwise direction : Lx")
  call esio_line_write_double(ohandle,"Ly",[DBLE(L(2)/sqSt(2))],0,"Length in wall normal direction : Ly")
  call esio_line_write_double(ohandle,"Lz",[DBLE(L(3)/sqSt(3))],0,"Length in spanwise direction : Lz")

  call esio_line_write_double(ohandle,"dt",[DBLE(0.0000001)],0,"Time advance interval")

  call esio_line_write_integer(ohandle,"nx",[INT(N(1)/sqSt(1))],0,"Number of x points")
  call esio_line_write_integer(ohandle,"ny",[INT(N(2)/sqSt(2))],0,"Number of y points")
  call esio_line_write_integer(ohandle,"nz",[INT(N(3)/sqSt(3))],0,"Number of z points")

  call esio_line_write_integer(ohandle,"forcing_on",[INT(forcing_on)],0,"Is Forcing On?")

  forcing_on = 0
  if (forcing_on.ne.0) then
      call esio_line_write_double(ohandle,"kfmin",[DBLE(kfmin)],0,"Minimum forced wavemode")
      call esio_line_write_double(ohandle,"kfmax",[DBLE(kfmax)],0,"Maximum forced wavemode")
      call esio_line_write_integer(ohandle,"nvIndex",[INT(nvIndex)],0,"Index of position in forcing filter")
      call esio_line_write_integer(ohandle,"nvCounter",[INT(nvCounter)],0,"How full the filter is")
      call esio_line_write_integer(ohandle,"nvFilterSize",[INT(nvFilterSize)],0,"Size of the filter")
      call esio_line_write_double(ohandle,"alphaNV",[DBLE(alphaNV)],0,"Current negative viscosity")
      call esio_line_write_double(ohandle,"targetDiss",[DBLE(targetDiss)],0,"Expected value of the dissipation rate")

      if (nvFilterSize.gt.0) then
           call esio_line_establish (ohandle, nvFilterSize, 1, nvFilterSize,ierr)
           call esio_line_write_double(ohandle,"nvFilter",[DBLE(energyInpStore(:))],1,"Filter for negative viscosity calculation")
           deallocate(energyInpStore)
      endif
  endif

  ! clean-up
  call esio_file_close(ohandle)
  call esio_handle_finalize(ohandle)

end subroutine exportField

!------------------------------------------------------------------------------
end module statistics

