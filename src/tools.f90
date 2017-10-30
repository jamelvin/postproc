

module tools

public check_div_wave, check_div_real


contains


!------------------------------------------------------------------------------
subroutine check_div_wave(N,L,u,v,w)
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

  integer  :: N(3)
  integer  :: ix, iy, iz
  real(8)  :: kx, ky, kz
  real(8)  :: constraint, L(3)
  complex(C_DOUBLE_COMPLEX), dimension(0:N(1)/2, 0:N(2)-1, 0:N(3)-1) :: u, v, w
 
  constraint = 0.0D0

  ! divergence in wave space
  do ix = 0, N(1)/2
     kx = dble(ix) * L(1)/2.0D0
     do iy = 0, N(2)-1
        ky = dble(iy) * L(2)/2.0D0
        if (ky.gt.N(2)/2) ky = dble(iy - N(2)) * L(2)/2.0D0
        do iz = 0, N(3)-1
           kz = dble(iz) * L(3)/2.0D0
           if (kz.gt.N(3)/2) kz = dble(iz - N(3)) * L(3)/2.0D0
           constraint = constraint + sqrt( & 
                        ( kx*real(u(ix,iy,iz)) + ky*real(v(ix,iy,iz)) + kz*real(w(ix,iy,iz)) )**2.0D0 + & 
                        ( kx*imag(u(ix,iy,iz)) + ky*imag(v(ix,iy,iz)) + kz*imag(w(ix,iy,iz)) )**2.0D0 )
        enddo
     enddo
  enddo
   
  print *, " -Spectral divergence (ave) =", constraint/((N(1)/2+1)*N(2)*N(3))

end subroutine check_div_wave


!------------------------------------------------------------------------------
subroutine check_div_real(N,u,v,w)
implicit none

  integer  :: N(3)
  integer  :: ix, iy, iz
  real(8)  :: pi = 3.14159265
  real(8)  :: constraint_c, constraint_s
  real(8), dimension(0:N(1)-1,0:N(2)-1,0:N(3)-1) :: u, v, w

  constraint_c = 0.0D0
  constraint_s = 0.0D0

  ! divergence in real space
  do ix = 0, N(1)-1
     do iy = 0, N(2)-1
        do iz = 0, N(3)-1
           constraint_c = constraint_c + abs( & 
                          u(mod(ix+1+N(1),N(1)),iy,iz) + v(ix,mod(iy+1+N(2),N(2)),iz) + w(ix,iy,mod(iz+1+N(3),N(3))) - &
                          u(mod(ix-1+N(1),N(1)),iy,iz) - v(ix,mod(iy-1+N(2),N(2)),iz) - w(ix,iy,mod(iz-1+N(3),N(3))) )
           constraint_s = constraint_s + abs( & 
                          u(mod(ix+1+N(1),N(1)),iy,iz) + v(ix,mod(iy+1+N(2),N(2)),iz) + w(ix,iy,mod(iz+1+N(3),N(3))) - &
                          u(mod(ix-0+N(1),N(1)),iy,iz) - v(ix,mod(iy-0+N(2),N(2)),iz) - w(ix,iy,mod(iz-0+N(3),N(3))) )
        enddo
     enddo
  enddo
   
  print *, " -Divergence collocated (ave) =", constraint_c /( (4.0D0*pi/N(1)) * N(1) * N(2) * N(3) )
  print *, " -Divergence staggered (ave)  =", constraint_s /( (2.0D0*pi/N(1)) * N(1) * N(2) * N(3) )

end subroutine check_div_real


end module
