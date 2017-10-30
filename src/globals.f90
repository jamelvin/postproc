

module globals
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

  save
  private
  public :: N, f3, delta, it_ini, it_fin, it_step, pad, &
            fieldnorm, fieldsize, sfieldsize, t, L,     &
            input_dir, spec, spec_x, spec_y, spec_z,    &
            strct, u, v, w, u_temp, s_u, s_v, s_w,      &
            plan_r2c, plan_c2r, type, type_fmt, nu,     &
            diss, diss_flag, pi, uo, s_uo, s_vo, s_wo,  &
            plan_r2c2, DS, FV, OutField, ftype, sqSt,   &
            output_prefix, stretch, ts, plan_c2r2


  integer        :: N(3), f3(3), fieldsize, sfieldsize, ts
  integer        :: it_ini, it_fin, it_step, pad, diss_flag, ftype
  real(8)        :: t, L(3), sqSt(3), diss, nu, delta(3), fieldnorm
  real(8)        :: pi = 3.141592653589
  real(8),   dimension(:,:), allocatable :: spec, spec_x, spec_y, spec_z
  real(8),   dimension(:),   allocatable :: strct
  real(C_DOUBLE),            dimension(:,:,:), allocatable :: u, v, w, u_temp, uo
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable :: s_u, s_v, s_w
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable :: s_uo, s_vo, s_wo
  character(100) :: input_dir, output_prefix
  character(3)   :: type
  character(7)   :: type_fmt
  logical :: DS, FV, OutField, stretch
  type(C_PTR) :: plan_r2c, plan_c2r, plan_r2c2, plan_c2r2


end module globals
