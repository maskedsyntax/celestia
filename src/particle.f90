module particle
  use kinds, only: dp
  implicit none
  private

  public :: body_t

  type :: body_t
     real(dp), dimension(3) :: pos = [0.0_dp, 0.0_dp, 0.0_dp]
     real(dp), dimension(3) :: vel = [0.0_dp, 0.0_dp, 0.0_dp]
     real(dp), dimension(3) :: acc = [0.0_dp, 0.0_dp, 0.0_dp]
     real(dp) :: mass = 1.0_dp
     real(dp) :: radius = 1.0_dp
     real(dp), dimension(3) :: color = [1.0_dp, 1.0_dp, 1.0_dp]
  end type body_t

end module particle
