module initial_conditions
  use kinds, only: dp
  use particle, only: body_t
  implicit none
  private

  public :: setup_solar_system
  public :: setup_galaxy

  contains

  subroutine setup_solar_system(bodies)
    type(body_t), dimension(:), allocatable, intent(out) :: bodies
    ! Placeholder for detailed solar system setup
    allocate(bodies(2))
    ! Sun
    bodies(1)%mass = 333000.0_dp
    bodies(1)%pos = [0.0_dp, 0.0_dp, 0.0_dp]
    bodies(1)%vel = [0.0_dp, 0.0_dp, 0.0_dp]
    bodies(1)%radius = 10.0_dp
    ! Earth
    bodies(2)%mass = 1.0_dp
    bodies(2)%pos = [149600.0_dp, 0.0_dp, 0.0_dp]
    bodies(2)%vel = [0.0_dp, 29.78_dp, 0.0_dp]
    bodies(2)%radius = 1.0_dp
  end subroutine setup_solar_system

  subroutine setup_galaxy(bodies, n, radius, total_mass, G)
    type(body_t), dimension(:), allocatable, intent(out) :: bodies
    integer, intent(in) :: n
    real(dp), intent(in) :: radius, total_mass, G
    integer :: i
    real(dp) :: r, theta, v_orbital, x, y
    real(dp) :: central_mass

    allocate(bodies(n))
    central_mass = total_mass * 0.1_dp ! 10% mass in central "black hole"
    
    ! Center mass
    bodies(1)%mass = central_mass
    bodies(1)%pos = [0.0_dp, 0.0_dp, 0.0_dp]
    bodies(1)%vel = [0.0_dp, 0.0_dp, 0.0_dp]
    bodies(1)%radius = radius * 0.01_dp

    do i = 2, n
       ! Random polar coordinates
       call random_number(r)
       r = r * radius + radius * 0.1_dp ! distributed from 0.1R to R
       call random_number(theta)
       theta = theta * 2.0_dp * 3.1415926535_dp
       
       x = r * cos(theta)
       y = r * sin(theta)
       
       bodies(i)%pos = [x, y, 0.0_dp]
       bodies(i)%mass = total_mass / n
       bodies(i)%radius = radius * 0.001_dp
       
       ! Circular orbital velocity v = sqrt(G*M_encl/r)
       v_orbital = sqrt(G * central_mass / r)
       bodies(i)%vel = [-v_orbital * sin(theta), v_orbital * cos(theta), 0.0_dp]
    end do
  end subroutine setup_galaxy

end module initial_conditions
