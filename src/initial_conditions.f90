module initial_conditions
  use kinds, only: dp
  use particle, only: body_t
  implicit none
  private

  public :: setup_solar_system
  public :: setup_galaxy
  public :: setup_collision

  contains

  subroutine setup_solar_system(bodies)
    type(body_t), dimension(:), allocatable, intent(out) :: bodies
    integer, parameter :: n_planets = 8
    real(dp), dimension(n_planets) :: distances, masses, v_orbit
    integer :: i

    ! Sun + 8 Planets
    allocate(bodies(n_planets + 1))

    ! Sun (scaled units)
    bodies(1)%mass = 1000.0_dp
    bodies(1)%pos = [0.0_dp, 0.0_dp, 0.0_dp]
    bodies(1)%vel = [0.0_dp, 0.0_dp, 0.0_dp]
    bodies(1)%radius = 5.0_dp

    ! Distances and Masses (approximate relative values)
    distances = [10.0, 15.0, 20.0, 30.0, 50.0, 70.0, 90.0, 110.0]
    masses = [0.01, 0.05, 0.05, 0.03, 0.5, 0.4, 0.1, 0.1]

    do i = 1, n_planets
       bodies(i+1)%mass = masses(i)
       bodies(i+1)%pos = [distances(i), 0.0_dp, 0.0_dp]
       ! Circular orbit: v = sqrt(G*M/r)
       v_orbit(i) = sqrt(1.0_dp * bodies(1)%mass / distances(i))
       bodies(i+1)%vel = [0.0_dp, v_orbit(i), 0.0_dp]
       bodies(i+1)%radius = 1.0_dp
    end do
  end subroutine setup_solar_system

  subroutine setup_collision(bodies, n, G)
    type(body_t), dimension(:), allocatable, intent(out) :: bodies
    integer, intent(in) :: n
    real(dp), intent(in) :: G
    type(body_t), dimension(:), allocatable :: gal1, gal2
    integer :: n2, i

    n2 = n / 2
    ! Setup two galaxies
    call setup_galaxy(gal1, n2, 30.0_dp, 500.0_dp, G)
    call setup_galaxy(gal2, n - n2, 30.0_dp, 500.0_dp, G)

    ! Offset positions and give them relative velocity
    do i = 1, n2
       gal1(i)%pos = gal1(i)%pos + [-50.0_dp, 0.0_dp, 0.0_dp]
       gal1(i)%vel = gal1(i)%vel + [10.0_dp, 0.0_dp, 0.0_dp]
    end do

    do i = 1, n - n2
       gal2(i)%pos = gal2(i)%pos + [50.0_dp, 0.0_dp, 0.0_dp]
       gal2(i)%vel = gal2(i)%vel + [-10.0_dp, 0.0_dp, 0.0_dp]
    end do

    allocate(bodies(n))
    bodies(1:n2) = gal1
    bodies(n2+1:n) = gal2
  end subroutine setup_collision

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
