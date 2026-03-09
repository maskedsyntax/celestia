module physics
  use kinds, only: dp
  use particle, only: body_t
  implicit none
  private

  public :: compute_forces_pairwise
  public :: update_positions
  public :: update_velocities

  contains

  subroutine compute_forces_pairwise(bodies, G)
    type(body_t), dimension(:), intent(inout) :: bodies
    real(dp), intent(in) :: G
    integer :: i, j, n
    real(dp), dimension(3) :: r_vec
    real(dp) :: r_mag_sq, r_mag_inv_cube, force_mag

    n = size(bodies)
    
    ! Reset accelerations
    do i = 1, n
       bodies(i)%acc = [0.0_dp, 0.0_dp, 0.0_dp]
    end do

    ! Pairwise computation
    do i = 1, n
       do j = i + 1, n
          r_vec = bodies(j)%pos - bodies(i)%pos
          r_mag_sq = sum(r_vec**2)
          
          ! Softening could be added here later (good-to-have)
          if (r_mag_sq > 0.0_dp) then
             r_mag_inv_cube = 1.0_dp / (sqrt(r_mag_sq) * r_mag_sq)
             
             ! Force on body i from body j: F = G * m_i * m_j * r_vec / |r|^3
             ! Acceleration on body i: a_i = F / m_i = G * m_j * r_vec / |r|^3
             force_mag = G * r_mag_inv_cube
             
             bodies(i)%acc = bodies(i)%acc + force_mag * bodies(j)%mass * r_vec
             bodies(j)%acc = bodies(j)%acc - force_mag * bodies(i)%mass * r_vec
          end if
       end do
    end do
  end subroutine compute_forces_pairwise

  subroutine update_positions(bodies, dt)
    type(body_t), dimension(:), intent(inout) :: bodies
    real(dp), intent(in) :: dt
    integer :: i, n

    n = size(bodies)
    do i = 1, n
       bodies(i)%pos = bodies(i)%pos + bodies(i)%vel * dt + 0.5_dp * bodies(i)%acc * dt**2
    end do
  end subroutine update_positions

  subroutine update_velocities(bodies, acc_old, dt)
    type(body_t), dimension(:), intent(inout) :: bodies
    real(dp), dimension(:, :), intent(in) :: acc_old ! acc at time t
    real(dp), intent(in) :: dt
    integer :: i, n

    n = size(bodies)
    do i = 1, n
       bodies(i)%vel = bodies(i)%vel + 0.5_dp * (acc_old(:, i) + bodies(i)%acc) * dt
    end do
  end subroutine update_velocities

end module physics
