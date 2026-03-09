module physics
  use kinds, only: dp
  use particle, only: body_t
  implicit none
  private

  public :: compute_forces_pairwise
  public :: compute_forces_bh
  public :: update_positions
  public :: update_velocities

  contains

  subroutine compute_forces_bh(bodies, root, G, theta)
    use tree, only: node_t
    type(body_t), dimension(:), intent(inout) :: bodies
    type(node_t), pointer, intent(in) :: root
    real(dp), intent(in) :: G, theta
    integer :: i, n

    n = size(bodies)
    ! Reset accelerations
    do i = 1, n
       bodies(i)%acc = [0.0_dp, 0.0_dp, 0.0_dp]
       call traverse_tree(bodies(i), root, G, theta)
    end do
  end subroutine compute_forces_bh

  recursive subroutine traverse_tree(body, node, G, theta)
    use tree, only: node_t
    type(body_t), intent(inout) :: body
    type(node_t), pointer, intent(in) :: node
    real(dp), intent(in) :: G, theta
    real(dp), dimension(3) :: r_vec
    real(dp) :: r_mag_sq, r_mag, force_mag
    integer :: i

    if (.not. associated(node)) return
    if (node%mass <= 0.0_dp) return

    ! Check if it's the same body (if leaf)
    if (node%is_leaf) then
       if (associated(node%body)) then
          if (loc(body) == loc(node%body)) return
       end if
    end if

    r_vec = node%com - body%pos
    r_mag_sq = sum(r_vec**2)
    r_mag = sqrt(r_mag_sq)

    if (node%is_leaf .or. (node%size / r_mag < theta)) then
       ! Use center of mass of the node
       if (r_mag_sq > 0.0_dp) then
          force_mag = G * node%mass / (r_mag_sq * r_mag)
          body%acc = body%acc + force_mag * r_vec
       end if
    else
       ! Recurse into children
       do i = 1, 8
          call traverse_tree(body, node%children(i), G, theta)
       end do
    end if
  end subroutine traverse_tree

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
