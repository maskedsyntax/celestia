module physics
  use kinds, only: dp
  use particle, only: body_t
  implicit none
  private

  public :: compute_forces_pairwise
  public :: compute_forces_bh
  public :: update_positions
  public :: update_velocities
  public :: integrate_rk4
  public :: compute_total_energy

  real(dp), parameter :: eps_sq = 0.5_dp**2 ! Increased from 0.1 to 0.5 for stability

  contains

  subroutine integrate_rk4(bodies, root_info, G, theta, dt)
    use tree, only: node_t, build_tree, delete_tree, update_node_mass
    type(body_t), dimension(:), intent(inout) :: bodies
    type(node_t), pointer, intent(inout) :: root_info
    real(dp), intent(in) :: G, theta, dt
    
    type(body_t), dimension(size(bodies)) :: b_tmp
    real(dp), dimension(3, size(bodies)) :: k1_v, k1_a, k2_v, k2_a, k3_v, k3_a, k4_v, k4_a
    real(dp), dimension(3) :: root_center
    real(dp) :: root_size
    integer :: i, n

    n = size(bodies)
    root_center = root_info%center
    root_size = root_info%size

    ! k1 = f(t, y)
    do i = 1, n
       k1_v(:, i) = bodies(i)%vel
       k1_a(:, i) = bodies(i)%acc
    end do

    ! k2 = f(t + dt/2, y + dt/2 * k1)
    b_tmp = bodies
    do i = 1, n
       b_tmp(i)%pos = bodies(i)%pos + 0.5_dp * dt * k1_v(:, i)
       b_tmp(i)%vel = bodies(i)%vel + 0.5_dp * dt * k1_a(:, i)
    end do
    call rebuild_tree_local(b_tmp, root_info, G, theta, root_center, root_size)
    do i = 1, n
       k2_v(:, i) = b_tmp(i)%vel
       k2_a(:, i) = b_tmp(i)%acc
    end do

    ! k3 = f(t + dt/2, y + dt/2 * k2)
    b_tmp = bodies
    do i = 1, n
       b_tmp(i)%pos = bodies(i)%pos + 0.5_dp * dt * k2_v(:, i)
       b_tmp(i)%vel = bodies(i)%vel + 0.5_dp * dt * k2_a(:, i)
    end do
    call rebuild_tree_local(b_tmp, root_info, G, theta, root_center, root_size)
    do i = 1, n
       k3_v(:, i) = b_tmp(i)%vel
       k3_a(:, i) = b_tmp(i)%acc
    end do

    ! k4 = f(t + dt, y + dt * k3)
    b_tmp = bodies
    do i = 1, n
       b_tmp(i)%pos = bodies(i)%pos + dt * k3_v(:, i)
       b_tmp(i)%vel = bodies(i)%vel + dt * k3_a(:, i)
    end do
    call rebuild_tree_local(b_tmp, root_info, G, theta, root_center, root_size)
    do i = 1, n
       k4_v(:, i) = b_tmp(i)%vel
       k4_a(:, i) = b_tmp(i)%acc
    end do

    ! Final update: y = y + dt/6 * (k1 + 2k2 + 2k3 + k4)
    do i = 1, n
       bodies(i)%pos = bodies(i)%pos + (dt/6.0_dp) * (k1_v(:, i) + 2.0_dp*k2_v(:, i) + 2.0_dp*k3_v(:, i) + k4_v(:, i))
       bodies(i)%vel = bodies(i)%vel + (dt/6.0_dp) * (k1_a(:, i) + 2.0_dp*k2_a(:, i) + 2.0_dp*k3_a(:, i) + k4_a(:, i))
    end do
    
    ! Update accelerations for next step
    call rebuild_tree_local(bodies, root_info, G, theta, root_center, root_size)

  end subroutine integrate_rk4

  subroutine rebuild_tree_local(b, root, G_in, theta_in, root_center, root_size)
    use tree, only: node_t, build_tree, delete_tree, update_node_mass
    type(body_t), dimension(:), target, intent(inout) :: b
    type(node_t), pointer, intent(inout) :: root
    real(dp), intent(in) :: G_in, theta_in, root_size
    real(dp), dimension(3), intent(in) :: root_center
    
    call delete_tree(root)
    allocate(root)
    root%center = root_center
    root%size = root_size
    call build_tree(root, b)
    call update_node_mass(root)
    call compute_forces_bh(b, root, G_in, theta_in)
  end subroutine rebuild_tree_local

  function compute_total_energy(bodies, G) result(total_energy)
    type(body_t), dimension(:), intent(in) :: bodies
    real(dp), intent(in) :: G
    real(dp) :: total_energy
    real(dp) :: kinetic, potential
    real(dp) :: r_mag
    integer :: i, j, n

    n = size(bodies)
    kinetic = 0.0_dp
    potential = 0.0_dp

    do i = 1, n
       kinetic = kinetic + 0.5_dp * bodies(i)%mass * sum(bodies(i)%vel**2)
       do j = i + 1, n
          r_mag = sqrt(sum((bodies(i)%pos - bodies(j)%pos)**2) + eps_sq)
          potential = potential - G * bodies(i)%mass * bodies(j)%mass / r_mag
       end do
    end do

    total_energy = kinetic + potential
  end function compute_total_energy

  subroutine compute_forces_bh(bodies, root, G, theta)
    use tree, only: node_t
    type(body_t), dimension(:), intent(inout) :: bodies
    type(node_t), pointer, intent(in) :: root
    real(dp), intent(in) :: G, theta
    integer :: i, n

    n = size(bodies)
    ! Reset accelerations and traverse tree in parallel
    !$omp parallel do private(i) shared(bodies, root, G, theta)
    do i = 1, n
       bodies(i)%acc = [0.0_dp, 0.0_dp, 0.0_dp]
       call traverse_tree(bodies(i), root, G, theta)
    end do
    !$omp end parallel do
  end subroutine compute_forces_bh

  recursive subroutine traverse_tree(body, node, G, theta)
    use tree, only: node_t
    type(body_t), intent(inout), target :: body
    type(node_t), pointer, intent(in) :: node
    real(dp), intent(in) :: G, theta
    real(dp), dimension(3) :: r_vec
    real(dp) :: r_mag_sq, r_mag, force_mag
    integer :: i
    type(body_t), pointer :: p_body

    if (.not. associated(node)) return
    if (node%mass <= 0.0_dp) return

    ! Check if it's the same body (if leaf)
    if (node%is_leaf) then
       if (associated(node%body)) then
          p_body => body
          if (associated(p_body, node%body)) return
       end if
    end if

    r_vec = node%com - body%pos
    r_mag_sq = sum(r_vec**2) + eps_sq
    r_mag = sqrt(r_mag_sq)

    if (node%is_leaf .or. (node%size / r_mag < theta)) then
       ! Use center of mass of the node
       force_mag = G * node%mass / (r_mag_sq * r_mag)
       body%acc = body%acc + force_mag * r_vec
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
          r_mag_sq = sum(r_vec**2) + eps_sq
          
          r_mag_inv_cube = 1.0_dp / (sqrt(r_mag_sq) * r_mag_sq)
          force_mag = G * r_mag_inv_cube
          
          bodies(i)%acc = bodies(i)%acc + force_mag * bodies(j)%mass * r_vec
          bodies(j)%acc = bodies(j)%acc - force_mag * bodies(i)%mass * r_vec
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
