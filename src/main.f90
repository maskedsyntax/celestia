program celestia
  use kinds, only: dp
  use particle, only: body_t
  use tree, only: node_t, build_tree, delete_tree, update_node_mass
  use physics, only: compute_forces_bh, update_positions, update_velocities, compute_total_energy
  use collisions, only: handle_collisions
  use initial_conditions, only: setup_galaxy
  use io_utils, only: export_csv
  implicit none

  type(body_t), dimension(:), allocatable, target :: bodies
  real(dp), dimension(:, :), allocatable :: acc_old
  type(node_t), pointer :: root => null()
  real(dp) :: G, dt, total_time, theta, energy_initial, energy_current
  integer :: i, n, steps, s, n_active

  ! Simulation parameters
  G = 1.0_dp
  dt = 0.0005_dp ! Half the time step
  total_time = 10.0_dp 
  steps = int(total_time / dt)
  n = 1000
  theta = 0.7_dp ! More stable theta

  ! Setup initial conditions (Galaxy)
  call setup_galaxy(bodies, n, 100.0_dp, 1000.0_dp, G)
  n_active = n
  allocate(acc_old(3, n))

  ! 1. Build tree and initial force
  allocate(root)
  root%center = [0.0_dp, 0.0_dp, 0.0_dp]
  root%size = 500.0_dp ! Large enough to contain all bodies
  call build_tree(root, bodies(1:n_active))
  call update_node_mass(root)
  call compute_forces_bh(bodies(1:n_active), root, G, theta)
  
  energy_initial = compute_total_energy(bodies(1:n_active), G)

  print *, "Starting simulation with ", n_active, " bodies..."
  print *, "Initial Energy: ", energy_initial

  ! Simulation loop
  do s = 1, steps
     ! a. Store old acceleration
     do i = 1, n_active
        acc_old(:, i) = bodies(i)%acc
     end do

     ! b. Update positions
     call update_positions(bodies(1:n_active), dt)

     ! c. Handle collisions
     call handle_collisions(bodies, n_active)

     ! d. Rebuild tree for new positions
     call delete_tree(root)
     allocate(root)
     root%center = [0.0_dp, 0.0_dp, 0.0_dp]
     root%size = 500.0_dp
     call build_tree(root, bodies(1:n_active))
     call update_node_mass(root)

     ! e. Compute new forces (Barnes-Hut)
     call compute_forces_bh(bodies(1:n_active), root, G, theta)

     ! f. Update velocities
     call update_velocities(bodies(1:n_active), acc_old(:, 1:n_active), dt)

     ! Print every 50 steps
     if (mod(s, 50) == 0) then
        energy_current = compute_total_energy(bodies(1:n_active), G)
        print '(A,I5,A,F8.3,A,F15.4,A,I5)', "Step: ", s, " Time: ", real(s, dp)*dt, &
              " Energy: ", energy_current, " Active: ", n_active
        call export_csv("galaxy_step.csv", bodies(1:n_active), real(s, dp)*dt)
     end if
  end do

  print *, "Final Energy: ", compute_total_energy(bodies(1:n_active), G)
  print *, "Energy Drift: ", (compute_total_energy(bodies(1:n_active), G) - energy_initial) / energy_initial

  call delete_tree(root)
  deallocate(bodies, acc_old)

end program celestia
