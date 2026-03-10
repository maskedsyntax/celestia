program celestia
  use kinds, only: dp
  use particle, only: body_t
  use tree, only: node_t, build_tree, delete_tree, update_node_mass
  use physics, only: compute_forces_bh, update_positions, update_velocities, integrate_rk4, compute_total_energy
  use collisions, only: handle_collisions
  use initial_conditions, only: setup_galaxy, setup_solar_system, setup_collision
  use io_utils, only: export_csv
  implicit none

  type(body_t), dimension(:), allocatable, target :: bodies
  real(dp), dimension(:, :), allocatable :: acc_old
  type(node_t), pointer :: root => null()
  real(dp) :: G, dt, total_time, theta, energy_initial, energy_current
  integer :: i, n, steps, s, n_active
  character(len=32) :: arg, val, scenario, integrator

  ! Default parameters
  G = 1.0_dp
  dt = 0.0005_dp
  total_time = 5.0_dp 
  n = 1000
  theta = 0.7_dp
  scenario = "galaxy"
  integrator = "verlet"

  ! Simple CLI Parser
  i = 1
  do while (i <= command_argument_count())
     call get_command_argument(i, arg)
     if (arg == "--n") then
        call get_command_argument(i+1, val)
        read(val, *) n
        i = i + 2
     else if (arg == "--time") then
        call get_command_argument(i+1, val)
        read(val, *) total_time
        i = i + 2
     else if (arg == "--dt") then
        call get_command_argument(i+1, val)
        read(val, *) dt
        i = i + 2
     else if (arg == "--theta") then
        call get_command_argument(i+1, val)
        read(val, *) theta
        i = i + 2
     else if (arg == "--scenario") then
        call get_command_argument(i+1, scenario)
        i = i + 2
     else if (arg == "--integrator") then
        call get_command_argument(i+1, integrator)
        i = i + 2
     else
        print *, "Unknown argument: ", trim(arg)
        print *, "Usage: ./celestia [--n 1000] [--time 10.0] [--dt 0.001] " // &
                 "[--theta 0.5] [--scenario galaxy] [--integrator verlet|rk4]"
        stop
     end if
  end do

  steps = int(total_time / dt)

  print "(A,I6,A,F8.3,A,F8.4,A,F5.2,A,A,A,A)", "Config: N=", n, ", Time=", total_time, &
        ", DT=", dt, ", Theta=", theta, ", Scenario=", trim(scenario), ", Integrator=", trim(integrator)

  ! Setup initial conditions
  if (trim(scenario) == "galaxy") then
     call setup_galaxy(bodies, n, 100.0_dp, 1000.0_dp, G)
  else if (trim(scenario) == "solar") then
     call setup_solar_system(bodies)
     n = size(bodies)
  else if (trim(scenario) == "collision") then
     call setup_collision(bodies, n, G)
  else
     print *, "Scenario not implemented: ", trim(scenario)
     stop
  end if

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
     if (trim(integrator) == "rk4") then
        call integrate_rk4(bodies(1:n_active), root, G, theta, dt)
     else
        ! Default: Velocity Verlet
        do i = 1, n_active
           acc_old(:, i) = bodies(i)%acc
        end do
        call update_positions(bodies(1:n_active), dt)
        call handle_collisions(bodies, n_active)
        call delete_tree(root)
        allocate(root)
        root%center = [0.0_dp, 0.0_dp, 0.0_dp]
        root%size = 500.0_dp
        call build_tree(root, bodies(1:n_active))
        call update_node_mass(root)
        call compute_forces_bh(bodies(1:n_active), root, G, theta)
        call update_velocities(bodies(1:n_active), acc_old(:, 1:n_active), dt)
     end if

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
