program celestia
  use kinds, only: dp
  use particle, only: body_t
  use physics, only: compute_forces_pairwise, update_positions, update_velocities
  implicit none

  type(body_t), dimension(:), allocatable :: bodies
  real(dp), dimension(:, :), allocatable :: acc_old
  real(dp) :: G, dt, total_time
  integer :: i, n, steps, s

  ! Simulation parameters
  G = 1.0_dp
  dt = 0.01_dp
  total_time = 1.0_dp
  steps = int(total_time / dt)
  n = 2

  allocate(bodies(n))
  allocate(acc_old(3, n))

  ! Initial conditions: Heavy central mass and a smaller orbiting body
  bodies(1)%mass = 100.0_dp
  bodies(1)%pos = [0.0_dp, 0.0_dp, 0.0_dp]
  bodies(1)%vel = [0.0_dp, 0.0_dp, 0.0_dp]
  bodies(1)%radius = 1.0_dp

  bodies(2)%mass = 1.0_dp
  bodies(2)%pos = [10.0_dp, 0.0_dp, 0.0_dp]
  ! Orbital velocity v = sqrt(G*M/r)
  bodies(2)%vel = [0.0_dp, sqrt(G * bodies(1)%mass / 10.0_dp), 0.0_dp]
  bodies(2)%radius = 0.1_dp

  ! Initial force calculation
  call compute_forces_pairwise(bodies, G)

  print *, "Starting simulation..."
  print *, "Step: 0, Time: 0.0, Body 2 Pos: ", bodies(2)%pos

  ! Simulation loop using Velocity Verlet
  do s = 1, steps
     ! 1. Store old acceleration (at time t)
     do i = 1, n
        acc_old(:, i) = bodies(i)%acc
     end do

     ! 2. Update positions: r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2
     call update_positions(bodies, dt)

     ! 3. Compute new forces/accelerations: a(t+dt)
     call compute_forces_pairwise(bodies, G)

     ! 4. Update velocities: v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
     call update_velocities(bodies, acc_old, dt)

     ! Print every 20 steps
     if (mod(s, 20) == 0) then
        print *, "Step: ", s, " Time: ", real(s, dp) * dt, " Body 2 Pos: ", bodies(2)%pos
     end if
  end do

  print *, "Simulation completed."
  deallocate(bodies, acc_old)

end program celestia
