module collisions
  use kinds, only: dp
  use particle, only: body_t
  implicit none
  private

  public :: handle_collisions

  contains

  subroutine handle_collisions(bodies, n_active)
    type(body_t), dimension(:), intent(inout) :: bodies
    integer, intent(inout) :: n_active
    integer :: i, j, k
    real(dp) :: dist_sq, r_sum_sq
    logical, dimension(size(bodies)) :: to_delete

    to_delete = .false.

    do i = 1, n_active
       if (to_delete(i)) cycle
       do j = i + 1, n_active
          if (to_delete(j)) cycle
          
          dist_sq = sum((bodies(i)%pos - bodies(j)%pos)**2)
          r_sum_sq = (bodies(i)%radius + bodies(j)%radius)**2
          
          if (dist_sq <= r_sum_sq) then
             ! Collision detected! Merge j into i
             call merge_bodies(bodies(i), bodies(j))
             to_delete(j) = .true.
          end if
       end do
    end do

    ! Compact the array (simple version)
    k = 0
    do i = 1, n_active
       if (.not. to_delete(i)) then
          k = k + 1
          if (k /= i) then
             bodies(k) = bodies(i)
          end if
       end if
    end do
    n_active = k

  end subroutine handle_collisions

  subroutine merge_bodies(b1, b2)
    type(body_t), intent(inout) :: b1
    type(body_t), intent(in) :: b2
    real(dp) :: total_mass

    total_mass = b1%mass + b2%mass
    
    ! Conservation of momentum: v_new = (m1*v1 + m2*v2) / (m1+m2)
    b1%vel = (b1%mass * b1%vel + b2%mass * b2%vel) / total_mass
    
    ! Center of mass for new position
    b1%pos = (b1%mass * b1%pos + b2%mass * b2%pos) / total_mass
    
    ! Update mass and radius (assuming volume conservation for radius)
    ! Volume V = 4/3 * pi * r^3 -> r_new = (r1^3 + r2^3)^(1/3)
    b1%radius = (b1%radius**3 + b2%radius**3)**(1.0_dp/3.0_dp)
    b1%mass = total_mass

  end subroutine merge_bodies

end module collisions
