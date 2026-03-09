module tree
  use kinds, only: dp
  use particle, only: body_t
  implicit none
  private

  public :: node_t
  public :: build_tree
  public :: delete_tree
  public :: update_node_mass

  type :: node_t
     real(dp), dimension(3) :: center = [0.0_dp, 0.0_dp, 0.0_dp]
     real(dp) :: size = 0.0_dp
     real(dp) :: mass = 0.0_dp
     real(dp), dimension(3) :: com = [0.0_dp, 0.0_dp, 0.0_dp] ! Center of Mass
     type(body_t), pointer :: body => null()
     type(node_t), pointer, dimension(:) :: children => null()
     logical :: is_leaf = .true.
  end type node_t

  contains

  recursive subroutine build_tree(node, bodies)
    type(node_t), pointer, intent(inout) :: node
    type(body_t), dimension(:), target, intent(in) :: bodies
    integer :: i, j, count, n
    integer, dimension(size(bodies)) :: octant_map
    integer, dimension(8) :: octant_counts
    real(dp), dimension(3) :: new_center
    real(dp) :: new_size
    type(body_t), dimension(:), allocatable, target :: sub_bodies
    type(node_t), pointer :: p_child

    n = size(bodies)
    if (n == 0) return

    if (n == 1) then
       node%body => bodies(1)
       node%mass = bodies(1)%mass
       node%com = bodies(1)%pos
       node%is_leaf = .true.
       return
    end if

    ! Internal node
    node%is_leaf = .false.
    allocate(node%children(8))
    new_size = node%size / 2.0_dp
    
    octant_counts = 0
    do i = 1, n
       j = 0
       if (bodies(i)%pos(1) > node%center(1)) j = ibset(j, 0)
       if (bodies(i)%pos(2) > node%center(2)) j = ibset(j, 1)
       if (bodies(i)%pos(3) > node%center(3)) j = ibset(j, 2)
       octant_map(i) = j + 1
       octant_counts(j+1) = octant_counts(j+1) + 1
    end do

    do i = 1, 8
       ! Boundaries
       new_center(1) = node%center(1) + merge(new_size/2.0_dp, -new_size/2.0_dp, btest(i-1, 0))
       new_center(2) = node%center(2) + merge(new_size/2.0_dp, -new_size/2.0_dp, btest(i-1, 1))
       new_center(3) = node%center(3) + merge(new_size/2.0_dp, -new_size/2.0_dp, btest(i-1, 2))
       
       node%children(i)%center = new_center
       node%children(i)%size = new_size
       node%children(i)%is_leaf = .true. ! Default

       if (octant_counts(i) > 0) then
          allocate(sub_bodies(octant_counts(i)))
          count = 0
          do j = 1, n
             if (octant_map(j) == i) then
                count = count + 1
                sub_bodies(count) = bodies(j)
             end if
          end do
          p_child => node%children(i)
          call build_tree(p_child, sub_bodies)
          deallocate(sub_bodies)
       end if
    end do
  end subroutine build_tree

  recursive subroutine delete_tree(node)
    type(node_t), pointer, intent(inout) :: node
    integer :: i
    if (.not. associated(node)) return
    if (associated(node%children)) then
       do i = 1, 8
          call delete_children(node%children(i))
       end do
       deallocate(node%children)
       nullify(node%children)
    end if
    deallocate(node)
    nullify(node)
  end subroutine delete_tree

  recursive subroutine delete_children(node)
    type(node_t), intent(inout) :: node
    integer :: i
    if (associated(node%children)) then
       do i = 1, 8
          call delete_children(node%children(i))
       end do
       deallocate(node%children)
       nullify(node%children)
    end if
  end subroutine delete_children

  recursive subroutine update_node_mass(node)
    type(node_t), pointer, intent(inout) :: node
    integer :: i
    type(node_t), pointer :: p_child
    if (.not. associated(node)) return
    if (node%is_leaf) then
       if (associated(node%body)) then
          node%mass = node%body%mass
          node%com = node%body%pos
       else
          node%mass = 0.0_dp
          node%com = [0.0_dp, 0.0_dp, 0.0_dp]
       end if
       return
    end if

    node%mass = 0.0_dp
    node%com = [0.0_dp, 0.0_dp, 0.0_dp]
    
    do i = 1, 8
       p_child => node%children(i)
       call update_node_mass(p_child)
       if (node%children(i)%mass > 0.0_dp) then
          node%mass = node%mass + node%children(i)%mass
          node%com = node%com + node%children(i)%mass * node%children(i)%com
       end if
    end do
    
    if (node%mass > 0.0_dp) then
       node%com = node%com / node%mass
    end if
  end subroutine update_node_mass

end module tree
