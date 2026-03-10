module profiler
  implicit none
  private

  public :: start_timer, stop_timer, print_profile_report

  real, save :: t_tree = 0.0, t_force = 0.0, t_total = 0.0
  integer(8), save :: c_tree, c_force, c_total, rate

  contains

  subroutine start_timer(category)
     character(len=*), intent(in) :: category
     integer(8) :: c
     call system_clock(c, rate)
     if (category == "tree") then
        c_tree = c
     else if (category == "force") then
        c_force = c
     else if (category == "total") then
        c_total = c
     end if
  end subroutine start_timer

  subroutine stop_timer(category)
     character(len=*), intent(in) :: category
     integer(8) :: c2
     call system_clock(c2)
     if (category == "tree") then
        t_tree = t_tree + real(c2 - c_tree) / real(rate)
     else if (category == "force") then
        t_force = t_force + real(c2 - c_force) / real(rate)
     else if (category == "total") then
        t_total = t_total + real(c2 - c_total) / real(rate)
     end if
  end subroutine stop_timer

  subroutine print_profile_report()
     print *, ""
     print *, "--- Performance Report ---"
     print "(A, F10.3, A)", " Tree Building:    ", t_tree, " s"
     print "(A, F10.3, A)", " Force Calculation: ", t_force, " s"
     print "(A, F10.3, A)", " Total Sim Time:    ", t_total, " s"
     if (t_total > 0) then
        print "(A, F10.1, A)", " Efficiency:        ", (t_tree + t_force)/t_total * 100.0, " % accounted"
     end if
     print *, "--------------------------"
  end subroutine print_profile_report

end module profiler
