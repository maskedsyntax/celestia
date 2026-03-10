module profiler
  implicit none
  private

  public :: start_timer, stop_timer, print_profile_report

  real, save :: t_tree = 0.0, t_force = 0.0, t_total = 0.0
  integer(8) :: c1, c2, rate

  contains

  subroutine start_timer()
     call system_clock(c1, rate)
  end subroutine start_timer

  subroutine stop_timer(category)
     character(len=*), intent(in) :: category
     call system_clock(c2)
     if (category == "tree") then
        t_tree = t_tree + real(c2 - c1) / real(rate)
     else if (category == "force") then
        t_force = t_force + real(c2 - c1) / real(rate)
     else if (category == "total") then
        t_total = t_total + real(c2 - c1) / real(rate)
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
