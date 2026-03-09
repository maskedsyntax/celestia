module io_utils
  use kinds, only: dp
  use particle, only: body_t
  implicit none
  private

  public :: export_csv
  public :: export_binary

  contains

  subroutine export_csv(filename, bodies, time)
    character(len=*), intent(in) :: filename
    type(body_t), dimension(:), intent(in) :: bodies
    real(dp), intent(in) :: time
    integer :: unit, i

    open(newunit=unit, file=filename, status='unknown', position='append')
    do i = 1, size(bodies)
       write(unit, '(F12.6, ",", I8, ",", F12.6, ",", F12.6, ",", F12.6, ",", F12.6, ",", F12.6, ",", F12.6, ",", F12.6)') &
            time, i, bodies(i)%mass, bodies(i)%pos(1), bodies(i)%pos(2), bodies(i)%pos(3), &
            bodies(i)%vel(1), bodies(i)%vel(2), bodies(i)%vel(3)
    end do
    close(unit)
  end subroutine export_csv

  subroutine export_binary(filename, bodies, time)
    character(len=*), intent(in) :: filename
    type(body_t), dimension(:), intent(in) :: bodies
    real(dp), intent(in) :: time
    integer :: unit

    open(newunit=unit, file=filename, status='unknown', position='append', form='unformatted', access='stream')
    write(unit) time
    write(unit) size(bodies)
    write(unit) bodies
    close(unit)
  end subroutine export_binary

end module io_utils
