module c

  use a, only: b

  implicit none

  private

  type :: d
  contains
    procedure :: e
  end type

  type(b) :: f

contains

  subroutine e(g)
    implicit none
    class(d) :: g
  end

end
