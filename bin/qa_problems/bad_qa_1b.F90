module c

  use a, only: b

  implicit none

  private

  type :: d
  contains
  end type

  type(b) :: f

end
