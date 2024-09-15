! Test for avalability of Fortran 2018 constructs

program F2018

   implicit none

   logical :: l
   integer, dimension(0) :: i
   real, dimension(0, 0, 0, 0, 0, 0) :: r

   write(*, '(3(a,3(i0," ")))')"[l,i,r]: ranks: ", rank(l), rank(i), rank(r), &
        &                          " is_vector: ", get_rank(l), get_rank(i), get_rank(r), &
        &                          " size: ", -1, size(i), size(r)

contains

   integer function get_rank(v) result(r)

      implicit none

      class(*), dimension(..), intent(in) :: v

      select rank (v)
         rank (0)  ! scalar
            r = 0
         rank (1)  ! vectors
            r = 1
         rank default  ! other arrays
            r = -1
      end select

   end function get_rank

end program F2018
