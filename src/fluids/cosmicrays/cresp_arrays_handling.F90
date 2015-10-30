module cresp_arrays_handling
! pulled by COSM_RAY_ELECTRONS
  implicit none
  
  interface allocate_with_index
     module procedure alloc_arr_ind_int4
     module procedure alloc_arr_ind_real8
     module procedure alloc_arr_ind_logical
  end interface allocate_with_index
  
  interface deallocate_with_index
     module procedure dealloc_arr_int4
     module procedure dealloc_arr_real8 !! <  to be replaced by procedures in diagnostics.F90
     module procedure dealloc_arr_logical
  end interface deallocate_with_index
  
!=================================================================
!
contains
!
!
  subroutine alloc_arr_ind_int4(array, array_first_ind, array_size)
!   integer(kind=4), parameter      :: name_length = 64
  integer(kind=4)                              :: array_size, array_first_ind
  integer(kind=4), allocatable, dimension(:), intent(out) :: array
    
     if(.not. allocated(array)) allocate(array(array_first_ind:array_size))
    
  end subroutine alloc_arr_ind_int4

!-----------------------------------------------------------------
  subroutine alloc_arr_ind_real8(array, array_first_ind, array_size)
!   integer(kind=4), parameter      :: name_length = 64
  integer(kind=4)                              :: array_size, array_first_ind
  real(kind=8), allocatable, dimension(:), intent(out) :: array
    
    if(.not. allocated(array)) allocate(array(array_first_ind:array_size))
    
  end subroutine alloc_arr_ind_real8
  
  subroutine alloc_arr_ind_logical(array,array_first_ind,array_size)
  integer(kind=4)                                 :: array_size, array_first_ind
  logical, allocatable, dimension(:), intent(out) :: array
  
    if (.not. allocated(array)) allocate(array(array_first_ind:array_size))
  
  end subroutine alloc_arr_ind_logical
!-----------------------------------------------------------------  
  subroutine dealloc_arr_int4(array)
  integer(kind=4), dimension(:), allocatable, intent(out) :: array
  
    if(allocated(array)) deallocate(array)
    
  end subroutine dealloc_arr_int4
!-----------------------------------------------------------------  

  subroutine dealloc_arr_real8(array)
  real(kind=8), dimension(:), allocatable, intent(out) :: array
  
    if(allocated(array)) deallocate(array)
  
  end subroutine dealloc_arr_real8
!-----------------------------------------------------------------  
  subroutine dealloc_arr_logical(array)
  logical, dimension(:), allocatable, intent(out) :: array
  
    if (allocated(array)) deallocate(array)
    
  end subroutine dealloc_arr_logical
end module cresp_arrays_handling
