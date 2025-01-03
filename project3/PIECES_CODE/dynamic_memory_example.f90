program molecular_dynamics_example
    implicit none

    ! Declare the dynamic 2D array and necessary variables
    double precision, allocatable :: a(:,:)  ! 2D array
    integer :: m, n                          ! Dimensions of the array
    integer :: i_stat                        ! Status code for allocation
    integer :: i, j                          ! Iteration variables for filling the array

    ! Example dimensions
    m = 3   ! Number of rows
    n = 4   ! Number of columns

    ! Allocate the 2D array
    allocate(a(m, n), stat=i_stat)

    ! Check if the allocation was successful
    if (i_stat /= 0) then
        print *, "Error: Memory allocation failed!"
        stop
    else
        print *, "Memory allocated successfully for array a(", m, ",", n, ")"
    end if

    ! Fill the 2D array with example values
    do i = 1, m
        do j = 1, n
            a(i, j) = i + j * 0.1d0  ! Example values based on row and column indices
        end do
    end do

    ! Print the 2D array to confirm its contents
    print *, "Array a:"
    do i = 1, m
        print *, a(i, :)
    end do

    ! Deallocate the 2D array
    deallocate(a)

    ! Confirm that the memory has been deallocated
    print *, "Memory deallocated successfully."
end program molecular_dynamics_example
