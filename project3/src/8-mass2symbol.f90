function get_symbol(mass) result(symbol)
    double precision, intent(in) :: mass  
    character(len=2) :: symbol 
    integer :: i
    double precision :: tolerance  ! Tolerance to mass values

    ! Predefined masses and their corresponding symbols
    double precision :: masses(1) = [39.948]  ! Mass array (can be expanded)
    character(len=2) :: symbols(1) = ['Ar']   ! Symbol array (can be expanded)

    symbol = 'Unknown'  ! Default symbol if no match is found

    tolerance = 1.0d-3  

    ! Finds the index of the mass and returns the corresponding symbol
    do i = 1, size(masses)
        if (abs(mass - masses(i)) < tolerance) then  ! To allow for small differences
            symbol = symbols(i)  ! Assigns the corresponding symbol
            return  
        end if
    end do
end function

