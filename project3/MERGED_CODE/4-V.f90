double precision function V(epsilon, sigma, Natoms, distance)
    implicit none

    ! Input arguments
    integer, intent(in) :: Natoms                ! Number of atoms
    double precision, intent(in) :: epsilon      ! Lennard-Jones epsilon parameter
    double precision, intent(in) :: sigma        ! Lennard-Jones sigma parameter
    double precision, intent(in) :: distance(Natoms, Natoms)  ! Pairwise distance matrix

    ! Local variables
    integer :: i, j                              ! Loop counters
    double precision :: r6, r12                  ! Intermediate calculations for r^6 and r^12

    ! Initialize the potential energy
    V = 0.0d0

    ! Loop over all unique pairs of atoms
    do i = 1, Natoms - 1
        do j = i + 1, Natoms
            if (distance(i, j) > 1.0d-10) then   ! Avoid division by zero or very small distances
                r6 = (sigma / distance(i, j))**6 ! Calculate (sigma / r)^6
                r12 = r6**2                      ! Calculate (sigma / r)^12
                V = V + 4.0d0 * epsilon * (r12 - r6)  ! Lennard-Jones potential formula
            end if
        end do
    end do
end function V

