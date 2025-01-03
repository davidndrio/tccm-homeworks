double precision function V(epsilon, sigma, Natoms, distance)
    implicit none
    integer, intent(in) :: Natoms
    double precision, intent(in) :: epsilon, sigma
    double precision, intent(in) :: distance(Natoms, Natoms)
    double precision :: r, r6, r12
    integer :: i, j

    V = 0.0d0
    do i = 1, Natoms
        do j = i + 1, Natoms
            r = distance(i, j)
            r6 = (sigma / r)**6
            r12 = r6**2
            V = V + 4.0d0 * epsilon * (r12 - r6)
        end do
    end do
end function V

