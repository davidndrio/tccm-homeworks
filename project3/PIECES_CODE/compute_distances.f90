subroutine compute_distances(Natoms, coord, distance)
    implicit none
    integer, intent(in) :: Natoms
    double precision, intent(in) :: coord(Natoms, 3)
    double precision, intent(out) :: distance(Natoms, Natoms)
    integer :: i, j
    double precision :: dx, dy, dz

    distance = 0.0
    do i = 1, Natoms
        do j = i + 1, Natoms
            dx = coord(i, 1) - coord(j, 1)
            dy = coord(i, 2) - coord(j, 2)
            dz = coord(i, 3) - coord(j, 3)
            distance(i, j) = sqrt(dx**2 + dy**2 + dz**2)
            distance(j, i) = distance(i, j)
        end do
    end do
end subroutine compute_distances
