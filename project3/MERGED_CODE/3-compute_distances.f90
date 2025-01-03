subroutine compute_distances(Natoms, coord, distance)
    implicit none
    integer, intent(in) :: Natoms
    double precision, intent(in) :: coord(Natoms, 3)
    double precision, intent(out) :: distance(Natoms, Natoms)
    integer :: i, j
    double precision :: dx, dy, dz

    distance = 0.0d0
    do i = 1, Natoms
        do j = i + 1, Natoms
            dx = coord(i, 1) - coord(j, 1)
            dy = coord(i, 2) - coord(j, 2)
            dz = coord(i, 3) - coord(j, 3)
            distance(i, j) = sqrt(dx**2 + dy**2 + dz**2)
            if (distance(i, j) < 1.0d-12) then
                print *, "Warning: Atoms too close. Setting minimum distance between atoms", i, "and", j
                distance(i, j) = 1.0d-12
            end if
            distance(j, i) = distance(i, j)
        end do
    end do
end subroutine compute_distances

