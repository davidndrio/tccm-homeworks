subroutine compute_distances(Natoms, coord, distance)
    implicit none

    ! Input arguments
    integer, intent(in) :: Natoms  ! Number of atoms
    double precision, intent(in) :: coord(Natoms, 3)  ! Atomic coordinates (x, y, z)

    ! Output arguments
    double precision, intent(out) :: distance(Natoms, Natoms)  ! Pairwise distance matrix

    ! Local variables
    integer :: i, j  ! Loop counters
    double precision :: dx, dy, dz  ! Differences in x, y, and z coordinates

    ! Initialize the distance matrix to zero
    distance = 0.0d0

    ! Compute the pairwise distances
    do i = 1, Natoms
        do j = i + 1, Natoms  ! Only compute for j > i to avoid redundant calculations
            dx = coord(i, 1) - coord(j, 1)  ! Difference in x-coordinates
            dy = coord(i, 2) - coord(j, 2)  ! Difference in y-coordinates
            dz = coord(i, 3) - coord(j, 3)  ! Difference in z-coordinates
            distance(i, j) = sqrt(dx**2 + dy**2 + dz**2)  ! Euclidean distance
            distance(j, i) = distance(i, j)  ! Symmetric matrix: distance(i, j) = distance(j, i)
        end do
    end do
end subroutine compute_distances

