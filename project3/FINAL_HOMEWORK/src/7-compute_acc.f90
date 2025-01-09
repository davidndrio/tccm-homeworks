subroutine compute_acc(Natoms, coord, mass, acceleration, epsilon, sigma)
    implicit none

    ! Input arguments
    integer, intent(in) :: Natoms                     ! Number of atoms
    double precision, intent(in) :: coord(Natoms, 3)  ! Coordinates of atoms in 3D space
    double precision, intent(in) :: mass(Natoms)      ! Mass of each atom
    double precision, intent(in) :: epsilon, sigma    ! Lennard-Jones potential parameters

    ! Output arguments
    double precision, intent(out) :: acceleration(Natoms, 3)  ! Acceleration of atoms in 3D space

    ! Local variables
    integer :: i, j                                   ! Loop indices
    double precision :: dx, dy, dz                   ! Components of distance vector between particles
    double precision :: r2, r                        ! Squared and actual distance between particles
    double precision :: r6, r12                      ! Distance terms for Lennard-Jones potential
    double precision :: f                            ! Force magnitude

    ! Initialize accelerations to zero
    acceleration = 0.0d0

    ! Calculate forces between particles
    do i = 1, Natoms
        do j = i + 1, Natoms
            ! Calculate distance vector components
            dx = coord(i, 1) - coord(j, 1)
            dy = coord(i, 2) - coord(j, 2)
            dz = coord(i, 3) - coord(j, 3)
            r2 = dx**2 + dy**2 + dz**2  ! Squared distance

            ! Avoid division by zero
            if (r2 > 1.0d-10) then
                r = sqrt(r2)               ! Actual distance
                r6 = (sigma / r)**6        ! r^6 term
                r12 = r6**2                ! r^12 term
                f = 24.0d0 * 1.0d9 * epsilon * (2.0d0 * r12 - r6) / r  ! Force magnitude

                ! Update accelerations for particle i
                acceleration(i, 1) = acceleration(i, 1) + f * (dx/r) * 1000d0 / mass(i)
                acceleration(i, 2) = acceleration(i, 2) + f * (dy/r) * 1000d0 / mass(i)
                acceleration(i, 3) = acceleration(i, 3) + f * (dz/r) * 1000d0 / mass(i)

                ! Update accelerations for particle j (opposite direction)
                acceleration(j, 1) = acceleration(j, 1) - f * (dx/r) * 1000d0 / mass(j)
                acceleration(j, 2) = acceleration(j, 2) - f * (dy/r) * 1000d0 / mass(j)
                acceleration(j, 3) = acceleration(j, 3) - f * (dz/r) * 1000d0 / mass(j)
            end if
        end do
    end do
end subroutine compute_acc

