double precision function T(Natoms, velocity, mass)
    implicit none

    ! Input arguments
    integer, intent(in) :: Natoms                ! Number of atoms
    double precision, intent(in) :: velocity(Natoms, 3)  ! Velocities of atoms in 3D
    double precision, intent(in) :: mass(Natoms)          ! Mass of each atom

    ! Local variables
    integer :: i                                 ! Loop counter
    double precision :: v_squared                ! Squared velocity magnitude

    ! Initialize the total kinetic energy
    T = 0.0d0

    ! Loop over all atoms to calculate the total kinetic energy
    do i = 1, Natoms
        v_squared = velocity(i, 1)**2 + velocity(i, 2)**2 + velocity(i, 3)**2  ! Velocity magnitude squared
        T = T + 0.5d0 * mass(i) * v_squared  ! Add the kinetic energy contribution of atom i
    end do
end function T

