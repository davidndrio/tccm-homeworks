! Function to calculate the total kinetic energy of the system
double precision function T(Natoms, velocity, mass)
    implicit none
    integer, intent(in) :: Natoms
    double precision, intent(in) :: velocity(Natoms,3)
    double precision, intent(in) :: mass(Natoms)
    integer :: i
    double precision :: v_squared

    T = 0.0d0
    do i = 1, Natoms
        v_squared = sum(velocity(i,:)**2)  ! Sum of squares of velocity components
        T = T + 0.5d0 * mass(i) * v_squared
    end do
end function T
