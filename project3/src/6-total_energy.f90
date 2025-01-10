double precision function E(Natoms, velocity, mass, distance, epsilon, sigma)
    implicit none

    ! Input arguments
    integer, intent(in) :: Natoms                     ! Number of atoms
    double precision, intent(in) :: velocity(Natoms, 3)  ! Velocities of atoms in 3D
    double precision, intent(in) :: mass(Natoms)         ! Mass of each atom
    double precision, intent(in) :: distance(Natoms, Natoms)  ! Pairwise distances between atoms
    double precision, intent(in) :: epsilon, sigma       ! Lennard-Jones potential parameters

    ! Interface definitions for functions
    interface
        double precision function V(epsilon, sigma, Natoms, distance)
            integer, intent(in) :: Natoms
            double precision, intent(in) :: epsilon, sigma, distance(Natoms, Natoms)
        end function V

        double precision function T(Natoms, velocity, mass)
            integer, intent(in) :: Natoms
            double precision, intent(in) :: velocity(Natoms, 3)
            double precision, intent(in) :: mass(Natoms)
        end function T
    end interface

    ! Total energy is the sum of kinetic and potential energy
    E = T(Natoms, velocity, mass) + V(epsilon, sigma, Natoms, distance)
end function E

