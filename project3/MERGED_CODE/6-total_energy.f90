double precision function E(Natoms, velocity, mass, distance, epsilon, sigma)
    implicit none
    integer, intent(in) :: Natoms
    double precision, intent(in) :: velocity(Natoms,3), mass(Natoms), distance(Natoms,Natoms)
    double precision, intent(in) :: epsilon, sigma
    double precision :: kinetic_energy, potential_energy

    interface
        double precision function V(epsilon, sigma, Natoms, distance)
            integer, intent(in) :: Natoms
            double precision, intent(in) :: epsilon, sigma, distance(Natoms,Natoms)
        end function V
    end interface

    interface
        double precision function T(Natoms, velocity, mass)
            integer, intent(in) :: Natoms
            double precision, intent(in) :: velocity(Natoms,3)
            double precision, intent(in) :: mass(Natoms)
        end function T
    end interface

    kinetic_energy = T(Natoms, velocity, mass)
    potential_energy = V(epsilon, sigma, Natoms, distance)
    
    E = kinetic_energy + potential_energy
end function E
