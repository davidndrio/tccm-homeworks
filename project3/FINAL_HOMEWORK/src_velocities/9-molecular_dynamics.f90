subroutine molecular_dynamics(Natoms, coord, velocity, acceleration, mass, epsilon, sigma, distance)
    implicit none

    ! Input arguments
    integer, intent(in) :: Natoms
    double precision, intent(inout) :: coord(Natoms, 3), velocity(Natoms, 3), acceleration(Natoms, 3)
    double precision, intent(in) :: mass(Natoms), epsilon, sigma
    character(len=2), dimension(:), allocatable :: symbol  
    double precision, intent(inout) :: distance(Natoms, Natoms)

    ! Parameters for the simulation
    double precision, parameter :: dt = 0.2d0        ! Time step in fs
    integer, parameter :: total_steps = 1000
    integer, parameter :: M = 10

    ! Local variables
    integer :: step, atom
    double precision :: kinetic_energy, potential_energy, total_energy
    double precision :: temp_velocity(Natoms, 3), new_acceleration(Natoms, 3)
    integer :: io_status


    ! Interfaces for external functions and subroutines
    interface
        double precision function T(Natoms, velocity, mass)
            integer, intent(in) :: Natoms
            double precision, intent(in) :: velocity(Natoms, 3)
            double precision, intent(in) :: mass(Natoms)
        end function T

        double precision function V(epsilon, sigma, Natoms, distance)
            integer, intent(in) :: Natoms
            double precision, intent(in) :: epsilon, sigma, distance(Natoms, Natoms)
        end function V

        subroutine compute_distances(Natoms, coord, distance)
            integer, intent(in) :: Natoms
            double precision, intent(in) :: coord(Natoms, 3)
            double precision, intent(out) :: distance(Natoms, Natoms)
        end subroutine compute_distances

        subroutine compute_acc(Natoms, coord, mass, acceleration, epsilon, sigma)
            integer, intent(in) :: Natoms
            double precision, intent(in) :: coord(Natoms, 3), mass(Natoms), epsilon, sigma
            double precision, intent(out) :: acceleration(Natoms, 3)
        end subroutine compute_acc
        
        function get_symbol(mass) result(symbols)
            double precision, intent(in) :: mass
            character(len=2) :: symbols
        end function get_symbol

    end interface
    
    ! Initialize velocities to zero
    do atom = 1, Natoms
        velocity(atom, 1) = 0.05d0
        velocity(atom, 2) = 0.05d0
        velocity(atom, 3) = 0.05d0
    end do


    ! Open output file for trajectories
    open(unit=20, file="trajectories.xyz", status="replace", action="write", iostat=io_status)
    if (io_status /= 0) then
        print *, "Error: Unable to open output file trajectories.xyz."
        stop
    end if

    ! Compute initial distances and accelerations

    call compute_distances(Natoms, coord, distance)

    call compute_acc(Natoms, coord, mass, acceleration, epsilon, sigma)

    ! Compute initial energies

    kinetic_energy = T(Natoms, velocity, mass)

    potential_energy = V(epsilon, sigma, Natoms, distance)

    total_energy = kinetic_energy + potential_energy

    ! Print energies before simulation

    write(*,*) "Initial energies:"

    write(*,*) "KE=", kinetic_energy, ", PE=", potential_energy, ", TE=", total_energy
        
    ! Main molecular dynamics loop
    do step = 1, total_steps
        ! Update positions
        coord = coord + (velocity * dt) / 1.0d6 + (0.5d0 * acceleration * dt**2) / 1.0d21

        ! Compute temporary velocities
        temp_velocity = velocity + (0.5d0 * acceleration * dt) / 1.0d15

        ! Update distances
        call compute_distances(Natoms, coord, distance)

        ! Compute new accelerations
        call compute_acc(Natoms, coord, mass, new_acceleration, epsilon, sigma)

        ! Update velocities
        velocity = temp_velocity + (0.5d0 * new_acceleration * dt) / 1.0d15

        ! Update accelerations
        acceleration = new_acceleration

        ! Compute energies
        kinetic_energy = T(Natoms, velocity, mass)
        potential_energy = V(epsilon, sigma, Natoms, distance)
        total_energy = kinetic_energy + potential_energy

        ! Write to output files
        if (mod(step, M) == 0) then
            write(20, '(I0)') Natoms    
            write(20, '(A,ES15.6,A,ES15.6,A,ES15.6,A, I6, A, F15.3, A)') &
                "Energies (in J/mol): KE=", kinetic_energy, ", PE=", potential_energy, &
                ", TE=", total_energy, ", Step=", step, ", Time=", step*dt, " fs"
            do atom = 1, Natoms
                write(20, '(A,3F15.8)') get_symbol(mass(atom)), coord(atom, 1), coord(atom, 2), coord(atom, 3)
            end do
            write(20, *)
        end if
    end do

    close(20)
end subroutine molecular_dynamics


