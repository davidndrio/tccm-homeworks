subroutine molecular_dynamics(Natoms, coord, velocity, acceleration, mass, epsilon, sigma, distance)
    implicit none

    ! Input arguments
    integer, intent(in) :: Natoms                     ! Number of atoms
    double precision, intent(inout) :: coord(Natoms, 3)  ! Atomic coordinates in 3D space
    double precision, intent(inout) :: velocity(Natoms, 3)  ! Velocities in 3D space
    double precision, intent(inout) :: acceleration(Natoms, 3)  ! Accelerations in 3D space
    double precision, intent(inout) :: distance(Natoms, Natoms)  ! Pairwise distances
    double precision, intent(in) :: mass(Natoms)          ! Mass of each atom
    double precision, intent(in) :: epsilon, sigma        ! Lennard-Jones parameters

    ! Parameters for the simulation
    double precision, parameter :: dt = 0.005d0          ! Time step
    integer, parameter :: total_steps = 1000            ! Total number of time steps
    integer, parameter :: M = 10                        ! Print frequency
    integer :: step, atom                               ! Loop variables

    ! Local variables
    double precision :: kinetic_energy, potential_energy, total_energy
    double precision :: temp_velocity(Natoms, 3), new_acceleration(Natoms, 3)
    integer :: io_status

    ! Interfaces for external subroutines and functions
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
    end interface

    ! Initialize velocities with a small random scale
    call random_seed()
    do atom = 1, Natoms
        call random_number(velocity(atom, 1))
        call random_number(velocity(atom, 2))
        call random_number(velocity(atom, 3))
        velocity(atom, :) = (velocity(atom, :) - 0.5d0) * 1.0d-8
    end do

    ! Open file for trajectory output
    open(unit=20, file="trajectories.xyz", status="replace", action="write", iostat=io_status)
    if (io_status /= 0) then
        print *, "Error: Unable to open output file."
        stop
    end if

    ! Compute initial distances and accelerations
    call compute_distances(Natoms, coord, distance)
    call compute_acc(Natoms, coord, mass, acceleration, epsilon, sigma)

    ! Molecular dynamics loop
    do step = 1, total_steps
        ! Update positions
        coord = coord + velocity * dt + 0.5d0 * acceleration * dt**2

        ! Compute temporary velocities
        temp_velocity = velocity + 0.5d0 * acceleration * dt

        ! Update distances
        call compute_distances(Natoms, coord, distance)

        ! Compute new accelerations
        call compute_acc(Natoms, coord, mass, new_acceleration, epsilon, sigma)

        ! Update velocities
        velocity = temp_velocity + 0.5d0 * new_acceleration * dt

        ! Update accelerations
        acceleration = new_acceleration

        ! Compute energies
        kinetic_energy = T(Natoms, velocity, mass)
        potential_energy = V(epsilon, sigma, Natoms, distance)
        total_energy = kinetic_energy + potential_energy

        ! Write to output file
        write(20, '(A,I6,A,F25.16,A,F25.16,A,F25.16)') "Step ", step, ": KE=", kinetic_energy, &
            ", PE=", potential_energy, ", TE=", total_energy
        write(20, '(I6)') Natoms
        write(20, '(A)') "Comment line"
        do atom = 1, Natoms
            write(20, '(A,3F15.8)') "Ar", coord(atom, 1), coord(atom, 2), coord(atom, 3)
        end do

        ! Print to console every M steps
        if (mod(step, M) == 0) then
            write(*, '(A,I6,A,F25.16,A,F25.16,A,F25.16)') "Step ", step, ": KE=", kinetic_energy, &
                ", PE=", potential_energy, ", TE=", total_energy
        end if
    end do

    ! Close output file
    close(20)
end subroutine molecular_dynamics

