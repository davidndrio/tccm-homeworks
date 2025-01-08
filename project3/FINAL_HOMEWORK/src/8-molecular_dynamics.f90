subroutine molecular_dynamics(Natoms, coord, velocity, acceleration, mass, epsilon, sigma, distance)
    implicit none

    ! Interfaces for external functions
    interface
        double precision function T(Natoms, velocity, mass)
            integer, intent(in) :: Natoms
            double precision, intent(in) :: velocity(Natoms,3)
            double precision, intent(in) :: mass(Natoms)
        end function T
    end interface

    interface
        double precision function V(epsilon, sigma, Natoms, distance)
            integer, intent(in) :: Natoms
            double precision, intent(in) :: epsilon, sigma, distance(Natoms,Natoms)
        end function V
    end interface

    ! Input variables
    integer, intent(in) :: Natoms
    double precision, intent(inout) :: coord(Natoms, 3), velocity(Natoms, 3), acceleration(Natoms, 3)
    double precision, intent(inout) :: distance(Natoms, Natoms)
    double precision, intent(in) :: mass(Natoms)
    double precision, intent(in) :: epsilon, sigma

    ! Parameters for the simulation
    double precision, parameter :: dt = 0.2d0
    integer, parameter :: total_steps = 1000
    integer, parameter :: M = 10

    ! Local variables
    integer :: step, i
    double precision :: kinetic_energy, potential_energy, total_energy
    double precision :: temp_velocity(Natoms, 3)
    character(len=100) :: output_file
    integer :: io_status

    ! Open the output file
    output_file = "trajectories.xyz"
    open(unit=10, file=output_file, status='replace', action='write', iostat=io_status)
    if (io_status /= 0) then
        print *, "Error: Could not open the output file."
        stop
    end if

    ! Molecular dynamics
    do step = 1, total_steps

        ! 1. Update coordinates (r(n+1) = r(n) + v(n)*dt + 0.5*a(n)*dt^2)
        coord = coord + velocity * dt + 0.5d0 * acceleration * dt**2

        ! 2. Compute velocities part dependent on a(n) (v(n+**) = v(n) + 0.5*a(n)*dt)
        temp_velocity = velocity + 0.5d0 * acceleration * dt

        ! 3. Update distances because positions changed
        call compute_distances(Natoms, coord, distance)

        ! 4. Compute new acceleration (a(n+1))
        call compute_acc(Natoms, coord, mass, distance, acceleration, epsilon, sigma)

        ! 5. Compute final velocities (v(n+1) = v(n+**) + 0.5*a(n+1)*dt)
        velocity = temp_velocity + 0.5d0 * acceleration * dt

        ! Compute energies
        kinetic_energy = T(Natoms, velocity, mass)
        potential_energy = V(epsilon, sigma, Natoms, distance)
        total_energy = kinetic_energy + potential_energy

        ! Write trajectory to file every M steps
        if (mod(step, M) == 0) then 
            write(10, '(A, I6, A, ES16.6, A, ES16.6, A, ES16.6)') &
                "Step: ", step, ", Energy (J/mol): KE:", kinetic_energy, ", PE:", &
                potential_energy, ", TE:", total_energy
            do i = 1, Natoms
                write(10, '(A, 3(1x, ES16.8))') "Ar", coord(i, 1), coord(i, 2), coord(i, 3)
            end do
            write(10,*)
        end if
    end do
    close(10)

end subroutine molecular_dynamics
