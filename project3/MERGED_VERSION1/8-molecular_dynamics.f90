subroutine molecular_dynamics(Natoms, coord, velocity, acceleration, mass, epsilon, sigma, distance)
    implicit none

    ! Input variables
    integer, intent(in) :: Natoms
    double precision, intent(inout) :: coord(Natoms, 3), velocity(Natoms, 3), acceleration(Natoms, 3)
    double precision, intent(in) ::  distance(Natoms, Natoms)
    double precision, intent(in) :: mass(Natoms)
    double precision, intent(in) :: epsilon, sigma

    ! Parameters for the simulation
    double precision, parameter :: dt = 0.2d0
    integer, parameter :: total_steps = 1000
    integer, parameter :: M = 10

    ! Local variables
    integer :: step, i
    double precision :: kinetic_energy, potential_energy, total_energy
    double precision :: new_acceleration(Natoms, 3), temp_velocity(Natoms, 3)
    character(len=100) :: output_file
    integer :: io_status

    ! Open the output file
    output_file = "trajectory.xyz"
    open(unit=10, file=output_file, status='replace', action='write', iostat=io_status)
    if (io_status /= 0) then ! IS THIS NEEDED HERE?
        print *, "Error: Could not open the output file."
        stop
    end if

    ! Molecular dynamics
    do step = 1, total_steps

        ! 1. Update coordinates (r(n+1) = r(n) + v(n)*dt + 0.5*a(n)*dt^2)
        coord = coord + velocity * dt + 0.5d0 * acceleration * dt**2

        ! 2. Compute velocities part dependent on a(n) (v(n+**) = v(n) + 0.5*a(n)*dt)
        temp_velocity = velocity + 0.5d0 * acceleration * dt

        ! 3. Compute new acceleration (a(n+1))
        call compute_acc(Natoms, coord, mass, distance, acceleration, epsilon, sigma) 

        ! 4. Compute final velocities (v(n+1) = v(n+**) + 0.5*a(n+1)*dt)
        velocity = temp_velocity + 0.5d0 * new_acceleration * dt

        ! Update acceleration for the next step
        acceleration = new_acceleration

        ! Compute energies
        call E(Natoms, velocity, mass, distance, epsilon, sigma)
        ! Write trajectory to file every M steps
        if (mod(step, M) == 0) then ! Ensures step is a multiple of M
            write(10, *) Natoms
            write(10, *) "Step ", step, ": KE=", kinetic_energy, ", PE=", potential_energy, ", TE=", total_energy 
            do i = 1, Natoms
                write(10, *) "Ar", coord(i, 1), coord(i, 2), coord(i, 3) ! WE ARE USING ONLY ARGON RIGHT?
            end do
        end if
    end do
    close(10)

end subroutine molecular_dynamics

