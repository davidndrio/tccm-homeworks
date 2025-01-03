subroutine molecular_dynamics(Natoms, coord, velocity, acceleration, mass, epsilon, sigma, distance)
    use angular_momentum_module
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    ! Interfaces para funciones externas
    interface
        double precision function T(Natoms, velocity, mass)
            integer, intent(in) :: Natoms
            double precision, intent(in) :: velocity(Natoms, 3)
            double precision, intent(in) :: mass(Natoms)
        end function T
    end interface

    interface
        double precision function V(epsilon, sigma, Natoms, distance)
            integer, intent(in) :: Natoms
            double precision, intent(in) :: epsilon, sigma, distance(Natoms, Natoms)
        end function V
    end interface

    ! Variables de entrada
    integer, intent(in) :: Natoms
    double precision, intent(inout) :: coord(Natoms, 3), velocity(Natoms, 3), acceleration(Natoms, 3)
    double precision, intent(inout) :: distance(Natoms, Natoms)
    double precision, intent(in) :: mass(Natoms)
    double precision, intent(in) :: epsilon, sigma

    ! Parámetros de simulación
    double precision, parameter :: dt = 0.00005d0
    integer, parameter :: total_steps = 1000
    integer, parameter :: M = 1

    ! Variables locales
    integer :: step, i
    double precision :: kinetic_energy, potential_energy, total_energy
    double precision :: prev_total_energy, energy_diff, scale_factor
    double precision :: temp_velocity(Natoms, 3)
    character(len=100) :: output_file, energy_file
    integer :: io_status

    ! Archivos de salida
    output_file = "trajectories.xyz"
    energy_file = "energies.txt"

    open(unit=10, file=output_file, status='replace', action='write', iostat=io_status)
    if (io_status /= 0) then
        print *, "Error: Could not open the trajectory output file."
        stop
    end if

    open(unit=11, file=energy_file, status='replace', action='write', iostat=io_status)
    if (io_status /= 0) then
        print *, "Error: Could not open the energy output file."
        stop
    end if

    ! Bucle principal
    prev_total_energy = 0.0d0  ! Inicializamos la energía anterior a 0.0
    do step = 1, total_steps
        coord = coord + velocity * dt + 0.5d0 * acceleration * dt**2
        temp_velocity = velocity + 0.5d0 * acceleration * dt

        call compute_distances(Natoms, coord, distance)
        call compute_acc(Natoms, coord, mass, distance, acceleration, epsilon, sigma, velocity)

        velocity = temp_velocity + 0.5d0 * acceleration * dt
        call correct_angular_momentum(Natoms, coord, velocity, mass)

        kinetic_energy = T(Natoms, velocity, mass)
        potential_energy = V(epsilon, sigma, Natoms, distance)
        total_energy = kinetic_energy + potential_energy

        ! Energy Conservation: Scale velocities if energy drifts
        if (prev_total_energy /= 0.0d0) then
            energy_diff = total_energy - prev_total_energy
            if (abs(energy_diff) > 1.0d-8) then
                ! Check if total energy is too small
                if (abs(total_energy) > 1.0d-10) then
                    scale_factor = prev_total_energy / total_energy
                    if (.not. ieee_is_finite(scale_factor) .or. scale_factor < 1.0d-6 .or. scale_factor > 1.0d4) then
                        print *, "Warning: Excessive energy correction. Skipping at step", step
                    else
                        velocity = velocity * sqrt(scale_factor)
                        print *, "Energy correction applied at step", step
                    end if
                else
                    print *, "Warning: Total energy too small to apply correction at step", step
                end if
            end if
        end if

        prev_total_energy = total_energy  ! Guardamos la energía total para la próxima comparación

        ! Verificar si alguna energía es NaN o infinita
        if (.not. ieee_is_finite(kinetic_energy) .or. .not. ieee_is_finite(potential_energy)) then
            print *, "Error: Non-finite energy detected at step", step
            print *, "Kinetic energy:", kinetic_energy, "Potential energy:", potential_energy
            print *, "Coordinates:", coord
            stop
        end if

        ! Escribir datos cada M pasos
        if (mod(step, M) == 0) then
            write(10, '(A, I5, A, E15.8, A, E15.8, A, E15.8)') "Step ", step, &
                ": KE=", kinetic_energy, ", PE=", potential_energy, ", TE=", total_energy
            do i = 1, Natoms
                write(10, '(A, 3(1x, ES16.8))') "Ar", coord(i, 1), coord(i, 2), coord(i, 3)
            end do
            write(11, '(I6, 3(ES16.8))') step, kinetic_energy, potential_energy, total_energy
        end if

        ! Impresiones para depuración
        if (mod(step, 100) == 0) then
            print *, "Step:", step
            do i = 1, Natoms
                print *, "Atom", i, "Velocity:", velocity(i, :), "Acceleration:", acceleration(i, :)
            end do
        end if
    end do

    close(10)
    close(11)

end subroutine molecular_dynamics

