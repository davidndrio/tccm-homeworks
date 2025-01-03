program main
    use angular_momentum_module  ! Incluye el módulo para corrección de momento angular
    implicit none

    ! Declaraciones
    integer :: Natoms
    double precision, allocatable :: coord(:,:), mass(:), distance(:,:), velocity(:,:), acceleration(:,:)
    double precision :: epsilon, sigma, total_potential_energy, kinetic_energy, total_energy
    character(len=100) :: input_file
    logical :: i_stat  ! Cambiado a logical
    integer :: i, j
    integer :: read_Natoms
    integer :: alloc_status

    ! Interfaces para funciones
    interface
        double precision function V(epsilon, sigma, Natoms, distance)
            integer, intent(in) :: Natoms
            double precision, intent(in) :: epsilon, sigma, distance(Natoms,Natoms)
        end function V

        double precision function T(Natoms, velocity, mass)
            integer, intent(in) :: Natoms
            double precision, intent(in) :: velocity(Natoms,3)
            double precision, intent(in) :: mass(Natoms)
        end function T
    end interface

    ! Solicitar el archivo de entrada
    write(*,*) "Enter the input file name:"
    read(*,*) input_file

    ! Verificar si el archivo existe
    inquire(file=input_file, exist=i_stat)
    if (.not. i_stat) then
        print *, "Error: Input file does not exist."
        stop
    end if

    ! Leer número de átomos
    Natoms = read_Natoms(input_file)
    if (Natoms <= 0) then
        print *, "Error: Number of atoms must be positive."
        stop
    end if
    write(*,*) "Number of atoms:", Natoms

    ! Asignar memoria
    allocate(coord(Natoms, 3), mass(Natoms), distance(Natoms, Natoms), &
             velocity(Natoms, 3), acceleration(Natoms, 3), stat=alloc_status)
    if (alloc_status /= 0) then
        print *, "Error: Memory allocation failed."
        stop
    end if

    ! Leer coordenadas y masas
    call read_molecule(input_file, Natoms, coord, mass)

    ! Calcular distancias internucleares
    call compute_distances(Natoms, coord, distance)

    ! Validar distancias demasiado pequeñas
    do i = 1, Natoms - 1
        do j = i + 1, Natoms
            if (distance(i, j) < 1.0d-6) then
                print *, "Error: Atoms too close. Distance between atoms", i, "and", j, "is", distance(i, j)
                stop
            end if
        end do
    end do

    ! Calcular energía potencial
    epsilon = 0.0661d0
    sigma = 0.3345d0
    total_potential_energy = V(epsilon, sigma, Natoms, distance)

    ! Inicializar velocidades
    velocity = 0.0d0

    ! Calcular energía cinética
    kinetic_energy = T(Natoms, velocity, mass)

    ! Calcular energía total
    total_energy = kinetic_energy + total_potential_energy

    ! Calcular aceleraciones
    call compute_acc(Natoms, coord, mass, distance, acceleration, epsilon, sigma)

    ! Corregir momento angular
    call correct_angular_momentum(Natoms, coord, velocity, mass)

    ! Imprimir energías iniciales
    write(*,*) "Initial energies:"
    write(*,*) "KE=", kinetic_energy, ", PE=", total_potential_energy, ", TE=", total_energy

    ! Ejecutar simulación
    call molecular_dynamics(Natoms, coord, velocity, acceleration, mass, epsilon, sigma, distance)

    ! Finalizar
    write(*,*) "Simulation completed. Output: trajectories.xyz"
    deallocate(coord, mass, distance, velocity, acceleration)
end program main

