program main
    implicit none

    ! Declarations
    integer :: Natoms
    double precision, allocatable :: coord(:,:), mass(:), distance(:,:), velocity(:,:), acceleration(:,:)
    double precision :: epsilon, sigma, total_potential_energy, kinetic_energy, total_energy
    character(len=100) :: input_file
    integer :: i_stat, i
    integer :: read_Natoms

    ! Interface for V function
    interface
        double precision function V(epsilon, sigma, Natoms, distance)
            integer, intent(in) :: Natoms
            double precision, intent(in) :: epsilon, sigma, distance(Natoms,Natoms)
        end function V
    end interface

       ! Interface for T function
    interface
        double precision function T(Natoms, velocity, mass)
            integer, intent(in) :: Natoms
            double precision, intent(in) :: velocity(Natoms,3)
            double precision, intent(in) :: mass(Natoms)
        end function T
    end interface

    ! Request the input file
    write(*,*) "Enter the input file name:"
    read(*,*) input_file

    ! Call the function read_Natoms
    Natoms = read_Natoms(input_file)
    write(*,*) "Number of atoms:", Natoms
    
    allocate(coord(Natoms, 3), mass(Natoms), distance(Natoms, Natoms), velocity(Natoms, 3), acceleration(Natoms, 3), stat=i_stat)
        if (i_stat /= 0) then
                print *, "Error: Memory allocation failed."
                stop
        end if

    ! Read coordinates and masses (assuming this function is defined elsewhere)
    call read_molecule(input_file, Natoms, coord, mass)

    ! Compute internuclear distances (assuming this function is defined elsewhere)
    call compute_distances(Natoms, coord, distance)

    ! Compute potential energy
    epsilon = 0.0661  ! J/mol
    sigma = 0.3345    ! nm
    total_potential_energy = V(epsilon, sigma, Natoms, distance)


    velocity = 0.0d0  ! All velocities are set to zero for this example

       ! Compute kinetic energy
    kinetic_energy = T(Natoms, velocity, mass)

    ! Compute total energy (sum of kinetic and potential)
    total_energy = kinetic_energy + total_potential_energy

    ! Compute acceleration for each atom
    call compute_acc(Natoms, coord, mass, distance, acceleration, epsilon, sigma)

    ! Print energies
    write(*,*) "The Total Potential Energy of the Ensemble is:", total_potential_energy, "J"
    write(*,*) "The Kinetic Energy of the Ensemble is:", kinetic_energy, "J"
    write(*,*) "The Total Energy of the Ensemble is:", total_energy, "J"
    
    ! Example: Print acceleration for the first atom as a check
    write(*,*) "Acceleration for atom 1:", acceleration(1,:)

    ! Free allocated memory
    deallocate(coord, mass, distance, acceleration)
end program main
