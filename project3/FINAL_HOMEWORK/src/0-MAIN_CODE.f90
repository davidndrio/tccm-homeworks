program main
    implicit none

    integer :: Natoms
    double precision, allocatable :: coord(:,:), mass(:), distance(:,:), velocity(:,:), acceleration(:,:)
    character(len=2), allocatable :: atom_symbols(:)
    character(len=100) :: input_file
    integer :: allocate_status

    interface
        function read_Natoms(input_file) result(Natoms)
            integer :: Natoms
            character(len=*) :: input_file
        end function read_Natoms

        subroutine read_molecule(input_file, Natoms, coord, mass)
            character(len=*), intent(in) :: input_file
            integer, intent(in) :: Natoms
            double precision, intent(out) :: coord(Natoms, 3), mass(Natoms)
        end subroutine read_molecule

        subroutine molecular_dynamics(Natoms, coord, velocity, acceleration, mass, epsilon, sigma, distance)
            integer, intent(in) :: Natoms
            double precision, intent(inout) :: coord(Natoms, 3), velocity(Natoms, 3), acceleration(Natoms, 3)
            double precision, intent(inout) :: distance(Natoms, Natoms)
            double precision, intent(in) :: mass(Natoms), epsilon, sigma
        end subroutine molecular_dynamics
    end interface

    write(*,*) "Enter the input file name:"
    read(*,*) input_file  ! Input file name

    ! Reads the number of atoms
    Natoms = read_Natoms(input_file)
    if (Natoms <= 0) then
        print *, "Error: Invalid number of atoms."
        stop
    end if

    ! To allocate the arrays
    allocate(coord(Natoms, 3), &
             mass(Natoms), & 
             distance(Natoms, Natoms), & 
             velocity(Natoms, 3), & 
             acceleration(Natoms, 3), & 
             atom_symbols(Natoms), & 
             stat=allocate_status)

    if (allocate_status /= 0) then
        print *, "Error: Memory allocation failed."
        stop
    end if


    ! Reads molecular data
    call read_molecule(input_file, Natoms, coord, mass)

    ! Executes molecular dynamics
    call molecular_dynamics(Natoms, coord, velocity, acceleration, mass, 0.0661d0, 0.3345d0, distance)

    ! Frees the memory
    deallocate(coord, mass, distance, velocity, acceleration, atom_symbols)
    
    write(*,*) "MD simulation succesfully finished!"
end program main

