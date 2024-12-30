program main
    implicit none

    ! Declarations
    integer :: Natoms
    double precision, allocatable :: coord(:,:), mass(:), distance(:,:)
    double precision :: epsilon, sigma, total_potential_energy
    character(len=100) :: input_file
    integer :: i_stat, i
    integer :: read_Natoms

    ! Request the input file
    write(*,*) "Enter the input file name:"
    read(*,*) input_file

    ! Call the function read_Natoms
    Natoms = read_Natoms(input_file)
    write(*,*) "Number of atoms:", Natoms

    ! Allocate memory
    allocate(coord(Natoms, 3), mass(Natoms), distance(Natoms, Natoms), stat=i_stat)
    if (i_stat /= 0) then
        print *, "Error: Memory allocation failed."
        stop
    end if

    ! Read coordinates and masses
    call read_molecule(input_file, Natoms, coord, mass)

    ! Compute internuclear distances
    call compute_distances(Natoms, coord, distance)

    ! Print distances
    write(*,*) "Internuclear distances:"
    do i = 1, Natoms
        write(*,*) distance(i, :)
    end do
    
    ! Free allocated memory
    deallocate(coord, mass, distance)
end program main
