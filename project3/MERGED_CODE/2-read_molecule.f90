subroutine read_molecule(input_file, Natoms, coord, mass)
    implicit none

    ! Input arguments
    character(len=*), intent(in) :: input_file  ! Input file name
    integer, intent(in) :: Natoms  ! Number of atoms

    ! Output arguments
    double precision, intent(out) :: coord(Natoms, 3)  ! Coordinates of atoms (x, y, z)
    double precision, intent(out) :: mass(Natoms)  ! Masses of atoms

    ! Local variables
    integer :: i, io_status  ! Loop counter and file I/O status

    ! Open the input file in read-only mode
    open(unit=20, file=input_file, status="old", action="read", iostat=io_status)
    if (io_status /= 0) then
        ! Stop execution if the file cannot be opened
        print *, "Error: Unable to open the input file."
        stop
    end if

    ! Skip the first line (contains the number of atoms)
    read(20, *, iostat=io_status)
    if (io_status /= 0) then
        ! Stop execution if unable to read the number of atoms
        print *, "Error: Unable to read the number of atoms."
        stop
    end if

    ! Read the coordinates (x, y, z) and masses of the atoms
    do i = 1, Natoms
        read(20, *, iostat=io_status) coord(i, 1), coord(i, 2), coord(i, 3), mass(i)
        if (io_status /= 0) then
            ! Stop execution if an error occurs while reading atom data
            print *, "Error reading data for atom ", i
            stop
        end if
    end do

    ! Close the file
    close(20)
end subroutine read_molecule

