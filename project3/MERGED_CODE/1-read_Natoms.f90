integer function read_Natoms(input_file) result(Natoms)
    implicit none
    character(len=*), intent(in) :: input_file  ! Input file name
    integer :: io_stat  ! Status of file operations

    ! Open the input file in read-only mode
    open(unit=10, file=input_file, status="old", action="read", iostat=io_stat)
    if (io_stat /= 0) then
        ! Stop the program if the file cannot be opened
        stop
    end if

    ! Read the number of atoms from the file
    read(10, *, iostat=io_stat) Natoms
    if (io_stat /= 0) then
        ! Print an error message if the read operation fails
        print *, "Error: Unable to read the number of atoms."
        Natoms = -1  ! Assign an invalid value to indicate failure
        stop
    end if

    ! Close the file
    close(10)
end function read_Natoms

