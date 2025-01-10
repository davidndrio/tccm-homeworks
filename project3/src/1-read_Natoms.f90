integer function read_Natoms(input_file) result(Natoms)
    implicit none
    character(len=*), intent(in) :: input_file
    integer :: io_stat

    open(unit=10, file=input_file, status="old", action="read", iostat=io_stat)
    if (io_stat /= 0) then
        print *, "Error: Cannot open input file ", trim(input_file)
        stop
    end if

    read(10, *, iostat=io_stat) Natoms
    if (io_stat /= 0) then
        print *, "Error: Failed to read the number of atoms."
        Natoms = -1
        stop
    end if

    close(10)
end function read_Natoms

