integer function read_Natoms(input_file) result(Natoms)
    implicit none
    character(len=*), intent(in) :: input_file
    integer :: io_stat
    integer :: unit_number

    unit_number = 10
    open(unit=unit_number, file=input_file, status="old", action="read", iostat=io_stat)
    if (io_stat /= 0) then
        print *, "Error: Unable to open the input file."
        stop
    end if

    read(unit_number, *) Natoms
    close(unit_number)
end function read_Natoms
