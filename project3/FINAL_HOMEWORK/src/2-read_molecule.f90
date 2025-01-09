subroutine read_molecule(input_file, Natoms, coord, mass)
    character(len=*), intent(in) :: input_file
    integer, intent(in) :: Natoms
    double precision, intent(out) :: coord(Natoms, 3), mass(Natoms)
    integer :: i, io_stat

    open(unit=10, file=input_file, status="old", action="read", iostat=io_stat)
    if (io_stat /= 0) then
        print *, "Error opening input file."
        stop
    end if

    ! Reads the number of atoms (already read in read_Natoms, so it skips that line)
    read(10, *, iostat=io_stat) i
    if (io_stat /= 0) then
        print *, "Error reading number of atoms."
        stop
    end if

    ! Reads molecule's data
    do i = 1, Natoms
        read(10, *, iostat=io_stat) coord(i, 1), coord(i, 2), coord(i, 3), mass(i)
        if (io_stat /= 0) then
            print *, "Error reading molecular data at atom index:", i
            stop
        end if
    end do

    close(10)
end subroutine read_molecule


