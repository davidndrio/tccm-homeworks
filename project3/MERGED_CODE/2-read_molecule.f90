subroutine read_molecule(input_file, Natoms, coord, mass)
    implicit none
    character(len=*), intent(in) :: input_file
    integer, intent(in) :: Natoms
    double precision, intent(out) :: coord(Natoms, 3)
    double precision, intent(out) :: mass(Natoms)
    integer :: i, io_stat
    integer :: unit_number

    unit_number = 20
    open(unit=unit_number, file=input_file, status="old", action="read", iostat=io_stat)
    if (io_stat /= 0) then
        print *, "Error: Unable to open the input file in read_molecule."
        stop
    end if

    read(unit_number, *)
    do i = 1, Natoms
        read(unit_number, *) coord(i, 1), coord(i, 2), coord(i, 3), mass(i)
    end do
    close(unit_number)
end subroutine read_molecule
