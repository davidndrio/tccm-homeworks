program lj
implicit none

double precision, allocatable :: coord(:,:), mass(:), dist(:,:) 
integer :: Natom, i, j, k, i_stat 
character(len=100) :: input_file
double precision :: epsilon, sigma, sumlj, V

epsilon = 0.0661  ! j/mol 
sigma = 0.3345    ! nm

write(6, *) "Enter the input file name:"
read(5, *) input_file

open(10, file=input_file, status="old", action="read", iostat=i_stat)
if (i_stat /= 0) then
    print *, "Error: Unable to open the file."
    stop
end if

read(10, *) Natom
close(10)

allocate(coord(Natom, 3), mass(Natom), dist(Natom, Natom), stat=i_stat)
if (i_stat /= 0) then
    print *, "Memory allocation failed!"
    stop
end if

call read_molecule(input_file, Natom, coord, mass) ! Subroutine used to read the coordinates and masses

call compute_distances(Natom, coord, dist) ! Subroutine used to calculate internuclear distances

! Function used to calculate the LJ potential for each pair
write(6, *) "The Total Potential Energy of the Ensemble is: ", V(epsilon, sigma, Natom, dist), "J"

deallocate(coord, mass, dist, stat=i_stat)

stop
end program lj


subroutine read_molecule(input_file, Natoms, coord, mass)
    implicit none
    character(len=*), intent(in) :: input_file
    integer, intent(in) :: Natoms
    double precision, intent(out) :: coord(Natoms,3)
    double precision, intent(out) :: mass(Natoms)
    integer :: i, io_stat

    open(20, file=input_file, status="old", action="read", iostat=io_stat)
    if (io_stat /= 0) then
        print *, "Error: Unable to open the file in subroutine."
        stop
    end if

    read(20, *) ! To skip the 1st line (already read)

    do i = 1, Natoms
        read(20, *) coord(i,1), coord(i,2), coord(i,3), mass(i)
    end do

    close(20)
end subroutine read_molecule


subroutine compute_distances(Natoms, coord, dist)
  implicit none
  integer, intent(in) :: Natoms
  double precision, intent(in) :: coord(Natoms, 3)
  double precision, intent(out) :: dist(Natoms, Natoms)
  integer :: i, j

  do i = 1, Natoms
      do j = i+1, Natoms
          dist(i, j) = sqrt((coord(i, 1) - coord(j, 1))**2 + (coord(i, 2) - coord(j, 2))**2 + &
                            (coord(i, 3) - coord(j, 3))**2)
          dist(j, i) = dist(i, j) ! Symmetric matrix
      end do
  end do
end subroutine compute_distances


! Function to calculate Lennard-Jones potential
double precision function V(epsilon, sigma, Natoms, distance)
  implicit none
  integer, intent(in) :: Natoms
  double precision, intent(in) :: epsilon, sigma, distance(Natoms, Natoms)
  integer :: i, j

  V = 0.0

  do i = 1, Natoms-1
   do j = i+1, Natoms
    V = V + 4.0 * epsilon * ((sigma / distance(i, j))**12 - (sigma / distance(i, j))**6)
   end do
  end do
end function V
