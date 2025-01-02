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

