subroutine compute_acc(Natoms, coord, mass, distance, acceleration, epsilon, sigma)
  implicit none
  integer, intent(in) :: Natoms
  double precision, intent(in) :: coord(Natoms,3), mass(Natoms), distance(Natoms,Natoms)
  double precision, intent(in) :: epsilon, sigma
  double precision, intent(out) :: acceleration(Natoms,3)
  double precision :: U, rij, dx, dy, dz, inv_mass
  integer :: i, j

  ! Initialize acceleration array to zero
  acceleration = 0.0d0

  do i = 1, Natoms
    inv_mass = 1.0d0 / mass(i)
    do j = 1, Natoms
      if (i /= j) then
        dx = coord(i,1) - coord(j,1)
        dy = coord(i,2) - coord(j,2)
        dz = coord(i,3) - coord(j,3)
        rij = distance(i,j)
        U = 24.0d0 * epsilon * ( (sigma/rij)**6 - 2.0d0*(sigma/rij)**12 ) / rij

        acceleration(i,1) = acceleration(i,1) - U * dx * inv_mass
        acceleration(i,2) = acceleration(i,2) - U * dy * inv_mass
        acceleration(i,3) = acceleration(i,3) - U * dz * inv_mass
      end if
    end do
  end do
end subroutine compute_acc
