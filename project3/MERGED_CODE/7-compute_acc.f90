subroutine compute_acc(Natoms, coord, mass, distance, acceleration, epsilon, sigma, velocity)
    implicit none
    integer, intent(in) :: Natoms
    double precision, intent(in) :: mass(Natoms), epsilon, sigma
    double precision, intent(inout) :: coord(Natoms, 3), distance(Natoms, Natoms), velocity(Natoms, 3)
    double precision, intent(out) :: acceleration(Natoms, 3)

    ! Variables locales
    integer :: i, j
    double precision :: r, r2, r6, r12, f, dx, dy, dz

    ! Inicializar aceleraciones a cero
    acceleration = 0.0d0

    ! Calcular las aceleraciones
    do i = 1, Natoms
        do j = i + 1, Natoms
            dx = coord(i, 1) - coord(j, 1)
            dy = coord(i, 2) - coord(j, 2)
            dz = coord(i, 3) - coord(j, 3)
            r2 = dx**2 + dy**2 + dz**2

            ! Validar distancia mínima para evitar divisiones por cero
            if (r2 < 1.0d-12) then
                print *, "Warning: Atoms too close. Skipping interaction between atoms", i, "and", j
                cycle
            end if

            r = sqrt(r2)
            r6 = r2**3
            r12 = r6**2

            ! Calcular la fuerza de Lennard-Jones
            f = 24.0d0 * epsilon * (2.0d0 * (sigma**12 / r12) - (sigma**6 / r6)) / r2

            ! Verificar si la fuerza es demasiado grande
            if (abs(f) > 1.0d5) then
                print *, "Warning: Large force detected between atoms", i, "and", j, "Force:", f
                cycle
            end if

            ! Actualizar aceleraciones
            acceleration(i, 1) = acceleration(i, 1) + f * dx / mass(i)
            acceleration(i, 2) = acceleration(i, 2) + f * dy / mass(i)
            acceleration(i, 3) = acceleration(i, 3) + f * dz / mass(i)
            acceleration(j, 1) = acceleration(j, 1) - f * dx / mass(j)
            acceleration(j, 2) = acceleration(j, 2) - f * dy / mass(j)
            acceleration(j, 3) = acceleration(j, 3) - f * dz / mass(j)
        end do
    end do

end subroutine compute_acc

