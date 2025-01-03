module angular_momentum_module
    implicit none

    contains

    subroutine correct_angular_momentum(Natoms, coord, velocity, mass)
        ! Subroutine to correct the angular momentum of the system
        ! Input/Output variables
        integer, intent(in) :: Natoms
        double precision, intent(inout) :: coord(Natoms, 3), velocity(Natoms, 3)
        double precision, intent(in) :: mass(Natoms)

        ! Local variables
        double precision :: angular_momentum(3), total_mass
        double precision :: correction(3), cross_product(3)
        integer :: i

        ! Initialize angular momentum to zero
        angular_momentum = 0.0d0

        ! Compute total angular momentum
        do i = 1, Natoms
            cross_product = cross_product_vector(coord(i, :), mass(i) * velocity(i, :))
            angular_momentum = angular_momentum + cross_product
        end do

        ! Compute the correction factor
        total_mass = sum(mass)
        correction = angular_momentum / total_mass

        ! Apply the correction to each particle
        do i = 1, Natoms
            velocity(i, :) = velocity(i, :) - cross_product_vector(coord(i, :), correction)
        end do
    end subroutine correct_angular_momentum

    ! Function to calculate the cross product of two 3D vectors
    function cross_product_vector(vec1, vec2) result(cross_product)
        double precision, intent(in) :: vec1(3), vec2(3)
        double precision :: cross_product(3)

        cross_product(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
        cross_product(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
        cross_product(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)
    end function cross_product_vector

end module angular_momentum_module

