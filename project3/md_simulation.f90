program molecular_dynamics
    implicit none

    ! Define double precision
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: num_atoms = 2  ! Adjust as needed
    real(dp), parameter :: epsilon = 0.0661_dp, sigma = 0.3345_dp
    real(dp), parameter :: dt = 0.01_dp  ! Time step
    integer :: i, j, step, num_steps, i_stat

    ! Variables
    real(dp), allocatable :: positions(:,:), velocities(:,:), forces(:,:)
    real(dp) :: r, v_lj

    ! Initialize matrices
    allocate(positions(3, num_atoms), velocities(3, num_atoms), forces(3, num_atoms), stat=i_stat)
    if (i_stat /= 0) then
        print *, "Error allocating memory for matrices."
        stop
    end if

    ! Read initial data from the input file
    call read_input(positions, velocities, num_atoms)

    ! Molecular dynamics simulation
    num_steps = 100  ! Number of simulation steps
    do step = 1, num_steps
        ! Initialize forces to zero
        forces = 0.0_dp

        ! Calculate forces between atom pairs
        do i = 1, num_atoms - 1
            do j = i + 1, num_atoms
                r = distance(positions(:, i), positions(:, j))
                v_lj = lennard_jones_potential(r, epsilon, sigma)
                call calculate_force(forces(:, i), forces(:, j), positions(:, i), positions(:, j), r, epsilon, sigma)
            end do
        end do

        ! Update positions and velocities
        do i = 1, num_atoms
            positions(:, i) = positions(:, i) + velocities(:, i) * dt
            velocities(:, i) = velocities(:, i) + forces(:, i) * dt
        end do

        ! Save updated positions
        call save_output(positions, step)
    end do

    ! Deallocate memory
    deallocate(positions, velocities, forces)

    print *, "Simulation completed successfully."
end program molecular_dynamics

!----------------------------------------
! Subroutines and auxiliary functions
!----------------------------------------

subroutine read_input(positions, velocities, num_atoms)
    implicit none
    integer, intent(in) :: num_atoms
    real(dp), intent(out) :: positions(3, num_atoms), velocities(3, num_atoms)
    real(dp) :: mass  ! Temporary variable to store mass
    integer :: i

    open(unit=10, file="file.inp", status="old")
    do i = 1, num_atoms
        read(10, *) positions(1, i), positions(2, i), positions(3, i), mass
    end do
    close(10)

    ! Initialize velocities to zero
    velocities = 0.0_dp
end subroutine read_input

!----------------------------------------

function lennard_jones_potential(r, epsilon, sigma) result(v_lj)
    implicit none
    real(dp), intent(in) :: r, epsilon, sigma
    real(dp) :: v_lj

    v_lj = 4.0_dp * epsilon * ((sigma / r)**12 - (sigma / r)**6)
end function lennard_jones_potential

!----------------------------------------

subroutine calculate_force(force1, force2, pos1, pos2, r, epsilon, sigma)
    implicit none
    real(dp), intent(inout) :: force1(3), force2(3)
    real(dp), intent(in) :: pos1(3), pos2(3), r, epsilon, sigma
    real(dp) :: f, dr(3)

    f = 24.0_dp * epsilon * (2.0_dp * (sigma / r)**12 - (sigma / r)**6) / r**2
    dr = pos2 - pos1

    force1 = force1 + f * dr / r
    force2 = force2 - f * dr / r
end subroutine calculate_force

!----------------------------------------

function distance(pos1, pos2) result(r)
    implicit none
    real(dp), intent(in) :: pos1(3), pos2(3)
    real(dp) :: r

    r = sqrt(sum((pos2 - pos1)**2))
end function distance

!----------------------------------------

subroutine save_output(positions, step)
    implicit none
    integer, intent(in) :: step
    real(dp), intent(in) :: positions(:,:)
    integer :: i

    open(unit=20, file="output.txt", status="unknown", position="append")
    write(20, '(a, i0)') "Step: ", step
    do i = 1, size(positions, 2)
        write(20, '(3f10.5)') positions(:, i)
    end do
    close(20)
end subroutine save_output
