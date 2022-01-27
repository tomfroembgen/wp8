program main
  use build_box, only: build_grid
  use print_matrix, only: write_matrix
  use algorithm, only: calc_forces, calc_e_pot_harm, calc_e_kin
  implicit none

  !> floating point accuracy
  integer, parameter :: wp = selected_real_kind(15)

  !> Number of atoms
  integer, parameter :: natom = 108

  !> length of box in Angstrom
  real(wp), parameter :: l = 17.158_wp

  !> atom coordinates of the system
  real(wp), dimension(3, natom) :: xyz

  !> forces of the system
  real(wp), dimension(3, natom) ::  f

  !> komischer parameter
  real(wp), dimension(3, natom) ::  a

  !> velocities
  real(wp), dimension(3, natom) :: v

  !> force constant
  real(wp), parameter :: k = 1.0_wp

  !> number of iterations
  integer, parameter :: itime = 200

  !> timestep
  real(wp), parameter :: delta = 5.E-3_wp

  !> mass of argon atoms
  real(wp), parameter :: mass = 39.948_wp

  !> harmonic potential energy
  real(wp):: e_pot_harm

  !> kinetic energy
  real(wp) :: e_kin

  !> total energy
  real(wp) :: e_tot

  !> temperatute to rescale velocities
  integer, parameter :: t_rescale = 8

  !> loop indices
  integer :: i, j, b

  a = 0.0_wp
  v = 0.0_wp

  !> call grid building subroutine, 1 for sc and 4 for fcc
  call build_grid(xyz, natom, l, 4)

  !> print the resulting coordinate matrix
  !call write_matrix(xyz, name="coordinate matrix")

  !> compute the forces
  call calc_forces(natom, xyz, f, a, k, l)

  !> print the resulting force matrix
  !call write_matrix(f, name="force matrix")

  !> compute potential energy
  call calc_e_pot_harm(natom, xyz, l, a, k, e_pot_harm)
  !write (*, *) "Harmonic potential energy =", e_pot_harm

  !> open file for energies and start printing
  open (unit=15, file="energy.csv")
  write (15, "(A)") "Step,Time,E_kin,E_pot,E_tot"

  !> open files for coordinates and start printing
  open (unit=16, file="traj.xyz")
  write (16, "(I4)") natom
  write (16, "(A,I5,A,F15.8)") "Step = ", i, ", Time = ", delta*i
  do b = 1, natom
    write (16, "(A,3F15.8)") "Ar", xyz(1, b), xyz(2, b), xyz(3, b)
  end do

  !> main loop
  main_loop: do i = 1, itime

    !> compute the new velocities (part 1) and positions
    do j = 1, natom
      v(:, j) = v(:, j) + 0.5_wp*(f(:, j)/mass)*delta
      xyz(:, j) = xyz(:, j) + (v(:, j)*delta)
    end do

    !> calculate new forces
    call calc_forces(natom, xyz, f, a, k, l)
    !> print the resulting force matrix
    !call write_matrix(f, name="force matrix")

    !> compute new velocities (part 2) and kinetic energy
    do j = 1, natom
      v(:, j) = v(:, j) + 0.5_wp*(f(:, j)/mass)*delta
    end do

    !> calculate energies
    call calc_e_kin(natom, mass, v, e_kin)
    call calc_e_pot_harm(natom, xyz, l, a, k, e_pot_harm)
    e_tot = e_kin + e_pot_harm

    !> print results to terminal
    !write (*, *) "Kinetic energy =           ", e_kin
    !write (*, *) "Harmonic potential energy =", e_pot_harm
    !write (*, *) "Total energy =             ", e_tot

    !> print results to energy file
    write (15, "(I3,A1,F15.8,A1,F15.8,A1,F15.8,A1,F15.8)") i, ",", delta*i, ",", e_kin, ",", e_pot_harm, ",", e_tot

    !> print results to trajectory file
    write (16, "(I4)") natom
    write (16, "(A,I5,A,F15.8)") "Step = ", i, ", Time = ", delta*i
    do b = 1, natom
      write (16, "(A,3F15.8)") "Ar", xyz(1, b), xyz(2, b), xyz(3, b)
    end do

  end do main_loop

end program main
