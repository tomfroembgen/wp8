program main
  use build_box, only: build_grid
  use print_matrix, only: write_matrix
  use algorithm

  implicit none
  intrinsic :: sqrt

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

  !> lennard-jones potential energy
  real(wp) :: e_pot_lj

  !> temperature
  real(wp) :: T

  !> temporary variable
  real(wp) :: temp

  !> temperatute to rescale velocities
  integer, parameter :: T_req = 8

  !> loop indices
  integer :: i, j, b

  !> initialze variables
  a = 0.0_wp
  v = 0.0_wp
  temp = 1._wp

  !> call grid building subroutine, 1 for sc and 4 for fcc
  call build_grid(xyz, natom, l, 4)

  !> print the resulting coordinate matrix
  !call write_matrix(xyz, name="coordinate matrix")

  !> compute the forces and potential energy
  call calc_forces(natom, xyz, f, l, e_pot_lj)

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

    if (i > 1) then
      !> temperature
      T = e_kin/(3._wp*natom)

      !> scaling factor for temperatur
      temp = sqrt(T_req/T)
    end if

    !> compute the new velocities (part 1) and positions
    do j = 1, natom
      v(:, j) = (v(:, j) + 0.5_wp*(f(:, j)/mass)*delta)*temp
      xyz(:, j) = xyz(:, j) + (v(:, j)*delta)
    end do

    !> calculate new forces and potential energy
    call calc_forces(natom, xyz, f, l, e_pot_lj)

    !> compute new velocities (part 2)
    do j = 1, natom
      v(:, j) = v(:, j) + 0.5_wp*(f(:, j)/mass)*delta
    end do

    !> calculate kinetic energy
    call calc_e_kin(natom, mass, v, e_kin)

    e_tot = e_kin + e_pot_lj

    !> print results to terminal
    !write (*, *) "Kinetic energy =           ", e_kin
    !write (*, *) "Harmonic potential energy =", e_pot_lj
    !write (*, *) "Total energy =             ", e_tot

    !> print results to energy file
    write (15, "(I3,A1,F15.8,A1,F15.8,A1,F15.8,A1,F15.8)") i, ",", delta*i, ",", e_kin, ",", e_pot_lj, ",", e_tot

    !> print results to trajectory file
    write (16, "(I4)") natom
    write (16, "(A,I5,A,F15.8)") "Step = ", i, ", Time = ", delta*i
    do b = 1, natom
      write (16, "(A,3F15.8)") "Ar", xyz(1, b), xyz(2, b), xyz(3, b)
    end do

  end do main_loop

end program main
