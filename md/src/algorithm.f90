module algorithm

  use print_matrix, only: write_matrix

  implicit none
  private

  !> floating point accuracy
  integer, parameter :: wp = selected_real_kind(15)

  public :: calc_forces, calc_e_pot_harm, calc_e_kin

contains

  subroutine calc_forces(natom, xyz, f, a, k, l)

    implicit none

    !> Number of atoms
    integer, intent(in) :: natom

    !> atom coordinates of the system
    real(wp), dimension(3, natom), intent(in) ::  xyz

    !> forces of the system
    real(wp), dimension(3, natom), intent(out) ::  f

    !> komischer parameter
    real(wp), dimension(3, natom), intent(in) ::  a

    !> length of box in Angstrom
    real(wp), intent(in) :: l

    !> force constant
    real(wp), intent(in) :: k

    !> loop indices
    integer :: i

    f = 0._wp
    do i = 1, natom
      f(:, i) = -k*(xyz(:, i) - a(:, i))
    end do

  end subroutine calc_forces

  subroutine calc_e_pot_harm(natom, xyz, l, a, k, e_pot_harm)

    implicit none

    !> Number of atoms
    integer, intent(in) :: natom

    !> atom coordinates of the system
    real(wp), dimension(3, natom), intent(in) ::  xyz

    !> length of box in Angstrom
    real(wp), intent(in) :: l

    !> komischer parameter
    real(wp), dimension(3, natom), intent(in) ::  a

    !> force constant
    real(wp), intent(in) :: k

    !> harmonic potential energy
    real(wp), intent(out) :: e_pot_harm

    !> loop indices
    integer :: i

    e_pot_harm = 0.0_wp
    do i = 1, natom
      e_pot_harm = e_pot_harm + 0.5_wp*k*sum((xyz(:,i) - a(:,i))**2)
    end do

  end subroutine calc_e_pot_harm

  subroutine calc_e_kin(natom, mass, v, e_kin)

    implicit none

    !> Number of atoms
    integer, intent(in) :: natom

    !> length of box in Angstrom
    real(wp), intent(in) :: mass

    !> komischer parameter
    real(wp), dimension(3, natom), intent(in) ::  v

    !> kinetic energy
    real(wp), intent(out) :: e_kin

    !> loop indices
    integer :: i

    e_kin = 0.0_wp
    do i = 1, natom
      e_kin = e_kin + 0.5_wp * mass * sum((v(:,i))**2)
    end do

  end subroutine calc_e_kin

end module algorithm
