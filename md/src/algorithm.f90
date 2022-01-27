module algorithm

  use print_matrix, only: write_matrix

  implicit none
  private

  !> floating point accuracy
  integer, parameter :: wp = selected_real_kind(15)

  public :: calc_forces_harm, calc_e_pot_harm, calc_e_kin, calc_e_pot_lj, calc_forces

contains

  subroutine calc_forces_harm(natom, xyz, f, a, k, l)

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

  end subroutine calc_forces_harm

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
      e_pot_harm = e_pot_harm + 0.5_wp*k*sum((xyz(:, i) - a(:, i))**2)
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
      e_kin = e_kin + 0.5_wp*mass*sum((v(:, i))**2)
    end do

  end subroutine calc_e_kin

  subroutine calc_e_pot_lj(natom, xyz, l, e_pot_lj)

    implicit none

    !> Number of atoms
    integer, intent(in) :: natom

    !> atom coordinates of the system
    real(wp), dimension(3, natom), intent(in) ::  xyz

    !> length of box in Angstrom
    real(wp), intent(in) :: l

    !> lennard-jones potential energy
    real(wp), intent(out) :: e_pot_lj

    !> lj parameter
    real(wp), parameter :: sigma = 3.405_wp

    !> lj parameter
    real(wp), parameter :: epsilon = 120.0_wp

    !> cutoff energy
    real(wp), parameter :: e_cutoff = 3.83738839608178386E-3_wp

    !> cutoff radii, squared, sigma squared
    real(wp) :: rcutoff, rcutoff2, sigma2

    !> roots of sigma
    real(wp) :: sigmar2, sigmar6, sigmar12

    !> squared distance between particles i,j
    real(wp) :: xyzij2

    !> loop indices
    integer :: i, j

    sigma2 = sigma*sigma
    rcutoff = 0.5_wp*l
    rcutoff2 = rcutoff*rcutoff

    e_pot_lj = 0.0_wp

    do i = 1, natom - 1
      do j = i + 1, natom

        xyzij2 = sum((xyz(:, i) - xyz(:, j))**2)
        sigmar2 = sigma2/xyzij2
        sigmar6 = sigmar2*sigmar2*sigmar2
        sigmar12 = sigmar6*sigmar6

        if (xyzij2 .lt. rcutoff2) then
          e_pot_lj = e_pot_lj + sigmar12 - sigmar6 + e_cutoff
        end if

      end do
    end do

    e_pot_lj = e_pot_lj*epsilon*4.0_wp

  end subroutine calc_e_pot_lj

  subroutine calc_forces(natom, xyz, f, l)

    implicit none

    !> Number of atoms
    integer, intent(in) :: natom

    !> atom coordinates of the system
    real(wp), dimension(3, natom), intent(in) ::  xyz

    !> forces of the system
    real(wp), dimension(3, natom), intent(out) ::  f

    !> length of box in Angstrom
    real(wp), intent(in) :: l

    !> lj parameter
    real(wp), parameter :: sigma = 3.405_wp

    !> lj parameter
    real(wp), parameter :: epsilon = 120.0_wp

    !> cutoff radii, squared, sigma squared
    real(wp) :: rcutoff, rcutoff2, sigma2

    !> roots of sigma
    real(wp) :: sigmar2, sigmar6, sigmar12

    !> squared distance between particles i,j
    real(wp) :: xyzij2

    !> distance vector of i,j
    real(wp), dimension(3) :: xyzij

    !> number needed for force update
    real(wp) :: ff

    !> loop indices
    integer :: i, j

    sigma2 = sigma*sigma
    rcutoff = 0.5_wp*l
    rcutoff2 = rcutoff*rcutoff

    f = 0._wp
    do i = 1, natom - 1
      do j = i + 1, natom

        xyzij = xyz(:, i) - xyz(:, j)
        xyzij = xyzij - l * (anint(xyzij/l))
        xyzij2 = sum(xyzij**2)

        if (xyzij2 .lt. rcutoff2) then
          sigmar2 = sigma2/xyzij2
          sigmar6 = sigmar2*sigmar2*sigmar2
          sigmar12 = sigmar6*sigmar6

          ff = 48.0_wp*epsilon*(sigmar12 - (0.5_wp*sigmar6))/xyzij2
          f(:, i) = f(:, i) + ff*xyzij
          f(:, j) = f(:, j) + ff*xyzij
        end if

      end do
    end do
  end subroutine calc_forces

end module algorithm
