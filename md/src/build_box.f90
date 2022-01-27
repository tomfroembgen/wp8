module build_box
  implicit none
  private

  !> floating point accuracy
  integer, parameter :: wp = selected_real_kind(15)

  !> export subroutines
  public :: build_grid

contains

  subroutine build_grid(xyz, natom, l, formula_units)

    implicit none

    !> number of formula units
    integer, intent(in) :: formula_units

    !> Number of atoms
    integer, intent(in) :: natom

    !> length of box in Angstrom
    real(wp), intent(in) :: l

    !> atom coordinates of the system
    real(wp), dimension(3, natom), intent(out) :: xyz

    !> number of grid points, lattice points along edge
    integer:: n, nl

    !> loop indices
    integer :: i, j, k

    !> half of the box length
    real(wp) :: hl

    !> distance between points along edge
    real(wp) :: dl

    !> coordinates of the atoms
    real(wp) :: x, y, z

    !> number of lattice points along edge of cube
    nl = int((natom/formula_units)**(1.0_wp/3.0_wp))

    !> adjust number according to number of atoms
    if ((nl*formula_units)**3 .lt. natom) then
      nl = nl + 1
    end if

    !> distance between atoms along edge of cube
    dl = l/nl
    hl = l/2.0_wp

    !> place atoms in box
    n = 1
    open (unit=14, file="box.xyz")
    write (14, "(I4)") natom
    write (14, *) ""

    loop1: do i = 0, nl - 1
      do j = 0, nl - 1
        do k = 0, nl - 1

          !> -hl results in centering of the box at 0,0,0
          x = i*dl - hl
          y = j*dl - hl
          z = k*dl - hl
          xyz(1, n) = x
          xyz(2, n) = y
          xyz(3, n) = z
          write (14, *) 'Ar', x, y, z! prints out “fort.14” file with coordinates

          if (formula_units == 4) then
            x = i*dl - hl
            y = (j + 0.5)*dl - hl
            z = (k + 0.5)*dl - hl
            xyz(1, n + 1) = x
            xyz(2, n + 1) = y
            xyz(3, n + 1) = z
            write (14, *) 'Ar', x, y, z

            x = (i + 0.5)*dl - hl
            y = j*dl - hl
            z = (k + 0.5)*dl - hl
            xyz(1, n + 2) = x
            xyz(2, n + 2) = y
            xyz(3, n + 2) = z
            write (14, *) 'Ar', x, y, z

            x = (i + 0.5)*dl - hl
            y = (j + 0.5)*dl - hl
            z = k*dl - hl
            xyz(1, n + 3) = x
            xyz(2, n + 3) = y
            xyz(3, n + 3) = z
            write (14, *) 'Ar', x, y, z
          end if

          n = n + formula_units
          if (n == (natom + 1)) then
            exit loop1
          end if

        end do
      end do
    end do loop1

  end subroutine build_grid

end module build_box

