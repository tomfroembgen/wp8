module build_box
  implicit none
  private

  !> floating point accuracy
  integer, parameter :: wp = selected_real_kind(15)

  !> export subroutines
  public :: build_grid

contains

  subroutine build_grid(formula_units)

    implicit none
    integer, parameter :: natom = 108         ! Number of atoms
    integer, parameter :: l = 18              ! length of box in Angstrom
    integer, intent(in) :: formula_units      ! numver of formula units
    integer:: n, nl, i, j, k
    real(wp), parameter :: hl = l/2
    real(wp) :: dl, x, y, z

    !> number of lattice points along edge of cube
    nl = int((natom/formula_units)**(1.0_wp/3.0_wp))

    !> adjust number according to number of atoms
    if ((nl*formula_units)**3 .lt. natom) then
      nl = nl + 1
    end if

    !> distance between atoms along edge of cube
    dl = l/nl

    !> place atoms in box
    n = 0
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
          write (14, *) 'Ar', x, y, z! prints out “fort.14” file with coordinates

          if (formula_units == 4) then
            x = i*dl - hl
            y = (j + 0.5)*dl - hl
            z = (k + 0.5)*dl - hl
            write (14, *) 'Ar', x, y, z

            x = (i + 0.5)*dl - hl
            y = j*dl - hl
            z = (k + 0.5)*dl - hl
            write (14, *) 'Ar', x, y, z

            x = (i + 0.5)*dl - hl
            y = (j + 0.5)*dl - hl
            z = k*dl - hl
            write (14, *) 'Ar', x, y, z
          end if

          n = n + formula_units
          if (n == natom) then
            exit loop1
          end if

        end do
      end do
    end do loop1

  end subroutine build_grid

end module build_box

