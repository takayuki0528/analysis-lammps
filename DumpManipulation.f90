module DumpManipulation
	use Globals
	implicit none
	private

	public :: ReadDump, CalcMolCoM, WrapAtom, WrapMolecule

contains

!===============================================================================
! Public Subroutines
!===============================================================================

subroutine ReadDump(nDump)

	integer, intent(in) :: nDump

  integer :: unitDump
  character(64) :: nameDump

  integer :: check_nDump, check_NumAtoms, tempID, tempTyp
  double precision :: tempPos(NUM_DIMS), tempVel(NUM_DIMS)

	integer :: iDim, counter

!-------------------------------------------------------------------------------

	unitDump = 100+nDump
	write(nameDump, "(A, i0, A)" )  trim(DumpPrefix), nDump, trim(DumpSuffix)
	open(unit = unitDump, file = nameDump, status = 'old')

	read(unitDump,'()')
  read(unitDump,*) check_nDump
  if (nDump .ne. check_nDump) then
    print *, "ERROR: Invalid dump number -> ", nDump, check_nDump
  end if
  read(unitDump,'()')
  read(unitDump,*) check_NumAtoms
  if (NumAtoms .ne. check_NumAtoms) then
    print *, "ERROR: Invalid number of atoms -> ", NumAtoms, check_NumAtoms
  end if
  read(unitDump,'()')
  do iDim = DIM_X, DIM_Z
    if (Periodicity(iDim)) then
      read(unitDump,*) c%minEdge(iDim), c%maxEdge(iDim)
      c%length(iDim) = c%maxEdge(iDim) - c%minEdge(iDim)
    else
      read(unitDump,'()')
			c%minEdge(iDim) = 0.0d0
			c%maxEdge(iDim) = 0.0d0
			c%length(iDim) = 0.0d0
    end if
  end do
  read(unitDump,'()')

	counter = 0
  do while (counter < NumAtoms)
    read(unitDump,*) tempID, tempTyp, &
      & tempPos(DIM_X), tempPos(DIM_Y), tempPos(DIM_Z), &
      & tempVel(DIM_X), tempVel(DIM_Y), tempVel(DIM_Z)
    a(tempID)%typ = tempTyp
    a(tempID)%pos(:) = tempPos(:)
    a(tempID)%vel(:) = tempVel(:)
		counter = counter+1
  end do

	close(unit = unitDump)

	print *, "Finish: Reading Dump - ", nDump

end subroutine ReadDump

!===============================================================================

subroutine CalcMolCoM

  double precision :: tempSumMass
  double precision :: tempSumPos(NUM_DIMS), tempSumVel(NUM_DIMS)

  integer :: iMol, iAtom

!-------------------------------------------------------------------------------

  do iMol = 1, NumMolecules
    tempSumMass = 0.0d0
    tempSumPos(:) = 0.0d0
    tempSumVel(:) = 0.0d0
    do iAtom = m(iMol)%minAtomID, m(iMol)%maxAtomID
      tempSumMass = tempSumMass + ATProp(1,a(iAtom)%typ)
      tempSumPos(:) = tempSumPos(:) + ATProp(1,a(iAtom)%typ)*a(iAtom)%pos(:)
      tempSumVel(:) = tempSumVel(:) + ATProp(1,a(iAtom)%typ)*a(iAtom)%vel(:)
    end do
		m(iMol)%mass = tempSumMass
    m(iMol)%pos(:) = tempSumPos(:)/tempSumMass
    m(iMol)%vel(:) = tempSumVel(:)/tempSumMass
  end do

	print *, "Finish: Calculating CoM of Molecules"

end subroutine CalcMolCoM

!===============================================================================

subroutine WrapAtom

  integer :: iMol, iAtom, iDim

!-------------------------------------------------------------------------------

  do iMol = 1, NumMolecules
    do iAtom = m(iMol)%minAtomID, m(iMol)%maxAtomID
      a(iAtom)%wpos(:) = a(iAtom)%pos(:)
      do iDim = DIM_X, DIM_Z
				if (.not.Periodicity(iDim)) cycle
        do while (a(iAtom)%wpos(iDim) < c%minEdge(iDim))
          a(iAtom)%wpos(iDim) = a(iAtom)%wpos(iDim) + c%length(iDim)
        end do
        do while (c%maxEdge(iDim) <= a(iAtom)%wpos(iDim))
          a(iAtom)%wpos(iDim) = a(iAtom)%wpos(iDim) - c%length(iDim)
        end do
      end do
    end do
  end do

	print *, "Finish: Wrapping Atoms"

end subroutine WrapAtom

!===============================================================================

subroutine WrapMolecule

  integer :: iMol, iDim

!-------------------------------------------------------------------------------

  do iMol = 1, NumMolecules
    m(iMol)%wpos(:) = m(iMol)%pos(:)
    do iDim = DIM_X, DIM_Z
			if (.not.Periodicity(iDim)) cycle
      do while (m(iMol)%wpos(iDim) < c%minEdge(iDim))
        m(iMol)%wpos(iDim) = m(iMol)%wpos(iDim) + c%length(iDim)
      end do
      do while (c%maxEdge(iDim) <= m(iMol)%wpos(iDim))
        m(iMol)%wpos(iDim) = m(iMol)%wpos(iDim) - c%length(iDim)
      end do
    end do
  end do

	print *, "Finish: Wrapping Molecules"

end subroutine WrapMolecule

end module
