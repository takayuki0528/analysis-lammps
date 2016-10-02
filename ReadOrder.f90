module ReadOrder
	use Globals
	use Analysis
	implicit none
	private

	public :: ReadOrderFile

contains

!===============================================================================
! Public Subroutines
!===============================================================================

subroutine ReadOrderFile

	character(64) :: firstChar, orderType

!-------------------------------------------------------------------------------

	open(unit = UnitOrder, file = "Order", status = 'old')

	do while (StatOrder == 0)
		read(UnitOrder,*,iostat=StatOrder) firstChar
		if (firstChar .eq. "$") then
			backspace(UnitOrder)
			read(UnitOrder,*) firstChar, orderType
			if (orderType .eq. "Molecule") then
				call readOrderMolecule
			else if (orderType .eq. "Atom") then
				call readOrderAtom
			else if (orderType .eq. "Bead") then
				call readOrderBead
			else if (orderType .eq. "DumpFile") then
				call readOrderDumpFile
			else if (orderType .eq. "Analysis") then
				call ReadOrderAnalysis
			else
				print *, "ERROR: Unknown order type -> ", orderType
				stop
			end if
		end if
	end do

	close(unit = UnitOrder)

	if (myrank == 0) then
		print *, "Finish: Reading Order"
	end if

	call Prepare

end subroutine ReadOrderFile

!===============================================================================
! Private Subroutines
!===============================================================================

subroutine readOrderMolecule

	character(64) :: curb, dummy																									! curb means "igeta" in Japanese
	integer :: rowCounter, i

!-------------------------------------------------------------------------------

	backspace(UnitOrder)
	read(UnitOrder,*) dummy, dummy, NumRowsOMM
	allocate(OrderMolMatrix(4,NumRowsOMM))

	rowCounter = 0
	do while (.true.)
		read(UnitOrder,*,iostat=StatOrder) curb
		if (StatOrder < 0) exit
		if (curb .eq. "$") then
			backspace(UnitOrder)
			exit
		else if (curb .ne. "#") then
			rowCounter = rowCounter+1
			backspace(UnitOrder)
			read(UnitOrder,*) (OrderMolMatrix(i,rowCounter), i=1,4)
		end if
	end do

	if (NumRowsOMM .ne. rowCounter) then
		print *, "ERROR: Invalid number of OrderMolMatrix's rows"
		stop
	end if

end subroutine readOrderMolecule

!===============================================================================

subroutine readOrderAtom

	character(64) :: 	curb, dummy																									! curb means "igeta" in Japanese
	integer :: 				rowCounter, dummy_i, i

!-------------------------------------------------------------------------------

	backspace(UnitOrder)
	read(UnitOrder,*) dummy, dummy, NumAtomTypes, dummy, dummy, NumColumnsATP
	allocate(ATProp(NumColumnsATP,NumAtomTypes))

	rowCounter = 0
	do while (.true.)
		read(UnitOrder,*,iostat=StatOrder) curb
		if (StatOrder < 0) exit
		if (curb .eq. "$") then
			backspace(UnitOrder)
			exit
		else if (curb .ne. "#") then
			rowCounter = rowCounter+1
			backspace(UnitOrder)
			read(UnitOrder,*) dummy_i, (ATProp(i,rowCounter), i = 1, NumColumnsATP)
		end if
	end do

	if (NumAtomTypes .ne. rowCounter) then
		print *, "ERROR: Invalid number of ATProp's rows"
	end if

end subroutine readOrderAtom

!===============================================================================

subroutine readOrderBead

	character(64) :: 	curb, dummy																									! curb means "igeta" in Japanese

!-------------------------------------------------------------------------------

	do while (.true.)
		read(UnitOrder,*,iostat=StatOrder) curb
		if (StatOrder < 0) exit
		if (curb .eq. "$") then
			backspace(UnitOrder)
			exit
		else if (curb .ne. "#") then

		end if
	end do

end subroutine readOrderBead

!===============================================================================

subroutine readOrderDumpFile

	character(64) :: 	curb, dummy																									! curb means "igeta" in Japanese

!-------------------------------------------------------------------------------

do while (.true.)
	read(UnitOrder,*,iostat=StatOrder) curb
	if (StatOrder < 0) exit
	if (curb .eq. "$") then
		backspace(UnitOrder)
		exit
	else if (curb .ne. "#") then
		backspace(UnitOrder)
		if (curb .eq. "Number_of_Files") then
			read(UnitOrder,*) dummy, NumDumps
		else if (curb .eq. "Initial_Dump") then
			read(UnitOrder,*) dummy, InitialDump
		else if (curb .eq. "Dump_Interval") then
			read(UnitOrder,*) dummy, DumpInterval
		else if (curb .eq. "Prefix_of_FileName") then
			read(UnitOrder,*) dummy, DumpPrefix
		else if (curb .eq. "Suffix_of_FileName") then
			read(UnitOrder,*) dummy, DumpSuffix
		else if (curb .eq. "Periodicity") then
			read(UnitOrder,*) dummy, &
				& Periodicity(DIM_X), Periodicity(DIM_Y), Periodicity(DIM_Z)
		else
			print *, "ERROR: Unknown character for DumpFile -> ", curb
		end if
	end if
end do

end subroutine readOrderDumpFile

!===============================================================================

subroutine Prepare

	integer :: numAtomsPerMol
	integer :: i, j ,iMol

!-------------------------------------------------------------------------------

	NumMoleculeTags = 0
	do i = 1, NumRowsOMM
		if (NumMoleculeTags < OrderMolMatrix(1,i))  then
			NumMoleculeTags = OrderMolMatrix(1,i)
		end if
	end do
	if (myrank == 0) then
		print *, "Number of Molecule Tags : ", NumMoleculeTags
	end if

	NumMolecules = 0
	do i = 1, NumRowsOMM
		if (OrderMolMatrix(1,i) == 0) cycle
		NumMolecules = NumMolecules + OrderMolMatrix(2,i)
	end do
	if (myrank == 0) then
		print *, "Number of Target Molecules : ", NumMolecules
	end if

	allocate(m(NumMolecules))
	do i = 1, NumMolecules
		m(i)%pos(:) = 0.0d0
		m(i)%wpos(:) = 0.0d0
		m(i)%vel(:) = 0.0d0
		m(i)%mass = 0.0d0
		m(i)%tag = 0
		m(i)%minAtomID = 0
		m(i)%maxAtomID = 0
	end do

	NumAtoms = 0
	do i = 1, NumRowsOMM
		NumAtoms = NumAtoms + OrderMolMatrix(4,i)-OrderMolMatrix(3,i)+1
	end do
	if (myrank == 0) then
		print *, "Number of Atoms : ", NumAtoms
	end if

	allocate(a(NumAtoms))
	do i = 1, NumAtoms
		a(i)%pos(:) = 0.0d0
		a(i)%wpos(:) = 0.0d0
		a(i)%vel(:) = 0.0d0
		a(i)%typ = 0
	end do

	iMol = 0
	do i = 1, NumRowsOMM
		if (OrderMolMatrix(1,i) == 0) cycle
		numAtomsPerMol &
			& = (OrderMolMatrix(4,i)-OrderMolMatrix(3,i)+1) / OrderMolMatrix(2,i)
		do j = 1, OrderMolMatrix(2,i)
			iMol = iMol+1
			m(iMol)%tag = OrderMolMatrix(1,i)
			m(iMol)%minAtomID = OrderMolMatrix(3,i) + (j-1)*numAtomsPerMol
			m(iMol)%maxAtomID = OrderMolMatrix(3,i) + j*numAtomsPerMol - 1
			if (myrank == 0) then
				print *, "Molecule ", iMol, " contains ", numAtomsPerMol, &
					& " atoms; from ", m(iMol)%minAtomID, " to ", m(iMol)%maxAtomID
			end if
		end do
	end do

end subroutine Prepare

end module
