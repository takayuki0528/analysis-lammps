program main
	use Globals
	use Analysis
	use ReadOrder
	use DumpManipulation
	implicit none
	include 'mpif.h'

	integer :: i

!-------------------------------------------------------------------------------

	call MPI_Init(ierr)
	call MPI_Comm_Size(MPI_COMM_WORLD,nproc,ierr)
	call MPI_Comm_Rank(MPI_COMM_WORLD,myrank,ierr)
	print *, "My MPI Rank : Total Process = ", myrank, " : ", nproc
	call MPI_Barrier(MPI_COMM_WORLD,ierr)
	if(myrank == 0) then
		print *, "Finish: Initializing MPI"
	end if

	call ReadOrderFile
	call AllocateAnalyses

	do i = myrank+1, NumDumps, nproc

		call ReadDump(InitialDump+(i-1)*DumpInterval)
		if (ToF_CalcMolCoM) call CalcMolCoM
		if (ToF_WrapAtom) call WrapAtom
		if (ToF_WrapMolecule) call WrapMolecule

		call CalcAnalyses(i)

	end do

	call ServeResults

	call MPI_Finalize(ierr)

end program
