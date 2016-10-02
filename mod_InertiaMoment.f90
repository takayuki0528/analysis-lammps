module mod_InertiaMoment
	use Globals
	use Mathematics
	implicit none
	include 'mpif.h'

	!--- From Order
	logical, public :: ToF_InertiaMoment = .false.
	integer :: aveDumps_im
	!--- Output
	integer :: unitIM = 20
	character(64) :: nameIM = "results_IM.dat"
	!--- Memory
	type type_im
		double precision :: matrix(NUM_DIMS,NUM_DIMS)
		double precision :: opVec(NUM_DIMS)																					! Order Parameter
		double precision :: rgVec(NUM_DIMS)																					! Gyration Radius
	end type type_im
	type(type_im), allocatable :: im(:,:)

	private

	public :: readInertiaMoment, allocateInertiaMoment, calcInertiaMoment, &
					& serveInertiaMoment

contains

!===============================================================================
! Public Subroutines
!===============================================================================

subroutine readInertiaMoment

	character(64) :: dummy

!-------------------------------------------------------------------------------

	ToF_InertiaMoment = .true.
	backspace(UnitOrder)
	read(UnitOrder,*) dummy, aveDumps_im

end subroutine readInertiaMoment

!===============================================================================

subroutine allocateInertiaMoment

		integer :: i, j

!-------------------------------------------------------------------------------

	allocate(im(NumMolecules,NumDumps))
	do i = 1, NumDumps
		do j = 1, NumMolecules
			im(j,i)%matrix(:,:) = 0.0d0
			im(j,i)%opVec(:) = 0.0d0
			im(j,i)%opVec(:) = 0.0d0
		end do
	end do

end subroutine allocateInertiaMoment

!===============================================================================

subroutine calcInertiaMoment(nDump)

	integer, intent(in) :: nDump

	double precision ::	tempVec(NUM_DIMS)
	double precision ::	eigens(NUM_DIMS)
	double precision ::	eigenVecs(NUM_DIMS,NUM_DIMS)
	double precision :: temp_cos2(NUM_DIMS)
	double precision :: mineigen, minvector(NUM_DIMS)

	double precision :: inv
	integer :: iMol, iAtom, iDim, jDim, i

!-------------------------------------------------------------------------------

	do iMol = 1, NumMolecules

		do iAtom = m(iMol)%minAtomID, m(iMol)%maxAtomID
			tempVec = a(iAtom)%pos(:)-m(iMol)%pos(:)
			do iDim = DIM_X, DIM_Z
				do jDim = DIM_X, DIM_Z
					im(iMol,nDump)%matrix(jDim,iDim) = im(iMol,nDump)%matrix(jDim,iDim) &
						& - ATProp(1,a(iAtom)%typ)*tempVec(jDim)*tempVec(iDim)
				end do
				im(iMol,nDump)%matrix(iDim,iDim) = im(iMol,nDump)%matrix(iDim,iDim) &
					& + ATProp(1,a(iAtom)%typ)*dot_product(tempVec,tempVec)
			end do
		end do

		call eigen33(im(iMol,nDump)%matrix,eigens(:),eigenVecs(:,:))
		mineigen = eigens(1)
		do iDim = 1, 3
			minvector(iDim) = eigenVecs(1,iDim)
		end do
		do i = 2, 3
			if (eigens(i) < mineigen) then
				mineigen = eigens(i)
				do iDim = 1, 3
					minvector(iDim) = eigenVecs(i,iDim)
				end do
			end if
		end do
		inv = 1.0d0/dot_product(minvector,minvector)
		temp_cos2(:) = inv*minvector(:)**2
		im(iMol,nDump)%opVec(:) = 0.5d0*(3.0d0*temp_cos2(:)-1.0d0)

		do iDim = DIM_X, DIM_Z
			im(iMol,nDump)%rgVec(iDim) &
				& = 0.01d0*im(iMol,nDump)%matrix(iDim,iDim)/m(iMol)%mass
		end do

	end do

	print *, "Finish: Calculating Inertia-Moment - ", nDump

end subroutine calcInertiaMoment

!===============================================================================

subroutine serveInertiaMoment

	double precision, allocatable :: data_im(:,:,:)
	double precision, allocatable :: stat_im(:,:)

	integer :: iDump, iMol, i

!-------------------------------------------------------------------------------

	allocate(data_im(7,NumMolecules,NumDumps))
	data_im(:,:,:) = 0.0d0

	call MPI_Barrier(MPI_COMM_WORLD,ierr)
	do iDump = 1, NumDumps
		do iMol = 1, NumMolecules
			call MPI_Reduce(im(iMol,iDump)%opVec(1),data_im(1,iMol,iDump), &
				& NUM_DIMS,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
			call MPI_Reduce(im(iMol,iDump)%rgVec(1),data_im(4,iMol,iDump), &
				& NUM_DIMS,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
		end do
	end do

	if (myrank == 0) then

		data_im(7,:,:) = 0.5d0*(data_im(4,:,:)+data_im(5,:,:)+data_im(6,:,:))
		call calcStatistics(data_im,aveDumps_im,stat_im)

		open(unit = unitIM, file = nameIM, position = 'append')
		write(unitIM, "('        time', a, &
			& '        OP_x', a, '        OP_y', a, '        OP_z', a, &
		 	& ' Rg_x [nm^2]', a, ' Rg_y [nm^2]', a, ' Rg_z [nm^2]', a, &
			& '   Rg [nm^2]')") &
			& TAB, TAB, TAB, TAB, TAB, TAB, TAB

		do i = 1, NumDumps/aveDumps_im

			write(unitIM, "(i12, a, f12.6, a, f12.6, a, f12.6, a, &
				& f12.6, a, f12.6, a, f12.6, a, f12.6)") i, TAB, &
				& stat_im(1,i), TAB, stat_im(2,i), TAB, stat_im(3,i), TAB, &
				& stat_im(4,i), TAB, stat_im(5,i), TAB, stat_im(6,i), TAB, stat_im(7,i)

		end do

		i = NumDumps/aveDumps_im + 1
		write(unitIM, "('Average    >', a, f12.6, a, f12.6, a, f12.6, a, &
			& f12.6, a, f12.6, a, f12.6, a, f12.6)") TAB, &
			& stat_im(1,i), TAB, stat_im(2,i), TAB, stat_im(3,i), TAB, &
			& stat_im(4,i), TAB, stat_im(5,i), TAB, stat_im(6,i), TAB, stat_im(7,i)

		i = NumDumps/aveDumps_im + 2
		write(unitIM, "('SD-Total   >', a, f12.6, a, f12.6, a, f12.6, a, &
			& f12.6, a, f12.6, a, f12.6, a, f12.6)") TAB, &
			& stat_im(1,i), TAB, stat_im(2,i), TAB, stat_im(3,i), TAB, &
			& stat_im(4,i), TAB, stat_im(5,i), TAB, stat_im(6,i), TAB, stat_im(7,i)

		i = NumDumps/aveDumps_im + 3
		write(unitIM, "('SysSD-Mol  >', a, f12.6, a, f12.6, a, f12.6, a, &
			& f12.6, a, f12.6, a, f12.6, a, f12.6)") TAB, &
			& stat_im(1,i), TAB, stat_im(2,i), TAB, stat_im(3,i), TAB, &
			& stat_im(4,i), TAB, stat_im(5,i), TAB, stat_im(6,i), TAB, stat_im(7,i)

		i = NumDumps/aveDumps_im + 4
		write(unitIM, "('AccSD-Mol  >', a, f12.6, a, f12.6, a, f12.6, a, &
			& f12.6, a, f12.6, a, f12.6, a, f12.6)") TAB, &
			& stat_im(1,i), TAB, stat_im(2,i), TAB, stat_im(3,i), TAB, &
			& stat_im(4,i), TAB, stat_im(5,i), TAB, stat_im(6,i), TAB, stat_im(7,i)

		i = NumDumps/aveDumps_im + 5
		write(unitIM, "('SysSD-Dump >', a, f12.6, a, f12.6, a, f12.6, a, &
			& f12.6, a, f12.6, a, f12.6, a, f12.6)") TAB, &
			& stat_im(1,i), TAB, stat_im(2,i), TAB, stat_im(3,i), TAB, &
			& stat_im(4,i), TAB, stat_im(5,i), TAB, stat_im(6,i), TAB, stat_im(7,i)

		i = NumDumps/aveDumps_im + 6
		write(unitIM, "('AccSD-Dump >', a, f12.6, a, f12.6, a, f12.6, a, &
			& f12.6, a, f12.6, a, f12.6, a, f12.6)") TAB, &
			& stat_im(1,i), TAB, stat_im(2,i), TAB, stat_im(3,i), TAB, &
			& stat_im(4,i), TAB, stat_im(5,i), TAB, stat_im(6,i), TAB, stat_im(7,i)

		close(unit = unitIM)

		print *, 'Finish: Writing Inertia-Moment'

	end if

end subroutine serveInertiaMoment

end module
