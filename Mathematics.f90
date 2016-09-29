module Mathematics
	use Globals
	implicit none
	private

  public :: eigen33, calcStatistics

contains

!===============================================================================
! Public Subroutines
!===============================================================================

! 3x3行列の固有値問題を解く
subroutine eigen33(matrix_original,eigenvalues,eigenvectors)

	double precision, intent(in) :: matrix_original(3,3)
	double precision, intent(out) :: eigenvalues(3)
	double precision, intent(out) :: eigenvectors(3,3)

	double precision ::	matrix(3,3), temp_matrix(3,3)
	double precision ::	temp_eigenvectors(3,3)

	double precision ::	diag_max, temp_max
	double precision ::	alpha, beta, gamma, sintheta, costheta, ss, cc, sc

	double precision :: threshold = 0.000001d0

	integer :: i, j, ro, co
!-------------------------------------------------------------------------------

	diag_max = 0.0d0
	do i = 1, 3
		if (diag_max < matrix_original(i,i)) then
			diag_max = matrix_original(i,i)
		end if
	end do
	matrix(:,:) = matrix_original(:,:)/diag_max

	eigenvalues(:) = 0.0d0
	eigenvectors(:,:) = 0.0d0
	do i = 1, 3
		eigenvectors(i,i) = 1.0d0
	end do

	temp_max = threshold

	do while (.true.)

		temp_max = dabs(matrix(2,1))
		ro = 1
		co = 2
		do i = 1, 2
			if (temp_max < dabs(matrix(3,i))) then
				temp_max = dabs(matrix(3,i))
				ro = i
				co = 3
			end if
		end do

		if (temp_max < threshold) exit

		alpha = (matrix(ro,ro)-matrix(co,co))*0.5d0
		beta = -matrix(co,ro)
		gamma = dabs(alpha)/dsqrt(alpha*alpha+beta*beta)
		sintheta = dsqrt((1.0d0-gamma)*0.5d0)
		costheta = dsqrt((1.0d0+gamma)*0.5d0)
		if (alpha*beta < 0.0d0) then
			sintheta = -sintheta
		end if
		ss = sintheta*sintheta
		cc = costheta*costheta
		sc = sintheta*costheta

		temp_matrix(:,:) = matrix(:,:)
		temp_matrix(ro,ro) &
			& = cc*matrix(ro,ro) + ss*matrix(co,co) - 2.0d0*sc*matrix(co,ro)
		temp_matrix(co,co) &
			& = ss*matrix(ro,ro) + cc*matrix(co,co) + 2.0d0*sc*matrix(co,ro)
		temp_matrix(co,ro) &
			& = sc*(matrix(ro,ro)-matrix(co,co)) + (cc-ss)*matrix(co,ro)
		temp_matrix(ro,co) = temp_matrix(co,ro)
		do i = 1, 3
			if (i /= ro .and. i /= co) then
				temp_matrix(i,ro) = costheta*matrix(i,ro) - sintheta*matrix(i,co)
				temp_matrix(i,co) = sintheta*matrix(i,ro) + costheta*matrix(i,co)
				temp_matrix(ro,i) = costheta*matrix(ro,i) - sintheta*matrix(co,i)
				temp_matrix(co,i) = sintheta*matrix(ro,i) + costheta*matrix(co,i)
			end if
		end do
		matrix(:,:) = temp_matrix(:,:)

		temp_eigenvectors(:,:) = eigenvectors(:,:)
		do i = 1, 3
				temp_eigenvectors(ro,i) &
					& = costheta*eigenvectors(ro,i) - sintheta*eigenvectors(co,i)
				temp_eigenvectors(co,i) &
					& = sintheta*eigenvectors(ro,i) + costheta*eigenvectors(co,i)
		end do
		eigenvectors(:,:) = temp_eigenvectors(:,:)

	end do

	do i = 1,3
		eigenvalues(i) = matrix(i,i)
	end do

end subroutine eigen33

!===============================================================================

! データ集団data(nDim,nMol,nDump)に対して統計量の計算を行う．
! nSubごとの平均値経時変化
! 全体の平均値
! 全体の不偏標準偏差
! Molに対しての系統不偏標準偏差
! Molに対しての偶然不偏標準偏差
! Dumpに対しての系統不偏標準偏差
! Dumpに対しての偶然不偏標準偏差
! 全体の標準偏差
subroutine calcStatistics(rawdata,nSub,statistics)

	double precision, allocatable, intent(in) :: rawdata(:,:,:)
	integer, intent(in) :: nSub
	double precision, allocatable, intent(out) :: statistics(:,:)

	double precision, allocatable :: subSum(:)
	double precision, allocatable :: totSum(:), molSum(:,:), dumpSum(:,:)
	double precision, allocatable :: totAve(:), molAve(:,:), dumpAve(:,:)
	double precision, allocatable :: totVar(:)
	double precision, allocatable :: molSysVar(:), molAccVar(:)
	double precision, allocatable :: dumpSysVar(:), dumpAccVar(:)

	integer :: dataShape(3), numTimeRows

	integer :: iDump, iMol, iDim, i

!-------------------------------------------------------------------------------

	dataShape = shape(rawdata)
	if (mod(dataShape(3),nSub) == 0) then
		numTimeRows = dataShape(3)/nSub
	else
		print *, "ERROR: Invalid nSub for statistics calculation"
		stop
	end if
	allocate(statistics(dataShape(1),numTimeRows+6))

	allocate(subSum(dataShape(1)))
	allocate(totSum(dataShape(1)))
	allocate(molSum(dataShape(1),dataShape(2)))
	allocate(dumpSum(dataShape(1),dataShape(3)))
	subSum(:) = 0.0d0
	totSum(:) = 0.0d0
	molSum(:,:) = 0.0d0
	dumpSum(:,:) = 0.0d0

	i = 0

	do iDump = 1, dataShape(3)

		do iMol = 1, dataShape(2)
			subSum(:) = subSum(:) + rawdata(:,iMol,iDump)
			molSum(:,iMol) = molSum(:,iMol) + rawdata(:,iMol,iDump)
			dumpSum(:,iDump) = dumpSum(:,iDump) + rawdata(:,iMol,iDump)
		end do

		if (mod(iDump,nSub) == 0) then
			i = i+1
			statistics(:,i) = subSum(:)/dble(dataShape(2)*nSub)
			totSum(:) = totSum(:) + subSum(:)
			subSum(:) = 0.0d0
		end if

	end do

	allocate(totAve(dataShape(1)))
	allocate(molAve(dataShape(1),dataShape(2)))
	allocate(dumpAve(dataShape(1),dataShape(3)))
	totAve(:) = totSum(:)/dble(dataShape(2)*dataShape(3))
	molAve(:,:) = molSum(:,:)/dble(dataShape(3))
	dumpAve(:,:) = dumpSum(:,:)/dble(dataShape(2))

	i = i+1
	statistics(:,i) = totAve(:)

	allocate(totVar(dataShape(1)))
	allocate(molSysVar(dataShape(1)))
	allocate(molAccVar(dataShape(1)))
	allocate(dumpSysVar(dataShape(1)))
	allocate(dumpAccVar(dataShape(1)))
	totVar(:) = 0.0d0
	molSysVar(:) = 0.0d0
	molAccVar(:) = 0.0d0
	dumpSysVar(:) = 0.0d0
	dumpAccVar(:) = 0.0d0

	do iDump = 1,dataShape(3)
		do iMol = 1, dataShape(2)
			totVar(:) = totVar(:) + (rawdata(:,iMol,iDump)-totAve(:))**2
			molAccVar(:) = molAccVar(:) &
				& + (rawdata(:,iMol,iDump)-molAve(:,iMol))**2
			dumpAccVar(:) = dumpAccVar(:) &
				& + (rawdata(:,iMol,iDump)-dumpAve(:,iDump))**2
		end do
	end do

	do iMol = 1, dataShape(2)
		molSysVar(:) = molSysVar(:) + (molAve(:,iMol)-totAve(:))**2
	end do

	do iDump = 1, dataShape(3)
		dumpSysVar(:) = dumpSysVar(:) + (dumpAve(:,iDump)-totAve(:))**2
	end do

	statistics(:,i+1) = dsqrt(totVar(:)/dble(dataShape(2)*dataShape(3)-1))
	statistics(:,i+2) = dsqrt(molSysVar(:)/dble(dataShape(2)-1))
	statistics(:,i+3) = dsqrt(molAccVar(:)/dble(dataShape(2)*(dataShape(3)-1)))
	statistics(:,i+4) = dsqrt(dumpSysVar(:)/dble(dataShape(3)-1))
	statistics(:,i+5) = dsqrt(dumpAccVar(:)/dble(dataShape(3)*(dataShape(2)-1)))

	if (i /= numTimeRows+4) then
		print *, "ERROR: Invalid i for statistics calculation"
	end if

end subroutine calcStatistics

end module
