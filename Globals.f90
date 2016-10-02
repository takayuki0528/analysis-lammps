module Globals
	implicit none

!===============================================================================
! Non Editabl Parameters
!===============================================================================

	integer, parameter :: NUM_DIMS = 3
	integer, parameter :: DIM_X = 1
	integer, parameter :: DIM_Y = 2
	integer, parameter :: DIM_Z = 3

	double precision, parameter :: PI		= 3.14159265358979323846264338327950288d0
	double precision, parameter :: R_gc = 8.3144621d0

	character, parameter :: TAB = '	'

!===============================================================================
! Parameters Read from Order
!===============================================================================

	!=== Order Fiel ===!
	integer :: UnitOrder = 10
	integer :: StatOrder = 0

	!=== Molecule & Atom ===!

	integer :: NumMoleculeTags
	integer :: NumMolecules

	integer :: NumRowsOMM
	integer, allocatable :: OrderMolMatrix(:,:)

	integer :: NumAtomTypes
	integer :: NumAtoms

	integer :: NumColumnsATP
	double precision, allocatable :: ATProp(:,:)																	! 1st Property must be Mass

	!=== Dump File ===!

	integer :: NumDumps
	integer :: InitialDump
	integer :: DumpInterval
	character(64) :: DumpPrefix
	character(64) :: DumpSuffix
	logical :: Periodicity(NUM_DIMS)

	!=== Dump File Manipulation ===!

	logical :: ToF_CalcMolCoM			= .false.
	logical :: ToF_WrapAtom				= .false.
	logical :: ToF_WrapMolecule		= .false.

!=======================================================================
! Non Editable Variables
!=======================================================================

	!=== MPI ===!
	integer:: nproc, ierr, myrank

	!=== Global Types ===!

	type cell
	  double precision :: minEdge(NUM_DIMS) = 0.0d0
		double precision :: maxEdge(NUM_DIMS) = 0.0d0
		double precision :: length(NUM_DIMS) = 0.0d0
	end type cell
	type(cell) :: c

	type molecule
	  double precision :: pos(NUM_DIMS), wpos(NUM_DIMS)
		double precision :: vel(NUM_DIMS)
		double precision :: mass
	  integer :: tag, minAtomID, maxAtomID
	end type molecule
	type(molecule), allocatable :: m(:)

	type atom
	  double precision :: pos(NUM_DIMS), wpos(NUM_DIMS)
		double precision :: vel(NUM_DIMS)
	  integer :: typ
	end type atom
	type(atom), allocatable :: a(:)

end module
