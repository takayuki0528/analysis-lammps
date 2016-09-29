module Analysis
	use Globals
	use mod_InertiaMoment
	implicit none
	private

	public :: ReadOrderAnalysis, AllocateAnalyses, CalcAnalyses, ServeResults

contains

!===============================================================================
! Public Subroutines
!===============================================================================

subroutine ReadOrderAnalysis

	character(64) :: curb																													! curb means "igeta" in Japanes

	integer :: i

!-------------------------------------------------------------------------------

	do while (.true.)
		read(UnitOrder,*,iostat=StatOrder) curb
		if (StatOrder < 0) exit
		if (curb .eq. "$") then
			backspace(UnitOrder)
			exit
		else if (curb .ne. "#") then
			if (curb .eq. "InertiaMoment") then
				call readInertiaMoment
			else
				print *, "ERROR: Unknown character for Analysis -> ", curb
			end if
		end if
	end do

	if (ToF_InertiaMoment) then
		ToF_CalcMolCoM = .true.
	end if

	if (.false.) then
		ToF_WrapAtom = .true.
	end if

	if (.false.) then
		ToF_WrapMolecule = .true.
	end if

end subroutine ReadOrderAnalysis

!===============================================================================

subroutine AllocateAnalyses

	if (ToF_InertiaMoment) call allocateInertiaMoment

end subroutine AllocateAnalyses

!===============================================================================

subroutine CalcAnalyses(nDump)

	integer, intent(in) :: nDump

!-------------------------------------------------------------------------------

	if (ToF_InertiaMoment) call calcInertiaMoment(nDump)

end subroutine CalcAnalyses

!===============================================================================

subroutine ServeResults

	if (ToF_InertiaMoment) call serveInertiaMoment

end subroutine ServeResults

end module
