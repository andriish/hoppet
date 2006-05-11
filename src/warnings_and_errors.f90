!======================================================================
!
! Routine for dealing with warnings that may crop up many times and 
! which should be output only a limited number of times
!
! In f95 could have used a special type which was preinitialised...
!
! $Id: warnings_and_errors.f90,v 1.4 2003/01/08 19:18:52 gsalam Exp $
!======================================================================
module warnings_and_errors
  use types
  implicit none
  private

  integer, parameter :: base = 10000
  integer            :: n_warn_sources   = 0
  integer, parameter, public :: warn_id_INIT=0
  integer, parameter, public :: default_max_warn = 5
  public :: wae_warn, wae_error

contains
  !---------------------------------------------------------------------
  ! On a call with warn_id <= 0 it allocates a new warn_id
  ! REMEMBER TO GIVE WARN_ID THE SAVE ATTRIBUTE
  !
  ! The warn_id should be saved: it indicates both the source of the 
  ! warning (so as to identify it as being from the same source on future 
  ! calls) and the number of warnings of this source written so far.
  subroutine wae_warn(max_warn, warn_id, text, text2, text3, intval, dbleval)
    integer, intent(in)          :: max_warn
    integer, intent(inout)       :: warn_id
    character(len=*), intent(in) :: text
    character(len=*), intent(in), optional :: text2
    character(len=*), intent(in), optional :: text3
    integer,          intent(in), optional :: intval
    real(kind(1d0)),  intent(in), optional :: dbleval
    !--------------------------------------
    integer :: warn_index, nwarn

    !-- generate a new warn_id
    if (warn_id <= 0) then
       n_warn_sources = n_warn_sources + 1
       warn_id = n_warn_sources * base
    end if
    
    warn_index = warn_id / base
    nwarn      = warn_id - warn_index*base

    if (nwarn < max_warn) then
       if (max_warn > base-2) call wae_error('wae_warn',&
            & 'max_warn exceeded maximum allowed value; message was', text)
       !-- does this make any sense at all (GPS 8/1/03)?
       if (warn_id > huge(warn_id)) call wae_error('wae_warn',&
            & 'exceeded max capicity for distinct warnings; message was', text)
    
       warn_id = warn_id + 1
       write(0,'(a)', advance='no') 'WARNING in '
       write(0,'(a)')  text
       if (present(text2)) write(0,'(a)') text2
       if (present(text3)) write(0,'(a)') text3
       if (present(intval))  write(0,*) intval
       if (present(dbleval)) write(0,*) dbleval
       !-- if there is only an 1 warning to be written then
       !   avoid cluttering screen with too many messages
       if (nwarn == max_warn - 1 .and. max_warn > 1) write(0,'(a)') &
            &'----- No more such warnings will be issued ------'
    end if
  end subroutine wae_warn


  !======================================================================
  ! Report an error and then crash if possible
  !======================================================================
  subroutine wae_error(text1, text2, text3,intval,dbleval)
    character(len=*), intent(in) :: text1
    character(len=*), intent(in), optional :: text2
    character(len=*), intent(in), optional :: text3
    integer,          intent(in), optional :: intval
    real(kind(1d0)),  intent(in), optional :: dbleval
    real :: a,b

    write(0,*)
    write(0,'(a)') '============================================================='
    write(0,'(a)', advance='no') 'FATAL ERROR in '
    write(0,'(a)') text1
    if (present(text2))   write(0,'(a)') text2
    if (present(text3))   write(0,'(a)') text3
    if (present(intval))  write(0,*) intval
    if (present(dbleval)) write(0,*) dbleval

    write(0,*)
    write(0,'(a)') &
         &'----- Error handler will now attempt to dump core and  stop'
    a = 1.0
    b = 1.0
    write(0,*) 1.0/sqrt(a-b)
    !-- in case division by zero didn't solve problem
    stop
  end subroutine wae_error
  
end module warnings_and_errors

