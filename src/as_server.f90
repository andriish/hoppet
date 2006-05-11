!---------------------------------------------------------------
! provides facilities related to alfas calculations
!
! there are two ways of accessing it: either with a handle which 
! contains the details of the particular alfas which is being dealt
! with (e.g. one loop etc...), or without a handle, in which case the
! global handle is used.
!
! GPS, with core based on BRW.
! Tested 25/04/00
!
! Innards may change one day... Make changes to utils version only,
! and then copy across to here.
!-------------------------------------------------------------------
module as_server
  use new_as; use types; use consts_dp
  use warnings_and_errors
  implicit none
  private
  !integer,  parameter :: dp = kind(1d0)
  !integer,  parameter :: sp = kind(1.0)
  !real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp

  real(dp), parameter :: rmass(6) = (/0D0,0D0,.15D0,1.7D0,5.3D0,175D0/)
  !real(dp), parameter :: rmass(6) = (/0D0,0D0,1D0,1.1D0,1.2D0,2005D0/)

  !-- inclusion of CN is a rather tricky way of introducing 
  !   support for 1 and 2-loops in a fairly similary way
  !   and B(3:6) is a way of introducing support for fixed or variable
  !   flavour numbers
  type as_handle 
     private
     real(dp) :: QCDL5
     real(dp) :: BN(3:6), CN(3:6), CN5(3:6)
     integer  :: nloop, nf
     type(na_handle) :: nah
     logical  :: use_nah
  end type as_handle
  type(as_handle), target, save :: ash_global
  public :: as_handle, ash_global

  interface as_Init
     module procedure  as_Init_ash, as_Init_noash
  end interface
  public :: as_Init
  interface as_Value
     module procedure  as_Value_ash, as_Value_noash
  end interface
  public :: as_Value
  interface as_Del
     !-- for digital f90 the second one causes problems
     module procedure  as_Del_ash!, as_Del_noash
  end interface
  public :: as_Del
  public :: as_nfRange, as_nfAtQ, as_QRangeAtNf

  real(dp), parameter :: CAFAC = 3.0_dp, CFFAC = 4.0_dp/3.0_dp
  real(dp), parameter :: PIFAC = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: B3=((11.0_dp*CAFAC)- 6.0_dp)/(12.0_dp*PIFAC) 
  real(dp), parameter :: B4=((11.0_dp*CAFAC)- 8.0_dp)/(12.0_dp*PIFAC) 
  real(dp), parameter :: B5=((11.0_dp*CAFAC)-10.0_dp)/(12.0_dp*PIFAC) 
  real(dp), parameter :: B6=((11.0_dp*CAFAC)-12.0_dp)/(12.0_dp*PIFAC) 
  !-- with absoft there are some difficulties in accessing this quant
  !   from outside 
  real(dp), parameter, public :: as_BN(3:6) = (/ B3, B4, B5, B6/)
  real(dp), parameter :: trueC3=((17.0_dp*CAFAC**2)-(5.0_dp*CAFAC+3.0_dp&
       &*CFFAC)*3.0_dp) /(24.0_dp*PIFAC**2)/B3**2 
  real(dp), parameter :: trueC4=((17.0_dp*CAFAC**2)-(5.0_dp*CAFAC+3.0_dp&
       &*CFFAC)*4.0_dp) /(24.0_dp*PIFAC**2)/B4**2 
  real(dp), parameter :: trueC5=((17.0_dp*CAFAC**2)-(5.0_dp*CAFAC+3.0_dp&
       &*CFFAC)*5.0_dp) /(24.0_dp*PIFAC**2)/B5**2 
  real(dp), parameter :: trueC6=((17.0_dp*CAFAC**2)-(5.0_dp*CAFAC+3.0_dp&
       &*CFFAC)*6.0_dp) /(24.0_dp*PIFAC**2)/B6**2 
  real(dp), parameter :: as_CN(3:6) = (/ trueC3, trueC4, trueC5, trueC6 /)


contains

  !-- alternative versions to allow access to a global alpha
  subroutine as_Init_noash(alfas, Q, qcdl5, nloop, fixnf,  &
       & muMatch_mQuark, use_nah)
    real(dp), intent(in), optional :: alfas, Q, qcdl5
    integer,  intent(in), optional :: nloop, fixnf
    real(dp), intent(in), optional :: muMatch_mQuark
    logical,  intent(in), optional :: use_nah
    type(as_handle), pointer :: ash
    ash => ash_global
    call as_Init_ash(ash, alfas, Q, qcdl5, nloop, fixnf,  &
       & muMatch_mQuark, use_nah)
  end subroutine as_Init_noash
  
  function as_Value_noash(Q,fixnf) result(as_Value)
    real(dp), intent(in) :: Q
    real(dp)             :: as_Value
    type(as_handle), pointer  :: ash
    integer, intent(in), optional :: fixnf
    ash => ash_global
    as_Value = as_Value_ash(ash, Q, fixnf)
  end function as_Value_noash

  subroutine as_Del_noash()
    type(as_handle), pointer :: ash
    ash => ash_global
    call as_Del_ash(ash)
  end subroutine as_Del_noash
  

  !----------------------------------------------------------------------
  ! The following are essentially interface routines to bryan's 
  ! code for alfas
  ! if present alfas then set ash such alfas at M_Z or Q is equal to alpha_s
  ! Alternatively set qcdl5
  ! Can also optionally modify the number of loops and fix the number
  ! of flavours
  !
  ! As of 30/06/2001, by sepcifying use_nah=.true. it should give
  ! fairly transparent access to the differential equation solution
  ! for alfas, using new_as. The only proviso is that initialisation 
  ! may well be quite a bit slower. Also if we want to reinitialise
  ! then it is necessary to call as_Del
  !
  ! As od 16/01/2001 use_nah=.true. will become the default
  subroutine as_Init_ash(ash, alfas, Q, qcdl5, nloop, fixnf, &
       & muMatch_mQuark, use_nah)
    use assertions
    type(as_handle) :: ash
    real(dp), intent(in), optional :: alfas, Q, qcdl5
    integer,  intent(in), optional :: nloop, fixnf
    real(dp), intent(in), optional :: muMatch_mQuark
    logical,  intent(in), optional :: use_nah
    real(dp) :: dummy, lower, upper, middle
    real(dp) :: Qloc
    real(dp), parameter :: mz = 91.2_dp, eps = 1e-8_dp

    ash%use_nah = default_or_opt(.true., use_nah)
    if (ash%use_nah) then
       call na_Init(ash%nah, alfas, Q, nloop, fixnf, muMatch_mQuark)
       return
    end if

    if (present(nloop)) then
       ash%nloop = nloop
    else
       ash%nloop = 2
    end if

    !-- a fudge to fix the number of flavours
    if (present(fixnf)) then
       ash%BN = as_BN(fixnf)
       ash%CN = as_CN(fixnf)
    else
       ash%BN = as_BN
       ash%CN = as_CN
    end if
    
    !-- a fudge to set the number of loops
    select case (ash%nloop)
    case(1)
       ash%CN   = zero
    case(2)
       ! keep ash%CN = as_CN as it was set before
    case default
       stop 'as_Init: nloop must be either 1 or 2'
    end select

    if (present(alfas)) then
       if (present(Q)) then
          Qloc = Q
       else
          Qloc = mz
       end if
       !-- set up search limits ----------------------------
       lower = 0.01_dp
       upper = 1.0_dp
       !-- we may need to revise lower and upper...
       do
          call priv_as_init_2l(ash,lower); dummy = as_Value(ash,Qloc)
          if (dummy > alfas) then
             !call hwwarn(100)
             upper = lower
             lower = lower**2
             cycle
          end if
          call priv_as_init_2l(ash,upper); dummy = as_Value(ash,Qloc)
          if (dummy < alfas) call hwwarn(101)
          exit
       end do
       do 
          !-- could be more efficient -- but too lazy for now
          !write(0,'(4f20.15)') upper, lower, dummy, alfas
          middle = sqrt(lower*upper)
          call priv_as_init_2l(ash,middle)
          dummy = as_Value(ash,Qloc)
          if (dummy >= alfas) upper = middle
          if (dummy < alfas) lower = middle
          if (upper/lower-1.0_dp < eps ) exit
       end do
    elseif (present(qcdl5)) then
       call priv_as_init_2l(ash,qcdl5)
    else
       !call priv_as_init_2l(ash,0.280_dp)  ! alfas(91) = 0.124
       call priv_as_init_2l(ash,0.214_dp)  ! alfas(91) = 0.117
       !call priv_as_init_2l(ash,0.158_dp)  ! alfas(91) = 0.112
    end if
  end subroutine as_Init_ash


  !------------------------------------------------------------------
  function as_Value_ash(ash, Q, fixnf) result(as_Value)
    real(dp), intent(in) :: Q
    real(dp)             :: as_Value
    type(as_handle)      :: ash
    integer, intent(in), optional :: fixnf
    if (ash%use_nah) then
       as_Value = na_Value(ash%nah, Q, fixnf)
    else
       if (present(fixnf)) call wae_error('as_Value_ash: &
            &fixnf not support for old alpha_s')
       as_Value = hwualf(ash, 1, Q)
    end if
  end function as_Value_ash

  !-----------------------------------------------------------------
  ! do any required cleaning up 
  subroutine as_Del_ash(ash)
    type(as_handle)      :: ash
    if (ash%use_nah) call na_Del(ash%nah)
  end subroutine as_Del_ash

  !-----------------------------------------------------------------
  ! Indicate the range of nf values supported by this handle.
  subroutine as_nfRange(ash, nflo, nfhi)
    type(as_handle), intent(in)  :: ash
    integer,         intent(out) :: nflo, nfhi
    if (ash%use_nah) then
       call na_nfRange(ash%nah,nflo,nfhi)
    else
       call wae_Error('as_nfRange: this routine&
            & is only supported with new alpha_s')
    end if
  end subroutine as_nfRange
  
  !-------------------------------------------------------------
  ! Returns the number of flavours relevant for scale Q
  ! Also returns the range of Qvalues (Qlo, Qhi) for which 
  ! this value of nf remains the same.
  subroutine as_nfAtQ(ash, Q, nfAtQ, Qlo, Qhi, muM_mQ)
    type(as_handle), intent(in)  :: ash
    real(dp),        intent(in)  :: Q
    integer,         intent(out) :: nfAtQ
    real(dp),        intent(out), optional :: Qlo, Qhi
    real(dp),        intent(in),  optional :: muM_mQ
    !-----------------------------
    if (ash%use_nah) then
       call na_nfAtQ(ash%nah, Q, nfAtQ, Qlo, Qhi, muM_mQ)
    else
       call wae_Error('as_nfAtQ: this routine&
            & is only supported with new alpha_s')
    end if
  end subroutine as_nfAtQ


  !-------------------------------------------------------------
  ! returns the Q range for a given value of nf. If supported.
  subroutine as_QRangeAtNf(ash, nflcl, Qlo, Qhi, muM_mQ)
    type(as_handle), intent(in)  :: ash
    integer,         intent(in)  :: nflcl
    real(dp),        intent(out) :: Qlo, Qhi
    real(dp),        intent(in),  optional :: muM_mQ
    if (ash%use_nah) then
       call na_QrangeAtNf(ash%nah, nflcl, Qlo, Qhi, muM_mQ)
    else
       call wae_Error('as_QRangeAtNf: this routine&
            & is only supported with new alpha_s')
    end if
  end subroutine as_QRangeAtNf
  
  !----------------------------------------------------------------------
  ! the routine that does the real initialisation of alfas
  subroutine priv_as_init_2l(ash, qcdl5)
    type(as_handle)      :: ash
    real(dp), intent(in) :: qcdl5
    real(dp) :: dummy
    ash%qcdl5 = qcdl5
    dummy = HWUALF(ash, 0, 1e2_dp)
  end subroutine priv_as_init_2l
  !-----------------------------------------------------------------------
  !     STRONG COUPLING CONSTANT
  !     IOPT == 0  INITIALIZES
  !          == 1  TWO-LOOP, FLAVOUR THRESHOLDS
  !
  ! BRW routine, with various bits and pieces shuffled around or
  ! removed by GPS on 24/04/00
  !---------------------------------------------------------------------
  FUNCTION HWUALF(ash, IOPT, SCALE) 
    REAL(DP) :: HWUALF 
    type(as_handle), target :: ash
    real(dp), intent(in) :: scale
    integer,  intent(in) :: iopt
    !--------------------------------------------------------------------
    REAL(DP)          :: RHO,RAT,RLF
    real(dp), pointer :: CN5(:), CN(:), QCDL5, BN(:)
    integer           :: nf

    CN5   => ash%CN5
    QCDL5 => ash%QCDL5
    CN    => ash%CN
    BN    => ash%BN

    !-- the rest is init related to this value of lambda5
    IF (IOPT == 0) THEN 
       !---QCDL5 IS 5-FLAVOUR LAMBDA-MS-BAR                                    
       !---COMPUTE THRESHOLD MATCHING                                          
       RHO=2.0_dp*LOG(RMASS(6)/QCDL5) 
       RAT=LOG(RHO)/RHO 
       CN5(6)=(ash%BN(5)/(1.0_dp-CN(5)*RAT)-ash%BN(6)/(1.0_dp-CN(6)*RAT))*RHO 
       RHO=2.0_dp*LOG(RMASS(5)/QCDL5) 
       RAT=LOG(RHO)/RHO 
       CN5(4)=(ash%BN(5)/(1.0_dp-CN(5)*RAT)-ash%BN(4)/(1.0_dp-CN(4)*RAT))*RHO 
       RHO=2.0_dp*LOG(RMASS(4)/QCDL5) 
       RAT=LOG(RHO)/RHO 
       CN5(3)=(ash%BN(4)/(1.0_dp-CN(4)*RAT)-&
            & ash%BN(3)/(1.0_dp-CN(3)*RAT))*RHO+CN5(4) 
       CN5(5) = zero
    ENDIF
    IF (SCALE <= QCDL5) then; CALL HWWARN(51); hwualf=0.0_dp; return
    end IF
    RHO=2.0_dp*LOG(SCALE/QCDL5) 
    RAT=LOG(RHO)/RHO 
    !-- this will allow us to fiddle nf later on, according to some
    !   flag in ash
    nf = priv_nf(ash, scale)
    RLF = ash%BN(nf)*RHO/(one - CN(nf)*RAT) + CN5(nf)
!!$    select case(nf)
!!$    case(6)
!!$       RLF=B6*RHO/(1.0_dp-CN(6)*RAT)+CN5(6) 
!!$    case(5)
!!$       RLF=B5*RHO/(1.0_dp-CN(5)*RAT) 
!!$    case(4)
!!$       RLF=B4*RHO/(1.0_dp-CN(4)*RAT)+CN5(4) 
!!$    case default
!!$       RLF=B3*RHO/(1.0_dp-CN(3)*RAT)+CN5(3) 
!!$    end select
    IF (RLF <= ZERO) then; CALL HWWARN(53); hwualf=0.0_dp; return
    end IF

    IF (IOPT == 1) THEN 
       HWUALF=1.0_dp/RLF 
    ENDIF
  end function hwualf

  !---------------------------------------------------------------
  ! returns the appropriate for the given scale
  function priv_nf(ash, scale) result(nf)
    type(as_handle), intent(in) :: ash
    real(dp),        intent(in) :: scale
    integer                     :: nf
    IF (SCALE > RMASS(6)) THEN 
       nf = 6
    ELSEIF (SCALE > RMASS(5)) THEN 
       nf = 5
    ELSEIF (SCALE > RMASS(4)) THEN 
       nf = 4
    ELSE 
       nf = 3
    ENDIF
  end function priv_nf
  

  !-- a little primitive? --------------------------------------------
  SUBROUTINE HWWARN(ICODE) 
    integer, intent(in) :: icode
    real(dp) :: a, b
    write(0,*) ' ALFAS WARNING CODE',ICODE
    a = sqrt(0.0_dp)
    b = a**2
    !write(0,*) 1.0_dp/(a-b)
    stop
  END SUBROUTINE HWWARN


end module as_server

!!$program astest
!!$  use types; use consts_dp
!!$  use as_server
!!$  implicit none
!!$  integer  :: i
!!$  real(dp) :: Q
!!$
!!$  call as_init(alfas=0.118_dp,nloop=2)
!!$  do i = 0,100
!!$     Q = 91.2_dp**(i/100.0_dp)
!!$     write(6,*) Q, as_Value(Q)
!!$  end do
!!$  
!!$end program astest
