! $Id: convolution.f90,v 1.26 2005/07/08 21:19:17 salam Exp $


!======================================================================
! Module exlcusively for communication between convolution and 
! any routine which might supply convolutions
!======================================================================
module convolution_communicator
  use types
  implicit none
  private

  integer, public , parameter :: cc_REAL=1,cc_VIRT=2,&
       & cc_REALVIRT=3,cc_DELTA =4 
  integer, public, save :: cc_piece
end module convolution_communicator








!======================================================================
! All the base types and routines for happy convolution
!
! Decide on the following policy: grid details will be copied rather
! than pointed to (avoids danger of object that is pointed to being 
! changed). Not clear if this is really the best policy, but we will
! give it a try...
!
! GPS 04/04/2000
!======================================================================
module convolution
  use types; use consts_dp; use assertions
  implicit none
  private

  integer, parameter :: conv_UndefinedInt = 1000000004

  !-------------------------------
  ! definition of a grid
  ! Includes possibility of one layer of subsiduary definitions
  type grid_def
     real(dp) :: dy, ymax
     integer  :: ny, order, nsub
     logical  :: locked
     integer,        pointer :: subiy(:) ! starting points of subsiduary grid
     type(grid_def), pointer :: subgd(:) ! subsiduary grid defs
  end type grid_def

  !--------------------------------------------------
  ! used for abbreviated access to conv_EvalGridQuant
  
  
  type grid_conv
     type(grid_def)           :: gd
     real(dp),        pointer :: conv(:,:)
     !-- support in construction...
     type(grid_conv), pointer :: subgc(:)
  end type grid_conv
  integer, parameter :: FULL=1, UPPR=2
  !-- for standard linear approach with proper end-point treatment.
  ! This must remain zero otherwise inconsistencies will arise
  integer, parameter :: LIN_ORDER=0
  

  public :: grid_def, grid_conv


  !-- public routines --
  interface conv_InitGridDef
     module procedure conv_InitGridDef_single, conv_InitGridDef_multi
  end interface
  public :: conv_InitGridDef_single, conv_InitGridDef_multi
  public :: conv_InitGridDef, conv_ValidateGD
  interface operator(==)
     module procedure conv_CmpGridDef
  end interface
  public :: operator(==)

  !-- quant routines -----------------------------------------------
  interface conv_AllocGridQuant
     module procedure conv_AllocGridQuant_1d,conv_AllocGridQuant_2d,&
          & conv_AllocGridQuant_3d
  end interface
  public :: conv_AllocGridQuant
  interface conv_InitGridQuant
     module procedure conv_InitGridQuant_func, conv_InitGridQuant_func2d,&
          & conv_InitGridQuant_func_a, conv_InitGridQuant_func2d_a,&
          & conv_InitGridQuant_func_ai, conv_InitGridQuant_func2d_ai
  end interface
  ! standard version gives problems with pgf90 5.1-2 (wrong answers)...
  interface conv_InitGridQuantSub
     module procedure conv_InitGridQuantSub_2d, &
          &conv_InitGridQuantSub_2d_a, conv_InitGridQuantSub_2d_ai
  end interface
  ! This is the workaround suggested by Giulia, across all variants
  ! of the subroutine... Testing with pgf90 however shows that we still 
  ! have the wrong answer (though whether because of this routine
  ! or others has not been established).
  !interface conv_InitGridQuantSub
  !   module procedure conv_InitGridQuantSub_2d_pgf, &
  !        &conv_InitGridQuantSub_2d_a_pgf, conv_InitGridQuantSub_2d_ai_pgf
  !end interface
  
  public :: conv_InitGridQuant, conv_InitGridQuantSub, conv_PrintGridQuant
  public :: conv_InitGridQuantLHAPDF
  interface conv_PrintGridQuant
     module procedure conv_PrintGridQuant_1, conv_PrintGridQuant_2,&
          & conv_PrintGridQuant_3, conv_PrintGridQuant_4
  end interface
  interface conv_EvalGridQuant
     module procedure conv_EvalGridQuant_0d, conv_EvalGridQuant_1d
  end interface
  public :: conv_MomGridQuant, conv_EvalGridQuant, conv_WgtGridQuant
  interface conv_DelGridQuant
     module procedure conv_DelGridQuant_1d, conv_DelGridQuant_2d
  end interface
  
  public :: conv_DelGridQuant

  !-- convolution routines -------------------------------------------
  interface conv_AllocGridConv
     module procedure conv_AllocGridConv_0d, conv_AllocGridConv_1d, &
          & conv_AllocGridConv_2d
  end interface
  interface conv_InitGridConv
     module procedure conv_InitGridConv_zero, conv_InitGridConv_zero_1d,&
          & conv_InitGridConv_zero_2d, conv_InitGridConv_func,&
          & conv_InitGridConv_gc, conv_InitGridConv_gc_1d, &
          & conv_InitGridConv_gc_2d, conv_InitGridConv_conv
  end interface
  !-- keep following for historical reasons? Not as of 28/12/01
  !interface conv_GridConvAdd
  !   module procedure conv_GridConvAdd_func, conv_GridConvAdd_gc
  !end interface
  interface conv_ZeroGridConv
     module procedure conv_ZeroGridConv_0d, conv_ZeroGridConv_1d,&
          & conv_ZeroGridConv_2d
  end interface
  interface conv_AddGridConv
     module procedure conv_AddGridConv_func, conv_AddGridConv_gc,&
          & conv_AddGridConv_gc_1d, conv_AddGridConv_gc_2d
  end interface
  interface conv_MultGridConv
     module procedure conv_MultGridConv_0d, conv_MultGridConv_1d, &
          & conv_MultGridConv_2d
  end interface
  interface conv_DelGridConv
     module procedure conv_DelGridConv_0d, conv_DelGridConv_1d, &
          & conv_DelGridConv_2d
  end interface
  interface conv_ConvGridConv
     module procedure conv_ConvGridConv_0d, conv_ConvGridConv_2dx2d
  end interface

  public :: conv_Seteps
  public :: conv_AllocGridConv, conv_InitGridConv, conv_DelGridConv
  public :: conv_AddGridConv, conv_MultGridConv, conv_ConvGridConv
  public :: conv_ZeroGridConv, conv_CommGridConv

  !-- REMEMBER: .conv. seems to have a precedence similar to .or., .and.
  !             and so is lower precednece than any arithmetic operation
  interface operator(.conv.)
     module procedure conv_ConvGridQuant_mat, conv_ConvGridQuant_scalar
  end interface
  public :: operator(.conv.)
  !-- put in this version as well because to ease the writing of 
  !   expressions: * should have the right precedence properties
  !   where .conv. seems to have a very low precedence
  interface operator(*)
     module procedure conv_ConvGridQuant_mat, conv_ConvGridQuant_scalar
  end interface
  public :: operator(*)

  !-- keep this for historical reasons (hopefully not too long)
  interface conv_ConvConv
     module procedure conv_InitGridConv_conv
  end interface
  public :: conv_ConvConv

  !--------------------------
  type gdval
     type(grid_def) :: gd
     real(dp)       :: val
  end type gdval
  public :: gdval
  interface operator(.with.)
     module procedure conv_gdval_gdv, conv_gdval_vgd
  end interface
  public :: operator(.with.)  
  interface operator(.atx.)
     module procedure conv_EvalGridQuant_atx, conv_EvalGridQuant_atx_1d
  end interface
  public :: operator(.atx.)

  !-- precision used for integration
  real(dp) :: eps=1e-7_dp
  !real(dp), parameter :: eps=1e-8_dp
  real(dp), parameter :: warn_tolerance = 1e-3_dp
  
  !-- used for generation of derived convoluters (e.g. exponentiations
  !   of existsing convoluters)
  logical :: override_grid_locking = .false.
  integer :: nconv_with_override_off = 0 ! 
  public :: conv_GetDerivedProbes, conv_SetDerivedConv
  public :: conv_SetDerivedConv_nodealloc

contains

  subroutine conv_WelcomeMessage
    write(0,*) '-----------------------------------------------------------'
    write(0,*) '             This is the PDFevln package'
    write(0,*) ' '
    write(0,*) '        Written by Gavin P. Salam (2001-2006)'
    write(0,*) ' '
    write(0,*) ' It is made available under the GNU public license,'
    write(0,*) ' with the additional restriction that if you use it or any'
    write(0,*) ' derivative of it in scientific work then you should cite:'
    write(0,*) ' M. Dasgupta and G.P. Salam, Eur.Phys.J.C24:213-236,2002.'
    write(0,*) ' '
    write(0,*) ' You are also encouraged to cite the original references,'
    write(0,*) ' for LO, NLO and NNLO splitting functions, the QCD'
    write(0,*) ' 1, 2 and 3 loop beta functions and the coupling and '
    write(0,*) ' PDF mass threshold matching functions.'
    write(0,*) '-----------------------------------------------------------'
    write(0,*) ' '
  end subroutine conv_WelcomeMessage
  
  
  !======================================================================
  ! Things just for Grid defs.
  ! updated for multi.
  !
  ! order = 0 -> standard linear
  ! order > 0 -> full (order) treatment
  ! order < 0 -> treatment correct to |order| except at end points 
  !              (this is similar to Ratcliffes proposal, hep-ph/0012376)
  subroutine conv_InitGridDef_single(gd,dy,ymax,order)
    type(grid_def), intent(out) :: gd
    real(dp),       intent(in)  :: dy, ymax
    integer,        intent(in), optional :: order
    logical, save :: first_call = .true.

    if (first_call) then
       first_call = .false.
       call conv_WelcomeMessage
    end if
    
    !-- this is a plain grid def ---------------------------
    gd%nsub = 0
    nullify(gd%subiy,gd%subgd)
    gd%locked = .false.
    !-------------------------------------------------------

    gd%ymax = ymax
    gd%ny   = nint(ymax / dy)
    if (gd%ny < 2) then
       write(0,*) 'conv_InitGridDef: requested too small a number of bins'
       write(0,*) '                 dy and ymax were',dy,ymax
       stop
    end if
    
    gd%dy   = ymax / gd%ny
    if (abs(gd%dy/dy - one) > 0.001_dp) then
       write(0,*) 'conv_InitGridDef: requested dy of', dy
       write(0,*) '                  provided  dy of', gd%dy
    end if

    if (present(order)) then
       !if (order < LIN_ORDER) then
       !   write(0,*) 'ERROR in conv_InitGridDef: order < 0&
       !        & is not valid in conv_InitGridDef'
       !   stop
       if (abs(order)+1 > gd%ny) then
          write(0,*) &
            &'ERROR in conv_InitGridDef: |order|+1 > number of grid points' 
       end if
       gd%order = order
    else
       gd%order = LIN_ORDER
    end if
    
  end subroutine conv_InitGridDef_single


  !--------------------------------------------------------------
  ! Create a multi grid def
  subroutine conv_InitGridDef_multi(gd,gdarray,locked)
    use sort
    type(grid_def), intent(out) :: gd
    type(grid_def), intent(in)  :: gdarray(:)
    logical,        intent(in), optional :: locked
    !----------------------------------
    integer :: i,j, indx(size(gdarray))
    logical :: used(size(gdarray))
    real(dp) :: dyratio, approx_ny, new_ymax
    type(grid_def), pointer :: subgd(:) ! shorthand
    ! temp needed for workaround on ifort 8.0.039
    real(dp) :: gdarraydy(size(gdarray))

    !-- enforce one layer only
    if (any(gdarray(:)%nsub /= 0)) then
       write(0,*) 'ERROR in conv_InitGridDef_multi:'
       write(0,*) 'One of grid defs in array was a compound grid def.'
       write(0,*) 'Only one layer of compounding is currently allowed.'
       stop
    end if
    
    gd%nsub = size(gdarray)
    allocate(gd%subiy(gd%nsub+1))
    allocate(gd%subgd(gd%nsub));   subgd => gd%subgd

    gd%locked = default_or_opt(.false.,locked)
    if (gd%locked) then
       ! this calls gives wrong results with ifort-8.0.039
       ! (dummy array gdarray%dy is corrupted in indexx). Issue
       ! submitted 27/01/2004: 225712 
       ! "Corrupted answer on call with array of derived type subcomponents"
       !call indexx(gdarray%dy, indx)

       ! workaround for ifort-8.0.039 
       gdarraydy = gdarray%dy
       call indexx(gdarraydy, indx)
       do i = 1, gd%nsub
          subgd(i) = gdarray(indx(i))
          if (i > 1 .and. subgd(i)%ymax <  subgd(i-1)%ymax) then
             write(0,*) 'ERROR in conv_InitGridDef_multi: for locking,'
             write(0,*) 'gdarray with smaller dy should&
                  & also have smaller gdarray%ymax'
             ! for testing ifort_8_0_039...
             !write(0,*) 'dy   values (i-1,i)', subgd(i-1:i)%dy
             !write(0,*) 'ymax values (i-1,i)', subgd(i-1:i)%ymax
             !write(0,*) indx
             !write(0,*) gdarray%dy
             !write(0,*) gdarray(indx(:))%dy
             !write(0,*) i, subgd(i)%ymax, subgd(i-1)%ymax
             stop
          end if
       end do
       
!!$       used = .false.
!!$       !-- sort into ascending order of dy
!!$       ! worlds least efficient sort algorithm -- but not critical...
!!$       do i = 1, gd%nsub
!!$          !-- intel compiler returns wrong value for minloc!
!!$          !   is that not fun?!
!!$          j = sum(minloc(gdarray%dy,mask=.not.used))
!!$          !write(0,*) j, used, real(gdarray%dy)
!!$          !write(0,*) '------', minval(gdarray%dy,mask=.not.used)
!!$          !used(j) = .true.
!!$          subgd(i) = gdarray(j)
!!$          if (i > 1 .and. subgd(i)%ymax <  subgd(i-1)%ymax) then
!!$             write(0,*) 'ERROR in conv_InitGridDef_multi: for locking,'
!!$             write(0,*) 'smallest gdarray%dy should&
!!$                  & correspond to smallest gdarray%ymax'
!!$             !write(0,*) i, subgd(i)%ymax, subgd(i-1)%ymax
!!$             stop
!!$          end if
!!$       end do
       !-- now ensure that there is proper matching between locked grids
       do i = gd%nsub-1, 1, -1
          ! dyratio must be an integer
          dyratio = subgd(i+1)%dy / subgd(i)%dy
          subgd(i)%dy = subgd(i+1)%dy / nint(dyratio)
          if (abs(dyratio-nint(dyratio)) > warn_tolerance*dyratio) then
             write(0,'(a,i2,a,f18.14)') ' conv_InitGridDef (locking):&
                  & redefined dy(', i,') to be ', subgd(i)%dy
          end if
          ! after fixing dy one must still have an integer number of bins
          approx_ny = subgd(i)%ymax/subgd(i)%dy
          subgd(i)%ny = ceiling(approx_ny - warn_tolerance)
          new_ymax = subgd(i)%ny * subgd(i)%dy
          if (abs(new_ymax-subgd(i)%ymax) > warn_tolerance*new_ymax) then
             write(0,'(a,i2,a,f18.14)') ' conv_InitGridDef (locking):&
                  & redefined ymax(', i,') to be ', new_ymax
          end if
          subgd(i)%ymax = new_ymax
          subgd(i)%ny   = nint(subgd(i)%ymax / subgd(i)%dy)
          ! condition on order must still hold
          if (abs(subgd(i)%order)+1 > subgd(i)%ny) then
             write(0,'(a)') 'Error in conv_InitGridDef (locking):'
             write(0,'(a,i2,a)') '       For grid def ',i,' |order|+1 > ny'
          end if
       end do
    else
       !-- no questions asked!
       subgd(:) = gdarray(:)
    end if
    
    gd%dy       = zero
    gd%ymax     = maxval(subgd(:)%ymax)
    gd%order    = conv_UndefinedInt
    !-- arrays will stretch from 0:gd%ny; so must get this right!
    gd%ny       = sum(subgd(:)%ny) + gd%nsub-1
    
    !-- indicate starting points of arrays
    gd%subiy(1) = 0
    do i = 2, gd%nsub+1
       gd%subiy(i) = gd%subiy(i-1) + subgd(i-1)%ny + 1
    end do
    
  end subroutine conv_InitGridDef_multi
  
  
  !----------------------------------------------------------------
  ! updated for multi; and have run checks to make sure that it
  ! behaves appropriately
  recursive function conv_CmpGridDef(gd1,gd2) result(equal)
    type(grid_def), intent(in) :: gd1, gd2
    logical :: equal
    integer :: i
    logical, parameter :: verbose = .false.

    if (gd1%nsub /= gd2%nsub) then
       equal = .false.
       return
    else if (gd1%nsub == 0) then
       if (verbose) write(0,*) gd1%dy,gd2%dy, gd1%ny, gd2%ny, gd1%ymax,gd2%ymax
       equal = (gd1%dy   == gd2%dy) .and. (gd1%ny == gd2%ny) .and.&
            &  (gd1%ymax == gd2%ymax)
    else
       equal = .true. 
       do i = 1, gd1%nsub
          equal = equal .and. conv_CmpGridDef(gd1%subgd(i),gd2%subgd(i))
       end do
    end if
    
  end function conv_CmpGridDef
  
  !-- useful to be able to check easily -------------------------
  ! no problem with multi
  subroutine conv_ValidateGD(gd1,gd2,source)
    type(grid_def), intent(in) :: gd1, gd2
    character(len=*), intent(in) :: source
    if (.not. (gd1 == gd2)) then
       write(0,*) 'Problem validating two grid defs in ',source
       stop
    end if
  end subroutine conv_ValidateGD


  !======================================================================
  ! Things just for grid quants
  ! multi makes no difference
  subroutine conv_AllocGridQuant_1d(gd,gq)
    type(grid_def), intent(in) :: gd
    real(dp),       pointer    :: gq(:)
    integer :: istat
    ! this form of deallocate ought to be OK (?), but on lahey,
    ! compaq+condor (and perhaps elsewhere) it causes problems
    !deallocate(gq,stat=istat)
    allocate(gq(0:gd%ny))
  end subroutine conv_AllocGridQuant_1d

  !----------------------------------------------------------------------
  ! multi makes no difference
  subroutine conv_AllocGridQuant_2d(gd,gq,nl,nh)
    type(grid_def), intent(in) :: gd
    real(dp),       pointer    :: gq(:,:)
    integer,        intent(in) :: nl,nh
    integer :: istat
    ! this form of deallocate ought to be OK (?), but on lahey,
    ! compaq+condor (and perhaps elsewhere) it causes problems
    !deallocate(gq,stat=istat)
    allocate(gq(0:gd%ny,nl:nh))
  end subroutine conv_AllocGridQuant_2d

  !----------------------------------------------------------------------
  ! multi makes no difference
  subroutine conv_AllocGridQuant_3d(gd,gq, nl2,nh2, nl3,nh3)
    type(grid_def), intent(in) :: gd
    real(dp),       pointer    :: gq(:,:,:)
    integer,        intent(in) :: nl2,nh2,nl3,nh3
    integer :: istat
    ! this form of deallocate ought to be OK (?), but on lahey,
    ! compaq+condor (and perhaps elsewhere) it causes problems
    !deallocate(gq,stat=istat)
    allocate(gq(0:gd%ny,nl2:nh2,nl3:nh3))
  end subroutine conv_AllocGridQuant_3d
  
  !----------------------------------------------------------
  ! multi makes no difference
  subroutine conv_DelGridQuant_1d(gq)
    real(dp),       pointer    :: gq(:)
    integer :: istat
    deallocate(gq,stat=istat)
  end subroutine conv_DelGridQuant_1d
  !----------------------------------------------------------
  ! multi makes no difference
  subroutine conv_DelGridQuant_2d(gq)
    real(dp),       pointer    :: gq(:,:)
    integer :: istat
    deallocate(gq,stat=istat)
  end subroutine conv_DelGridQuant_2d
  

  
  !----------------------------------------------------------------------
  ! updated for multi
  recursive subroutine conv_InitGridQuant_func(gd, gq, func)
    real(dp),         intent(inout) :: gq(0:)
    type(grid_def),   intent(in)    :: gd
    interface
       function func(x)
         use types; implicit none
         real(dp), intent(in) :: x
         real(dp)             :: func
       end function func
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny

    ny = assert_eq(gd%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (gd%nsub /= 0) then
       do isub = 1, gd%nsub
          call conv_InitGridQuant_func(gd%subgd(isub), &
               &gq(gd%subiy(isub):gd%subiy(isub+1)-1), func)
       end do
    else
       do iy = 0, ny
          gq(iy) = func(iy*gd%dy)
       end do
    end if
    
  end subroutine conv_InitGridQuant_func


  !----------------------------------------------------------------------
  ! updated for multi
  recursive subroutine conv_InitGridQuant_func_a(gd, gq, func, axtra)
    real(dp),         intent(inout) :: gq(0:)
    type(grid_def),   intent(in)    :: gd
    real(dp),         intent(in)    :: axtra
    interface
       function func(x,axtra)
         use types; implicit none
         real(dp), intent(in) :: x, axtra
         real(dp)             :: func
       end function func
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny

    ny = assert_eq(gd%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (gd%nsub /= 0) then
       do isub = 1, gd%nsub
          call conv_InitGridQuant_func_a(gd%subgd(isub), &
               &gq(gd%subiy(isub):gd%subiy(isub+1)-1), func, axtra)
       end do
    else
       do iy = 0, ny
          gq(iy) = func(iy*gd%dy, axtra)
       end do
    end if
    
  end subroutine conv_InitGridQuant_func_a


  !----------------------------------------------------------------------
  ! updated for multi
  recursive subroutine conv_InitGridQuant_func_ai(gd, gq, func, axtra, ixtra)
    real(dp),         intent(inout) :: gq(0:)
    type(grid_def),   intent(in)    :: gd
    real(dp),         intent(in)    :: axtra
    integer,          intent(in)    :: ixtra
    interface
       function func(x,axtra,ixtra)
         use types; implicit none
         real(dp), intent(in) :: x, axtra
         integer,  intent(in) :: ixtra
         real(dp)             :: func
       end function func
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny

    ny = assert_eq(gd%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (gd%nsub /= 0) then
       do isub = 1, gd%nsub
          call conv_InitGridQuant_func_ai(gd%subgd(isub), &
               &gq(gd%subiy(isub):gd%subiy(isub+1)-1), func, axtra, ixtra)
       end do
    else
       do iy = 0, ny
          gq(iy) = func(iy*gd%dy, axtra, ixtra)
       end do
    end if
    
  end subroutine conv_InitGridQuant_func_ai


  !----------------------------------------------------------------------
  ! updated for multi
  recursive subroutine conv_InitGridQuant_func2d(gd, gq, func)
    real(dp),         intent(inout) :: gq(0:,:)
    type(grid_def),   intent(in)    :: gd
    interface
       function func(x,n)
         use types; implicit none
         real(dp), intent(in) :: x
         integer , intent(in) :: n
         real(dp)             :: func(n)
       end function func
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny, n

    ny = assert_eq(gd%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (gd%nsub /= 0) then
       do isub = 1, gd%nsub
          call conv_InitGridQuant_func2d(gd%subgd(isub), &
               &gq(gd%subiy(isub):gd%subiy(isub+1)-1,:), func)
       end do
    else
       n  = ubound(gq, dim=2)
       do iy = 0, ny
          gq(iy,:) = func(iy*gd%dy, n)
       end do
    end if
  end subroutine conv_InitGridQuant_func2d

  !----------------------------------------------------------------------
  ! version intended for use when there is an extra argument whose
  ! value is fixed and needs to be passed to func
  !
  ! updated for multi
  recursive subroutine conv_InitGridQuant_func2d_a(gd, gq, func, axtra)
    real(dp),         intent(inout) :: gq(0:,:)
    type(grid_def),   intent(in)    :: gd
    real(dp),         intent(in)    :: axtra
    interface
       function func(x,axtra,n)
         use types; implicit none
         real(dp), intent(in) :: x,axtra
         integer , intent(in) :: n
         real(dp)             :: func(n)
       end function func
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny, n

    ny = assert_eq(gd%ny,ubound(gq,dim=1),"conv_InitGridQuant_func2d_a")
    if (gd%nsub /= 0) then
       do isub = 1, gd%nsub
          call conv_InitGridQuant_func2d_a(gd%subgd(isub), &
               &gq(gd%subiy(isub):gd%subiy(isub+1)-1,:), func, axtra)
       end do
    else
       n  = ubound(gq, dim=2)
       do iy = 0, ny
          gq(iy,:) = func(iy*gd%dy, axtra, n)
       end do
    end if
  end subroutine conv_InitGridQuant_func2d_a

  !----------------------------------------------------------------------
  ! version intended for use when there is an extra argument whose
  ! value is fixed and needs to be passed to func
  !
  ! updated for multi
  recursive subroutine conv_InitGridQuant_func2d_ai(gd, gq, func, axtra, ixtra)
    real(dp),         intent(inout) :: gq(0:,:)
    type(grid_def),   intent(in)    :: gd
    real(dp),         intent(in)    :: axtra
    integer,          intent(in)    :: ixtra
    interface
       function func(x,axtra,ixtra,n)
         use types; implicit none
         real(dp), intent(in) :: x,axtra
         integer , intent(in) :: ixtra,n
         real(dp)             :: func(n)
       end function func
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny, n

    ny = assert_eq(gd%ny,ubound(gq,dim=1),"conv_InitGridQuant_func2d_a")
    if (gd%nsub /= 0) then
       do isub = 1, gd%nsub
          call conv_InitGridQuant_func2d_ai(gd%subgd(isub), &
               &gq(gd%subiy(isub):gd%subiy(isub+1)-1,:), func, axtra, ixtra)
       end do
    else
       n  = ubound(gq, dim=2)
       do iy = 0, ny
          gq(iy,:) = func(iy*gd%dy, axtra, ixtra, n)
       end do
    end if
  end subroutine conv_InitGridQuant_func2d_ai


  !----------------------------------------------------------------------
  !! Version added specially for initialising a PDF from a LHAPDF style 
  !! interface.
  recursive subroutine conv_InitGridQuantLHAPDF(gd, gq, LHAsub, Q)
    real(dp),         intent(inout) :: gq(0:,:)
    type(grid_def),   intent(in)    :: gd
    real(dp),         intent(in)    :: Q
    interface
       subroutine LHAsub(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine LHAsub
    end interface
    !-----------------------------------------
    integer  :: iy, isub, ny
    real(dp) :: f(size(gq,dim=2))

    ny = assert_eq(gd%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (gd%nsub /= 0) then
       do isub = 1, gd%nsub
          call conv_InitGridQuantLHAPDF(gd%subgd(isub), &
               &gq(gd%subiy(isub):gd%subiy(isub+1)-1,:), LHAsub, Q)
       end do
    else
       do iy = 0, ny
          call LHAsub(exp(-iy*gd%dy), Q, f)
          gq(iy,:) = f
       end do
    end if
  end subroutine conv_InitGridQuantLHAPDF

  !----------------------------------------------------------------------
  ! updated for multi
  recursive subroutine conv_InitGridQuantSub_2d(gd, gq, sub)
    real(dp),         intent(inout) :: gq(0:,:)
    type(grid_def),   intent(in)    :: gd
    interface
       subroutine sub(y,res)
         use types; implicit none
         real(dp), intent(in)  :: y
         real(dp), intent(out) :: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny

    ny = assert_eq(gd%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (gd%nsub /= 0) then
       do isub = 1, gd%nsub
          call conv_InitGridQuantSub_2d(gd%subgd(isub), &
               &gq(gd%subiy(isub):gd%subiy(isub+1)-1,:), sub)
       end do
    else
       do iy = 0, ny
          call sub(iy*gd%dy, gq(iy,:))
       end do
    end if
  end subroutine conv_InitGridQuantSub_2d

  recursive subroutine conv_InitGridQuantSub_2d_a(gd, gq, sub, axtra)
    real(dp),         intent(inout) :: gq(0:,:)
    type(grid_def),   intent(in)    :: gd
    real(dp),         intent(in)    :: axtra
    interface
       subroutine sub(y, axtra, res)
         use types; implicit none
         real(dp), intent(in)  :: y, axtra
         real(dp), intent(out) :: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny

    ny = assert_eq(gd%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (gd%nsub /= 0) then
       do isub = 1, gd%nsub
          call conv_InitGridQuantSub_2d_a(gd%subgd(isub), &
               &gq(gd%subiy(isub):gd%subiy(isub+1)-1,:), sub, axtra)
       end do
    else
       do iy = 0, ny
          call sub(iy*gd%dy, axtra, gq(iy,:))
       end do
    end if
  end subroutine conv_InitGridQuantSub_2d_a

  recursive subroutine conv_InitGridQuantSub_2d_ai(gd, gq, sub, axtra, ixtra)
    real(dp),         intent(inout) :: gq(0:,:)
    type(grid_def),   intent(in)    :: gd
    real(dp),         intent(in)    :: axtra
    integer,          intent(in)    :: ixtra
    interface
       subroutine sub(y, axtra, ixtra, res)
         use types; implicit none
         real(dp), intent(in)  :: y, axtra
         integer,  intent(in)  :: ixtra
         real(dp), intent(out) :: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny

    ny = assert_eq(gd%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (gd%nsub /= 0) then
       do isub = 1, gd%nsub
          call conv_InitGridQuantSub_2d_ai(gd%subgd(isub), &
               &gq(gd%subiy(isub):gd%subiy(isub+1)-1,:), sub, axtra, ixtra)
       end do
    else
       do iy = 0, ny
          call sub(iy*gd%dy, axtra, ixtra, gq(iy,:))
       end do
    end if
  end subroutine conv_InitGridQuantSub_2d_ai


  !======================================================================
  ! Variants for buggy pgf90-5.1.2, using workaround suggested by
  ! Giulia
  !----------------------------------------------------------------------
  ! updated for multi
  recursive subroutine conv_InitGridQuantSub_2d_pgf(gd, gq, sub)
    real(dp),         intent(out)   :: gq(0:,:)
    type(grid_def),   intent(in)    :: gd
    interface
       subroutine sub(y,res)
         use types; implicit none
         real(dp), intent(in)  :: y
         real(dp), intent(out) :: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny
    ! -- GZ:make a copy shifting indices to make it compatible with pgf90 
    real(dp) :: gq_cpy(1:size(gq,dim=1),1:size(gq,dim=2))

    ny = assert_eq(gd%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (gd%nsub /= 0) then
       do isub = 1, gd%nsub
          call conv_InitGridQuantSub_2d(gd%subgd(isub), &
               &gq_cpy(gd%subiy(isub)+1:gd%subiy(isub+1),:), sub)
       end do
    else
       do iy = 0, ny
          call sub(iy*gd%dy, gq_cpy(iy+1,:))
       end do
    end if
    do iy=0,ny
       gq(iy,:) = gq_cpy(iy+1,:)
    end do
  end subroutine conv_InitGridQuantSub_2d_pgf

  recursive subroutine conv_InitGridQuantSub_2d_a_pgf(gd, gq, sub, axtra)
    real(dp),         intent(out)   :: gq(0:,:)
    type(grid_def),   intent(in)    :: gd
    real(dp),         intent(in)    :: axtra
    interface
       subroutine sub(y, axtra, res)
         use types; implicit none
         real(dp), intent(in)  :: y, axtra
         real(dp), intent(out) :: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny
    real(dp) :: gq_cpy(1:size(gq,dim=1),1:size(gq,dim=2))

    ny = assert_eq(gd%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (gd%nsub /= 0) then
       do isub = 1, gd%nsub
          call conv_InitGridQuantSub_2d_a(gd%subgd(isub), &
               &gq_cpy(gd%subiy(isub)+1:gd%subiy(isub+1),:), sub, axtra)
       end do
    else
       do iy = 0, ny
          call sub(iy*gd%dy, axtra, gq_cpy(iy+1,:))
       end do
    end if
    do iy=0,ny
       gq(iy,:) = gq_cpy(iy+1,:)
    end do
  end subroutine conv_InitGridQuantSub_2d_a_pgf

  recursive subroutine conv_InitGridQuantSub_2d_ai_pgf(gd, gq, sub, axtra, ixtra)
    real(dp),         intent(out)   :: gq(0:,:)
    type(grid_def),   intent(in)    :: gd
    real(dp),         intent(in)    :: axtra
    integer,          intent(in)    :: ixtra
    interface
       subroutine sub(y, axtra, ixtra, res)
         use types; implicit none
         real(dp), intent(in)  :: y, axtra
         integer,  intent(in)  :: ixtra
         real(dp), intent(out) :: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny
    real(dp) :: gq_cpy(1:size(gq,dim=1),1:size(gq,dim=2))

    ny = assert_eq(gd%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (gd%nsub /= 0) then
       do isub = 1, gd%nsub
          call conv_InitGridQuantSub_2d_ai(gd%subgd(isub), &
               &gq_cpy(gd%subiy(isub)+1:gd%subiy(isub+1),:), sub, axtra, ixtra)
       end do
    else
       do iy = 0, ny
          call sub(iy*gd%dy, axtra, ixtra, gq_cpy(iy+1,:))
       end do
    end if
    do iy=0,ny
       gq(iy,:) = gq_cpy(iy+1,:)
    end do
  end subroutine conv_InitGridQuantSub_2d_ai_pgf


  !--------------------------------------------------------------------
  ! multi-grid done
  recursive function conv_EvalGridQuant_0d(gd, gq, y) result(f)
    use interpolation
    type(grid_def), intent(in) :: gd
    real(dp), intent(in) :: gq(0:)
    real(dp), intent(in) :: y
    real(dp) :: f
    !-----------------------------------------
    integer, parameter :: npnt_min = 4, npnt_max = 10
    integer :: i, ny, npnt, isub
    real(dp) :: ey, df
    real(dp), parameter :: resc_yvals(npnt_max) = (/ (i,i=0,npnt_max-1) /)
    real(dp) :: wgts(npnt_max)

    !write(0,*) y,gd%ny, ubound(gq,dim=1)
    ny = assert_eq(gd%ny,ubound(gq,dim=1),"conv_EvalGridQuant")
    if (y > gd%ymax*(one+warn_tolerance)) then
       write(0,*) 'conv_EvalGridQuant: &
            &requested function value beyond maximum'
       write(0,*) 'y = ', y, 'ymax=',gd%ymax
       stop
    end if
    if (gd%nsub /= 0) then
       isub = conv_BestIsub(gd,y)
       f = conv_EvalGridQuant(gd%subgd(isub), &
            & gq(gd%subiy(isub):gd%subiy(isub+1)-1), y)
    else
       npnt = min(npnt_max, max(npnt_min, abs(gd%order)))
       
       i = min(max(floor(y / gd%dy)-(npnt-1)/2,0),ny-npnt+1)
       call uniform_interpolation_weights(y/gd%dy - i, wgts(1:npnt))
       f = sum(wgts(1:npnt)*gq(i:i+npnt-1))
       !-- this was less efficient...
       !call polint(resc_yvals(1:npnt),gq(i:i+npnt-1),y/gd%dy-i,f,df)
!!$    i = min(gd%ny - 1, floor(y / gd%dy))
!!$    ey = y/gd%dy - i
!!$    f  = (one-ey)*gq(i) + ey*gq(i+1)
    end if
    
  end function conv_EvalGridQuant_0d

  !--------------------------------------------------------------------
  !! 1-D version of grid evaluation.
  !! 
  !! NB: we rewrite everything above so as to avoid unnecessary looping.
  !!     It would be better to have a common routine that gives the
  !!     required set of points and weights?
  recursive function conv_EvalGridQuant_1d(gd, gq, y) result(f)
    use interpolation
    type(grid_def), intent(in) :: gd
    real(dp), intent(in) :: gq(0:,1:)
    real(dp), intent(in) :: y
    real(dp) :: f(size(gq,dim=2))
    !-----------------------------------------
    integer, parameter :: npnt_min = 4, npnt_max = 10
    integer :: i, j, ny, npnt, isub
    real(dp) :: ey, df
    real(dp), parameter :: resc_yvals(npnt_max) = (/ (i,i=0,npnt_max-1) /)
    real(dp) :: wgts(npnt_max)

    !write(0,*) y,gd%ny, ubound(gq,dim=1)
    ny = assert_eq(gd%ny,ubound(gq,dim=1),"conv_EvalGridQuant")
    if (y > gd%ymax*(one+warn_tolerance)) then
       write(0,*) 'conv_EvalGridQuant: &
            &requested function value beyond maximum'
       write(0,*) 'y = ', y, 'ymax=',gd%ymax
       stop
    end if
    if (gd%nsub /= 0) then
       isub = conv_BestIsub(gd,y)
       f = conv_EvalGridQuant_1d(gd%subgd(isub), &
            & gq(gd%subiy(isub):gd%subiy(isub+1)-1,:), y)
    else
       npnt = min(npnt_max, max(npnt_min, abs(gd%order)))
       
       i = min(max(floor(y / gd%dy)-(npnt-1)/2,0),ny-npnt+1)
       call uniform_interpolation_weights(y/gd%dy - i, wgts(1:npnt))
       do j = 1, size(f)
          f(j) = sum(wgts(1:npnt)*gq(i:i+npnt-1,j))
       end do
       
       !-- this was less efficient...
       !call polint(resc_yvals(1:npnt),gq(i:i+npnt-1),y/gd%dy-i,f,df)
!!$    i = min(gd%ny - 1, floor(y / gd%dy))
!!$    ey = y/gd%dy - i
!!$    f  = (one-ey)*gq(i) + ey*gq(i+1)
    end if
  end function conv_EvalGridQuant_1d

!O  !--------------------------------------------------------------------
!O  !! 1-D version of grid evaluation.
!O  !! 
!O  !! NB: this is currently VERY slow because interpolation is recalculated
!O  !!     from scratch for each component of the grid quantity
!O  function conv_EvalfGridQuant_1d(gd, gq, y) result(f)
!O    type(grid_def), intent(in) :: gd
!O    real(dp), intent(in) :: gq(0:,1:)
!O    real(dp), intent(in) :: y
!O    real(dp) :: f(size(gq,dim=2))
!O    integer :: i
!O    do i = 1, size(gq,dim=2)
!O       f(i) = conv_EvalGridQuant_0d(gd,gq(:,i),y)
!O    end do
!O  end function conv_EvalfGridQuant_1d
  

  !--------------------------------------------------------------------
  ! Returns starting iymin point and a set of weights in order to calculate
  ! the value of the function at y -- one day we might introduce some
  ! option of setting the number of points; but not for now...
  !
  ! Qu: is the relation between number of points and order correct? It seems
  !     like we ought to have abs(gd%order)+1...
  recursive subroutine conv_WgtGridQuant(gd, y, iylo, wgts)
    use interpolation
    type(grid_def), intent(in) :: gd
    real(dp), intent(in)  :: y
    integer,  intent(out) :: iylo
    real(dp), pointer     :: wgts(:)
    !-----------------------------------------
    integer, parameter :: npnt_min = 4, npnt_max = 10
    integer :: ny, npnt, isub

    ny = gd%ny
    if (gd%nsub /= 0) then
       isub = conv_BestIsub(gd,y)
       call conv_WgtGridQuant(gd%subgd(isub), y, iylo, wgts)
       iylo = iylo + gd%subiy(isub)
    else
       if (y > gd%ymax*(one+warn_tolerance) .or. y < -warn_tolerance) then
          write(0,*) 'conv_WgtGridQuant: &
               &requested function value outside y range'
          write(0,*) 'y = ', y, 'ymax=',gd%ymax
          stop
       end if
       
       npnt = min(npnt_max, max(npnt_min, abs(gd%order)))
       allocate(wgts(0:npnt-1))
       
       iylo = min(max(floor(y / gd%dy)-(npnt-1)/2,0),ny-npnt+1)
       call uniform_interpolation_weights(y/gd%dy-iylo, wgts)
    end if
  end subroutine conv_WgtGridQuant
  

  

  !-- for internal use only
  function conv_BestIsub(gd,y) result(isub)
    type(grid_def), intent(in) :: gd
    real(dp),       intent(in) :: y
    integer                    :: isub
       !-- this will probably slow things down, but
       !   for the time being accept this penalty
       ! find the grid with the smallest ymax > y
       if (y>gd%ymax) then
          isub = sum(maxloc(gd%subgd%ymax))
       else
          isub = sum(minloc(gd%subgd%ymax,mask=(gd%subgd%ymax>=y)))
       end if
  end function conv_BestIsub
  
  !----------------------------------------------------------------------
  ! not yet updated for multi-grid purposes, because would require
  ! a little bit of thought as to how to treat different grid regions
  function conv_MomGridQuant(gd,gq,omega) result(res)
    type(grid_def), intent(in) :: gd
    real(dp), intent(in) :: gq(0:)
    real(dp), intent(in) :: omega
    real(dp) :: res
    !-----------------------------------------
    real(dp) :: weight, weightprod, dy
    integer :: i, ny

    if (gd%nsub /= 0) then
       write(0,*) 'ERROR in conv_MomGridQuant:&
            & multiple grids not yet supported'
    end if
    
    ny = assert_eq(gd%ny,ubound(gq,dim=1),"conv_MomGridQuant")
    dy = gd%dy
    if (omega == zero) then
       weight = gd%dy
       weightprod = one
       res = half*weight*gq(0)
    else
       weightprod = exp(-dy*omega)
       weight = (exp(-dy*omega) - one + dy*omega)/(dy*omega**2)
       res = weight*gq(0)
       weight = weight + &
            & (exp(+dy*omega) - one - dy*omega)/(dy*omega**2)
    end if
    do i = 1, ny
       weight = weight * weightprod
       if (i == ny) weight = weight * half
       res = res + weight*gq(i)
    end do
  end function conv_MomGridQuant
  

  !----------------------------------------------------------------------
  ! Use a non-sophisticated fix for the multi-grid option, just take the
  ! largest of the dy values if none is specified (so as to avoid information
  ! overload).
  subroutine conv_PrintGridQuant_1(gd,gq,dy,dev)
    type(grid_def), intent(in) :: gd
    real(dp), intent(in) :: gq(0:)
    real(dp), intent(in), optional :: dy
    integer,  intent(in), optional :: dev
    real(dp) :: dy_local, y, x, q
    integer :: ny, i, dev_local
    
    ny = assert_eq(gd%ny,ubound(gq,dim=1),'conv_PrintGridQuant')
    if (gd%nsub /= 0) then
       dy_local = default_or_opt(maxval(gd%subgd%dy),dy)
    else
       dy_local = default_or_opt(gd%dy,dy)
    end if
    dev_local = default_or_opt(6,dev)

    ny = floor(gd%ymax / dy_local)
    do i = 0, ny
       y = i*dy_local
       q = conv_EvalGridQuant(gd, gq, y)
       write(dev_local,*) y, exp(-y),q
    end do
    
  end subroutine conv_PrintGridQuant_1

  !----------------------------------------------------------------------
  ! See conv_PrintGridQuant_1 re multigrid
  subroutine conv_PrintGridQuant_2(gd,gq,gq2,dy,dev)
    type(grid_def), intent(in) :: gd
    real(dp), intent(in) :: gq(0:),gq2(0:)
    real(dp), intent(in), optional :: dy
    integer,  intent(in), optional :: dev
    real(dp) :: dy_local, y, x, q,q2
    integer :: ny, i, dev_local
    
    ny = assert_eq(gd%ny,ubound(gq,dim=1),&
         & ubound(gq2,dim=1),'conv_PrintGridQuant')
    if (gd%nsub /= 0) then
       dy_local = default_or_opt(maxval(gd%subgd%dy),dy)
    else
       dy_local = default_or_opt(gd%dy,dy)
    end if
    dev_local = default_or_opt(6,dev)

    ny = floor(gd%ymax / dy_local)
    do i = 0, ny
       y = i*dy_local
       q  = conv_EvalGridQuant(gd, gq, y)
       q2 = conv_EvalGridQuant(gd, gq2, y)
       write(dev_local,'(25es25.16)') y, exp(-y),q, q2
    end do
    
  end subroutine conv_PrintGridQuant_2

  !----------------------------------------------------------------------
  ! See conv_PrintGridQuant_1 re multigrid
  subroutine conv_PrintGridQuant_3(gd,gq,gq2,gq3,dy,dev)
    type(grid_def), intent(in) :: gd
    real(dp), intent(in) :: gq(0:),gq2(0:),gq3(0:)
    real(dp), intent(in), optional :: dy
    integer,  intent(in), optional :: dev
    real(dp) :: dy_local, y, x, q,q2,q3
    integer :: ny, i, dev_local
    
    ny = assert_eq(gd%ny,ubound(gq,dim=1),&
         & ubound(gq2,dim=1),ubound(gq3,dim=1),&
         & 'conv_PrintGridQuant: distributions must be same size')
    if (gd%nsub /= 0) then
       dy_local = default_or_opt(maxval(gd%subgd%dy),dy)
    else
       dy_local = default_or_opt(gd%dy,dy)
    end if
    dev_local = default_or_opt(6,dev)

    ny = floor(gd%ymax / dy_local)
    do i = 0, ny
       y = i*dy_local
       q  = conv_EvalGridQuant(gd, gq, y)
       q2 = conv_EvalGridQuant(gd, gq2, y)
       q3 = conv_EvalGridQuant(gd, gq3, y)
       write(dev_local,'(25es25.16)') y, exp(-y),q, q2, q3
    end do
    
  end subroutine conv_PrintGridQuant_3
  
  !----------------------------------------------------------------------
  ! See conv_PrintGridQuant_1 re multigrid
  subroutine conv_PrintGridQuant_4(gd,gq,gq2,gq3,gq4,dy,dev)
    type(grid_def), intent(in) :: gd
    real(dp), intent(in) :: gq(0:),gq2(0:),gq3(0:),gq4(0:)
    real(dp), intent(in), optional :: dy
    integer,  intent(in), optional :: dev
    real(dp) :: dy_local, y, x, q,q2,q3, q4
    integer :: ny, i, dev_local
    
    ny = assert_eq(gd%ny,ubound(gq,dim=1),&
         & ubound(gq2,dim=1),ubound(gq3,dim=1),ubound(gq4,dim=1),&
         & 'conv_PrintGridQuant: distributions must be same size')
    if (gd%nsub /= 0) then
       dy_local = default_or_opt(maxval(gd%subgd%dy),dy)
    else
       dy_local = default_or_opt(gd%dy,dy)
    end if
    dev_local = default_or_opt(6,dev)

    ny = floor(gd%ymax / dy_local)
    do i = 0, ny
       y = i*dy_local
       q  = conv_EvalGridQuant(gd, gq, y)
       q2 = conv_EvalGridQuant(gd, gq2, y)
       q3 = conv_EvalGridQuant(gd, gq3, y)
       q4 = conv_EvalGridQuant(gd, gq4, y)
       write(dev_local,'(25es25.16)') y, exp(-y),q, q2, q3, q4
    end do
    
  end subroutine conv_PrintGridQuant_4


  !======================================================================
  ! Routines related to convolutions
  !======================================================================

  !----------------------------------------------------------------------
  !integration precision
  subroutine conv_Seteps(eps_in)
    real(dp), intent(in) :: eps_in
    eps = eps_in
  end subroutine conv_Seteps


  !------------------------------------------------------------------
  ! Just does the memory allocation
  recursive subroutine conv_AllocGridConv_0d(gd,gc)
    type(grid_def),  intent(in)    :: gd
    type(grid_conv), intent(inout) :: gc
    !-------------------------------------------
    integer :: isub
    
    gc%gd = gd
    if (gc%gd%nsub /= 0) then
       nullify(gc%conv) ! avoid doing anything too stupid...
       allocate(gc%subgc(gc%gd%nsub))
       do isub = 1, gc%gd%nsub
          call conv_AllocGridConv_0d(gd%subgd(isub),gc%subgc(isub))
       end do
    else
       nullify(gc%subgc)
       ! this form of deallocate ought to be OK (?), but on lahey,
       ! compaq+condor (and perhaps elsewhere) it causes problems
       !deallocate(gc%conv,stat=istat)
       if (gd%order == LIN_ORDER) then
          allocate(gc%conv(0:gd%ny,2))
       else if (gd%order > LIN_ORDER) then
          allocate(gc%conv(0:gd%ny,0:gd%order+1))
       else
          allocate(gc%conv(0:gd%ny,1))
       end if
    end if
  end subroutine conv_AllocGridConv_0d

  !-----------------------------------------------------------------
  ! Handy to have multi-dimensional versions of the above
  subroutine conv_AllocGridConv_1d(gd,gc)
    type(grid_def),  intent(in)    :: gd
    type(grid_conv), intent(inout) :: gc(:)
    !-------------------------------------------
    integer :: i
    do i = 1, size(gc)
       call conv_AllocGridConv_0d(gd,gc(i))
    end do
  end subroutine conv_AllocGridConv_1d
  subroutine conv_AllocGridConv_2d(gd,gc)
    type(grid_def),  intent(in)    :: gd
    type(grid_conv), intent(inout) :: gc(:,:)
    !-------------------------------------------
    integer :: i,j
    do i = 1, size(gc,dim=2)
       do j = 1, size(gc,dim=1)
          call conv_AllocGridConv_0d(gd,gc(j,i))
       end do
    end do
  end subroutine conv_AllocGridConv_2d


  !--------------------------------------------------------------------
  ! initialise a grid convolution with zero
  ! Default for alloc is .true.; writing this this way allows
  ! in principle the possibility of a more levels of recursion
  ! in the grid def, should one ever want to have the option...
  recursive subroutine conv_InitGridConv_zero(gd,gc,alloc)
    type(grid_def),  intent(in)    :: gd
    type(grid_conv), intent(inout) :: gc
    logical,         intent(in), optional :: alloc
    !-------------------------------------------
    integer :: isub

    if (default_or_opt(.true.,alloc)) then
       !write(0,*) 'it was true', gd%nsub
       call conv_AllocGridConv(gd,gc)
    else
       !write(0,*) 'it was false'
    end if
    if (gd%nsub /= 0) then
       do isub = 1, gd%nsub
          ! remember not to do double allocate
          call conv_InitGridConv_zero(gd%subgd(isub),&
               &gc%subgc(isub),alloc=.false.)
       end do
    else
       gc%conv = zero
    end if
    
  end subroutine conv_InitGridConv_zero

  !-----------------------------------------------------------------
  ! Handy to have multi-dimensional versions of the above
  subroutine conv_InitGridConv_zero_1d(gd,gc,alloc)
    type(grid_def),  intent(in)    :: gd
    type(grid_conv), intent(inout) :: gc(:)
    logical,         intent(in), optional :: alloc
    !-------------------------------------------
    integer :: i
    do i = 1, size(gc)
       call conv_InitGridConv_zero(gd,gc(i),alloc)
    end do
  end subroutine conv_InitGridConv_zero_1d
  subroutine conv_InitGridConv_zero_2d(gd,gc,alloc)
    type(grid_def),  intent(in)    :: gd
    type(grid_conv), intent(inout) :: gc(:,:)
    logical,         intent(in), optional :: alloc
    !-------------------------------------------
    integer :: i,j
    do i = 1, size(gc,dim=2)
       do j = 1, size(gc,dim=1)
          call conv_InitGridConv_zero(gd,gc(j,i),alloc)
       end do
    end do
  end subroutine conv_InitGridConv_zero_2d
  

  !--------------------------------------------------------------------
  ! initialise a grid convolution with another gc, potentially multiplied
  ! by a factor. Qu: is alloc supposed to be used from outside?
  recursive subroutine conv_InitGridConv_gc(gc,gcin,fact,alloc)
    type(grid_conv), intent(inout) :: gc
    type(grid_conv), intent(in)    :: gcin
    real(dp),        intent(in), optional :: fact
    logical,         intent(in), optional :: alloc
    !-------------------------------------------
    integer :: isub

    if (default_or_opt(.true.,alloc)) then
       call conv_AllocGridConv(gcin%gd,gc)
    else
       call conv_ValidateGD(gcin%gd,gc%gd,'conv_InitGridConv_gc')
    end if
    if (gcin%gd%nsub /= 0) then
       do isub = 1, gcin%gd%nsub
          ! remember not to do double allocate
          call conv_InitGridConv_gc(gc%subgc(isub),&
               &gcin%subgc(isub),fact=fact,alloc=.false.)
       end do
    else
       if (present(fact)) then
          gc%conv = gcin%conv * fact
       else
          gc%conv = gcin%conv
       end if
    end if
  end subroutine conv_InitGridConv_gc

  !-----------------------------------------------------------------
  ! Handy to have multi-dimensional versions of the above
  subroutine conv_InitGridConv_gc_1d(gc,gcin,fact,alloc)
    type(grid_conv), intent(inout) :: gc(:)
    type(grid_conv), intent(in)    :: gcin(:)
    real(dp),        intent(in), optional :: fact
    logical,         intent(in), optional :: alloc
    !-------------------------------------------
    integer :: i, ni
    ni = assert_eq(size(gc),size(gcin),'conv_InitGridConv_gc_1d')
    do i = 1, ni
       call conv_InitGridConv_gc(gc(i),gcin(i),fact,alloc)
    end do
  end subroutine conv_InitGridConv_gc_1d
  subroutine conv_InitGridConv_gc_2d(gc,gcin,fact,alloc)
    type(grid_conv), intent(inout) :: gc(:,:)
    type(grid_conv), intent(in)    :: gcin(:,:)
    real(dp),        intent(in), optional :: fact
    logical,         intent(in), optional :: alloc
    !-------------------------------------------
    integer :: i,j,ni,nj
    ni = assert_eq(size(gc,dim=2),size(gcin,dim=2),'conv_InitGridConv_gc_1d')
    nj = assert_eq(size(gc,dim=1),size(gcin,dim=1),'conv_InitGridConv_gc_1d')
    do i = 1, ni
       do j = 1, nj
          call conv_InitGridConv_gc(gc(j,i),gcin(j,i),fact,alloc)
       end do
    end do
  end subroutine conv_InitGridConv_gc_2d


  !----------------------------------------------------------------------
  ! Initialise a convoluter with the function to use in the 
  ! convolution. 
  subroutine conv_InitGridConv_func(gd,gc,func,alloc)
    type(grid_def),  intent(in)    :: gd
    type(grid_conv), intent(inout) :: gc
    interface
       function func(x)
         use types; implicit none
         real(dp), intent(in) :: x
         real(dp)             :: func
       end function func
    end interface
    logical,         intent(in), optional :: alloc
    !-------------------------------------------
    integer :: isub

    if (default_or_opt(.true.,alloc)) call conv_InitGridConv_zero(gd,gc)
    
    call conv_AddGridConv(gc,func)
  end subroutine conv_InitGridConv_func


  !----------------------------------------------------------------------
  ! Initialise a convoluter with the function to use in the 
  ! convolution. 
  subroutine conv_InitGridConv_conv(gc,gca,gcb,alloc)
    type(grid_conv),  intent(inout) :: gc
    type(grid_conv),  intent(in)    :: gca, gcb
    logical,          intent(in), optional :: alloc
    !-------------------------------------------
    integer :: isub

    if (default_or_opt(.true.,alloc)) then
       call conv_AllocGridConv(gca%gd,gc)
    else
       call conv_ValidateGD(gc%gd, gca%gd,&
            & 'conv_InitGridConv_conv: gc and gca')
    end if
    
    call conv_ConvGridConv(gc,gca, gcb)
  end subroutine conv_InitGridConv_conv

  !------------------------------------------------------------------
  ! Zero the contents of the grid convolution
  !
  recursive subroutine conv_ZeroGridConv_0d(gc)
    type(grid_conv), intent(inout) :: gc
    integer :: isub
    if (gc%gd%nsub /= 0) then
       do isub = 1, gc%gd%nsub
          call conv_ZeroGridConv(gc%subgc(isub))
       end do
    else
       gc%conv = zero
    end if
  end subroutine conv_ZeroGridConv_0d
  ! Handy to have multi-dimensional versions of the above
  subroutine conv_ZeroGridConv_1d(gc)
    type(grid_conv), intent(inout) :: gc(:)
    !-------------------------------------------
    integer :: i
    do i = 1, size(gc)
       call conv_ZeroGridConv_0d(gc(i))
    end do
  end subroutine conv_ZeroGridConv_1d
  subroutine conv_ZeroGridConv_2d(gc)
    type(grid_conv), intent(inout) :: gc(:,:)
    !-------------------------------------------
    integer :: i,j
    do i = 1, size(gc,dim=2)
       do j = 1, size(gc,dim=1)
          call conv_ZeroGridConv_0d(gc(j,i))
       end do
    end do
  end subroutine conv_ZeroGridConv_2d
  
  !--------------------------------------------------------
  ! Remove memory associated with a given gc.
  recursive subroutine conv_DelGridConv_0d(gc)
    type(grid_conv), intent(inout) :: gc
    !-------------------------------------------
    integer :: istat, isub
    
    if (gc%gd%nsub /= 0) then
       do isub = 1, gc%gd%nsub
          call conv_DelGridConv(gc%subgc(isub))
       end do
       deallocate(gc%subgc,stat=istat)
    else
       deallocate(gc%conv,stat=istat)
       if (istat /= 0) then
          write(0,*) 'ERROR: problems with deallocation in conv_DelGridConv'
          !write(0,*) one/sqrt(sin(zero))
          stop
       end if
    end if

    if (istat /= 0) then
       write(0,*) 'ERROR: problems with deallocation in conv_DelGridConv'
       !write(0,*) one/sqrt(sin(zero))
       stop
    end if
    
  end subroutine conv_DelGridConv_0d
  !-----------------------------------------------------------------
  ! Handy to have multi-dimensional versions of the above
  subroutine conv_DelGridConv_1d(gc)
    type(grid_conv), intent(inout) :: gc(:)
    !-------------------------------------------
    integer :: i
    do i = 1, size(gc)
       call conv_DelGridConv_0d(gc(i))
    end do
  end subroutine conv_DelGridConv_1d
  subroutine conv_DelGridConv_2d(gc)
    type(grid_conv), intent(inout) :: gc(:,:)
    !-------------------------------------------
    integer :: i,j
    do i = 1, size(gc,dim=2)
       do j = 1, size(gc,dim=1)
          call conv_DelGridConv_0d(gc(j,i))
       end do
    end do
  end subroutine conv_DelGridConv_2d


  !--------------------------------------------------------------------
  ! Multiply the given GridConv by the relevant factor
  recursive subroutine conv_MultGridConv_0d(gc,fact)
    type(grid_conv), intent(inout) :: gc
    real(dp),        intent(in)    :: fact
    !-------------------------------------------
    integer :: istat, isub

    if (gc%gd%nsub /= 0) then
       do isub = 1, gc%gd%nsub
          call conv_MultGridConv(gc%subgc(isub),fact)
       end do
    else
       gc%conv = gc%conv * fact
    end if
    
  end subroutine conv_MultGridConv_0d
  !-----------------------------------------------------------------
  ! Handy to have multi-dimensional versions of the above
  subroutine conv_MultGridConv_1d(gc,fact)
    type(grid_conv), intent(inout) :: gc(:)
    real(dp),        intent(in)    :: fact
    !-------------------------------------------
    integer :: i
    do i = 1, size(gc)
       call conv_MultGridConv_0d(gc(i),fact)
    end do
  end subroutine conv_MultGridConv_1d
  subroutine conv_MultGridConv_2d(gc,fact)
    type(grid_conv), intent(inout) :: gc(:,:)
    real(dp),        intent(in)    :: fact
    !-------------------------------------------
    integer :: i,j
    do i = 1, size(gc,dim=2)
       do j = 1, size(gc,dim=1)
          call conv_MultGridConv_0d(gc(j,i),fact)
       end do
    end do
  end subroutine conv_MultGridConv_2d

  

  !-------------------------------------------------------------
  ! To gc add a function for convolution
  recursive subroutine conv_AddGridConv_func(gc,func)
    use integrator; use convolution_communicator
    type(grid_conv), intent(inout), target :: gc
    interface
       function func(x)
         use types; implicit none
         real(dp), intent(in) :: x
         real(dp)             :: func
       end function func
    end interface
    !----------------------------------------------------
    real(dp), pointer :: dy
    integer,  pointer :: ny
    integer  :: order
    real(dp) :: upper,lower, yl,ym,yh
    integer  :: i,j, k
    !------------------------------------------------------
    real(dp) :: res!, nodes(gc%gd%order+1)
    integer  :: inode_one, il, ih, iy, jy, inode
    real(dp), allocatable :: nodes(:)
    !------------------------------------------------------
    integer :: isub

    if (gc%gd%nsub /= 0) then
       do isub = 1, gc%gd%nsub
          call conv_AddGridConv_func(gc%subgc(isub),func)
       end do
       return
    end if
    
    dy => gc%gd%dy
    ny => gc%gd%ny
    order = gc%gd%order

    !-- this used to be an automatic, but ran into problems 
    !   with multi-grid splitting functions because compound holder
    !   had an undefined value for gc%gd%order
    allocate(nodes(abs(gc%gd%order)+1))

    cc_piece = cc_REAL
    if (gc%gd%order == LIN_ORDER) then
       do i = 1, ny
          ym = i*dy
          yl = ym - dy
          yh = ym + dy
          lower = ig_LinWeight(func,yl,ym,zero,one,eps) ! specify 0 at yl
          upper = ig_LinWeight(func,ym,yh,one,zero,eps) ! specify 0 at yh
          gc%conv(i,FULL) = gc%conv(i,FULL) + lower + upper
          gc%conv(i,UPPR) = gc%conv(i,UPPR) + upper
       end do
       !-- see CCN17-62 ----
       ym = zero
       yh = dy
       cc_piece = cc_REAL
       upper    = ig_LinWeight(func,ym,yh,zero,-one,eps)
       cc_piece = cc_REALVIRT
       upper    = upper + ig_LinWeight(func,ym,yh,one,one,eps) 
       cc_piece = cc_VIRT
       upper    = upper + ig_LinWeight(func,yh,-two*log(eps),one,one,eps)
       cc_piece = cc_DELTA
       upper    = upper + func(zero)
       gc%conv(0,FULL) = gc%conv(0,FULL) + upper
       gc%conv(0,UPPR) = gc%conv(0,UPPR) + upper
    else if (gc%gd%order < LIN_ORDER) then
       ! should be similar to Ratcliffes proposal (and the sort of thing
       ! regularly done in BFKL). NB Not documented in any CCN -- hopefully
       ! straightforward enough that it can be done in one's head?
       order = -gc%gd%order
       do i = 1, ny
          !-- this is the range of interest
          yl = (i-1) * dy
          yh =  i    * dy
          !-- this is range of nodes for that interval
          il = i-1
          ih = il + order
          nodes = (/ (j,j=il,ih) /) * dy
          do iy = il, min(ih,ny)
             res = conv_GCAf_Helper(func,yl,yh,il,iy,nodes)
             gc%conv(iy,1) = gc%conv(iy,1) + res
          end do
       end do
       !write(0,*) 'conv_AddGridConv_func: Negative orders not yet supported'
       !stop
    else
       !-- CCN19 p.4 -------------
       ! index 0:          central points
       ! index 1..order+1: for dealing with last order+1 points
       !-- first do central points. Do an interval at a time
       !   i = 1, means interval from 0 to 1
       !-- the first loop is missing something:: it should also be filling up
       !   some pieces in the second loop, at least in some cases
       do i = 1, ny
          !-- this is the range of interest
          yl = (i-1) * dy
          yh =  i    * dy
          !-- these are the interpolation points to be used
          il = max(0,i - (order+2)/2)
          ih = min(ny, il+order)
          il = ih - order
          !-- fill things up?
          nodes = (/ (j,j=il,ih) /) * dy
          do iy = il, ih
             res = conv_GCAf_Helper(func,yl,yh,il,iy,nodes)
             gc%conv(iy,0) = gc%conv(iy,0) + res
             !-- deal properly with special end array --
             !   try to be clever with limits, but there are no guarantees.
             !   It does however seem to work!
             do jy = max(i+order,iy), iy+order
                if (jy <= ny) then
                   gc%conv(jy,order+1-(jy-iy)) &
                        &= gc%conv(jy,order+1-(jy-iy)) + res
                end if
             end do
          end do
       end do
       !-- now deal with integrations close to end points
       ! k represents the end points
       do k = 1, ny
          if (k <= order) then
             yl = max(zero, (k-order)*dy)
             yh = k*dy
             !-- the interpolation points & yes, it is normal for them 
             !   to be negative
             !   it is related to the way the convolution is done for the case 
             !   i < order, i.e. effectively conv(0:order)*gq(order:0:-1)
             !   (actually conv(i, 1:order+1))
             ih = k
             il = k - order
             !-- fill things up
             nodes = (/ (j,j=il,ih) /) * dy
             do iy = il, ih
                res = conv_GCAf_Helper(func,yl,yh,il,iy,nodes)
                gc%conv(k,iy-il+1) = gc%conv(k,iy-il+1) + res
             end do
             cycle
          end if
          ! now do things region by region
          !-- i represents the region that we will study
          do i = max(1,k-order+1), k
             !-- this is the range of interest
             yl = (i-1) * dy
             yh =  i    * dy
             !-- these are the interpolation points to be used
             !   note extra min(k,...) compared to above formula
             il = max(0,i - (order+2)/2)
             ih = min(k,min(ny, il+order)) !-- could simplify expression
             il = ih - order
             nodes = (/ (j,j=il,ih) /) * dy
             do iy = max(il,k-order), ih
                res = conv_GCAf_Helper(func,yl,yh,il,iy,nodes)
                gc%conv(k,order+1-(k-iy)) = &
                     &gc%conv(k,order+1-(k-iy)) + res
             end do
          end do
       end do
    end if
    
    deallocate(nodes)

  end subroutine conv_AddGridConv_func


  !---------------------------------------------------------------------
  ! look after some repetitive work...
  !
  ! Guess that this function may do the following:
  !   Work out the weight corresponding to the integral, between yl and yh
  !   of the function func(y) mutiplied by a polynomial which is zero at 
  !   all the nodes (which start from il) except that indicated by iy.
  !
  function conv_GCAf_Helper(func,yl,yh,il,iy,nodes) result(res)
    use integrator; use convolution_communicator
    interface
       function func(x)
         use types; implicit none
         real(dp), intent(in) :: x
         real(dp)             :: func
       end function func
    end interface
    real(dp), intent(in) :: yl, yh, nodes(:)
    integer,  intent(in) :: il, iy
    real(dp)             :: res
    integer :: inode
    inode = iy - il + 1
    res = zero
    if (yl == zero .and. iy == 0) then
       cc_piece = cc_VIRT
       res = res + ig_LinWeight(func, yh,-two*log(eps), one,one, eps)
       cc_piece = cc_DELTA
       res = res + func(zero)
       cc_piece = cc_REALVIRT
       res = res + ig_LinWeight(func, yl, yh, one, one, eps)
       cc_piece = cc_REAL
       res = res + ig_PolyWeight(func, yl, yh, nodes, inode, eps,wgtadd=-one)
    else
       cc_piece = cc_REAL
       res = ig_PolyWeight(func, yl, yh, nodes, inode, eps)
    end if
  end function conv_GCAf_Helper
  


  !-------------------------------------------------------------
  ! To gc add another grid convolution.
  !
  ! Perhaps there are some inefficiences here (validation may sometimes
  ! be carried out twice), but for the time being, leave it as is.
  recursive subroutine conv_AddGridConv_gc(gc,gcadd,fact)
    type(grid_conv), intent(inout) :: gc
    type(grid_conv), intent(in)  :: gcadd
    real(dp),        intent(in), optional  :: fact
    !----------------------------------------------
    integer :: isub

    call conv_ValidateGD(gc%gd, gcadd%gd, 'conv_AddGridConv_gc')
    if (gc%gd%nsub /= 0) then
       do isub = 1, gc%gd%nsub
          call conv_AddGridConv_gc(gc%subgc(isub),gcadd%subgc(isub),fact)
       end do
    else
       if (present(fact)) then
          gc%conv = gc%conv + gcadd%conv * fact
       else
          gc%conv = gc%conv + gcadd%conv
       end if
    end if
    
  end subroutine conv_AddGridConv_gc
  !-----------------------------------------------------------------
  ! Handy to have multi-dimensional versions of the above
  subroutine conv_AddGridConv_gc_1d(gc,gcadd,fact)
    type(grid_conv), intent(inout) :: gc(:)
    type(grid_conv), intent(in)    :: gcadd(:)
    real(dp),        intent(in), optional :: fact
    !-------------------------------------------
    integer :: i
    do i = 1, size(gc)
       call conv_AddGridConv_gc(gc(i),gcadd(i),fact)
    end do
  end subroutine conv_AddGridConv_gc_1d
  subroutine conv_AddGridConv_gc_2d(gc,gcadd,fact)
    type(grid_conv), intent(inout) :: gc(:,:)
    type(grid_conv), intent(in)    :: gcadd(:,:)
    real(dp),        intent(in), optional    :: fact
    !-------------------------------------------
    integer :: i,j
    do i = 1, size(gc,dim=2)
       do j = 1, size(gc,dim=1)
          call conv_AddGridConv_gc(gc(j,i),gcadd(j,i),fact)
       end do
    end do
  end subroutine conv_AddGridConv_gc_2d


  !--------------------------------------------------------------
  ! Carry out the convolution of gc on gq
  recursive function conv_ConvGridQuant_scalar(gc,gq) result(gqout)
    type(grid_conv),  intent(in),target    :: gc
    real(dp),         intent(in)    :: gq(0:)
    real(dp)                        :: gqout(0:ubound(gq,dim=1))
    !---------------------------------------------
    integer :: i, ny, j
    integer :: order
    integer :: isub, iy, dy_ratio

    if (gc%gd%nsub /= 0) then
       do isub = 1, gc%gd%nsub
          gqout(gc%gd%subiy(isub):gc%gd%subiy(isub+1)-1) = &
               & conv_ConvGridQuant_scalar(gc%subgc(isub),&
               &        gq(gc%gd%subiy(isub):gc%gd%subiy(isub+1)-1))
       end do
        
       if (gc%gd%locked .and. .not.override_grid_locking) then
          !-- decant information from finer grids into coarser grids
          ! (remember: finest grid has lowest isub)
          do isub = 2, gc%gd%nsub
             ! the ratio should be an exact integer, but use
             ! nint() to avoid the dangers of rounding errors
             dy_ratio = nint(gc%gd%subgd(isub)%dy / gc%gd%subgd(isub-1)%dy)
             do iy = 0, gc%gd%subgd(isub-1)%ny / dy_ratio
                gqout(gc%gd%subiy(isub)+iy) = &
                     &gqout(gc%gd%subiy(isub-1)+iy*dy_ratio)
             end do
          end do
          nconv_with_override_off = nconv_with_override_off + 1
       end if
       return
    end if

    ny = assert_eq(gc%gd%ny,ubound(gq,dim=1),"conv_ConvGridQuant")
    order = gc%gd%order

    !-- Hopefully this will avoid some wasted convolutions?
    if (all(gq == zero)) then
       gqout = zero
       return
    end if
    

    if (order == LIN_ORDER) then
       !-- a test to avoid N^2 operations...
       if (all(gc%conv(:,FULL) == zero)) then
          gqout = zero
          return
       end if
       
       gqout(0) = zero
       do i = 1, ny
          !-- following is legal fortran but it is very slow sometimes
          gqout(i) =  sum(gq(0:i)*gc%conv(i:0:-1,FULL))
          !-- maybe the version that follows will be faster on some compilers
          !   on absoft it definitely is
!!$       gqout(i) = zero
!!$       do j = 0, i
!!$          gqout(i) = gqout(i) + gq(j)*gc%conv(i-j,FULL)
!!$       end do
          !write(0,*) i*gc%gd%dy,gqout(i), gc%conva(1)
          gqout(i) = (gqout(i)-gq(0)*gc%conv(i,UPPR))
       end do
    else if (order < 0) then
       ! Ratcliffes proposal
       !gqout(0) = zero
       ! NB: explicitly do also the i=0 case. Under normal circumstances
       !     it should give zero; however we still let the user enjoy
       !     their folly if they really want to have gq(0)/=0 --- this
       !     is something that comes in handy when generating "derived"
       !     convolutors...
       do i = 0, ny
          gqout(i) = sum(gc%conv(0:i,1)*gq(i:0:-1))
       end do
    else
       gqout = zero
       !-- a test to avoid N^2 operations...
       if (all(gc%conv(:,0) == zero)) return
       !-- current writing is designed to reduce thrashing of cache. 
       !   It is not clear that it actually helps in any way
       !   Commented version theoretically has higher cache thrashing.
       do i = 1, ny
          if (i > order) then
             gqout(i) = sum(gc%conv(0:i-order-1,0)*gq(i:order+1:-1))
          end if
          !-- HIGH CACHE THRASH
          !gqout(i) = gqout(i) + sum(gc%conv(i,1:order+1)*gq(order:0:-1))
       end do
       !-- LOW CACHE THRASH
       do i = 1, order+1
          gqout(:) = gqout(:) + gc%conv(:,i)*gq(order+1-i)
       end do
    end if
    
  end function conv_ConvGridQuant_scalar



  !--------------------------------------------------------------
  ! Carry out the convolution of two convolution functions
  ! It is recommended that if one of them is singular at x=1, then this be 
  ! gca: this helps minimise the discretisation error (in some cases).
  !
  ! NB: subroutine name is not quite consistent with things like
  !     Add or Mult since the operands are completely distinct from the
  !     result.
  recursive subroutine conv_ConvGridConv_0d(gc,gca, gcb, reorder)
    type(grid_conv),  intent(inout) :: gc
    type(grid_conv),  intent(in), target   :: gca, gcb
    logical,          intent(in), optional :: reorder
    !---------------------------------------------
    type(grid_conv), pointer :: gcap, gcbp
    integer :: i, ny, j, order, ix
    !-- these might be better if on the fly, given recursive nature?
    real(dp) :: deltafn(0:gca%gd%ny), res(0:gca%gd%ny)
    integer :: isub

    call conv_ValidateGD(gca%gd, gcb%gd, 'conv_ConvGridConv: gca and gcb')
    call conv_ValidateGD(gc%gd, gca%gd, 'conv_ConvGridConv: gc and gca')

    if (gc%gd%nsub /= 0) then
       do isub = 1, gc%gd%nsub
          call conv_ConvGridConv_0d(gc%subgc(isub),&
               &  gca%subgc(isub), gcb%subgc(isub),reorder)
       end do
       return
    end if

    if (default_or_opt(.true.,reorder)) then
       !-- decide which part of conv matrix to use
       select case (gca%gd%order)
       case(LIN_ORDER)
          ix = FULL
       case(LIN_ORDER+1:)
          ix = 0
       case(:LIN_ORDER-1)
          ix = 1
       end select
       
       !-- use as only criterion of singularity that there should be
       !   opposite signs for the first two grid points? This is not
       !   by any means foolproof, but in basic cases it seems to work.
       !   Put the less singular convolutor to the right.
       if (product(gca%conv(0:1,ix)) >= zero .and. &
            & product(gcb%conv(0:1,ix)) < zero) then
          call conv_ConvGridConv_0d(gc,gcb, gca, reorder=.false.)
          return
       end if
    end if


    if (gca%gd%order == LIN_ORDER) then
       !-- this is just a guess. Now we will see if it is right
       gc%conv(:,FULL) = (gca .conv. gcb%conv(:,FULL)) +&
            &gca%conv(:,UPPR) * gcb%conv(0,FULL)
       
       !-- next 3 lines actually calculate LOWER. Then correct it to upper.
       !   Do the convolution of a deltafn by hand in order to
       !   be more efficient: tested that it gave identical answers
       deltafn = gcb%conv(:,FULL) - gcb%conv(:,UPPR)
       deltafn(0) = zero
       gc%conv(:,UPPR) = gca .conv. deltafn
       !gc%conv(:,UPPR) = gca .conv. (gcb .conv. deltafn)
       gc%conv(:,UPPR) = gc%conv(:,FULL) - gc%conv(:,UPPR)
    else if (gca%gd%order < LIN_ORDER) then
       ! the ratcliffe approach
       do i = 0, gc%gd%ny
          gc%conv(i,1) = sum(gca%conv(0:i,1)*gcb%conv(i:0:-1,1))
       end do
       !write(0,*) 'conv_ConvGridConv: Negative orders not yet supported'
       !stop
    else
       order = gca%gd%order
       !-- let us hope that the following works?
       !   in this first incarnation it will be roughly half the optimum
       !   speed because the inner convolution will be done "stupidly",
       !   using O(N^2) steps rather than the possible O(N)
       do j = 0, order+1
          !-- set up a delta function and its convolution
          deltafn = zero; deltafn(j) = one
          !-- try to avoid too much rubbish in the result?
          res = zero
          res = gca .conv. (gcb .conv. deltafn)
          !-- redistribute the result
          if (j /= order+1) then
             gc%conv(:,order+1-j) = res(:)
          else
             gc%conv(0:gca%gd%ny-(order+1),0) = res(order+1:gca%gd%ny)
             !-- The region conv(gca%gd%ny-order:) remains unfilled. This 
             !   is OK since in practice (current incarnation) it never 
             !   gets used.
          end if
       end do
       !-- now remember to test it!
    end if


  end subroutine conv_ConvGridConv_0d

  !--------------------------------------------------------------------
  ! version for automatically doing convolution of two 2x2 arrays
  ! with array multiplication in the usual fortran sense
  subroutine conv_ConvGridConv_2dx2d(gc,gca, gcb)
    type(grid_conv),  intent(inout) :: gc(:,:)
    type(grid_conv),  intent(in)    :: gca(:,:), gcb(:,:)
    !----------------------------------------------------
    integer :: ir, im, ic, nr, nm, nc
    type(grid_conv) :: conv
    
    nr = assert_eq(size(gc,dim=1),size(gca,dim=1) ,'conv_ConvGridConv_2dx2d')
    nc = assert_eq(size(gc,dim=2),size(gcb,dim=2) ,'conv_ConvGridConv_2dx2d')
    nm = assert_eq(size(gca,dim=2),size(gcb,dim=1),'conv_ConvGridConv_2dx2d')
    
    call conv_ZeroGridConv(gc)
    call conv_AllocGridConv(gc(1,1)%gd,conv)
    do ic = 1, nc
       do ir = 1, nr
          do im = 1, nm
             call conv_ConvGridConv_0d(conv,gca(ir,im),gcb(im,ic))
             call conv_AddGridConv(gc(ir,ic),conv)
          end do
       end do
    end do
    call conv_DelGridConv(conv)
    
  end subroutine conv_ConvGridConv_2dx2d
  
  !--------------------------------------------------------------------------
  ! Returns in gc (assumed preallocated) the commutator of
  ! gca and gcb. Does nothing clever about ordering of
  ! the arrays...
  !
  ! This is a very inefficient version because it does many more convolutions
  ! than are needed. Only for 2x2 arrays does it do a reasonably sensible job.
  !
  subroutine conv_CommGridConv(gc,gca, gcb)
    type(grid_conv),  intent(inout) :: gc(:,:)
    type(grid_conv),  intent(in)    :: gca(:,:), gcb(:,:)
    !----------------------------------------------------
    type(grid_conv) :: cc(size(gc,dim=1),size(gc,dim=2))
    type(grid_conv) :: prod

    call conv_AllocGridConv(gc(1,1)%gd,cc)
    call conv_ZeroGridConv(gc)

    !-- use a shortcut for 2-dim arrays; do things long-hand
    !   for others
    if (size(gc,dim=1) == 2 .and. size(gc,dim=2) == 2) then
    !-- for 2x2 arrays we will want to do something along the following lines
    ! (with a similar definition for B).
!> A:=matrix(2,2,[a11,a12,a21,a22]);
!                                                     [a11    a12]
!                                                A := [          ]
!                                                     [a21    a22]
!----------------------------------------------------------------------- 
!> matadd(multiply(A,B),-multiply(B,A));
!  [         a12 b21 - a21 b12           a11 b12 + a12 b22 - b11 a12 - b12 a22]
!  [                                                                          ]
!  [a21 b11 + a22 b21 - b21 a11 - b22 a21          a21 b12 - a12 b21          ]
       call conv_AllocGridConv(gc(1,1)%gd,prod)
       
       ! res11 = a12 b21 - a21 b12
       call conv_ConvGridConv(prod,gca(1,2),gcb(2,1))
       call conv_AddGridConv(gc(1,1),prod)
       call conv_ConvGridConv(prod,gca(2,1),gcb(1,2))
       call conv_AddGridConv(gc(1,1),prod,-one)
       ! res22 = -res 11
       call conv_AddGridConv(gc(2,2),gc(1,1),-one)
       ! res12 = a11 b12 + a12 b22 - b11 a12 - b12 a22
       call conv_ConvGridConv(prod,gca(1,1),gcb(1,2))
       call conv_AddGridConv(gc(1,2),prod)
       call conv_ConvGridConv(prod,gca(1,2),gcb(2,2))
       call conv_AddGridConv(gc(1,2),prod)
       call conv_ConvGridConv(prod,gca(1,2),gcb(1,1))
       call conv_AddGridConv(gc(1,2),prod,-one)
       call conv_ConvGridConv(prod,gca(2,2),gcb(1,2))
       call conv_AddGridConv(gc(1,2),prod,-one)
       ! res21 = a21 b11 + a22 b21 - b21 a11 - b22 a21
       call conv_ConvGridConv(prod,gca(2,1),gcb(1,1))
       call conv_AddGridConv(gc(2,1),prod)
       call conv_ConvGridConv(prod,gca(2,2),gcb(2,1))
       call conv_AddGridConv(gc(2,1),prod)
       call conv_ConvGridConv(prod,gca(1,1),gcb(2,1))
       call conv_AddGridConv(gc(2,1),prod,-one)
       call conv_ConvGridConv(prod,gca(2,1),gcb(2,2))
       call conv_AddGridConv(gc(2,1),prod,-one)
       
       call conv_DelGridConv(prod)
    else
       call conv_ConvGridConv(cc,gca,gcb)
       call conv_AddGridConv(gc,cc)
       call conv_ConvGridConv(cc,gcb,gca)
       call conv_AddGridConv(gc,cc,-one)
    end if
    
    call conv_DelGridConv(cc)
  end subroutine conv_CommGridConv
  


  !--------------------------------------------------------------
  ! Carry out the convolution of gc on gq: but not always the
  ! most efficient way of dealing with the problem...
  function conv_ConvGridQuant_mat(gc,gq) result(gqout)
    type(grid_conv),  intent(in)    :: gc(:,:)
    real(dp),         intent(in)    :: gq(0:,:)
    real(dp)                        :: gqout(0:ubound(gq,dim=1),size(gc,dim=1))
    !---------------------------------------------
    integer :: i, ny, ic, ir, ncol, nrow

    ny = assert_eq(gc(1,1)%gd%ny,ubound(gq,dim=1),"conv_ConvGridQuant")
    ncol = assert_eq(size(gc,dim=2),size(gq,dim=2),"conv_ConvGridQuant")
    nrow = size(gc,dim=1)

    gqout = zero
    do ir = 1, nrow
       do ic = 1, ncol
          !gqout(:,ir) = gqout(:,ir) + gc(ir,ic) .conv. gq(:,ic)
          gqout(:,ir) = gqout(:,ir) +&
               & conv_ConvGridQuant_scalar(gc(ir,ic), gq(:,ic))
       end do
    end do
  end function conv_ConvGridQuant_mat

  

  !======================================================================
  ! things for easier access to a given grid point
  ! (but may not be as fast as just calling the routine)
  function conv_gdval_gdv(gd,val) result(gdv)
    type(grid_def), intent(in) :: gd
    real(dp),       intent(in) :: val
    type(gdval) :: gdv
    gdv%gd  = gd
    gdv%val = val
  end function conv_gdval_gdv
  function conv_gdval_vgd(val,gd) result(gdv)
    type(grid_def), intent(in) :: gd
    real(dp),       intent(in) :: val
    type(gdval) :: gdv
    gdv%gd  = gd
    gdv%val = val
  end function conv_gdval_vgd
  function conv_EvalGridQuant_atx(gq, gdv) result(res)
    real(dp),    intent(in) :: gq(:)
    type(gdval), intent(in) :: gdv
    real(dp) :: res
    res = conv_EvalGridQuant(gdv%gd, gq, -log(gdv%val))
  end function conv_EvalGridQuant_atx
  function conv_EvalGridQuant_atx_1d(gq, gdv) result(res)
    real(dp),    intent(in) :: gq(:,:)
    type(gdval), intent(in) :: gdv
    real(dp) :: res(size(gq,dim=2))
    res = conv_EvalGridQuant(gdv%gd, gq, -log(gdv%val))
  end function conv_EvalGridQuant_atx_1d
  

  !======================================================================
  ! Routines for getting derived convolution operators
  !
  ! See also notes on CCN25-96
  !======================================================================
  ! Allocates and returns an array of probes, each of which must 
  ! be operated on with the convolver, after which the user should
  ! call conv_SetDerivedConv, which formalises the results
  subroutine conv_GetDerivedProbes(gd,probes)
    type(grid_def), intent(in) :: gd
    real(dp),       pointer    :: probes(:,:)
    !---------------------
    integer :: nprobes

    nprobes = 0
    call conv_GetNProbes(gd,nprobes)
    call conv_AllocGridQuant(gd,probes,1,nprobes)
    call conv_SetProbes(gd,probes)

    ! safety measures -- this is a NASTY business, but cannot
    ! thing of another way to ensure & check that grid locking is
    ! off...
    override_grid_locking = .true.
    nconv_with_override_off = 0
  end subroutine conv_GetDerivedProbes

  !----------------------------------------------------------------------
  ! Sets up the "derived" convolution operator.
  !
  ! This routine actually just does the housekeeping -- the actual "hard"
  ! work is done by conv_SetDerivedConv_rec
  !
  ! NB: gc must previously have been allocated
  subroutine conv_SetDerivedConv(gc,probes)
    real(dp),        pointer       :: probes(:,:)
    type(grid_conv), intent(inout) :: gc

    call conv_SetDerivedConv_nodealloc(gc,probes)
    deallocate(probes)
  end subroutine conv_SetDerivedConv

  !----------------------------------------------------------------------
  ! runs through each of the elements of gd and establishes
  ! the number of probes needed. Currently only supports cases
  ! in which there is just a single probe...
  recursive subroutine conv_GetNProbes(gd,nprobes)
    use warnings_and_errors
    type(grid_def), intent(in)  :: gd
    integer,        intent(out) :: nprobes
    !---------
    integer :: nprobes_tmp, isub
    if (gd%nsub /= 0) then
       nprobes = 0
       do isub = 1, gd%nsub
          call conv_GetNProbes(gd%subgd(isub),nprobes_tmp)
          nprobes = max(nprobes, nprobes_tmp)
       end do
    else
       select case(gd%order)
       case(:LIN_ORDER-1)
          nprobes = 1
       case(LIN_ORDER)
          nprobes = 2
       case(LIN_ORDER+1:)
          nprobes = gd%order+2
       end select
    end if
  end subroutine conv_GetNProbes

  !----------------------------------------------------------------------
  ! Sets the exact form of the probes... Only works for order<0 for now.
  recursive subroutine conv_SetProbes(gd,probes)
    type(grid_def), intent(in)  :: gd
    real(dp),       intent(out) :: probes(0:,:)
    integer :: isub, iprobe
    if (gd%nsub /= 0) then
       do isub = 1, gd%nsub
          call conv_SetProbes(gd%subgd(isub),&
               &probes(gd%subiy(isub):gd%subiy(isub+1)-1,:))
       end do
    else
       select case(gd%order)
       case(:LIN_ORDER-1)
          probes(0,1) = one
          probes(1:,:) = zero
       case(LIN_ORDER)
          probes(:,:) = zero
          probes(0,1) = one
          probes(1,2) = one
       case(LIN_ORDER+1:)
          probes = zero
          do iprobe = 1, gd%order+2
             probes(iprobe-1,iprobe) = one
          end do
       end select
    end if
  end subroutine conv_SetProbes
  
  
  !----------------------------------------------------------------------
  ! Does the work in the setup of the "derived" convolution operator.
  recursive subroutine conv_SetDerivedConv_nodealloc(gc,probes)
    use warnings_and_errors
    real(dp),        intent(in)    :: probes(0:,:)
    type(grid_conv), intent(inout) :: gc
    integer :: isub, order, iprobe

    override_grid_locking = .false.
    if (nconv_with_override_off /= 0) call wae_error(&
         &'conv_SetDerivedConv_nodealloc',&
         &'Detected convolutions while lock overried off')

    if (gc%gd%nsub /= 0) then
       do isub = 1, gc%gd%nsub
          call conv_SetDerivedConv_nodealloc(gc%subgc(isub),&
               &probes(gc%gd%subiy(isub):gc%gd%subiy(isub+1)-1,:))
       end do
    else
       select case(gc%gd%order)
       case(:LIN_ORDER-1)
          gc%conv = probes
       case(LIN_ORDER)
          gc%conv(0:gc%gd%ny-1,FULL) = probes(1:,2)
          ! fake value here -- but since we will actually only be 
          ! interested in lower for this point, it does not
          ! matter!
          gc%conv(gc%gd%ny,FULL) = zero
          ! actually it is lower that has been calculated here
          gc%conv(0:,UPPR) = probes(0:,1)
          ! now convert it to upper
          gc%conv(0:,UPPR) = gc%conv(0:,FULL) - gc%conv(0:,UPPR)
       case(LIN_ORDER+1:)
          order = gc%gd%order
          ! a safety check -- if I forget anything it should
          ! jump out...
          gc%conv = 1e90_dp ! dummy value to detect uninitialized pieces...
          do iprobe = 1, order+1
             gc%conv(:,iprobe) = probes(:,order+2-iprobe)
          end do
          gc%conv(0:gc%gd%ny-order-1,0) = probes(order+1:,order+2)
          ! originally known, but unknowable with this method (and
          ! irrelevant).
          gc%conv(gc%gd%ny-order:,0) = zero
       end select
    end if
  end subroutine conv_SetDerivedConv_nodealloc
  
end module convolution
