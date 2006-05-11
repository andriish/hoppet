!-----------------------------------------------------------
!! Routines related to allocation and initialisation of parton 
!! distributions.
module pdf_general
  use types; use consts_dp
  !-- do some recycling (actually only the 2d versions will
  !   be of any use, but hopefully that will not cause problems)
  use convolution
  use pdf_representation
  implicit none
  private

  interface pdfgen_AllocPDF
     module procedure pdfgen_AllocPDF_0d, pdfgen_AllocPDF_1d
  end interface
  public :: pdfgen_AllocPDF


  interface pdfgen_InitPDF
     module procedure pdfgen_InitPDF_, pdfgen_InitPDF_a,&
          & pdfgen_InitPDF_ai
  end interface
  public :: pdfgen_InitPDF

  interface pdfgen_InitPDFSub
     module procedure pdfgen_InitPDFSub_, pdfgen_InitPDFSub_a,&
          & pdfgen_InitPDFSub_ai
  end interface
  public :: pdfgen_InitPDFSub

  public :: pdfgen_InitPDF_LHAPDF

  interface pdfgen_AllocInitPDF
     module procedure pdfgen_AllocInitPDF_, pdfgen_AllocInitPDF_a,&
          & pdfgen_AllocInitPDF_ai
  end interface
  public :: pdfgen_AllocInitPDF

  interface pdfgen_AllocInitPDFSub
     module procedure pdfgen_AllocInitPDFSub_, &
          & pdfgen_AllocInitPDFSub_a, pdfgen_AllocInitPDFSub_ai
  end interface
  public :: pdfgen_AllocInitPDFSub

  interface operator(.anti.)
     module procedure pdfgen_anti_0d, pdfgen_anti_1d
  end interface
  public :: operator(.anti.)
contains
  !-------------------------------------------------------
  ! allocate a parton distribution; dimensionality refers
  ! to extra directions over and above (x,nf).
  subroutine pdfgen_AllocPDF_0d(gd,q)
    type(grid_def), intent(in) :: gd
    real(dp),       pointer    :: q(:,:)

    call conv_AllocGridQuant(gd,q,ncompmin, ncompmax)
  end subroutine pdfgen_AllocPDF_0d

  !-------------------------------------------------------
  subroutine pdfgen_AllocPDF_1d(gd,q, nl, nh)
    type(grid_def), intent(in) :: gd
    real(dp),       pointer    :: q(:,:,:)
    integer ,       intent(in) :: nl, nh
    
    call conv_AllocGridQuant(gd,q,ncompmin, ncompmax, nl, nh)
  end subroutine pdfgen_AllocPDF_1d



  !----------------------------------------------------------------------
  ! All of these routine just redirect to the official routine, 
  ! with the added features that they redirect only the components 
  ! id_min:id_max, and that they label the result as being in the
  ! "Human" representation
  recursive subroutine pdfgen_InitPDF_(gd, gq, func)
    real(dp),         intent(inout) :: gq(0:,ncompmin:)
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
    call conv_InitGridQuant(gd, gq(:,id_min:id_max), func)
    call pdfr_LabelRep(gq,pdfr_Human)
  end subroutine pdfgen_InitPDF_

  !----------------------------------------------------------------------
  ! version intended for use when there is an extra argument whose
  ! value is fixed and needs to be passed to func
  !
  ! updated for multi
  recursive subroutine pdfgen_InitPDF_a(gd, gq, func, axtra)
    real(dp),         intent(inout) :: gq(0:,ncompmin:)
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
    call conv_InitGridQuant(gd, gq(:,id_min:id_max), func, axtra)
    call pdfr_LabelRep(gq,pdfr_Human)
  end subroutine pdfgen_InitPDF_a

  !----------------------------------------------------------------------
  ! version intended for use when there is an extra argument whose
  ! value is fixed and needs to be passed to func
  !
  ! updated for multi
  recursive subroutine pdfgen_InitPDF_ai(gd, gq, func, axtra, ixtra)
    real(dp),         intent(inout) :: gq(0:,ncompmin:)
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
    call conv_InitGridQuant(gd, gq(:,id_min:id_max), func, axtra, ixtra)
    call pdfr_LabelRep(gq,pdfr_Human)
  end subroutine pdfgen_InitPDF_ai


  !----------------------------------------------------------------------
  ! updated for multi
  recursive subroutine pdfgen_InitPDFSub_(gd, gq, sub)
    real(dp),         intent(inout) :: gq(0:,ncompmin:)
    type(grid_def),   intent(in)    :: gd
    interface
       subroutine sub(y,res)
         use types; implicit none
         real(dp), intent(in)  :: y
         real(dp), intent(out) :: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------
    call conv_InitGridQuantSub(gd, gq(:,id_min:id_max), sub)
    call pdfr_LabelRep(gq,pdfr_Human)
  end subroutine pdfgen_InitPDFSub_

  recursive subroutine pdfgen_InitPDFSub_a(gd, gq, sub, axtra)
    real(dp),         intent(inout) :: gq(0:,ncompmin:)
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
    call conv_InitGridQuantSub(gd, gq(:,id_min:id_max), sub, axtra)
    call pdfr_LabelRep(gq,pdfr_Human)
  end subroutine pdfgen_InitPDFSub_a

  recursive subroutine pdfgen_InitPDFSub_ai(gd, gq, sub, axtra, ixtra)
    real(dp),         intent(inout) :: gq(0:,ncompmin:)
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
    call conv_InitGridQuantSub(gd, gq(:,id_min:id_max), sub, axtra, ixtra)
    call pdfr_LabelRep(gq,pdfr_Human)
  end subroutine pdfgen_InitPDFSub_ai



  !======================================================================
  !! Initialise the subroutine using an LHAPDF style subroutine
  subroutine pdfgen_InitPDF_LHAPDF(gd, gq, LHAsub, Q)
    real(dp),         intent(inout) :: gq(0:,ncompmin:)
    type(grid_def),   intent(in)    :: gd
    real(dp),         intent(in)    :: Q
    interface
       subroutine LHAsub(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine LHAsub
    end interface
    !-------------------------------------------
    call conv_InitGridQuantLHAPDF(gd, gq(:,id_min:id_max), LHAsub, Q)
    call pdfr_LabelRep(gq,pdfr_Human)
  end subroutine pdfgen_InitPDF_LHAPDF
  



  !------------------------------------------------------
  ! allocate and initialise in various forms!
  subroutine pdfgen_AllocInitPDF_(gd,q, func)
    type(grid_def), intent(in) :: gd
    real(dp),       pointer    :: q(:,:)
    interface
       function func(y,n)
         use types; implicit none
         real(dp), intent(in) :: y
         integer , intent(in) :: n
         real(dp)             :: func(n)
       end function func
    end interface
    !-----------------------------------------

    call pdfgen_AllocPDF(gd,q)
    call pdfgen_InitPDF(gd, q, func)
  end subroutine pdfgen_AllocInitPDF_

  !--------------------------------------------------
  subroutine pdfgen_AllocInitPDF_a(gd,q, func, axtra)
    type(grid_def), intent(in) :: gd
    real(dp),       pointer    :: q(:,:)
    real(dp),       intent(in) :: axtra
    interface
       function func(y,axtra,n)
         use types; implicit none
         real(dp), intent(in) :: y,axtra
         integer , intent(in) :: n
         real(dp)             :: func(n)
       end function func
    end interface
    !-----------------------------------------

    call pdfgen_AllocPDF(gd,q)
    call pdfgen_InitPDF(gd, q, func, axtra)
  end subroutine pdfgen_AllocInitPDF_a

  !-------------------------------------------------------------
  subroutine pdfgen_AllocInitPDF_ai(gd,q, func, axtra, ixtra)
    type(grid_def), intent(in) :: gd
    real(dp),       pointer    :: q(:,:)
    real(dp),       intent(in) :: axtra
    integer,        intent(in) :: ixtra
    interface
       function func(y,axtra,ixtra,n)
         use types; implicit none
         real(dp), intent(in) :: y,axtra
         integer , intent(in) :: ixtra,n
         real(dp)             :: func(n)
       end function func
    end interface
    !-----------------------------------------

    call pdfgen_AllocPDF(gd,q)
    call pdfgen_InitPDF(gd, q, func, axtra, ixtra)
  end subroutine pdfgen_AllocInitPDF_ai

  !------------------------------------------------------
  ! allocate and initialise in various forms!
  subroutine pdfgen_AllocInitPDFSub_(gd,q, sub)
    type(grid_def), intent(in) :: gd
    real(dp),       pointer    :: q(:,:)
    interface
       subroutine sub(y,res)
         use types; implicit none
         real(dp), intent(in) :: y
         real(dp), intent(out):: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------

    call pdfgen_AllocPDF(gd,q)
    call pdfgen_InitPDFSub(gd, q, sub)
  end subroutine pdfgen_AllocInitPDFSub_

  !--------------------------------------------------
  subroutine pdfgen_AllocInitPDFSub_a(gd,q, sub, axtra)
    type(grid_def), intent(in) :: gd
    real(dp),       pointer    :: q(:,:)
    real(dp),       intent(in) :: axtra
    interface
       subroutine sub(y,axtra,res)
         use types; implicit none
         real(dp), intent(in) :: y,axtra
         real(dp), intent(out):: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------

    call pdfgen_AllocPDF(gd,q)
    call pdfgen_InitPDFSub(gd, q, sub, axtra)
  end subroutine pdfgen_AllocInitPDFSub_a

  !-------------------------------------------------------------
  subroutine pdfgen_AllocInitPDFSub_ai(gd,q, sub, axtra, ixtra)
    type(grid_def), intent(in) :: gd
    real(dp),       pointer    :: q(:,:)
    real(dp),       intent(in) :: axtra
    integer,        intent(in) :: ixtra
    interface
       subroutine sub(y,axtra,ixtra,res)
         use types; implicit none
         real(dp), intent(in) :: y,axtra
         integer , intent(in) :: ixtra
         real(dp), intent(out):: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------

    call pdfgen_AllocPDF(gd,q)
    call pdfgen_InitPDFSub(gd, q, sub, axtra, ixtra)
  end subroutine pdfgen_AllocInitPDFSub_ai

  !---------------------------------------------------
  ! write it out in full rather than using f90 array
  ! notation shortcuts to reduce need for assumptions about
  ! relation between id_min and ncompmin
  function pdfgen_anti_0d(q) result(antiq)
    real(dp), intent(in) :: q(:,ncompmin:)
    real(dp)             :: antiq(size(q,1),lbound(q,2):ubound(q,2))
    integer :: i

    !-- try to catch stupid and illegal uses...
    if (ubound(q,2) /= ncompmax) call wae_error('pdfgen_anti_0d',&
         &'ubound of second dimension of q should be ncompmax, instead it is',&
         &intval=ubound(q,2))

    do i = ncompmin, ncompmax
       if (i >= id_min .and. i <= id_max) then
          antiq(:,i) = q(:,-i)
       else
          antiq(:,i) = q(:,i)
       end if
    end do
  end function pdfgen_anti_0d
  
  !---------------------------------------------------
  function pdfgen_anti_1d(q) result(antiq)
    real(dp), intent(in) :: q(:,ncompmin:,:)
    real(dp)             :: antiq(size(q,1),lbound(q,2):ubound(q,2),size(q,3))
    integer :: i
    do i = 1, size(q,dim=3)
       antiq(:,:,i) = pdfgen_anti_0d(q(:,:,i))
    end do
  end function pdfgen_anti_1d
  
end module pdf_general
