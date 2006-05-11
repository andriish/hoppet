! $Id: conv_objects.f90,v 1.21 2004/09/18 14:39:15 salam Exp $




!======================================================================
! 
! Here we adopt a representation suitable for arbitrary order of evolution,
! along the lines of what is outlined by Vogt et al in hep-ph/9907472 and 
! hep-ph/0006154, and CCN19-93
!
! Currently the representation need not be defined 100%. The basic structure
! however is:
!
!             -6..-2    q_minus_NS_ik == (q_i - qbar_i) - (q_k - qbar_k) 
!                 -1    q_V_NS        == sum_i (q_i - qbar_i)
!                  0    g             == gluon
!                  1    Sigma         == sum_i (q_i + qbar_i)
!               2..6    q_plus_NS_ik  == (q_i + qbar_i) - (q_k + qbar_k) 
!
! Exact allocation of the ik is something which is not necessary to know
! here.
!
!======================================================================
module conv_objects_hidden
  use types; use consts_dp; use splitting_functions
  !use qcd
  use dglap_choices
  use warnings_and_errors
  use pdf_representation
  use convolution
  use assertions
  implicit none
  private

  !-----------------------------------------------------------------
  ! holds all the components of the splitting functions
  type PMat
     !private
     !-- These are the singlet evolution pieces
     !   qg is defined including a 2nf factor
     type(grid_conv)           :: singlet(iflv_g:iflv_sigma,iflv_g:iflv_sigma)
     type(grid_conv), pointer  :: gg, qq, gq, qg
     !-- These are the non-singlet pieces
     type(grid_conv)    :: NS_plus, NS_minus, NS_V
     !-- LO -> loops=1, etc...
     integer            :: loops, nf_int
  end type PMat
  public             :: PMat

  integer, parameter :: default_ibase = 1

  !---------------------------------------------------------------
  ! for going from nf to nf+1. Currently contains pieces
  ! required for O(as^2), but not the full general structure.
  ! Initially, support for this might be a bit limited?
  type MassThresholdMat
     type(grid_conv) :: PShq, PShg
     type(grid_conv) :: NSqq_H, Sgg_H, Sgq_H
     ! LOOPS == 1+POWER OF AS2PI, NF_INT = nf including heavy quark.
     ! This is potentially adaptable.
     integer         :: loops, nf_int
  end type MassThresholdMat
  public :: MassThresholdMat

  !-----------------------------------------------------------------
  ! holds the components of a coefficient function
  type CMat
     !private
     !-- need a grid def. Sometimes no g or q will be defined
     !   but still want the coefficient function associated with a grid
     !   def
     !   delta is the magnitude of a delta function piece
     !   quark charge contributions are included by the routines here.
     !   for gluons the implicitness is \sum_{q,qbar} e_q^2
     type(grid_def)     :: gd
     type(grid_conv)    :: g, q
     real(dp)           :: delta
     logical            :: HO
  end type CMat
  public             :: CMat


  !---------- avoid need for main program to access convolution
  !  give it both possible names. Maybe not necessary?
  public :: grid_def
  interface cobj_InitGridDef
     ! perhaps illegal, so do it explicitly?
     !module procedure InitGridDef
     module procedure conv_InitGridDef_single, conv_InitGridDef_multi
  end interface
  public :: cobj_InitGridDef
  public :: InitGridDef

  public :: InitSplitMatLO, InitSplitMatNLO
  public :: InitSplitMatNNLO
  public :: InitSplitMatPolLO, InitSplitMatPolNLO
  public :: cobj_InitSplit, cobj_AddSplit, cobj_DelSplit
  public :: cobj_PConv, cobj_PConv_1d

  public :: InitMTMNNLO, cobj_SetNfMTM, cobj_ConvMTM, cobj_DelMTM

  !-------- things for splitting functions --------------------
  interface cobj_InitCoeff
     module procedure cobj_InitCoeffLO, cobj_InitCoeffHO, cobj_InitCoeff_cf
  end interface
  public :: cobj_InitCoeff, cobj_InitCoeffLO, cobj_InitCoeffHO
  public :: cobj_CConv
  public :: cobj_AddCoeff, cobj_DelCoeff
  public :: cobj_DefaultPrep
  public :: cobj_GetDerivedProbes, cobj_AllocSplit, cobj_SetDerivedSplit

contains

  function cobj_DefaultPrep(nf_lcl) 
    integer, intent(in) :: nf_lcl
    type(pdf_rep)       :: cobj_DefaultPrep
    cobj_DefaultPrep%nf    = nf_lcl
    cobj_DefaultPrep%ibase = default_ibase
  end function cobj_DefaultPrep

  !======================================================================
  ! make sure all required links are set up
  subroutine cobj_InitSplitLinks(P)
    type(PMat), target,  intent(inout) :: P
    
    if ((iflv_sigma - iflv_g) /= 1) call wae_error(&
         &'cobj_InitSplitLinks:', 'local gluon and sigma ids not consistent')
    
    !-- NB singlet matrix is not written with usual convention of 
    !   (sigma g), but rather (g sigma)
    P%gg => P%singlet(iflv_g,iflv_g)
    P%gq => P%singlet(iflv_g,iflv_sigma)
    P%qg => P%singlet(iflv_sigma,iflv_g)
    P%qq => P%singlet(iflv_sigma,iflv_sigma)
  end subroutine cobj_InitSplitLinks
  
  !======================================================================
  subroutine InitSplitMatLO(gd, P)
    use qcd
    type(grid_def), intent(in)    :: gd
    type(PMat),     intent(inout) :: P

    !-- info to describe the splitting function
    P%loops  = 1
    P%nf_int = nf_int

    call cobj_InitSplitLinks(P)

    call InitGridConv(gd, P%gg, sf_Pgg)
    call InitGridConv(gd, P%qq, sf_Pqq)
    call InitGridConv(gd, P%gq, sf_Pgq)
    call InitGridConv(gd, P%qg, sf_Pqg)

    !-- now fix up pieces so that they can be directly used as a matrix
    call Multiply(P%qg, 2*nf)

    !-- PqqV +- PqqbarV
    call InitGridConv(P%NS_plus,  P%qq)
    call InitGridConv(P%NS_minus, P%qq)

    !-- PNSminus + nf * (PqqS - PqqbarS)
    call InitGridConv(P%NS_V, P%NS_minus)
  end subroutine InitSplitMatLO

  !======================================================================
  subroutine InitSplitMatNLO(gd, P, factscheme)
    use qcd
    type(grid_def), intent(in)    :: gd
    type(PMat),     intent(inout) :: P
    integer, optional, intent(in) :: factscheme
    integer :: factscheme_local
    !-----------------------------------------
    type(grid_conv) :: P1qqV, P1qqbarV, P1qqS
    !-- needed for DIS schemes
    type(grid_conv) :: Cq,Cg
    type(grid_conv) :: CqPqq, CqPqg, CqPgq, CqPgg
    type(grid_conv) :: CgPgg, CgPqg, CgPgq, CgPqq


    factscheme_local = default_or_opt(factscheme_default, factscheme)
    if (factscheme_local /= factscheme_MSbar) then
       ! NB: do not support DIS scheme here since it involves 
       !     determination of Pqq etc at LO (already done once, so do 
       !     not want to repeat) -- rather do this stuff in
       !     holders, where Pqq(LO) will in any case be available.
       !     (Is this "chickening out"?)
       write(0,*) 'InitSplitMatNLO: unsupported fact scheme', factscheme
       call wae_error('InitSplitMatNLO: stopping')
    end if
    
    !-- info to describe the splitting function
    P%loops = 2
    P%nf_int = nf_int
    
    call cobj_InitSplitLinks(P)

    !-- these are the building blocks
    write(0,*) associated(P1qqV%conv), associated(P1qqV%subgc)
    call InitGridConv(gd, P1qqV, sf_P1qqV)
    write(0,*) associated(P1qqV%conv), associated(P1qqV%subgc)
    write(0,*) associated(P1qqbarV%conv), associated(P1qqbarV%subgc)
    call InitGridConv(gd, P1qqbarV, sf_P1qqbarV)
    write(0,*) associated(P1qqbarV%conv), associated(P1qqbarV%subgc)
    call Delete(P1qqV%subgc); write(0,*) 'deleted with success';stop
    call InitGridConv(gd, P1qqS, sf_P1qqS)

    !-- PqqV + PqqbarV
    call InitGridConv(P%NS_plus, P1qqV)
    call AddWithCoeff(P%NS_plus, P1qqbarV, one)
    !-- PqqV - PqqbarV
    call InitGridConv(P%NS_minus, P1qqV)
    call AddWithCoeff(P%NS_minus, P1qqbarV, -one)

    !-- PNSminus + nf * (PqqS - PqqbarS) 
    !   [NB at NLO, PqqS = PqqbarS] 
    call InitGridConv(P%NS_V, P%NS_minus)
    
    !-- Pqq in matrix:  PNS_plus + nf*(PqqS + PqqbarS)
    !   [NB at NLO, PqqS = PqqbarS] 
    call InitGridConv(P%qq, P%NS_plus)
    call AddWithCoeff(P%qq, P1qqS, two*nf)

    !-- rest of singlet matrix
    call InitGridConv(gd, P%gq, sf_P1gq)
    call InitGridConv(gd, P%gg, sf_P1gg)
    call InitGridConv(gd, P%qg, sf_P1qg)
    !-- recall that the way it is defined it needs a factor 2nf
    call Multiply(P%qg, two*nf)

    !-- tidy up 
    call Delete(P1qqV)
    write(0,*) 'Got here'; stop
    call Delete(P1qqbarV)
    call Delete(P1qqS)

  end subroutine InitSplitMatNLO



  !======================================================================
  subroutine InitSplitMatNNLO(gd, P, factscheme)
    use qcd; use convolution_communicator
    type(grid_def), intent(in)    :: gd
    type(PMat),     intent(inout) :: P
    integer, optional, intent(in) :: factscheme
    integer :: factscheme_local
    !-----------------------------------------
    type(grid_conv) :: P2NSS
    real(dp) :: dummy

    !call wae_error('InitSplitMatNNLO: NNLO not yet implemented')
    factscheme_local = default_or_opt(factscheme_default, factscheme)
    if (factscheme_local /= factscheme_MSbar) then
       write(0,*) 'InitSplitMatNNLO: unsupported fact scheme', factscheme
       call wae_error('InitSplitMatNNLO: stopping')
    end if
    
    !-- info to describe the splitting function
    P%loops = 3
    P%nf_int = nf_int

    ! NO LONGER NECESSARY
!!$    !-- dummy to initialize Vogt routines (needed in exact cases 
!!$    !   for the A3 piece to be set up). Do it better later on if it works?
!!$    cc_piece = cc_real
!!$    dummy = sf_P2NSMinus(0.5_dp)
!!$    dummy = sf_P2gg(0.5_dp)

    call cobj_InitSplitLinks(P)

    call InitGridConv(gd, P%NS_plus, sf_P2NSPlus)
    call InitGridConv(gd, P%NS_minus, sf_P2NSMinus)
    
    !-- if understanding of convention is right then P_V = P_- + P_S
    call InitGridConv(P%NS_V, P%NS_minus)
    call InitGridConv(gd, P2NSS, sf_P2NSS)
    call AddWithCoeff(P%NS_V, P2NSS)    
    call Delete(P2NSS)

    !-- now the singlet functions
    call InitGridConv(gd, P%qg, sf_P2qg2nf)
    call InitGridConv(gd, P%gg, sf_P2gg)
    call InitGridConv(gd, P%gq, sf_P2gq)
    !-- qq = "pure-singlet" + P+
    call InitGridConv(gd, P%qq, sf_P2PS)
    call AddWithCoeff(P%qq, P%NS_plus)
  end subroutine InitSplitMatNNLO
  

  !======================================================================
  subroutine InitSplitMatPolLO(gd, P)
    use qcd
    type(grid_def), intent(in)    :: gd
    type(PMat),     intent(inout) :: P

    !-- info to describe the splitting function
    P%loops  = 1
    P%nf_int = nf_int

    call cobj_InitSplitLinks(P)

    call InitGridConv(gd, P%gg, sf_DPgg)
    call InitGridConv(gd, P%qq, sf_DPqq)
    call InitGridConv(gd, P%gq, sf_DPgq)
    call InitGridConv(gd, P%qg, sf_DPqg)

    !-- now fix up pieces so that they can be directly used as a matrix
    call Multiply(P%qg, 2*nf)

    !-- PqqV +- PqqbarV
    call InitGridConv(P%NS_plus,  P%qq)
    call InitGridConv(P%NS_minus, P%qq)

    !-- PNSminus + nf * (PqqS - PqqbarS)
    call InitGridConv(P%NS_V, P%NS_minus)
  end subroutine InitSplitMatPolLO

  !======================================================================
  subroutine InitSplitMatPolNLO(gd, P, factscheme)
    use qcd
    type(grid_def), intent(in)    :: gd
    type(PMat),     intent(inout) :: P
    integer, optional, intent(in) :: factscheme
    integer :: factscheme_local
    !-----------------------------------------
    type(grid_conv) :: DP1qqV, DP1qqbarV, DP1qqS

    factscheme_local = default_or_opt(factscheme_Poldefault, factscheme)
    if (factscheme_local /= factscheme_PolMSbar) then
       ! NB: do not support DIS scheme here since it involves 
       !     determination of Pqq etc at LO (already done once, so do 
       !     not want to repeat) -- rather do this stuff in
       !     holders, where Pqq(LO) will in any case be available.
       !     (Is this "chickening out"?)
       write(0,*) 'InitSplitMatPolNLO: unsupported fact scheme', factscheme
       call wae_error('InitSplitMatPolNLO: stopping')
    end if
    
    !-- info to describe the splitting function
    P%loops = 2
    P%nf_int = nf_int
    
    call cobj_InitSplitLinks(P)

    !-- these are the building blocks
    call InitGridConv(gd, DP1qqV, sf_DP1qqV)
    call InitGridConv(gd, DP1qqbarV, sf_DP1qqbarV)
    call InitGridConv(gd, DP1qqS, sf_DP1qqS)

    !-- PqqV + PqqbarV
    call InitGridConv(P%NS_plus, DP1qqV)
    call AddWithCoeff(P%NS_plus, DP1qqbarV, one)
    !-- PqqV - PqqbarV
    call InitGridConv(P%NS_minus, DP1qqV)
    call AddWithCoeff(P%NS_minus, DP1qqbarV, -one)

    !-- PNSminus + nf * (PqqS - PqqbarS) 
    !   [NB at NLO, PqqS = PqqbarS] 
    call InitGridConv(P%NS_V, P%NS_minus)
    
    !-- Pqq in matrix:  PNS_plus + nf*(PqqS + PqqbarS)
    !   [NB at NLO, PqqS = PqqbarS] 
    call InitGridConv(P%qq, P%NS_plus)
    call AddWithCoeff(P%qq, DP1qqS, two*nf)

    !-- rest of singlet matrix
    call InitGridConv(gd, P%gq, sf_DP1gq)
    call InitGridConv(gd, P%gg, sf_DP1gg)
    call InitGridConv(gd, P%qg, sf_DP1qg)
    !-- recall that the way it is defined it needs a factor 2nf
    call Multiply(P%qg, two*nf)

    !-- tidy up 
    call Delete(DP1qqV)
    call Delete(DP1qqbarV)
    call Delete(DP1qqS)

  end subroutine InitSplitMatPolNLO

  !----------------------------------------------------------------------
  ! init a splitting function set with another one (potentially 
  ! multiplied by some factor)
  subroutine cobj_InitSplit(P, Pin, fact)
    type(PMat),     intent(inout) :: P
    type(PMat),     intent(in)    :: Pin
    real(dp),       intent(in), optional :: fact
    !if (nf_d /= (nf_int+1)/2) write(0,*) 'WARNING: non-standard nf_d'

    P%loops  = Pin%loops
    P%nf_int = Pin%nf_int
    call cobj_InitSplitLinks(P)

    call InitGridConv(P%gg,    Pin%gg,     fact)
    call InitGridConv(P%qq,    Pin%qq,     fact)
    call InitGridConv(P%gq,    Pin%gq,     fact)
    call InitGridConv(P%qg,    Pin%qg,     fact)

    !-- from here on, one could exploit info in loops so as to 
    !   reduce number of multiplications... But I am too lazy in this first
    !   go
    call InitGridConv(P%NS_plus,  Pin%NS_plus,  fact)
    call InitGridConv(P%NS_minus, Pin%NS_minus, fact)
    call InitGridConv(P%NS_V,     Pin%NS_V,     fact)
    
  end subroutine cobj_InitSplit

  !----------------------------------------------------------------------
  ! init a splitting function set with another one (potentially 
  ! multiplied by some factor)
  subroutine cobj_AddSplit(P, Padd, fact)
    type(PMat),     intent(inout) :: P
    type(PMat),     intent(in)    :: Padd
    real(dp),       intent(in), optional :: fact
    !if (nf_d /= (nf_int+1)/2) write(0,*) 'WARNING: non-standard nf_d'

    P%loops = max(P%loops, Padd%loops)
    P%nf_int = assert_eq(P%nf_int, Padd%nf_int, &
         &'cobj_AddSplit: nf must be the same')
    call AddWithCoeff(P%gg,       Padd%gg,       fact)
    call AddWithCoeff(P%qq,       Padd%qq,       fact)
    call AddWithCoeff(P%gq,       Padd%gq,       fact)
    call AddWithCoeff(P%qg,       Padd%qg,       fact)
    call AddWithCoeff(P%NS_plus,  Padd%NS_plus,  fact)
    call AddWithCoeff(P%NS_minus, Padd%NS_minus, fact)
    call AddWithCoeff(P%NS_V,     Padd%NS_V,     fact)

  end subroutine cobj_AddSplit

  !----------------------------------------------------------------------
  ! Clear a splitting function (dealloc memory)
  subroutine cobj_DelSplit(P)
    type(PMat),     intent(inout) :: P

    call Delete(P%gg); nullify(P%gg)
    call Delete(P%qq); nullify(P%qq)
    call Delete(P%gq); nullify(P%gq)
    call Delete(P%qg); nullify(P%qg)
    call Delete(P%NS_plus)
    call Delete(P%NS_minus)
    call Delete(P%NS_V)
  end subroutine cobj_DelSplit

  !----------------------------------------------------------------------
  ! This will be overloaded...
  function cobj_PConv(P,q_in) result(Pxq)
    type(PMat), intent(in)         :: P
    real(dp),   intent(in), target :: q_in(0:,ncompmin:)
    real(dp)                     :: Pxq(0:ubound(q_in,dim=1),ncompmin:ncompmax)
    !-- for when we have to change rep
    integer                        :: pdfr
    type(pdf_rep)                  :: prep
    real(dp), allocatable, target  :: q_ev(:,:)
    real(dp), pointer              :: q(:,:)
    integer :: i

    if (ncomponents < P%nf_int) then
       call wae_error('cobj_Pconv:',&
            &'ncomponents in representation is < nf in P.')
    end if

    pdfr = pdfr_GetRep(q_in)
    if (pdfr == pdfr_Evln) then
       q => q_in
    else
       allocate(q_ev(0:ubound(q_in,dim=1),ncompmin:ncompmax))
       prep = cobj_DefaultPrep(p%nf_int)
       call pdfr_HumanToEvln(prep, q_in, q_ev)
       q => q_ev
    end if
    
    PxQ(:, iflv_V)     = P%NS_V * q(:,iflv_V)
    !PxQ(:, iflv_g)     = P%gg * q(:,iflv_g) + P%gq * q(:,iflv_sigma)
    !PxQ(:, iflv_sigma) = P%qg * q(:,iflv_g) + P%qq * q(:,iflv_sigma)
    PxQ(:,iflv_g:iflv_sigma) = P%singlet .conv. q(:,iflv_g:iflv_sigma)

    !-- we need nf-1 NS (+ and -) pieces
    do i = 2, P%nf_int
       PxQ(:,+i) = P%NS_plus  * q(:,+i)
       PxQ(:,-i) = P%NS_minus * q(:,-i)
    end do

    !-- everything else should be set to zero
    do i = P%nf_int + 1, ncomponents
       PxQ(:,+i) = zero
       PxQ(:,-i) = zero
    end do
!!$    !---- TMP TMP ------------------
!!$    write(0,*) '-------------------------------------'
!!$    write(0,*) all(q(:,-1)==zero), all(Pxq(:,-1)==zero)
!!$    write(0,*) real(Pxq(:,-1))

    call pdfr_LabelRep(Pxq,pdfr_Evln)
    if (pdfr == pdfr_Human) then
       q_ev = Pxq ! avoid an extra temporary... (real cleverness might
                  ! have been to avoid this copy...?)
       call pdfr_EvlnToHuman(prep, q_ev, Pxq) ! does labelling as well
       deallocate(q_ev)
    end if
    
  end function cobj_PConv

  !-----------------------------------------------------------------
  ! In some situations the 1d overloaded version is useful.
  ! Do just for PConv, but think also about doing it for others
  ! later on...
  function cobj_PConv_1d(P,q_in) result(Pxq)
    type(PMat), intent(in)         :: P
    real(dp),   intent(in), target :: q_in(0:,ncompmin:,:)
    real(dp)                       :: &
         &Pxq(0:ubound(q_in,dim=1),ncompmin:ncompmax,size(q_in,dim=3))
    integer :: i
    do i = 1, size(q_in,dim=3)
       Pxq(:,:,i) = cobj_PConv(P,q_in(:,:,i))
    end do
  end function cobj_PConv_1d
  

  
  !-----------------------------------------------------------------
  !! Returns the set of probes needed to establish a matrix of
  !! "derived" effective splitting functions.
  subroutine cobj_GetDerivedProbes(gd,probes)
    type(grid_def), intent(in) :: gd
    real(dp),       pointer    :: probes(:,:,:)
    !-----------
    real(dp), pointer :: probes_1d(:,:)
    integer :: nprobes, nprobes_1d, iprobe

    call GetDerivedProbes(gd,probes_1d)

    nprobes_1d = ubound(probes_1d,2)
    nprobes    = 2*nprobes_1d

    ! need more probes because 
    allocate(probes(0:ubound(probes_1d,1),ncompmin:ncompmax, 1:nprobes))

    probes = zero
    do iprobe = 1, nprobes
       ! make sure representation is correct...
       call pdfr_LabelRep(probes(:,:,iprobe),pdfr_Evln)
    end do

    probes(:,iflv_V, 1:nprobes_1d) = probes_1d       ! NS_V
    probes(:,   2, 1:nprobes_1d) = probes_1d       ! NS+
    probes(:,  -2, 1:nprobes_1d) = probes_1d       ! NS-
    probes(:,iflv_sigma, 1:nprobes_1d) = probes_1d   ! Quark column of singlet

    ! gluon column of singlet
    probes(:,iflv_g, nprobes_1d+1:nprobes) = probes_1d       

    ! normally this would be part of housekeeping job of convolution
    ! module (SetDerivedConv), but since probes_1d is lost from 
    ! view at this point, we have to do it ourselves...
    deallocate(probes_1d)
  end subroutine cobj_GetDerivedProbes


  !------------------------------------------------------------
  !! Allocate space for a set of splitting functions, labelling
  !! also with the appropriate nf and nloops...
  subroutine cobj_AllocSplit(gd,P,nf_in, nloops)
    type(grid_def),   intent(in)   :: gd
    type(Pmat),       intent(out)  :: P
    integer,          intent(in)   :: nf_in, nloops
    
    P%nf_int = nf_in
    P%loops  = nloops
    call cobj_InitSplitLinks(P)
    call AllocGridConv(gd,P%gg      )
    call AllocGridConv(gd,P%qq      )
    call AllocGridConv(gd,P%gq      )
    call AllocGridConv(gd,P%qg      )
    call AllocGridConv(gd,P%NS_plus )
    call AllocGridConv(gd,P%NS_minus)
    call AllocGridConv(gd,P%NS_V    )
  end subroutine cobj_AllocSplit
  
  !---------------------------------------------------------------------
  !! Given an allocated Pmat and the results of operating on the probes
  !! determine the resulting "derived" splitting function.
  subroutine cobj_SetDerivedSplit(P,probes)
    type(Pmat),      intent(inout) :: P
    real(dp),        pointer       :: probes(:,:,:)
    !-----------------------------------------
    integer :: nprobes_1d,nprobes, il, ih

    nprobes   = size(probes,dim=3)
    nprobes_1d = nprobes/2

    il = 1; ih = nprobes_1d
    call SetDerivedConv_nodealloc(P%NS_V,     probes(:,iflv_V,il:ih))
    call SetDerivedConv_nodealloc(P%NS_plus,  probes(:,2,il:ih))
    call SetDerivedConv_nodealloc(P%NS_minus, probes(:,-2,il:ih))
    call SetDerivedConv_nodealloc(P%gq,       probes(:,iflv_g,il:ih))
    call SetDerivedConv_nodealloc(P%qq,       probes(:,iflv_sigma,il:ih))
    
    il = nprobes_1d+1; ih = nprobes
    call SetDerivedConv_nodealloc(P%gg,       probes(:,iflv_g,il:ih))
    call SetDerivedConv_nodealloc(P%qg,       probes(:,iflv_sigma,il:ih))

    deallocate(probes)
  end subroutine cobj_SetDerivedSplit
  
  !======================================================================
  !         MASS THRESHOLDS
  !======================================================================
  ! Here we keep all things to do with mass thresholds

  subroutine InitMTMNNLO(gd,MTM)
    type(grid_def),         intent(in)  :: gd
    type(MassThresholdMat), intent(out) :: MTM
    !logical, parameter :: vogt_A2PShg = .false.
    !logical, parameter :: vogt_A2PShg = .true.

    call InitGridConv(gd, MTM%PSHq, sf_A2PShq)
    select case (nnlo_nfthreshold_variant)
    case(nnlo_nfthreshold_param)
       call InitGridConv(gd, MTM%PSHg, sf_A2PShg_vogt)
       write(0,*) 'WARNING in InitMTMNNLO:&
            & using less accuracte vogt parametrisation for A2PShg'
    case(nnlo_nfthreshold_exact)
       call InitGridConv(gd, MTM%PSHg, sf_A2PShg)
    case default
       call wae_error('InitMTMNNLO', 'Unknown nnlo_threshold_variant',&
            &intval=nnlo_nfthreshold_variant)
    end select
   
    call InitGridConv(gd, MTM%NSqq_H, sf_A2NSqq_H)
    call InitGridConv(gd, MTM%Sgg_H, sf_A2Sgg_H)
    call InitGridConv(gd, MTM%Sgq_H, sf_A2Sgq_H)
    ! just store info that it is NNLO. For now it is obvious, but
    ! one day when we have NNNLO it may be useful, to indicate
    ! which structures exist and which do not. 
    ! (Mind you inexistent structures have yet to be implemented here...)
    MTM%loops = 3
    !-- no default value
    MTM%nf_int = 0
  end subroutine InitMTMNNLO

  !---------------------------------------------------------------------
  ! want to be able to set nf, defined as number of flavours
  ! after heavy matching
  subroutine cobj_SetNfMTM(MTM,nf_lcl)
    type(MassThresholdMat), intent(inout) :: MTM
    integer,                intent(in)    :: nf_lcl
    if (MTM%loops > 3) then
       call wae_Error('cobj_SetNfMTM: MTM had loops > 3; nf is probably fixed')
    end if
    MTM%nf_int = nf_lcl
  end subroutine cobj_SetNfMTM
  

  !----------------------------------------------------------------------
  ! Returns the amount to be added to go from nf-1 to nf flavours.
  ! Will try some tests to make sure that nf is consistent?
  !
  ! PDFs are assumed to be in "HUMAN" REPRESENTATION.
  function cobj_ConvMTM(MTM,q) result(Pxq)
    type(MassThresholdMat), intent(in) :: MTM
    real(dp),   intent(in) :: q(0:,ncompmin:)
    real(dp)               :: Pxq(0:ubound(q,dim=1),ncompmin:ncompmax)
    real(dp) :: singlet(0:ubound(q,dim=1))
    integer :: i, nf_light, nf_heavy

    !-- general sanity checks
    if (MTM%loops <=0 .or. MTM%nf_int <=0) call wae_error('cobj_ConvMTM:',&
         &'Mass threshold matrix is undefined')

    if (pdfr_GetRep(q) /= pdfr_Human) call wae_error('cobj_ConvMTM',&
         &'q is not in Human representation')

    nf_heavy = MTM%nf_int
    nf_light = nf_heavy - 1
    !write(0,*) 'Doing a MT convolution with nf_heavy =',nf_heavy

    if (ncomponents < nf_heavy) call wae_error('cobj_ConvMTM:',&
            &'ncomponents in representation is < nf in MTM.')
    if (iflv_g /= 0) call wae_Error('cobj_ConvMTM:','iflv_g/=0')
    
    !-- sanity check on nf
    if (any(q(:,-nf_heavy)/=zero) .or. any(q(:,nf_heavy)/=zero)) &
         & call wae_error('cobj_ConvMTM:',&
         &'Distribution already has non-zero components at nf_heavy')
    
    singlet = sum(q(:,-nf_light:-1),dim=2) + sum(q(:,1:nf_light),dim=2)

    Pxq(:,nf_heavy) = half*(&
         &(MTM%PShq .conv. singlet) + (MTM%PShg .conv. q(:,iflv_g)) )
    Pxq(:,-nf_heavy) = Pxq(:,nf_heavy)
    Pxq(:,iflv_g)     = (MTM%Sgq_H.conv. singlet) + (MTM%Sgg_H.conv. q(:,iflv_g))
    
    do i = -ncomponents, ncomponents
       if (abs(i) > nf_heavy) then
          Pxq(:,i) = zero
       else if (i == iflv_g .or. abs(i) == nf_heavy) then
          cycle
       else
          Pxq(:,i) = MTM%NSqq_H .conv. q(:,i)
       end if
    end do

    call pdfr_LabelRep(Pxq,pdfr_Human)
  end function cobj_ConvMTM
  

  !---- usual cleaning up --------------------------------------------
  subroutine cobj_DelMTM(MTM)
    type(MassThresholdMat), intent(inout) :: MTM
    call Delete(MTM%PSHq)
    call Delete(MTM%PSHg)
    call Delete(MTM%NSqq_H)
    call Delete(MTM%Sgg_H)
    call Delete(MTM%Sgq_H)
    MTM%loops = -1
    MTM%nf_int = -1
  end subroutine cobj_DelMTM


  !======================================================================
  ! From here onwards we have things to do with coefficient functions
  subroutine cobj_InitCoeffLO(gd,C, fact)
    type(grid_def), intent(in)    :: gd
    type(Cmat),     intent(inout) :: C
    real(dp),       intent(in), optional :: fact

    C%gd = gd
    C%HO = .false.
    if (present(fact)) then
       C%delta = fact
    else
       C%delta = one
    end if
  end subroutine cobj_InitCoeffLO
  
  !-- HO. coeff function is not quite so simple. Note that we do simply 
  !   call it NLO, because it could also be relevant at O(as^2).
  !   NOTE UNUSUAL ORDER...
  subroutine cobj_InitCoeffHO(gd, C, coeff_g, coeff_q)
    use convolution_communicator
    type(grid_def), intent(in)    :: gd
    type(Cmat),     intent(inout) :: C
    interface
       function coeff_g(x)
         use types; implicit none
         real(dp), intent(in) :: x
         real(dp)             :: coeff_g
       end function coeff_g
    end interface
    interface
       function coeff_q(x)
         use types; implicit none
         real(dp), intent(in) :: x
         real(dp)             :: coeff_q
       end function coeff_q
    end interface
    real(dp) :: sanity1, sanity2

    C%gd = gd
    C%HO = .true.
    C%delta = zero
    call InitGridConv(gd,C%g, coeff_g)
    call InitGridConv(gd,C%q, coeff_q)
    !-- a sanity check: quite often the arguments for the
    !   quark and gluon pieces get exchanged.
    !   Cannot always ensure correctness (e.g. FL), but can 
    !   detect some cases of misuse
    cc_piece = cc_VIRT
    sanity1 = coeff_g(one)
    cc_piece = cc_DELTA
    sanity2 = coeff_g(zero)
    if (sanity1 /= zero .or. sanity2 /= zero) then
       write(0,*) 'WARNING in cobj_InitCoeffHO **********************'
       write(0,*) 'gluon coefficient function has virtual corrections'
       write(0,*) 'this could be a sign that quark and gluon cf fns have&
            & been exchanged'
       write(0,*) '**************************************************'
    end if
  end subroutine cobj_InitCoeffHO

  !-- initialise C as being fact*Cin ----------------------
  subroutine cobj_InitCoeff_cf(C, Cin, fact)
    type(Cmat),     intent(inout) :: C
    type(Cmat),     intent(in)    :: Cin
    real(dp),       intent(in), optional :: fact

    C%gd = Cin%gd
    C%HO = Cin%HO
    C%delta = Cin%delta
    if (present(fact)) C%delta = C%delta * fact
    if (C%HO) then
       call InitGridConv(C%q, Cin%q, fact)
       call InitGridConv(C%g, Cin%g, fact)
    end if
  end subroutine cobj_InitCoeff_cf

  
  !-- initialise C as being fact*Cin ----------------------
  subroutine cobj_AddCoeff(C, Cin, fact)
    type(Cmat),     intent(inout) :: C
    type(Cmat),     intent(in)    :: Cin
    real(dp),       intent(in), optional :: fact

    call ValidateGD(C%gd,Cin%gd,'cobj_AddCoeff')
    if (present(fact)) then 
       C%delta = C%delta  + Cin%delta * fact
    else
       C%delta = C%delta  + Cin%delta
    end if
    
    if (Cin%HO) then
       if (C%HO) then
          call AddWithCoeff(C%q, Cin%q, fact)
          call AddWithCoeff(C%g, Cin%g, fact)
       else
          call InitGridConv(C%q, Cin%q, fact)
          call InitGridConv(C%g, Cin%g, fact)
       end if
    end if
    C%HO = C%HO .or. Cin%HO
  end subroutine cobj_AddCoeff

  !----------------------------------------------------------------------
  ! Clear a coefficient function (dealloc memory)
  subroutine cobj_DelCoeff(C)
    type(CMat),     intent(inout) :: C

    if (C%HO) then
       call Delete(C%g)
       call Delete(C%q)
    end if
  end subroutine cobj_DelCoeff



  !-------------------------------------------------------------------
  function cobj_CConv(C,q) result(Cxq)
    type(CMat), intent(in) :: C
    real(dp),   intent(in) :: q(0:,0:)
    real(dp)               :: Cxq(0:ubound(q,dim=1))
    !--------------------------------------------------
    real(dp) :: ud(0:ubound(q,dim=1))
    real(dp), parameter :: uq2 = (4.0_dp/9.0_dp)
    real(dp), parameter :: dq2 = (1.0_dp/9.0_dp)

    call wae_error('cobj_CConv:', &
         &'coefficient functions not yet suppoprted with new representation')
    Cxq = zero
    !REPLACE ud = q(:,iflv_u)*uq2 + q(:,iflv_d)*dq2
    !REPLACE if (C%HO) then
    !REPLACE    !-- CHECK: factor of two in front of nf_u and nf_d is there both 
    !REPLACE    !   in ESW and the ADS paper.
    !REPLACE    Cxq = (two*(nf_u*uq2 + nf_d*dq2)) * (C%g .conv. q(:,iflv_g))&
    !REPLACE         & + (C%q .conv. ud)
    !REPLACE    if (C%delta /= zero) Cxq = Cxq + C%delta * ud
    !REPLACE else
    !REPLACE    Cxq = C%delta * ud
    !REPLACE end if
  end function cobj_CConv
    
end module conv_objects_hidden


!======================================================================
! In the above module it was not possible to overload * and .conv. for
! ev_PConv. The reason for this is not immediately clear, but it may be
! connected to the fact that these already existed (and had been included
! through a use specification) and the public statement did not make clear
! whether it referred to the included operators or the newly defined ones. 
!
! A solution is to not use evolution_hidden, but evolution, which includes
! everything that is public in evolution_hidden, but not those things that
! were accessible only there (notably operators included from
! convolution). Then here we define the new additional operators.
!
! NOTE: exceptionally this module is not private.
module conv_objects
  use conv_objects_hidden
  implicit none


  interface operator(*)
     module procedure cobj_PConv, cobj_PConv_1d, cobj_CConv, cobj_ConvMTM
  end interface
  public :: operator(*)
  interface operator(.conv.)
     module procedure cobj_PConv, cobj_PConv_1d, cobj_CConv, cobj_ConvMTM
  end interface
  public  :: operator(.conv.)
  public  :: cobj_Eval2LConv
  !-- avoid access through here
  private :: cobj_PConv, cobj_PConv_1d, cobj_CConv


contains
  !-------------------------------------------------------------------
  ! Calculates [(C2 - y^2/(1+(1-y)^2) CL) .conv. pdf] .atx. xbj
  ! put it here to have access to .conv. for coefficient functions
  function cobj_Eval2LConv(C2, CL, pdf, xbj, ybj) result(res)
    use types; use consts_dp; use convolution
    type(CMat), intent(in)  :: C2, CL
    real(dp),   intent(in)  :: pdf(0:,0:), xbj, ybj
    real(dp)                :: res
    !----------------------------------------------
    type(gdval) :: gdx
    type(Cmat)  :: C2L

    call cobj_InitCoeff(C2L, C2)
    call cobj_AddCoeff (C2L, CL, - ybj**2/(1+(1-ybj)**2) )
    gdx = xbj .with. C2%gd
    res = EvalGridQuant(gdx%gd, (C2L .conv. pdf), -log(gdx%val))
    !res = ((C2L.conv.pdf) .atx. gdx ) 
    call cobj_DelCoeff (C2L)
  end function cobj_Eval2LConv

end module conv_objects


