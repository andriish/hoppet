! $Id: evolution.f90,v 1.23 2004/09/29 17:27:13 salam Exp $

module evolution
  use types; use consts_dp
  use conv_objects; use as_server
  use assertions; use qcd
  use warnings_and_errors
  implicit none
  private

  public :: ev_MSBar2DIS
  !public :: ev_evolveLO, ev_evolveNLO, ev_evolveNNLO
  public :: ev_evolve, ev_evolve_varnf, ev_evolve_varnf_gen
  public :: ev_Setdt
  
  !!
  !! A type that allows one to store the result of the evolution as an
  !! operator that can be "applied" to any parton distribution, eliminating 
  !! the need to repeat the whole Runge-Kutta evolution for each new PDF
  !!
  type evln_operator
     type(Pmat)              :: P
     type(MassThresholdMat)  :: MTM ! assume we have just one of these...
     real(dp)                :: MTM_coeff, Q_init, Q_end
     logical                 :: cross_mass_threshold
     type(evln_operator), pointer :: next
  end type evln_operator
  public :: evln_operator, ev_DelEvOp
  !!
  !! related operation and operator
  !!
  public :: ev_conv_evop
  interface operator(*)
     module procedure ev_conv_evop
  end interface
  public :: operator(*)
  interface operator(.conv.)
     module procedure ev_conv_evop
  end interface
  public  :: operator(.conv.)

  !!
  !! fact that alpha_s/2pi is small means that we can quite comfortably
  !! take large steps. Have tested with 0.4 down to 1.5 GeV
  !!   
  !! Originally steps were uniform in dlnQ^2 (dt). As of September 2004,
  !! steps are now uniform in d(alpha_s * lnQ^2) = du. They are taken
  !! "matched" for alpha_s = 0.25
  !!
  real(dp), parameter :: du_dt = 0.25_dp, dt_ballpark_def = 0.4_dp
  real(dp) :: dt_ballpark = dt_ballpark_def
  real(dp) :: du_ballpark = du_dt * dt_ballpark_def
  !real(dp), parameter :: dt_ballpark = 0.2

  type(PMat),      pointer :: ev_PLO, ev_PNLO, ev_PNNLO
  type(as_handle), pointer :: ev_ash
  integer                  :: ev_nloop
  real(dp) :: ev_muR_Q, fourpibeta0_lnmuR_Q
  logical  :: ev_untie_nf = .false.

  ! this was introduced for testing --- should now never be
  ! necessary to set it to .true.
  logical,  parameter :: ev_force_old_dt = .false.
  real(dp), parameter :: ev_minalphadiff     = 0.02_dp
  integer,  parameter :: ev_du_is_dt         = 1
  integer,  parameter :: ev_du_is_dtas_fixed = 2
  integer,  parameter :: ev_du_is_dtas_run   = 3
  integer             :: ev_du_type
  real(dp)            :: ev_u_asref, ev_u_tref, ev_u_bval, ev_du_dt

contains

  !-----------------------------------------------------
  ! At least give the user the option of setting this...
  subroutine ev_Setdt(dt)
    real(dp), intent(in) :: dt
    dt_ballpark = dt
    du_ballpark = dt * du_dt
  end subroutine ev_Setdt
  
  !-------------------------------------------------------------------
  ! This now serves as a simplified entry point to the general routine
  ! which carries out both evolution of a given pdf and the calculation
  ! of the evln_operator.
  subroutine ev_evolve_varnf(sh, pdfset, nloop, ash, Q_init, Q_end, muR_Q,&
       & untie_nf)
    use holders; use pdf_representation
    type(sigma_holder), intent(in), target :: sh
    real(dp),        intent(inout)         :: pdfset(:,:)
    integer,         intent(in)            :: nloop
    type(as_handle), intent(in), target    :: ash
    real(dp),        intent(in)            :: Q_init, Q_end
    real(dp),        intent(in), optional  :: muR_Q
    logical,         intent(in), optional  :: untie_nf
    call ev_evolve_varnf_gen(sh, nloop, ash, Q_init, Q_end, pdfset=pdfset,  &
         &muR_Q = muR_Q, untie_nf = untie_nf)
  end subroutine ev_evolve_varnf
  
  !----------------------------------------------------------------
  ! Give intent(inout) to sh because we will be playing games
  ! with holder_SetNf so as to vary the value of nf.
  !
  ! untie_nf option allows alpha_s to have its natural value 
  ! for nf, even if range of available nf in splitting functions
  ! does not allow this...
  subroutine ev_evolve_varnf_gen(sh, nloop, ash, Q_init, Q_end, &
       & pdfset, evop, muR_Q, untie_nf)
    use holders; use pdf_representation
    type(sigma_holder), intent(in), target   :: sh
    integer,         intent(in)              :: nloop
    type(as_handle), intent(in),    target   :: ash
    real(dp),        intent(in)              :: Q_init, Q_end
    real(dp),        intent(inout), optional :: pdfset(:,:)
    type(evln_operator),intent(inout),  optional, target :: evop
    real(dp),        intent(in),    optional :: muR_Q
    logical,         intent(in),    optional :: untie_nf
    !-----------------------------------------------------
    !real(dp)           :: pdf_ev(size(pdfset,dim=1), size(pdfset,dim=2))
    type(sigma_holder) :: shcopy
    integer            :: nfstore, nf_init, nf_end, nflcl, nfuse
    integer            :: shnf_init, shnf_end, direction, i
    type(evln_operator), pointer :: this_evop
    real(dp),          pointer :: probes(:,:,:)
    !-- ranges
    real(dp) :: QRlo, QRhi, lcl_Q_init, lcl_Q_end
    
    !-- this will just copy pointers, so that we don't
    !   alter info held in original
    shcopy = sh
    !   We want to be able to restore things at the end...
    nfstore = nf_int


    call ev_SetModuleConsts(nloop, muR_Q, untie_nf)

    call as_nfAtQ(ash,Q_init,nf_init, muM_mQ=one)
    call as_nfAtQ(ash,Q_end,nf_end, muM_mQ=one)

    !write(0,*) real(Q_init), real(Q_end), nf_init, nf_end
    
    !-- could use sign, but cannot remember subtleties of signed
    !   zeroes in f95, and do not have manual at hand...
    if (Q_end >= Q_init) then
       direction = 1
    else
       direction = -1
    end if

    call ev_limit_nf(shcopy, nf_init, shnf_init,'initial nf')
    call ev_limit_nf(shcopy, nf_end , shnf_end ,'final nf')

    ! in loop we may want to refer to subsequent elements in chain, 
    ! so make use of a pointer which is kept updated
    if (present(evop)) then
       this_evop => evop
    else
       nullify(this_evop)
    end if

    do nflcl = shnf_init, shnf_end, direction
       call as_QRangeAtNf(ash,nflcl,QRlo,QRhi, muM_mQ=one)
       if (direction == 1) then
          lcl_Q_init = max(Q_init, QRlo)
          lcl_Q_end  = min(Q_end,  QRhi)
       else
          lcl_Q_init = min(Q_init, QRhi)
          lcl_Q_end  = max(Q_end,  QRlo)
       end if
       ! make sure that end points are correct...
       if (nflcl == shnf_init) lcl_Q_init = Q_init
       if (nflcl == shnf_end)  lcl_Q_end  = Q_end
       !-- this will also set the global (qcd) nf_int which will
       !   be used elsewhere in this module
       call holder_SetNf(shcopy, nflcl)

       ! convention: cross mass thresholds before a step in the evolution
       if (nflcl /= shnf_init) then
          ! act directly
          if (present(pdfset))&
               & call ev_CrossMassThreshold(shcopy,ash,direction,pdfset)
          ! or store the information that will be needed to
          ! act subsequently.
          if (present(evop)) &
               & call ev_CrossMassThreshold(shcopy,ash,direction,evop=this_evop)
       else if (present(evop)) then
          this_evop%cross_mass_threshold = .false.
       end if
       
       ! do evolution
       if (present(pdfset)) call ev_evolve(shcopy, &
            &pdfset, nloop, ash, lcl_Q_init, lcl_Q_end, muR_Q)

       ! do the fake evolutions that allow us to do an accelerated 
       ! evolution later on.
       if (present(evop)) then
          ! recall: memory management of probes is done automatically
          call cobj_GetDerivedProbes(sh%gd,probes)
          do i = 1, size(probes,3)
             call ev_evolve(shcopy, &
                  &probes(:,:,i), nloop, ash, lcl_Q_init, lcl_Q_end, muR_Q)
          end do
          call cobj_AllocSplit(sh%gd, this_evop%P, nflcl, nloop)
          call cobj_SetDerivedSplit(this_evop%P,probes)
          if (nflcl == shnf_end) then
             nullify(this_evop%next)
          else
             allocate(this_evop%next)
             this_evop => this_evop%next
          end if
       end if
    end do

    !-- clean up
    call holder_SetNf(shcopy, nfstore)
  end subroutine ev_evolve_varnf_gen


  !-----------------------------------------------------------------
  ! Currently only supports mass thresholds at muF = m_H.
  ! And only when evolving up in scale.
  subroutine ev_CrossMassThreshold(sh,ash,direction,pdfset,evop)
    use holders; use pdf_representation; use dglap_choices
    type(sigma_holder), intent(in)   :: sh
    type(as_handle),   intent(in)    :: ash
    integer,           intent(in)    :: direction
    real(dp),          intent(inout), optional :: pdfset(:,:)
    type(evln_operator), intent(inout), optional :: evop
    !----------------------------------------------
    integer, parameter :: max_warn = 1
    integer            :: warn_id_DIS = warn_id_INIT
    integer            :: warn_id_Direction = warn_id_INIT
    real(dp) :: as2pi, muR
    
    !-- CHANGE THIS IF HAVE MATCHING AT MUF/=MH
    if (ev_nloop < 3) return
    if (.not. mass_steps_on) return
    if (sh%factscheme /= factscheme_MSBar) then
       call wae_Warn(max_warn,warn_id_DIS,&
            &'ev_CrossMassThreshold',&
            &'Factscheme is not MSBar;&
            & mass thresholds requested but not implemented')
       return
    end if
    
    if (direction /= 1) then
       call wae_Warn(max_warn,warn_id_Direction,&
            &'ev_CrossMassThreshold',&
            &'Direction/=1; mass thresholds requested but not implemented')
       return
    end if
    
    !-- now actually do something!
    muR   = quark_masses(nf_int) * ev_MuR_Q
    !write(0,*) 'evolution crossing threshold ', nf_int, muR
    !-- fix nf so as to be sure of getting alpha value corresponding
    !   to the desired nf value, despite proximity to threshold.
    as2pi = as_value(ash, muR, fixnf=nf_int) / twopi
    if (present(pdfset)) pdfset = pdfset + as2pi**2 * (sh%MTM2 .conv. pdfset)
    if (present(evop)) then
       evop%cross_mass_threshold = .true.
       evop%MTM = sh%MTM2
       evop%MTM_coeff = as2pi**2
    end if
  end subroutine ev_CrossMassThreshold
  

  !-----------------------------------------------------
  !! Return the action of the evln_operator on the pdfdist 
  function ev_conv_evop(evop,pdfdist) result(res)
    type(evln_operator), intent(in), target :: evop
    real(dp),          intent(in)         :: pdfdist(:,:)
    real(dp) :: res(size(pdfdist,dim=1),size(pdfdist,dim=2))
    !------------
    type(evln_operator), pointer :: this_evop
    
    this_evop => evop
    res = pdfdist
    do
       if (this_evop%cross_mass_threshold) then
          ! NB: this never eccurs on first pass
          res = res + this_evop%MTM_coeff * (this_evop%MTM .conv. res)
       end if
       
       res = this_evop%P .conv. res

       if (associated(this_evop%next)) then
          this_evop => this_evop%next
       else
          return
       end if
    end do
  end function ev_conv_evop
  

  !-----------------------------------------------------------
  ! Return shnflcl, which is nflcl modified to as to be
  ! within the supported limits of sh.
  !
  ! Thus we can evolve with 5 flavours even into the 
  ! 6 flavour region, without the whole house crashing down
  !
  ! Beware that this sort of thing is dangerous, because
  ! alphas will have one value for nf, beta0 a different value, 
  ! and overall ca sera la pagaille.
  subroutine ev_limit_nf(sh, nflcl, shnflcl, nfname)
    use holders
    type(sigma_holder), intent(in)  :: sh
    integer,            intent(in)  :: nflcl
    integer,            intent(out) :: shnflcl
    character(len=*),   intent(in)  :: nfname
    !----------------------------------------
    integer, parameter :: max_warn = 4
    integer            :: warn_id = warn_id_INIT
    integer            :: nflo, nfhi
    character(len=80)  :: nfstring
    
    nflo = lbound(sh%allP,dim=2)
    nfhi = ubound(sh%allP,dim=2)

    shnflcl = max(min(nflcl,nfhi),nflo)
    
    if (nflcl /= shnflcl) then
       write(nfstring,'(a,i1,a,i1,a)') " changed from ",nflcl," to ",&
            &shnflcl,"."
       call wae_warn(max_warn, warn_id, 'ev_limit_nf: '//&
            &nfname//trim(nfstring))
    end if
  end subroutine ev_limit_nf
  


  !======================================================================
  ! General entry routine for evolution. Takes a human distribution
  ! and deals with the necessary conversion to and from "evolution" format. 
  ! It works for LO and NLO evolution: the order being set by "nloop"
  subroutine ev_evolve(sh, pdfset, nloop, ash, Q_init, Q_end, muR_Q,&
       & untie_nf)
    use holders; use pdf_representation
    type(sigma_holder), intent(in), target :: sh
    real(dp),        intent(inout)         :: pdfset(iflv_min:,0:)
    integer,         intent(in)            :: nloop
    type(as_handle), intent(in), target    :: ash
    real(dp),        intent(in)            :: Q_init, Q_end
    real(dp),        intent(in), optional  :: muR_Q
    logical,         intent(in), optional  :: untie_nf
    !-----------------------------------------------------
    real(dp) :: pdf_ev(size(pdfset,dim=1), size(pdfset,dim=2))
    integer  :: pdfrep

    !-- make sure that number of flavours is correct?
    if (nf_int /= sh%prep%nf) call wae_error('ev_evolve: &
         &global nf and representation nf are not equal.')

    !-- put things into the right format just once (it would
    !   be done automatically by the conv_object routines, but
    !   that would be wasteful)
    pdfrep = pdfr_GetRep(pdfset)
    if (pdfrep == pdfr_Human) then
       call pdfr_HumanToEvln(sh%prep, pdfset, pdf_ev)
    else
       pdf_ev = pdfset
    end if
    

    if (nloop > sh%nloop) &
         &call wae_error('ev_evolve: sh%nloop must be >= nloop')

    ev_PLO  => sh%P
    if (nloop >= 2) ev_PNLO => sh%P1
    if (nloop >= 3) ev_PNNLO => sh%P2
    ev_ash => ash
    call ev_SetModuleConsts(nloop, muR_Q, untie_nf)
!!$    ev_nloop = nloop
!!$    ev_muR_Q = default_or_opt(one,muR_Q)

    call ev_evolveLocal(pdf_ev, Q_init, Q_end)

    !-- put things back into a "human" format
    if (pdfrep == pdfr_Human) then
       call pdfr_EvlnToHuman(sh%prep, pdf_ev, pdfset)
    else
       pdfset = pdf_ev
    end if
  end subroutine ev_evolve


  !---------------------------------------------------------
  ! A shortcut for setting up copies of information
  ! which may be useful module-wide
  !
  ! Is it really needed?
  subroutine ev_SetModuleConsts(nloop, muR_Q, untie_nf)
    integer,         intent(in)            :: nloop
    real(dp),        intent(in), optional  :: muR_Q
    logical,         intent(in), optional  :: untie_nf

    ev_nloop = nloop
    ev_muR_Q = default_or_opt(one,muR_Q)
    ev_untie_nf = default_or_opt(.false., untie_nf)
  end subroutine ev_SetModuleConsts
  


  !======================================================================
  ! Takes pdfset in the MSBar scheme and converts it into the DIS scheme,
  ! hopefully correctly!
  subroutine ev_MSBar2DIS(sh,pdfset,nloop,ash,Q)
    use holders; use pdf_representation; use convolution
    type(sigma_holder), intent(in)         :: sh
    real(dp),        intent(inout)         :: pdfset(:,:)
    integer,         intent(in)            :: nloop
    type(as_handle), intent(in), target    :: ash
    real(dp),        intent(in)            :: Q
    !-----------------------------------------------------
    real(dp) :: pdf_ev(size(pdfset,dim=1), -ncomponents:ncomponents )
    real(dp) :: Cq_x_q(size(pdfset,dim=1)), Cg_x_g(size(pdfset,dim=1))
    real(dp) :: as2pi
    integer  :: id

    !-- put things into the right format
    call pdfr_HumanToEvln(sh%prep, pdfset, pdf_ev)

    as2pi = as_value(ash,Q) / twopi
    
    Cq_x_q = as2pi * (sh%C2_1%q .conv. pdf_ev(:,iflv_sigma))
    Cg_x_g = (2*sh%prep%nf * as2pi) * (sh%C2_1%g .conv. pdf_ev(:,iflv_g))

    pdf_ev(:,iflv_sigma) = pdf_ev(:,iflv_sigma) + Cq_x_q + Cg_x_g
    pdf_ev(:,iflv_g)     = pdf_ev(:,iflv_g)     - Cq_x_q - Cg_x_g

    do id = -sh%prep%nf, sh%prep%nf
       if (id == iflv_sigma .or. id == iflv_g) cycle
       pdf_ev(:,id) = pdf_ev(:,id) + as2pi * (sh%C2_1%q .conv. pdf_ev(:,id))
    end do

    !-- put things back into a "human" format
    call pdfr_EvlnToHuman(sh%prep, pdf_ev, pdfset)
    !write(0,*) 'hello', lbound(pdfset), lbound

  end subroutine ev_MSBar2DIS
  

  

  !----------------------------------------------------------------------
  ! this bit here does the actual evolution. Could make it more
  ! sophisticated later, if necessary (e.g. variable step size).
  subroutine ev_evolveLocal(pdfset, Q_init, Q_end)
    use runge_kutta
    real(dp),        intent(inout)         :: pdfset(:,:)
    real(dp),        intent(in)            :: Q_init, Q_end
    !------------------------------------------------------------
    real(dp) :: t1, t2, dt, t, u1, u2, du, u, as1, as2
    integer  :: n, i
    integer :: ntot=0

    fourpibeta0_lnmuR_Q = four*pi*beta0*log(ev_muR_Q)
    !write(0,*) fourpibeta0_lnmuR_Q
    t1 = two*log(Q_init)
    t2 = two*log(Q_end)

    if (t2 == t1) return

    !-- now work out jacobians...
    as1 = ev_asval(Q_init)
    as2 = ev_asval(Q_end)

    ! allow for both signs of coupling (and sign changes)
    ! See CCN25-95 for formulae (NB: 0->1 and 1->2)
    if (ev_force_old_dt) then
       ev_du_type = ev_du_is_dt
       u1 = t1; u2 = t2
       n  = ceiling(abs(t2-t1)/dt_ballpark)
    else if (abs(as1 - as2)/max(as1,as2) < ev_minalphadiff &
         &.or. as1*as2 <= zero) then
    !else if (.true.) then
       ev_du_type = ev_du_is_dtas_fixed
       ev_du_dt = max(abs(as1),abs(as2))
       u1 = t1 * ev_du_dt; u2 = t2 * ev_du_dt
       n  = ceiling(abs(u2-u1)/du_ballpark)
    else
       ev_du_type = ev_du_is_dtas_run
       ev_u_asref = as1
       ev_u_tref  = t1
       ev_u_bval  = (as1/as2-1)/(as1*(t2-t1))
       u1 = zero
       u2 = log(1+ev_u_bval*ev_u_asref*(t2-t1)) / ev_u_bval
       n  = ceiling(abs(u2-u1)/du_ballpark)
       !write(0,*) ev_u_asref, ev_u_bval
    end if
    
    !n  = ceiling(abs(t2-t1)/dt_ballpark)
    !dt = (t2 - t1)/n
    du = (u2 - u1)/n

    ntot = ntot + n
    !write(0,*) 'Qinit,end, nsteps', Q_init, Q_end, n, ntot
    u = u1
    do i = 1, n
       call rkstp(du, u, pdfset, ev_conv)
    end do
  end subroutine ev_evolveLocal


  !----------------------------------------------------------------------
  subroutine ev_conv(u, pdfset, dpdfset)
    use pdf_representation
    ! The following fails with absoft compiler: see ABSOFT_BUG.txt
    real(dp), intent(in)  :: u, pdfset(:,:)
    real(dp), intent(out) :: dpdfset(:,:)
    ! The following fails with the intel compiler! See INTEL_BUG.txt
    !real(dp), intent(in)  :: t, pdfset(0:,-ncomponents:)
    !real(dp), intent(out) :: dpdfset(0:,-ncomponents:)
    !--------------------------------------
    real(dp) :: as2pi, Q, t, jacobian
    type(PMat) :: Pfull

    ! for analysing Intel bug
    !write(0,*) 'X',lbound(pdfset),lbound(pdfset,dim=1),size(pdfset,dim=1)

    select case(ev_du_type)
    case(ev_du_is_dt)
       t = u
       jacobian = one
    case(ev_du_is_dtas_fixed)
       t = u / ev_du_dt
       jacobian = one / ev_du_dt
    case(ev_du_is_dtas_run)
       t = ev_u_tref + (exp(ev_u_bval*u)-one)/(ev_u_bval*ev_u_asref)
       jacobian = (1+ev_u_bval*ev_u_asref*(t-ev_u_tref))/ev_u_asref
    case default
       call wae_error('evconv: unknown ev_du_type',intval=ev_du_type)
    end select
    
    Q     = exp(half*t)
    as2pi = ev_asval(Q)/twopi
    
    select case (ev_nloop)
    case(1)
       dpdfset = (jacobian * as2pi) * (ev_PLO .conv. pdfset)
    case(2)
       if (fourpibeta0_lnmuR_Q /= zero) then
          call cobj_InitSplit(Pfull, ev_PLO, one + as2pi*fourpibeta0_lnmuR_Q)
       else
          call cobj_InitSplit(Pfull, ev_PLO)
       end if
       call cobj_AddSplit(Pfull, ev_PNLO, as2pi)
       dpdfset = (jacobian * as2pi) * (Pfull .conv. pdfset)
       call cobj_DelSplit(Pfull)
!!$       !---- TMP TMP ------------------
!!$       write(0,*) real(dpdfset(:,-1)/as2pi)
    case(3)
       if (fourpibeta0_lnmuR_Q /= zero) then
          call cobj_InitSplit(Pfull, ev_PLO, one + as2pi*fourpibeta0_lnmuR_Q&
               & + (as2pi*fourpibeta0_lnmuR_Q)**2&
               & + as2pi**2*(twopi**2*beta1)*two*log(ev_muR_Q))
          call cobj_AddSplit(Pfull, ev_PNLO, &
               &as2pi*(one + two*as2pi*fourpibeta0_lnmuR_Q))
          call cobj_AddSplit(Pfull, ev_PNNLO, as2pi**2)
          !call wae_error('ev_conv: NNL evolution not supported with muR_Q/=1')
       else
          call cobj_InitSplit(Pfull, ev_PLO)
          call cobj_AddSplit(Pfull, ev_PNLO, as2pi)
          call cobj_AddSplit(Pfull, ev_PNNLO, as2pi**2)
       end if
       dpdfset = (jacobian * as2pi) * (Pfull .conv. pdfset)
       call cobj_DelSplit(Pfull)
    end select
  end subroutine ev_conv
  
  
  !--------------------------------------------------------
  ! alphas will be needed potentially in several locations,
  ! so make it common here....
  function ev_asval(Q) result(res)
    real(dp), intent(in) :: Q
    real(dp)             :: res
    if (.not. ev_untie_nf) then
       !-- fixnf option here will be quite restrictive. It means that
       !   if we have ev_muR_Q/=1 and variable numbers of flavours,
       !   then either ev_ash supports "extrapolation" with the same nf
       !   beyond the strictly legal region, or else it has been defined 
       !   so as to have flavour thresholds at ev_muR_Q*masses
       res   = as_value(ev_ash, Q*ev_muR_Q, fixnf=nf_int)
    else
       ! sometimes (e.g. comparisons with others) it is useful to
       ! allow nf in alpha_s to have a value different from the nf
       ! being used in the splitting functions...
       res   = as_value(ev_ash, Q*ev_muR_Q)
    end if

  end function ev_asval
  

  !---------------------------------------------------------
  !! Delete the objects allocated in evop, including any
  !! subsiduary objects.
  recursive subroutine ev_DelEvOp(evop)
    type(evln_operator), intent(inout) :: evop

    if (associated(evop%next)) then
       call ev_DelEvOp(evop%next)
       deallocate(evop%next)
    end if
    
    ! do nothing here since MTM has not actually been allocated
    ! but rather just set "equal" to sh%MTM2 -- i.e. it just points
    ! to the contents of sh%MTM2
    !if (evop%cross_mass_threshold) then ! do nothing

    call cobj_DelSplit(evop%P)
    return
  end subroutine ev_DelEvOp
  

end module evolution

