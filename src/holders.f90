!======================================================================
! Module containing types which hold properties of QCD.
!
! E.g. splitting and coefficient functions.
!
! Details of event shapes
!
module holders
  use types; use consts_dp
  use conv_objects; use pdf_representation
  use assertions; use warnings_and_errors
  implicit none
  private

  !-- Everything needed to calculate a cross section, given structure
  !   functions and alpha_s.
  !
  !   nomenclature is: objects are leading order, unless there is a suffix n, 
  !   in which case they are NLO.
  type sigma_holder
     type(grid_def) :: gd
     !-----------------------  nloop, nf  --------------
     type(Pmat), pointer :: allP(:,   :)
     type(PMat), pointer :: P, P1, P2
     type(Cmat), pointer :: allC(:,:)
     type(CMat), pointer :: C2, C2_1, CL_1
     !-- indep of nf for the time being ----------------
     type(MassThresholdMat) :: MTM2
     logical    :: MTM2_exists=.false.
     integer    :: factscheme, nloop
     !--------------------------------  nf  ------------
     type(pdf_rep), pointer :: all_prep(:)
     type(pdf_rep), pointer :: prep
     integer                :: nf
  end type sigma_holder
  
  !-- this is potentially makeshift?
  integer, parameter :: nC = 3

  public :: sigma_holder, holder_InitSigma, holder_SetNf

contains
  
  !-------------------------------------------------------
  ! Sets up eveything needed to calculate cross sections
  ! i.e. splitting functions and coefficient functions
  subroutine holder_InitSigma(gd, sh, factscheme, nloop, nflo, nfhi)
    use coefficient_functions; use convolution
    use qcd; use dglap_choices
    type(grid_def),     intent(in)    :: gd
    type(sigma_holder), intent(inout), target :: sh
    integer, optional,  intent(in)    :: factscheme
    integer, optional,  intent(in)    :: nloop
    integer, optional,  intent(in)    :: nflo, nfhi
    !-- holds temporary results
    type(grid_conv) :: dconv
    !-- holds all possible combinations of coefficient and splitting functions
    !   needed for DIS schemes
    type(grid_conv) :: Cq,Cg
    type(grid_conv) :: CqPqq, CqPqg, CqPgq, CqPgg
    type(grid_conv) :: CgPqq, CgPqg, CgPgq, CgPgg
    !-- more compactly written version of DIS scheme
    type(grid_conv) :: CC(id_g:id_sigma,id_g:id_sigma)
    type(grid_conv) :: tmp2d(id_g:id_sigma,id_g:id_sigma)
    logical :: newDIS = .true.
    !logical :: newDIS = .false.
    integer :: nfstore, nflcl

    sh%factscheme = default_or_opt(factscheme_default, factscheme)
    sh%nloop     = default_or_opt(2, nloop)
    if (sh%nloop > 3 .or. sh%nloop < 1) then
       call wae_error('holder_InitSigma: nloop must be between 1 and 3')
    end if
    sh%gd = gd

    if (present(nflo) .and. present(nfhi)) then
       allocate(sh%allP(sh%nloop,nflo:nfhi))
       allocate(sh%all_prep(nflo:nfhi))
       allocate(sh%allC(nC,nflo:nfhi))
    else
       !-- otherwise, use whatever nf is currently set
       allocate(sh%allP(sh%nloop,nf_int:nf_int))
       allocate(sh%all_prep(nf_int:nf_int))
       allocate(sh%allC(nC,nf_int:nf_int))
    end if
    
    !-- want to reset it at end
    nfstore = nf_int
    
    do nflcl = lbound(sh%allP,dim=2), ubound(sh%allP,dim=2)
       !-- this sets up all the pointers to make it look like a fixed-nf
       !   sh!
       call holder_SetNf(sh, nflcl)

       !----- this is the remainder, essentially unchanged ----
       select case (sh%factscheme)
       case (factscheme_MSbar)
          call cobj_InitSplitLO (gd, sh%P)
          if (sh%nloop >= 2) call cobj_InitSplitNLO(gd, sh%P1, sh%factscheme)
          if (sh%nloop >= 3) then
             call cobj_InitSplitNNLO(gd, sh%P2, sh%factscheme)
             !-- do this once, and only if really needed
             if (lbound(sh%allP,dim=2) /= ubound(sh%allP,dim=2) &
                  &.and. mass_steps_on &
                  &.and. nflcl == lbound(sh%allP,dim=2)) then
                call cobj_InitMTMNNLO(gd,sh%MTM2)
                sh%MTM2_exists = .true.
             end if
          end if

          call cobj_InitCoeff(gd, sh%C2)
          call cobj_InitCoeff(gd, sh%C2_1, cf_CgF2MSbar, cf_CqF2MSbar)
          call cobj_InitCoeff(gd, sh%CL_1,  cf_CgFL, cf_CqFL)
       case (factscheme_DIS)
          call cobj_InitSplitLO (gd, sh%P)
          if (sh%nloop >= 2) &
               &call cobj_InitSplitNLO(gd, sh%P1, factscheme_MSbar)
          if (sh%nloop >= 3) write(0,*) &
               &'DIS factorisation scheme not supported for 3 loops or more'
          call cobj_InitCoeff(gd, sh%C2)
          call cobj_InitCoeff(sh%C2_1, sh%C2, zero)
          call cobj_InitCoeff(gd, sh%CL_1,  cf_CgFL, cf_CqFL)
          
          !-- now convert MSbar splitting functions into DIS scheme -------
          ! See CCN21-6 (and also CCN17-61)
          
          if (newDIS) then
             !
             ! NB THIS VERSION OF THE DIS SCHEME IS NOT YET IN AGREEMENT
             !    WITH THE OLDER VERSION...
             !
             !-- create the matrix C for use in 
             !   P_matrix -> P_matrix + [C,P]
             !           /  Cq   Cg \
             !     C ==  |          |   where Cg includes 2nf factor
             !           \ -Cq  -Cg /
             call conv_InitGridConv(gd,CC(id_sigma,id_sigma),cf_CqF2MSbar)
             call conv_InitGridConv(gd,CC(id_sigma,id_g),    cf_CgF2MSbar)
             call conv_MultGridConv(CC(id_sigma,id_g), two*nf)
             call conv_InitGridConv(CC(id_g,:),CC(id_sigma,:),-one)
             !-- now get a temporary to hold the commutator
             call conv_AllocGridConv(gd,tmp2d)
             !-- work out the commutator: tmp2d=[C,P]
             !write(0,*) 'lb',lbound(sh%P%singlet), ubound(sh%P%singlet)
             !write(0,*) 'sh%P%singlet(0,0)%nsub',sh%P%singlet(0,0)%gd%nsub
             !write(0,*) 'sh%P%singlet(0,1)%nsub',sh%P%singlet(0,1)%gd%nsub
             !-----------------------------------------------
             ! putting the explicit singlet bounds somehow eliminates
             ! an ifc memory error, which was associated with 
             ! gcb in conv_CommGridConv obtaining ubound=(/220,1/)
             ! No understanding of origin of error and locally
             ! singlet has right bounds
             !
             ! result however seems to be wrong.
             !
             ! Message: be wary of intel here
             !call conv_CommGridConv(tmp2d,CC,sh%P%singlet(id_g:id_sigma,id_g:id_sigma))
             call conv_CommGridConv(tmp2d,CC,sh%P%singlet)
             !-- add it to P1
             call conv_AddGridConv(sh%P1%singlet,tmp2d)
             !-- add the beta function pieces as well
             call conv_AddGridConv(sh%P1%singlet,CC, -twopi_beta0)
             call conv_AddGridConv(sh%P1%NS_plus, CC(id_sigma,id_sigma),&
                  & -twopi_beta0)
             !-- quark number conservation remains OK because
             !   Cq has the om=0 moment equal to zero.
             call conv_AddGridConv(sh%P1%NS_minus, CC(id_sigma,id_sigma), &
                  &-twopi_beta0)
             call conv_AddGridConv(sh%P1%NS_V, CC(id_sigma,id_sigma), &
                  &-twopi_beta0)
             
             !-- clean up
             call conv_DelGridConv(CC)
             call conv_DelGridConv(tmp2d)
          else
             call conv_InitGridConv(gd, Cq, cf_CqF2MSbar)
             call conv_InitGridConv(gd, Cg, cf_CgF2MSbar)
             call conv_MultGridConv(Cg, two*nf)
             !   where possible put the smoother distribution to the right
             !   (only makes a difference when pdfs are non zero at x=1?)
             call conv_ConvConv(CqPqq, Cq, sh%P%qq)
             call conv_ConvConv(CqPqg, Cq, sh%P%qg)
             call conv_ConvConv(CqPgq, Cq, sh%P%gq)
             call conv_ConvConv(CqPgg, Cq, sh%P%gg)
             call conv_ConvConv(CgPqq, sh%P%qq, Cg)
             call conv_ConvConv(CgPqg, sh%P%qg, Cg)
             call conv_ConvConv(CgPgq, sh%P%gq, Cg)
             call conv_ConvConv(CgPgg, sh%P%gg, Cg)
             
             !
             !   First deal with P_matrix -> P_matrix + [C,P]
             !           /  Cq   Cg \
             !     C ==  |          |   where Cg includes 2nf factor
             !           \ -Cq  -Cg /
             !
             !   Pqq -> Pqq + (Cg Pgq + Cq Pqg)
             call conv_AddGridConv(sh%P1%qq, CgPgq, one)
             call conv_AddGridConv(sh%P1%qq, CqPqg, one)
             !   Pqg -> Pqg + (Cq Pqg + Cg Pgg - Cg Pqq + Cg Pqg)
             call conv_AddGridConv(sh%P1%qg, CqPqg, one)
             call conv_AddGridConv(sh%P1%qg, CgPgg, one)
             call conv_AddGridConv(sh%P1%qg, CgPqq,-one)
             call conv_AddGridConv(sh%P1%qg, CgPqg, one)
             !   Pgq -> Pqg + (-Cq Pqq - Cg Pgq - Cq Pgq + Cq Pgg)
             call conv_AddGridConv(sh%P1%gq, CqPqq,-one)
             call conv_AddGridConv(sh%P1%gq, CgPgq,-one)
             call conv_AddGridConv(sh%P1%gq, CqPgq,-one)
             call conv_AddGridConv(sh%P1%gq, CqPgg, one)
             !   Pgg -> Pgg + (-Cq Pqg - Cg Pgq)
             call conv_AddGridConv(sh%P1%gg, CqPqg,-one)
             call conv_AddGridConv(sh%P1%gg, CgPgq,-one)
             !
             !   Now deal with P_matrix -> P_matrix - beta0 * C
             !                 P_+      -> P_+      - beta0 * C_q
             call conv_AddGridConv(sh%P1%qq, Cq, -twopi_beta0)
             call conv_AddGridConv(sh%P1%qg, Cg, -twopi_beta0)
             call conv_AddGridConv(sh%P1%gq, Cq, +twopi_beta0)
             call conv_AddGridConv(sh%P1%gg, Cg, +twopi_beta0)
             call conv_AddGridConv(sh%P1%NS_plus, Cq, -twopi_beta0)
             !-- quark number conservation remains OK because
             !   Cq has the om=0 moment equal to zero.
             call conv_AddGridConv(sh%P1%NS_minus, Cq, -twopi_beta0)
             call conv_AddGridConv(sh%P1%NS_V, Cq, -twopi_beta0)
             !
             !   tidy up
             !write(0,*) 'Hey:',sh%P%singlet(0,0)%subgc(1)%conv(0:3,1)
             !write(0,*) 'Hey:',sh%P%singlet(0,1)%subgc(1)%conv(0:3,1)
             !write(0,*) 'Hey:',sh%P%singlet(1,0)%subgc(1)%conv(0:3,1)
             !write(0,*) 'Hey:',sh%P%singlet(1,1)%subgc(1)%conv(0:3,1)
             !write(0,*) 'Hey:',               Cq%subgc(1)%conv(0:3,1)
             !write(0,*) 'Hey:',               Cg%subgc(1)%conv(0:3,1)
             
             call conv_DelGridConv(Cq)
             call conv_DelGridConv(Cg)
             call conv_DelGridConv(CqPqq)
             call conv_DelGridConv(CqPqg)
             call conv_DelGridConv(CqPgq)
             call conv_DelGridConv(CqPgg)
             call conv_DelGridConv(CgPqq)
             call conv_DelGridConv(CgPqg)
             call conv_DelGridConv(CgPgq)
             call conv_DelGridConv(CgPgg)
             write(0,*) 'result:',sh%P1%singlet(id_sigma,id_sigma)%conv(10,1)
          end if
          
       case (factscheme_PolMSbar)
          write(0,*) "SETTING UP POLARIZED EVOLUTION"
          call cobj_InitSplitPolLO (gd, sh%P)
          if (sh%nloop >= 2) call cobj_InitSplitPolNLO(gd, &
               &sh%P1, sh%factscheme)
          if (sh%nloop >= 3) call wae_error('holder_InitSigma',&
               &'nloop >= 3 not supported for polarized case')
          !!if (sh%nloop >= 3) then
          !!   call cobj_InitSplitNNLO(gd, sh%P2, sh%factscheme)
          !!   !-- do this once, and only if really needed
          !!   if (lbound(sh%allP,dim=2) /= ubound(sh%allP,dim=2) &
          !!        &.and. mass_steps_on &
          !!        &.and. nflcl == lbound(sh%allP,dim=2))&
          !!        & call cobj_InitMTMNNLO(gd,sh%MTM2)
          !!end if
          !!
          !!call cobj_InitCoeff(gd, sh%C2)
          !!call cobj_InitCoeff(gd, sh%C2_1, cf_CgF2MSbar, cf_CqF2MSbar)
          !!call cobj_InitCoeff(gd, sh%CL_1,  cf_CgFL, cf_CqFL)
       case default
          write(0,*) 'factorisation scheme ',sh%factscheme,&
               &' is not currently supported' 
          stop
       end select
       
       !-- used for converting to and from the evolution representation
       sh%prep = cobj_DefaultPrep(nf_int)
!!$       sh%prep%nf    = nf_int
!!$       !sh%prep%ibase = nf_int
!!$       sh%prep%ibase = 1
    end do

    

    !-- be "clean" ----------
    call holder_SetNf(sh,nfstore)
  end subroutine holder_InitSigma


  !--------------------------------------------------------------
  ! set up all pointers so that it looks like an old holder in the
  ! good-old fixed-nf days; it also sets the global nf.
  ! (NB: perhaps that is not too good?)
  subroutine holder_SetNf(sh, nflcl)
    use qcd
    type(sigma_holder), intent(inout) :: sh
    integer, intent(in) :: nflcl

    if (nflcl < lbound(sh%allP,dim=2) .or. nflcl > ubound(sh%allP,dim=2)) then
       call wae_Error('holder_SetNf: tried to set unsupported nf')
    end if
    
    !-- want general nf to be consistent. Not really all that nice a 
    !   way of doing things, but cannot think of a better solution 
    !   given the current structure...
    call qcd_SetNf(nflcl)
    if (sh%mtm2_exists) call cobj_SetNfMTM(sh%MTM2, nflcl)

    !-- set up links so that the remainder of the routine
    !   can stay unchanged
    sh%P  => sh%allP(1,nflcl)
    if (sh%nloop >= 2) then
       sh%P1 => sh%allP(2,nflcl)
    else
       nullify(sh%P1)
    end if
    
    if (sh%nloop >= 3) then
       sh%P2 => sh%allP(3,nflcl)
    else
       nullify(sh%P2)
    end if
    
    sh%C2 => sh%allC(1,nflcl)
    sh%C2_1 => sh%allC(2,nflcl)
    sh%CL_1 => sh%allC(3,nflcl)
    sh%prep => sh%all_prep(nflcl)
    sh%nf = nflcl
  end subroutine holder_SetNf
  
end module holders
