
! on tycho the following combination runs in 4ms per pdf and has
! a relative accuracy (on g,d,u,s) that is typically \lesssim 1e-4, though
! one might want to check the accuracy more systematically.
! 
! ./small_fast_tab -nrep 1 -nxQ 500 -output -dy 0.25 -dlnlnQ 0.07 -dt 0.4 -olnlnQ 4 -order 6 -order2 -6
!
! roughly speaking: evolution takes about 2.5ms and evaluations 1.5ms
!
! Funny thing is that earlier it seemed that runs were much faster
! (factor of two for evolution) and I can't figure out what has
! changed...

module initpdf
  use pdfevln_all_public
  implicit none
  
contains
  subroutine VogtInitSub(y,res)
    real(dp), intent(in)  :: y
    real(dp), intent(out) :: res(ncompmin:)
    real(dp) :: x
    real(dp) :: uv, dv, ubar, dbar
    !---------------------
    real(dp), parameter :: N_g = 1.7_dp, N_ls = 0.387975_dp
    real(dp), parameter :: N_uv=5.107200_dp, N_dv = 3.064320_dp
    real(dp), parameter :: N_db = half*N_ls
    
    res = zero
    x = exp(-y)
    !-- remember that my definitions that these are all x*q(x)
    uv = N_uv * x**0.8_dp * (1-x)**3
    dv = N_dv * x**0.8_dp * (1-x)**4
    dbar = N_db * x**(-0.1_dp) * (1-x)**6
    ubar = dbar * (1-x)
    res(id_g) = N_g * x**(-0.1_dp) * (1-x)**5

    res(-id_s) = 0.2_dp*(dbar + ubar)
    res( id_s) = res(-id_s)
    res( id_u) = uv + ubar
    res(-id_u) = ubar
    res( id_d) = dv + dbar
    res(-id_d) = dbar
  end subroutine VogtInitSub
end module initpdf

!======================================================================
program small_fast_tab
  use pdfevln_all_public
  use initpdf
  use sub_defs_io
  implicit none
  !-------------------------------
  type(grid_def)     :: gd, gdarray(3)
  type(sigma_holder) :: sh
  type(as_handle)    :: ash
  integer            :: order, order2, nloop, i, j, nrep, nxQ,nn, olnlnQ
  real(dp)           :: dy, Qinit, Qmax, y, Q, pdfval(-6:6), dt, dlnlnQ
  real(dp)           :: ymax
  real(dp), pointer  :: vogt_init(:,:)
  type(pdftab)       :: table
  logical            :: output

  ! set the details of the y=ln1/x grid
  dy    = dble_val_opt('-dy',0.25_dp)
  ymax  = dble_val_opt('-ymax',5.0_dp)
  order = int_val_opt('-order',5)
  order2 = int_val_opt('-order2',order)
  if (log_val_opt('-alt')) then
     call conv_InitGridDef(gdarray(3),dy/7.0_dp, 1.0_dp, order=order2)
     call conv_InitGridDef(gdarray(2),dy/2.0_dp, 3.0_dp, order=order2)
     call conv_InitGridDef(gdarray(1),dy,        ymax, order=order)
  else
     call conv_InitGridDef(gdarray(3),dy/9.0_dp, 0.5_dp, order=order2)
     call conv_InitGridDef(gdarray(2),dy/3.0_dp, 2.0_dp, order=order2)
     call conv_InitGridDef(gdarray(1),dy,        ymax, order=order)
  end if
  call conv_InitGridDef(gd,gdarray(1:3),locked=.true.)
  ! set parameter for evolution step in Q
  dt = dble_val_opt('-dt',0.4_dp)
  call ev_Setdt(dt)

  ! set up the splitting functions
  nloop = 3
  if (log_val_opt('-exactth')) &
       &   call dglap_Set_nnlo_nfthreshold(nnlo_nfthreshold_exact)
  if (log_val_opt('-exactsp')) &
       &   call dglap_Set_nnlo_splitting(nnlo_splitting_exact)
  call holder_InitSigma(gd, sh, factscheme=factscheme_MSbar, &
       &                              nloop=nloop, nflo=3, nfhi=6)
  
  ! get the initial distribution
  call pdfgen_AllocInitPDFSub(gd,vogt_init,VogtInitSub)

  ! set up the coupling
  Qinit = sqrt(two); Qmax = dble_val_opt('-Qmax',50.0_dp)
  call as_init(ash, alfas=0.35_dp, Q=Qinit, nloop=nloop, use_nah=.true.)

  ! set up the table
  dlnlnQ = dble_val_opt('-dlnlnQ',0.07_dp)
  olnlnQ = int_val_opt('-olnlnQ',4)
  call pdftab_AllocTab(gd,table,Qinit,Qmax,dlnlnQ,lnlnQ_order=olnlnQ)
  call pdftab_AssocNfInfo(table,ash)
  call pdftab_PreEvolve(table,Qinit,sh,nloop,ash)

  nrep  = int_val_opt('-nrep',1)
  nxQ = int_val_opt('-nxQ',0); y = 3.14_dp; Q = 13.354_dp
  output = log_val_opt('-output')
  if (output) call output_info
 
  do i = 1, nrep
     call pdftab_InitTabEvolve(table,vogt_init)
     
     nn = nxQ/3
     do j = 1, nn
        y = j*5.0_dp/nn
        Q = Qmax - j*(Qmax-Qinit)/nn
        call pdftab_ValTab_yQ(table,y,Q,pdfval)
        if (output .and. i==1) write(6,'(20es20.8)') y,Q,pdfval(0:3)
        !if (output .and. i==1) write(6,'(20es20.8)') y,Q,vogt_init(:,0:3).atx.(gd.with.exp(-y))
     end do

     if (output) write(6,*)
     if (output) write(6,*)
     do j = nn,1,-1
        y = j*5.0_dp/nn
        Q = 4.0_dp + j*5.0_dp/nn
        call pdftab_ValTab_yQ(table,y,Q,pdfval) 
        if (output .and. i==1) write(6,'(20es20.8)') y,Q,pdfval(0:3)
     end do
     
     if (output) write(6,*)
     if (output) write(6,*)
     do j = nn,1,-1
        y = j*5.0_dp/nn
        Q = Qmax*(1-j*0.2_dp/nn)
        call pdftab_ValTab_yQ(table,y,Q,pdfval) 
        if (output .and. i==1) write(6,'(20es20.8)') y,Q,pdfval(0:3)
     end do
  end do
  

contains
  subroutine output_info
    write(6,'(a)') '# '//trim(command_line())
    write(6,'(a,f10.5,a,f10.5)') '# dy = ',dy,      ';    ymax = ',ymax
    write(6,'(a,f10.5,a,f10.5)') '# dt = ',dt,      ';    Qmax = ',Qmax
    write(6,'(a,i5,a,i5)')    '# order  = ',order,  ';    order2 = ',order2
    write(6,'(a,f10.5,a,i5)') '# dlnlnQ = ',dlnlnQ, ';    olnlnQ = ',olnlnQ
  end subroutine output_info
  
end program small_fast_tab
