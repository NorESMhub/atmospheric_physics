module dme_adjust_camnor

  use shr_kind_mod,  only: r8 => shr_kind_r8

  implicit none
  private          ! Make default type private to the module

  public :: dme_adjust_camnor_run

  logical :: levels_are_moist=.true.
  logical :: hydrostatic = .true.
  real(r8), parameter :: rtiny = 1e-14_r8    ! a small number (relative to total q change)

contains

  subroutine dme_adjust_camnor_run(lchnk, ncol, &
       state_psetcols, state_pint, state_pmid, &
       state_pdel, state_rpdel, state_pdeldry, state_lnpint, state_lnpmid, &
       state_ps, state_phis, state_zm, state_zi, &
       state_t, state_u, state_v, state_q, state_s, &
       tend_dudt, tend_dvdt, tend_dtdt, &
       qini, liqini, iceini, dt, &
       ntrnprd, ntsnprd, tevap, tprec, mflx, eflx, eflx_out, mflx_out, &
       ent_tnd, pdel_rf)

    !-----------------------------------------------------------------------
    !
    ! Purpose: Adjust the dry mass in each layer back to the value of physics input state
    !          Adjust air specific enthalpy accordingly. Diagnose boundary enthalpy flux.
    !
    ! Method
    !     Revised adjustment towards consistency with local energy conservation.
    !     Hydrostatic pressure work, de = alpha * dp, where alpha is the specific volume
    !     pressure adjustment, is added locally as an source of enthalpy. An enthalpy of
    !     mass (water) exchange with the surface is also defined, which should be passed
    !     to the surface model components (ocean/land/ice etc).
    !     If moist thermodynamics where handled correctly in CAM, the two terms would
    !     match, guaranteeing local energy conservation.
    !     With the present CAM formulation (constant dry heat capacity, constant latent
    !     heat of condensation valid for 0 degree C), consistency demands one of these
    !     choices:
    !        1. no pressure work and no boundary enthalpy flux (CESM)
    !        2. correct local pressure work and boundary enthalpy flux equal to (S dp/g)
    !          where S=local *dry* static energy of air
    !     The boundary enthalpy flux is at present not passed to other model components,
    !     so it is treated as internal CAM non-conservation and folded into fix_energy.
    !
    ! Author: Thomas Toniazzo (17.07.21)
    !
    !-----------------------------------------------------------------------

    use constituents,    only: pcnst, qmin
    use shr_const_mod,   only: shr_const_rwv
    use ppgrid,          only: pcols, pver
    use geopotential,    only: geopotential_t
    use phys_control,    only: waccmx_is
    use air_composition, only: dry_air_species_num
    use air_composition, only: thermodynamic_active_species_num
    use air_composItion, only: thermodynamic_active_species_idx
    use air_composition, only: cpairv, rairv, cp_or_cv_dycore
    use constituents,    only: cnst_get_ind, cnst_type
    use cam_thermo,      only: inv_conserved_energy
    use cam_thermo,      only: get_conserved_energy
    use cam_thermo,      only: cam_thermo_water_update
    use dyn_tests_utils, only: vc_dycore
    use qneg_module,     only: qneg3
    use cam_history,     only: outfld
    use physconst,       only: cpair, cpwv, cpliq, cpice, gravit, zvir
    !
    ! Arguments
    !
    integer,          intent(in)    :: lchnk
    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: state_psetcols
    real(r8),         intent(inout) :: state_pint(:,:)
    real(r8),         intent(out)   :: state_pmid(:,:)
    real(r8),         intent(inout) :: state_pdel(:,:)
    real(r8),         intent(out)   :: state_rpdel(:,:)
    real(r8),         intent(in)    :: state_pdeldry(:,:)
    real(r8),         intent(out)   :: state_lnpint(:,:)
    real(r8),         intent(out)   :: state_lnpmid(:,:)
    real(r8),         intent(in)    :: state_phis(:)
    real(r8),         intent(inout) :: state_ps(:)
    real(r8),         intent(inout) :: state_zm(:,:)
    real(r8),         intent(inout) :: state_zi(:,:)
    real(r8),         intent(inout) :: state_t(:,:)
    real(r8),         intent(inout) :: state_u(:,:)
    real(r8),         intent(inout) :: state_v(:,:)
    real(r8),         intent(inout) :: state_q(:,:,:)
    real(r8),         intent(inout) :: state_s(:,:)
    real(r8),         intent(inout) :: tend_dudt(:,:)
    real(r8),         intent(inout) :: tend_dvdt(:,:)
    real(r8),         intent(inout) :: tend_dtdt(:,:)
    real(r8),         intent(in)    :: qini(pcols,pver)     ! initial specific humidity
    real(r8),         intent(in)    :: liqini(pcols,pver)   ! initial total liquid
    real(r8),         intent(in)    :: iceini(pcols,pver)   ! initial total ice
    real(r8),         intent(in)    :: dt
    real(r8),         intent(in)    :: ntrnprd(pcols,pver)  ! net precip (liq+ice) production in layer
    real(r8),         intent(in)    :: ntsnprd(pcols,pver)  ! net snow production in layer
    real(r8),         intent(in)    :: tevap(pcols)         ! temperature of surface evaporation
    real(r8),         intent(in)    :: tprec(pcols)         ! temperature of surface precipitation
    real(r8),         intent(in)    :: mflx(pcols)          ! mass   flux for use in check_energy
    real(r8),         intent(in)    :: eflx(pcols)          ! energy flux for use in check_energy
    real(r8),         intent(out)   :: mflx_out(pcols)      ! column (surfce) enthalpy flux from bflx (sanity check)
    real(r8),         intent(out)   :: eflx_out(pcols)      ! column (surfce) enthalpy flux from bflx (sanity check)
    real(r8),         intent(out)   :: ent_tnd (pcols)      ! column-integrated enthalpy tendency
    real(r8),         intent(out)   :: pdel_rf (pcols,pver) ! ratio  old pdel / new pdel
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: klev                 ! Level index
    integer  :: m_cnst               ! Constituent index
    integer  :: m_thermo             ! Thermodynamic constituent indes
    real(r8) :: fdq(pcols)           ! mass adjustment factor
    real(r8) :: utmp(pcols)          ! temp variable for recalculating the initial u values
    real(r8) :: vtmp(pcols)          ! temp variable for recalculating the initial v values
    real(r8) :: te(pcols,pver)       ! conserved energy in layer
    real(r8) :: emce(pcols,pver)     ! total enthalpy - conserved energy in layer
    real(r8) :: cpm(pcols,pver)      ! moist air heat capacity
    real(r8) :: ttsc(pcols,pver)     ! moist air heat capacity
    integer  :: vcoord
    real(r8) :: zvirv(pcols,pver)    ! Local zvir array pointer
    real(r8) :: tot_water(pcols  )   ! total water (initial, present)
    real(r8) :: ps_old(pcols)        ! old surface pressure
    real(r8) :: pdel_new(pcols,pver) ! Layer thickness (pint(k+1) - pint(k))
    real(r8) :: pdot(pcols)          ! total(lagrangian) pressure adjustment
    real(r8) :: pdzp(pcols)          ! pressure work term in press adjustment
    real(r8) :: edot(pcols)          ! advective pressure adjustment
    real(r8) :: uf(pcols), vf(pcols) ! work arrays
    real(r8) :: tp(pcols,pver)       ! work array for T/Tv
    real(r8) :: latent(pcols,pver)   ! work array for Lq
    integer  :: ixnumice, ixnumliq
    integer  :: ixnumsnow, ixnumrain
    real(r8) :: htx_cond(pcols,pver) ! enthalpy tendency due to heat exchange with "condensates"
    real(r8) :: mdq(pcols,pver)      ! total water tendency
    !-----------------------------------------------------------------------

    ! Diagnose boundary enthalpy flux and local heating rates associated to
    ! atmospheric moisture change

    call dme_bflx(lchnk, ncol, &
         state_ps, state_pint, state_zm, state_q, state_pdel, state_phis, state_t, &
         qini, liqini, iceini, tevap, tprec, dt, &
         ntrnprd=ntrnprd, ntsnprd=ntsnprd, &
         mflx=mflx, eflx=eflx, mflx_out=mflx_out, eflx_out=eflx_out, htx_cond=htx_cond, mdq=mdq )

    ! Ajust the dry mass in each layer back to the value of physics input state
    ! Adjust air specific enthalpy accordingly
    ! Diagnose boundary enthalpy flux

    call cnst_get_ind('NUMICE', ixnumice, abort=.false.)
    call cnst_get_ind('NUMLIQ', ixnumliq, abort=.false.)
    call cnst_get_ind('NUMRAI', ixnumrain, abort=.false.)
    call cnst_get_ind('NUMSNO', ixnumsnow, abort=.false.)

    !------------------------------------
    ! initialise adjustment loop
    !------------------------------------

    ! old surface pressure
    ps_old  (:ncol) = state_ps(:ncol)
    state_ps(:ncol) = state_pint(:ncol,1)

    vcoord=vc_dycore
    cpm(:ncol,:) = cp_or_cv_dycore(:ncol,:,lchnk)

    do klev = 1, pver
       tp(:ncol,klev) = state_t(:ncol,klev)  ! TODO - remoe and use state_t instead below
    end do

    call get_conserved_energy(levels_are_moist, &
         1, pver, &
         cpm(:ncol,:), &
         state_t(:ncol,:), state_q(:ncol,:,:) ,state_pdel(:ncol,:), &
         pdel_new(:ncol,:), state_s(:ncol,:), &
         qini=qini(:ncol,:), liqini=liqini(:ncol,:), iceini=iceini(:ncol,:), &
         phis=state_phis(:ncol), gph=state_zm(:ncol,:), &
         U=state_u(:ncol,:), V=state_v(:ncol,:), rairv=rairv(:ncol,:,lchnk), &
         vcoord=vcoord, refstate='liq', &
         flatent=latent(:ncol,:), temce=emce(:ncol,:))

    do klev = 1, pver
       ! Dp'/Dp
       tot_water(:ncol) = 0.0_r8
       do m_cnst=dry_air_species_num+1,thermodynamic_active_species_num
          m_thermo = thermodynamic_active_species_idx(m_cnst)
          tot_water(:ncol) = tot_water(:ncol)+state_q(:ncol,klev,m_thermo)
       end do

       ! new surface pressure
       state_ps(:ncol) = state_ps(:ncol) + state_pdel(:ncol,klev)*(1._r8 + mdq(:ncol,klev))

       ! make all tracers wet
       do m_cnst=1,pcnst
          if (cnst_type(m_cnst) == 'dry') then
             state_q(:ncol,klev,m_cnst) = state_q(:ncol,klev,m_cnst)*(1._r8-tot_water(:ncol))
          end if
       end do
    end do

    ! lagrangian & advective pressure change at top interface
    pdot(:ncol) = 0._r8
    pdzp(:ncol) = 0._r8
    edot(:ncol) = 0._r8

    ! store old enthalpy integral
    ent_tnd(:ncol)=0._r8
    do klev = 1,pver
       ent_tnd(:ncol) = ent_tnd(:ncol) - state_pdel(:ncol,klev)*state_s(:ncol,klev)
    end do

    !------------------------------------
    ! start adjustment loop
    !------------------------------------
    do klev = 1, pver

       ! new Dp (=:Dp")
       pdel_new(:ncol,klev) = state_pdel(:ncol,klev)*(1._r8 + mdq(:ncol,klev))


       ! compute Dp"/Dp
       fdq(:ncol) = pdel_new(:ncol,klev)/state_pdel(:ncol,klev)

       ! wind adjustment increments
       uf(:ncol) = 0.
       vf(:ncol) = 0.

       ! set utmp and vtmp pre-physics u,v from the updated values and the tendencies
       utmp(:ncol) = state_u(:ncol,klev) - dt * tend_dudt(:ncol,klev)
       vtmp(:ncol) = state_v(:ncol,klev) - dt * tend_dvdt(:ncol,klev)

       ! adjust specific enthalpy
       te(:ncol,klev) = 0._r8

       ! lagrangian pressure change *zi at upper interfac
       pdzp(:ncol) = pdot(:ncol)*gravit*state_zi(:ncol,klev)

       ! lagrangian pressure change at next interface
       if (hydrostatic) then
          pdot(:ncol) = pdot(:ncol) + state_pdel(:ncol,klev)*mdq(:ncol,klev)
       end if

       ! layer increment = work (~alpha*dp)
       pdzp(:ncol) = (pdot(:ncol)*gravit*state_zi(:ncol,klev+1)-pdzp(:ncol))/pdel_new(:ncol,klev)

       ! enthalpy change due to mass loss and to hydrost. pressure work in full adjustment
       te(:ncol,klev) = te(:ncol,klev) &
            + state_s(:ncol,klev)/(fdq(:ncol)/(1._r8+mdq(:ncol,klev)))  & ! te *(Dp'/Dp")
            + emce(:ncol,klev)*mdq(:ncol,klev)/fdq(:ncol)               & ! (phi-phis)*dq*(Dp/Dp")
            - pdzp(:ncol)                                         & ! del(g*zm*dp)
            + htx_cond(:ncol,klev)                                     ! EFLX

       ! momentum
       uf(:ncol) = uf(:ncol) +state_u(:ncol,klev)/(fdq(:ncol)/(1._r8+mdq(:ncol,klev)))
       vf(:ncol) = vf(:ncol) +state_v(:ncol,klev)/(fdq(:ncol)/(1._r8+mdq(:ncol,klev)))

       ! adjust constituents to conserve mass in each layer
       do m_cnst = 1, pcnst
          ! store unadjusted q for use in next k
          state_q(:ncol,klev,m_cnst) = state_q(:ncol,klev,m_cnst) / fdq(:ncol)
       end do

       ! adjust L-dependent part of local total enthalpy accordingly
       latent(:ncol,klev) = latent(:ncol,klev)/fdq(:ncol)

       ! adjusted u,v tendencies
       tend_dudt(:ncol,klev) = (uf(:ncol) - utmp(:ncol)) / dt
       tend_dvdt(:ncol,klev) = (vf(:ncol) - vtmp(:ncol)) / dt

       ! store unadjusted u,v for use in next k
       utmp(:ncol) = state_u(:ncol,klev)
       vtmp(:ncol) = state_v(:ncol,klev)

       ! write adjusted u,v
       state_u(:ncol,klev) = uf(:ncol)
       state_v(:ncol,klev) = vf(:ncol)

       ! compute new total pressure variables
       state_pint  (:ncol,klev+1) = state_pint(:ncol,klev  ) + pdel_new(:ncol,klev)
       state_lnpint(:ncol,klev+1) = log(state_pint(:ncol,klev+1))

       ! also update pmid for geopotential
       state_pmid  (:ncol,klev) = .5_r8*(state_pint(:ncol,klev)+state_pint(:ncol,klev+1))
       state_lnpmid(:ncol,klev) = log(state_pmid(:ncol,klev  ))

       pdel_rf(:ncol,klev)=state_pdel(:ncol,klev)/pdel_new(:ncol,klev)
       state_pdel  (:ncol,klev  ) = pdel_new(:ncol,klev)
       state_rpdel (:ncol,klev  ) = 1._r8/state_pdel(:ncol,klev)

    end do

    !------------------------------------
    ! end adjustment loop
    !------------------------------------

    ! make dry tracers dry again
    do klev = 1, pver
       tot_water(:ncol) = 0.0_r8
       do m_cnst=dry_air_species_num+1,thermodynamic_active_species_num
          m_thermo = thermodynamic_active_species_idx(m_cnst)
          tot_water(:ncol) = tot_water(:ncol)+state_q(:ncol,klev,m_thermo)
       end do
       do m_cnst=1,pcnst
          if (cnst_type(m_cnst) == 'dry') then
             state_q(:ncol,klev,m_cnst) = state_q(:ncol,klev,m_cnst)/(1._r8-tot_water(:ncol))
          end if
       end do
    end do

    ! call QNEG3 (cf physics_update)
    do m_cnst = 1, pcnst
       if (m_cnst /= ixnumice  .and.  m_cnst /= ixnumliq .and. &
           m_cnst /= ixnumrain .and.  m_cnst /= ixnumsnow ) then
          call qneg3('dme_adjust', lchnk, ncol, state_psetcols, pver, m_cnst, m_cnst, &
               qmin(m_cnst:m_cnst), state_q(:,1:pver,m_cnst:m_cnst))
       else
          do klev = 1,pver
             state_q(:ncol,klev,m_cnst) = max(1.e-12_r8,state_q(:ncol,klev,m_cnst))
             state_q(:ncol,klev,m_cnst) = min(1.e10_r8,state_q(:ncol,klev,m_cnst))
          end do
       end if
    end do

    call cam_thermo_water_update(state_q(:ncol,:,:), lchnk, ncol, vc_dycore, &
         to_dry_factor=state_pdel(:ncol,:)/state_pdeldry(:ncol,:))

    ttsc(:ncol,:) = cpm(:ncol,:)/cp_or_cv_dycore(:ncol,:,lchnk)
    cpm(:ncol,:) = cp_or_cv_dycore(:ncol,:,lchnk)

    call inv_conserved_energy(levels_are_moist, &
         1, pver, &
         te(:ncol,:), &
         cpm(:ncol,:), &
         state_q(:ncol,:,:), state_pdel(:ncol,:), &
         pdel_new(:ncol,:), tp(:ncol,:), &
         flatent=latent(:ncol,:)*0._r8, &
         phis=state_phis(:ncol), gph=state_zm(:ncol,:), &
         vcoord=vcoord, refstate='liq', &
         U=state_u(:ncol,:), V=state_v(:ncol,:))

    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
       zvirv(:,:) = shr_const_rwv / rairv(:,:,lchnk) - 1._r8
    else
       zvirv(:,:) = zvir
    endif

    ! diagnostics: dme T tendency
    ttsc(:ncol,:) = (tp(:ncol,:) - state_t(:ncol,:))/dt ! &

    call outfld('PTTEND_DME', ttsc, pcols, lchnk)

    ! update ttend and T (cf physics_update)
    tend_dtdt(:ncol,:) = tend_dtdt(:ncol,:) + (tp(:ncol,:) - state_t(:ncol,:))/dt
    state_t(:ncol,:) = tp(:ncol,:)

    ! diagnose total internal enthalpy change
    do klev=1,pver
       ent_tnd(:ncol) = ent_tnd(:ncol) + state_pdel(:ncol,klev)*te(:ncol,klev)
    end do
    ent_tnd(:ncol) = ent_tnd(:ncol)/dt/gravit
    call geopotential_t  (                                                                    &
         state_lnpint, state_lnpmid, state_pint  , state_pmid  , state_pdel  , state_rpdel  , &
         state_t     , state_q(:,:,:), rairv(:,:,lchnk), gravit  , zvirv              , &
         state_zi    , state_zm      , ncol         )

    ! update original dry static energy
    do klev = 1, pver
       state_s(:ncol,klev) = state_t(:ncol,klev  )*cpairv(:ncol,klev,lchnk) &
                        + gravit*state_zm(:ncol,klev) + state_phis(:ncol)
    end do

  contains

    !===============================================================================

    subroutine dme_bflx(lchnk, ncol, &
         state_ps, state_pint, state_zm, state_q, state_pdel, state_phis, state_t, &
         qini, liqini, iceini, tevp, tprc, dt, &
         ntrnprd, ntsnprd, &
         mflx, eflx, mflx_out, eflx_out, htx_cond, mdq)

      !-----------------------------------------------------------------------
      !
      ! Purpose: Diagnose boundary enthalpy flux and local heating rates associated to
      ! atmospheric moisture change
      !
      ! Method
      !   boundary enthalpy flux is *local* total enthalpy (\epsilon dp/g)
      !
      ! Author: Thomas Toniazzo (17.07.21)
      !
      !-----------------------------------------------------------------------

      use air_composition, only: thermodynamic_active_species_idx
      use air_composition, only: thermodynamic_active_species_liq_idx
      use air_composition, only: thermodynamic_active_species_ice_idx
      use air_composition, only: thermodynamic_active_species_num
      use air_composition, only: thermodynamic_active_species_liq_num
      use air_composition, only: thermodynamic_active_species_ice_num
      use air_composition, only: dry_air_species_num
      use air_composition, only: t00a, h00a
      !
      ! Arguments
      !
      integer,          intent(in)    :: lchnk
      integer,          intent(in)    :: ncol
      real(r8),         intent(inout) :: state_ps(:)
      real(r8),         intent(inout) :: state_pint(:,:)
      real(r8),         intent(in)    :: state_zm(:,:)
      real(r8),         intent(in)    :: state_q(:,:,:)
      real(r8),         intent(in)    :: state_pdel(:,:)
      real(r8),         intent(in)    :: state_phis(:)
      real(r8),         intent(in)    :: state_t(:,:)
      real(r8),         intent(in)    :: qini(pcols,pver)     ! initial specific humidity
      real(r8),         intent(in)    :: liqini(pcols,pver)   ! initial total liquid
      real(r8),         intent(in)    :: iceini(pcols,pver)   ! initial total ice
      real(r8),         intent(in)    :: tevp(pcols)          ! temperature of evaporation at bottom of atmo
      real(r8),         intent(in)    :: tprc(pcols)          ! temperature of precipitation at bottom of atmo
      real(r8),         intent(in)    :: dt                   ! model physics timestep
      real(r8),         intent(in)    :: ntrnprd(pcols,pver)  ! net precip (liq+ice) production in layer
      real(r8),         intent(in)    :: ntsnprd(pcols,pver)  ! net snow production in layer
      real(r8),         intent(in)    :: eflx(pcols)          ! boundary enthalpy flux
      real(r8),         intent(in)    :: mflx(pcols)          ! boundary mass     flux
      real(r8),         intent(out)   :: eflx_out(pcols)      ! diagnostic: boundary enthalpy flux
      real(r8),         intent(out)   :: mflx_out(pcols)      ! diagnostic: boundary mass flux
      real(r8),         intent(out)   :: htx_cond(pcols,pver) ! exchange enthalpy increment for dme_adjust
      real(r8),         intent(out)   :: mdq(pcols,pver)      ! total water       increment for dme_adjust

      !---------------------------Local workspace-----------------------------

      integer  :: icol,klev, ixq          ! column, level, water vapor indices
      integer  :: m_cnst                  ! constituent index
      integer  :: m_thermo                ! thermodynamic constituent index
      integer  :: ierr                    ! error flag
      real(r8) :: fdq   (pcols)           ! mass adjustment factor
      real(r8) :: dcvap(pcols)            ! total column vapour change
      real(r8) :: dcliq(pcols)            ! total column liquid change
      real(r8) :: dcice(pcols)            ! total column ice    change
      real(r8) :: dcwat(pcols)            ! total column water  change
      real(r8) :: dcwatr(pcols)           ! residual column water change (in excess of surface flux)
      real(r8) :: tot_water(pcols,2)      ! work array: total water (initial, present)
      real(r8) :: ps_old(pcols)           ! old surface pressure
      real(r8) :: pdel_new(pcols,pver)    ! Layer thickness (pint(k+1) - pint(k))
      real(r8) :: dvap(pcols,pver)        ! wv  mass adjustment
      real(r8) :: dliq(pcols,pver)        ! liq mass adjustment
      real(r8) :: dice(pcols,pver)        ! ice mass adjustment
      real(r8) :: mdqr(pcols,pver)        ! residual mass change (work array)
      real(r8) :: dcqm(pcols)             ! fraction of total/absolute mass change
      real(r8) :: te(pcols,pver)          ! conserved energy in layer
      real(r8) :: emce(pcols,pver)        ! total enthalpy - conserved energy in layer
      real(r8) :: zm(pcols,pver)          ! (phi-phis)/g
      real(r8) :: condeps_ref(pcols,pver) ! local specific enthalpy of "condensates" (mass source)
      real(r8) :: condepss(pcols,pver)    ! specific enthalpy of source reservoir for q changes
      real(r8) :: condepsf(pcols,pver)    ! specific enthalpy of final reservoir for q changes
      real(r8) :: condcp(pcols,pver)      ! species-increment-weighted cp
      real(r8) :: pint_old(pcols,pver+1)  ! work array
      real(r8) :: dummy(pcols,pver)       ! work array
      logical  :: has_dcwat(pcols)
      !-----------------------------------------------------------------------

      ! store old pressure
      ps_old  (:ncol)   = state_ps(:ncol)
      pint_old(:ncol,:) = state_pint(:ncol,:)

      zm(:ncol,:) = state_zm(:ncol,:)

      ! get local specific enthalpy, excluding latent heats
      call get_conserved_energy(levels_are_moist, &
           1, pver, &
           cp_or_cv_dycore(:ncol,:,lchnk) , &
           state_t(:ncol,:) ,state_q(:ncol,:,:) ,state_pdel(:ncol,:), &
           pdel_new(:ncol,:) ,te(:ncol,:) , &
           qini=qini(:ncol,:),liqini=liqini(:ncol,:),iceini=iceini(:ncol,:), &
           phis=state_phis(:ncol) ,gph=zm(:ncol,:), &
           U=state_u(:ncol,:) ,V=state_v(:ncol,:), &
           vcoord=vc_dycore ,refstate='liq', &
           flatent=dummy, temce=emce, rairv=rairv(:ncol,:,lchnk))

      call cnst_get_ind('Q', ixq)

      ! change in water
      dcvap(:ncol)=0._r8
      dcliq(:ncol)=0._r8
      dcice(:ncol)=0._r8
      dcwat(:ncol)=0._r8

      ! heat associated with cp change
      do klev = 1, pver
         ! mass increments Dp'/Dp
         tot_water(:ncol,1) = qini(:ncol,klev)+liqini(:ncol,klev)+iceini(:ncol,klev) !initial total  H2O
         tot_water(:ncol,2) = 0.0_r8
         do m_cnst=dry_air_species_num+1,thermodynamic_active_species_num
            m_thermo = thermodynamic_active_species_idx(m_cnst)
            tot_water(:ncol,2) = tot_water(:ncol,2)+state_q(:ncol,klev,m_thermo)
         end do
         mdq(:ncol,klev)=(tot_water(:ncol,2)-tot_water(:ncol,1))

         dvap(:ncol,klev) = state_q(:ncol,klev,ixq) - qini(:ncol,klev)
         dliq(:ncol,klev) = -liqini(:ncol,klev)
         do m_cnst=1,thermodynamic_active_species_liq_num
            m_thermo = thermodynamic_active_species_liq_idx(m_cnst)
            dliq(:ncol,klev) = dliq(:ncol,klev)+state_q(:ncol,klev,m_thermo)
         end do
         dice(:ncol,klev) = -iceini(:ncol,klev)
         do m_cnst=1,thermodynamic_active_species_ice_num
            m_thermo = thermodynamic_active_species_ice_idx(m_cnst)
            dice(:ncol,klev) = dice(:ncol,klev)+state_q(:ncol,klev,m_thermo)
         end do

         dcvap(:ncol)=dcvap(:ncol)+dvap(:ncol,klev)*state_pdel(:ncol,klev)/gravit
         dcliq(:ncol)=dcliq(:ncol)+dliq(:ncol,klev)*state_pdel(:ncol,klev)/gravit
         dcice(:ncol)=dcice(:ncol)+dice(:ncol,klev)*state_pdel(:ncol,klev)/gravit
         dcwat(:ncol)=dcwat(:ncol)+ mdq(:ncol,klev)*state_pdel(:ncol,klev)/gravit
      end do

      do icol = 1, ncol
         if (dcwat(icol)*mflx(icol) > 0._r8) then
            has_dcwat(icol) = .true.
         else
            has_dcwat(icol) = .false.
         end if
      end do

      ! local specific enthalpy
      do klev = 1, pver
         condeps_ref(:ncol,klev) = te(:ncol,klev) +emce(:ncol,klev)
      end do

      ! exchange specific enthalpies, incremental
      ! we can partition between source and destination
      dcwatr(:ncol) = 0._r8
      do klev=1,pver
         mdqr(:ncol,klev)=mdq(:ncol,klev)+ntrnprd(:ncol,klev)+ntsnprd(:ncol,klev) ! residual: integrates to vapour change

         condcp  (:ncol,klev) = dvap(:ncol,klev)*cpwv + dliq(:ncol,klev)*cpliq+dice (:ncol,klev)*cpice
         condepss(:ncol,klev) = condcp(:ncol,klev)*(state_t(:ncol,klev)-t00a) &
              + (zm(:ncol,klev)*gravit + state_phis(:ncol))*mdq (:ncol,klev)
         condepss(:ncol,klev) = condepss(:ncol,klev)+(cpliq*t00a+h00a)*mdq(:ncol,klev)

         condepsf(:ncol,klev) =-(cpliq*(tprc(:ncol)-t00a  )+state_phis(:ncol))*ntrnprd(:ncol,klev) &
              -(cpice*(tprc(:ncol)-t00a  )+state_phis(:ncol))*ntsnprd(:ncol,klev)
         condepsf(:ncol,klev) = condepsf(:ncol,klev)-(ntrnprd(:ncol,klev)+ntsnprd(:ncol,klev))*(cpliq*t00a+h00a)
         condepsf(:ncol,klev) = condepsf(:ncol,klev) &
              + mdqr(:ncol,klev)*(cpwv*(tevp(:ncol)-t00a)+state_phis(:ncol)+(cpliq*t00a+h00a))

         ! residual column water change: integrates to surface evaporation
         dcwatr(:ncol) = dcwatr(:ncol)  + mdqr(:ncol,klev)*state_pdel(:ncol,klev)/gravit
      end do

      ! partition arbitrarily based on sign match
      ! EFLX_OUT here: work array for part of input EFLX not accounted for by NTSN/RNPR
      eflx_out(:ncol) = eflx(:ncol)*dt
      do klev = 1, pver
         where(.not. has_dcwat(:ncol))
            eflx_out(:ncol) = eflx_out(:ncol) - state_pdel(:ncol,klev)/gravit*condepsf(:ncol,klev)
         elsewhere
            eflx_out(:ncol) = 0._r8
         endwhere
      end do
      dcqm(:ncol)=0._r8
      do klev=1,pver
         where(mdqr(:ncol,klev)*dcwatr(:ncol) > 0._r8)
            dcqm(:ncol)=dcqm(:ncol)+mdqr(:ncol,klev)*state_pdel(:ncol,klev)/gravit
         endwhere
      end do
      where(abs(dcwatr(:ncol)) > rtiny)
         dcqm(:ncol)=dcwatr(:ncol)/dcqm(:ncol)
      elsewhere
         dcqm(:ncol)=0._r8
      endwhere
      do klev=1,pver
         where(mdqr(:ncol,klev)*dcwatr(:ncol) > 0._r8)
            condepsf(:ncol,klev) = condepsf(:ncol,klev)+eflx_out(:ncol)/dcwatr(:ncol)*mdqr(:ncol,klev)*dcqm(:ncol)
         endwhere
         where(has_dcwat(:ncol))
            condepsf(:ncol,klev) = 0._r8
         endwhere
      end do

      ! boundary flux of energy due to mass sources (diagnostic)
      mflx_out(:ncol) = 0._r8
      do klev = 1, pver
         where(.not. has_dcwat(:ncol))
            ! boundary-flux diagnostic associated with water exchanged (column water gained/lost)
            mflx_out(:ncol) = mflx_out(:ncol) + state_pdel(:ncol,klev)/gravit*mdq(:ncol,klev)/dt
         endwhere
      end do

      ! boundary flux of energy due to mass sources (diagnostic)
      eflx_out(:ncol  ) = 0._r8
      do klev = 1, pver
         where(.not. has_dcwat(:ncol))
            ! boundary-flux diagnostic associated with water exchanged (column water gained/lost)
            eflx_out(:ncol) = eflx_out(:ncol) + state_pdel(:ncol,klev)/gravit*condepsf(:ncol,klev)/dt
         endwhere
      end do

      ! make local specific enthalpy incremental
      do klev = 1, pver
         condeps_ref(:ncol,klev) = condeps_ref(:ncol,klev)*mdq(:ncol,klev)
      end do

      ! new surface pressure
      state_ps(:ncol) = state_pint(:ncol,1)
      do klev = 1, pver
         state_ps(:ncol) = state_ps(:ncol) + state_pdel(:ncol,klev)*(1._r8 + mdq(:ncol,klev))
      end do

      ! heat exchange with condensates
      htx_cond(:ncol,:) = 0._r8
      do klev = 1, pver
         do icol=1,ncol
            ! diff. between destination enthalpy and LOCAL     enthalpy (or zero) is distributed in column below
            if (klev == 1) then
               condepsf(icol,klev)=(condepsf(icol,klev)-condepss(icol,klev)) &
                    *state_pdel(icol,klev)/(state_ps(icol)-state_pint(icol,klev))
            else
               condepsf(icol,klev)=(condepsf(icol,klev)-condepss(icol,klev)) &
                    *state_pdel(icol,klev)/(state_ps(icol)-state_pint(icol,klev))   &
                    + condepsf(icol,klev-1)
            endif
            htx_cond(icol,klev) = condepsf(icol,klev) &
                 ! diff. between LOCAL  enthalpy and reference enthalpy is applied locally
                 +(condepss(icol,klev)-condeps_ref(icol,klev))/(1._r8 + mdq(icol,klev))
         end do

         pdel_new(:ncol,klev) = state_pdel(:ncol,klev)*(1._r8 + mdq(:ncol,klev))

         ! compute new total pressure variables
         state_pint(:ncol,klev+1) = state_pint(:ncol,klev  ) + pdel_new(:ncol,klev)

      end do

      ! original pressure
      state_ps  (:ncol)   = ps_old  (:ncol)
      state_pint(:ncol,:) = pint_old(:ncol,:)

    end subroutine dme_bflx

  end subroutine dme_adjust_camnor_run

end module dme_adjust_camnor
