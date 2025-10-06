module dme_adjust_camnor

  use shr_kind_mod,  only: r8 => shr_kind_r8

  implicit none
  private          ! Make default type private to the module

  public :: dme_adjust_camnor_run

  logical :: levels_are_moist=.true. ! TODO: put in namelist?

  ! 5 possibilities (-> = currently reccommended):
  !    1) conserve_dycore=.false. , conserve_physics=.false.  (no conservation = current CAM)
  !    2) conserve_dycore=.true.  , bndry_flx_surface=.true.  (full conservation, bad climatology)
  ! -> 3) conserve_dycore=.true.  , bndry_flx_local=.true.    (requires fixer to match correct surface fluxes)
  !    4) conserve_physics=.true. , bndry_flx_local=.true.    (as 3., plus fixer for atmo energy)
  !    5) conserve_physics=.true. , bndry_flx_surface=.true.  (no advantage wrt option 2)

  ! N.B. old case CONDEPSF=CONDEPS_REF (with CONDEPSS consistent with dycore) not allowed here, since its
  !      rationale isn't clear. For FV, only three of these options (e.g. 1,2,3) are distinct.

  logical, parameter :: conserve_dycore   = .true.
  logical, parameter :: bndry_flx_surface = .true.
  logical, parameter :: conserve_physics  = .not. conserve_dycore
  logical, parameter :: bndry_flx_local   = .not. bndry_flx_surface
  logical, parameter :: conserve = conserve_dycore .or. conserve_physics

  real(r8), parameter :: rtiny = 1e-14_r8    ! a small number (relative to total q change)

  ! set to T to use distribute implied heating over column section to the surface
  logical, parameter  :: l_nolocdcpttend=.true.

  logical :: hydrostatic = .true.

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
    use dyn_tests_utils, only: vc_dycore, vc_physics
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
    integer  :: i,k,m                ! Longitude, level indices
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
    integer  :: m_cnst
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

    if (conserve_dycore) then
       vcoord=vc_dycore
       cpm(:ncol,:) = cp_or_cv_dycore(:ncol,:,lchnk)
    else
       vcoord=vc_physics
       cpm(:ncol,:) = cpairv(:ncol,:,lchnk)
    endif

    do k = 1, pver
       tp(:ncol,k) = state_t(:ncol,k)  ! TODO - remoe and use state_t instead below
    enddo

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

    do k = 1, pver
       ! Dp'/Dp
       tot_water(:ncol) = 0.0_r8
       do m_cnst=dry_air_species_num+1,thermodynamic_active_species_num
          m = thermodynamic_active_species_idx(m_cnst)
          tot_water(:ncol) = tot_water(:ncol)+state_q(:ncol,k,m)
       enddo
       ! new surface pressure
       state_ps(:ncol) = state_ps(:ncol) + state_pdel(:ncol,k)*(1._r8 + mdq(:ncol,k))
       ! make all tracers wet
       do m=1,pcnst
          if (cnst_type(m).eq.'dry') then
             state_q(:ncol,k,m) = state_q(:ncol,k,m)*(1._r8-tot_water(:ncol))
          end if
       enddo
    enddo

    ! lagrangian & advective pressure change at top interface
    pdot(:ncol) = 0._r8
    pdzp(:ncol) = 0._r8
    edot(:ncol) = 0._r8

    ! store old enthalpy integral
    ent_tnd(:ncol)=0._r8
    do k = 1,pver
       ent_tnd(:ncol) = ent_tnd(:ncol) - state_pdel(:ncol,k)*state_s(:ncol,k)
    enddo

    !------------------------------------
    ! start adjustment loop
    !------------------------------------
    do k = 1, pver

       ! new Dp (=:Dp")
       pdel_new(:ncol,k) = state_pdel(:ncol,k)*(1._r8 + mdq(:ncol,k))


       ! compute Dp"/Dp
       fdq(:ncol) = pdel_new(:ncol,k)/state_pdel(:ncol,k)

       ! wind adjustment increments
       uf(:ncol) = 0.
       vf(:ncol) = 0.

       ! set utmp and vtmp pre-physics u,v from the updated values and the tendencies
       utmp(:ncol) = state_u(:ncol,k) - dt * tend_dudt(:ncol,k)
       vtmp(:ncol) = state_v(:ncol,k) - dt * tend_dvdt(:ncol,k)

       ! adjust specific enthalpy
       te(:ncol,k) = 0._r8

       ! lagrangian pressure change *zi at upper interfac
       pdzp(:ncol) = pdot(:ncol)*gravit*state_zi(:ncol,k)

       ! lagrangian pressure change at next interface
       if (hydrostatic) then
          pdot(:ncol) = pdot(:ncol) + state_pdel(:ncol,k)*mdq(:ncol,k)
       end if

       ! layer increment = work (~alpha*dp)
       pdzp(:ncol) = (pdot(:ncol)*gravit*state_zi(:ncol,k+1)-pdzp(:ncol))/pdel_new(:ncol,k)

       ! enthalpy change due to mass loss and to hydrost. pressure work in full adjustment
       te(:ncol,k) = te(:ncol,k) &
            + state_s(:ncol,k)/(fdq(:ncol)/(1._r8+mdq(:ncol,k)))  & ! te *(Dp'/Dp")
            + emce(:ncol,k)*mdq(:ncol,k)/fdq(:ncol)               & ! (phi-phis)*dq*(Dp/Dp")
            - pdzp(:ncol)                                         & ! del(g*zm*dp)
            + htx_cond(:ncol,k)                                     ! EFLX

       ! momentum
       uf(:ncol) = uf(:ncol) +state_u(:ncol,k)/(fdq(:ncol)/(1._r8+mdq(:ncol,k)))
       vf(:ncol) = vf(:ncol) +state_v(:ncol,k)/(fdq(:ncol)/(1._r8+mdq(:ncol,k)))

       ! adjust constituents to conserve mass in each layer
       do m = 1, pcnst
          ! store unadjusted q for use in next k
          state_q(:ncol,k,m) = state_q(:ncol,k,m) / fdq(:ncol)
       end do

       ! adjust L-dependent part of local total enthalpy accordingly
       latent(:ncol,k) = latent(:ncol,k)/fdq(:ncol)

       ! adjusted u,v tendencies
       tend_dudt(:ncol,k) = (uf(:ncol) - utmp(:ncol)) / dt
       tend_dvdt(:ncol,k) = (vf(:ncol) - vtmp(:ncol)) / dt

       ! store unadjusted u,v for use in next k
       utmp(:ncol) = state_u(:ncol,k)
       vtmp(:ncol) = state_v(:ncol,k)

       ! write adjusted u,v
       state_u(:ncol,k) = uf(:ncol)
       state_v(:ncol,k) = vf(:ncol)

       ! compute new total pressure variables
       state_pint  (:ncol,k+1) = state_pint(:ncol,k  ) + pdel_new(:ncol,k)
       state_lnpint(:ncol,k+1) = log(state_pint(:ncol,k+1))

       ! also update pmid for geopotential
       state_pmid  (:ncol,k) = .5_r8*(state_pint(:ncol,k)+state_pint(:ncol,k+1))
       state_lnpmid(:ncol,k) = log(state_pmid(:ncol,k  ))

       pdel_rf(:ncol,k)=state_pdel(:ncol,k)/pdel_new(:ncol,k)
       state_pdel  (:ncol,k  ) = pdel_new(:ncol,k)
       state_rpdel (:ncol,k  ) = 1._r8/state_pdel(:ncol,k)

    end do

    !------------------------------------
    ! end adjustment loop
    !------------------------------------

    ! make dry tracers dry again
    do k = 1, pver
       tot_water(:ncol) = 0.0_r8
       do m_cnst=dry_air_species_num+1,thermodynamic_active_species_num
          m = thermodynamic_active_species_idx(m_cnst)
          tot_water(:ncol) = tot_water(:ncol)+state_q(:ncol,k,m)
       enddo
       do m=1,pcnst
          if (cnst_type(m).eq.'dry') then
             state_q(:ncol,k,m) = state_q(:ncol,k,m)/(1._r8-tot_water(:ncol))
          end if
       enddo
    enddo

    ! call QNEG3 (cf physics_update)
    do m = 1, pcnst
       if (m /= ixnumice  .and.  m /= ixnumliq .and. &
           m /= ixnumrain .and.  m /= ixnumsnow ) then
          call qneg3('dme_adjust', lchnk, ncol, state_psetcols, pver, m, m, qmin(m:m), state_q(:,1:pver,m:m))
       else
          do k = 1,pver
             state_q(:ncol,k,m) = max(1.e-12_r8,state_q(:ncol,k,m))
             state_q(:ncol,k,m) = min(1.e10_r8,state_q(:ncol,k,m))
          end do
       end if
    enddo

    if (conserve_dycore) then
       call cam_thermo_water_update(state_q(:ncol,:,:), lchnk, ncol, vc_dycore, &
            to_dry_factor=state_pdel(:ncol,:)/state_pdeldry(:ncol,:))
       ttsc(:ncol,:)=cpm(:ncol,:)/cp_or_cv_dycore(:ncol,:,lchnk)
       cpm(:ncol,:)=cp_or_cv_dycore(:ncol,:,lchnk)
    endif

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
    do k=1,pver
       ent_tnd(:ncol) = ent_tnd(:ncol) + state_pdel(:ncol,k)*te(:ncol,k)
    enddo
    ent_tnd(:ncol) = ent_tnd(:ncol)/dt/gravit
    call geopotential_t  (                                                                    &
         state_lnpint, state_lnpmid, state_pint  , state_pmid  , state_pdel  , state_rpdel  , &
         state_t     , state_q(:,:,:), rairv(:,:,lchnk), gravit  , zvirv              , &
         state_zi    , state_zm      , ncol         )

    ! update original dry static energy
    do k = 1, pver
       state_s(:ncol,k) = state_t(:ncol,k  )*cpairv(:ncol,k,lchnk) &
                        + gravit*state_zm(:ncol,k) + state_phis(:ncol)
    enddo

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
      !        1. boundary enthalpy flux is *local* total enthalpy (\epsilon dp/g)
      !        2. same as 1., but with different specific enthalpy of boundary mass exchange,
      !          CONDEPS, and a matching heat exchange betweeen air and condensated
      !          = (\epsilon - CONDEPS) dp/g (sign is for a heat source for air).
      !     Choice 2. is taken with dme_ ohf_adjust=.true. For CONDEPS then the following
      !     choice is made: CONDEPS = cpcond *ocnfrac *SST + cpcond *(1-ocnfrac) *TS
      !     cpcond is a parameter representing the heat capacity of the condensate phase.
      !     The heating rates and enthalpy boundary fluxes are not applied here,
      !     they are intended to be passed to dme_adjust.
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
      real(r8),         intent(out)   :: mflx_out(pcols)      ! diagnostic: boundary enthalpy flux
      real(r8),         intent(out)   :: htx_cond(pcols,pver) ! exchange enthalpy increment for dme_adjust
      real(r8),         intent(out)   :: mdq(pcols,pver)      ! total water       increment for dme_adjust

      !---------------------------Local workspace-----------------------------

      integer  :: i,k,m, ixq              ! Longitude, level indices
      integer  :: ierr                    ! error flag
      real(r8) :: fdq   (pcols)           ! mass adjustment factor
      real(r8) :: dcvap(pcols)            ! total column vapour change
      real(r8) :: dcliq(pcols)            ! total column liquid change
      real(r8) :: dcice(pcols)            ! total column ice    change
      real(r8) :: dcwat(pcols)            ! total column water  change
      real(r8) :: dcwatr(pcols)           ! residual column water change (in excess of surface flux)
      real(r8) :: tot_water(pcols,2)      ! work array: total water (initial, present)
      integer  :: m_cnst                  ! index
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
      integer  :: is_invalid(pcols)
      !-----------------------------------------------------------------------

      ! store old pressure
      ps_old  (:ncol)   = state_ps(:ncol)
      pint_old(:ncol,:) = state_pint(:ncol,:)

      zm(:ncol,:) = state_zm(:ncol,:)

      ! get local specific enthalpy, excluding latent heats
      if (conserve_dycore) then
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
      else
         call get_conserved_energy(levels_are_moist, &
              1, pver, &
              cpairv(:ncol,:,lchnk) , &
              state_t(:ncol,:) ,state_q(:ncol,:,:) ,state_pdel(:ncol,:), &
              pdel_new(:ncol,:) ,te(:ncol,:), &
              qini=qini(:ncol,:),liqini=liqini(:ncol,:),iceini=iceini(:ncol,:), &
              phis=state_phis(:ncol), gph=zm(:ncol,:), &
              U=state_u(:ncol,:) ,V=state_v(:ncol,:), &
              refstate='liq', &
              flatent=dummy, temce=emce, rairv=rairv(:ncol,:,lchnk))
      endif

      call cnst_get_ind('Q', ixq)

      ! change in water
      dcvap(:ncol)=0._r8
      dcliq(:ncol)=0._r8
      dcice(:ncol)=0._r8
      dcwat(:ncol)=0._r8

      ! heat associated with cp change
      do k = 1, pver
         ! mass increments Dp'/Dp
         tot_water(:ncol,1) = qini(:ncol,k)+liqini(:ncol,k)+iceini(:ncol,k) !initial total  H2O
         tot_water(:ncol,2) = 0.0_r8
         do m_cnst=dry_air_species_num+1,thermodynamic_active_species_num
            m = thermodynamic_active_species_idx(m_cnst)
            tot_water(:ncol,2) = tot_water(:ncol,2)+state_q(:ncol,k,m)
         end do
         mdq(:ncol,k)=(tot_water(:ncol,2)-tot_water(:ncol,1))

         dvap(:ncol,k) = state_q(:ncol,k,ixq) - qini(:ncol,k)
         dliq(:ncol,k) = -liqini(:ncol,k)
         do m_cnst=1,thermodynamic_active_species_liq_num
            m = thermodynamic_active_species_liq_idx(m_cnst)
            dliq(:ncol,k) = dliq(:ncol,k)+state_q(:ncol,k,m)
         end do
         dice(:ncol,k) = -iceini(:ncol,k)
         do m_cnst=1,thermodynamic_active_species_ice_num
            m = thermodynamic_active_species_ice_idx(m_cnst)
            dice(:ncol,k) = dice(:ncol,k)+state_q(:ncol,k,m)
         end do

         dcvap(:ncol)=dcvap(:ncol)+dvap(:ncol,k)*state_pdel(:ncol,k)/gravit
         dcliq(:ncol)=dcliq(:ncol)+dliq(:ncol,k)*state_pdel(:ncol,k)/gravit
         dcice(:ncol)=dcice(:ncol)+dice(:ncol,k)*state_pdel(:ncol,k)/gravit
         dcwat(:ncol)=dcwat(:ncol)+ mdq(:ncol,k)*state_pdel(:ncol,k)/gravit
      end do

      is_invalid(:ncol)=0
      where(dcwat(:ncol)*mflx(:ncol) .gt. 0._r8)
         is_invalid(:ncol) = 1
      endwhere

      ! local specific enthalpy
      if (conserve)  then
         do k = 1, pver
            condeps_ref(:ncol,k) = te(:ncol,k) +emce(:ncol,k)
         enddo
      else
         condeps_ref(:ncol,:) = 0._r8
      endif

      ! exchange specific enthalpies, incremental
      if (conserve) then ! we can partition between source and destination
         dcwatr(:ncol) = 0._r8
         do k=1,pver
            mdqr(:ncol,k)=mdq(:ncol,k)+ntrnprd(:ncol,k)+ntsnprd(:ncol,k) ! residual: integrates to vapour change

            if (conserve_physics .or. .not. l_nolocdcpttend)  then
               condepss(:ncol,k) = condeps_ref(:ncol,k)*mdq (:ncol,k)
            else if (conserve_dycore) then
               condcp  (:ncol,k) = dvap  (:ncol,k)*cpwv +dliq (:ncol,k)*cpliq+dice (:ncol,k)*cpice
               condepss(:ncol,k) = condcp(:ncol,k)*(state_t(:ncol,k)-t00a) &
                    +(zm(:ncol,k)*gravit+state_phis(:ncol))*mdq (:ncol,k)
               condepss(:ncol,k) = condepss(:ncol,k)+(cpliq*t00a+h00a)*mdq (:ncol,k)
            endif

            if (bndry_flx_surface) then
               condepsf(:ncol,k) =-(cpliq*(tprc(:ncol)-t00a  )+state_phis(:ncol))*ntrnprd(:ncol,k) &
                    -(cpice*(tprc(:ncol)-t00a  )+state_phis(:ncol))*ntsnprd(:ncol,k)
               condepsf(:ncol,k) = condepsf(:ncol,k)-(ntrnprd(:ncol,k)+ntsnprd(:ncol,k))*(cpliq*t00a+h00a)
               condepsf(:ncol,k) = condepsf(:ncol,k)+mdqr(:ncol,k)*(cpwv*(tevp(:ncol)-t00a)+state_phis(:ncol)+(cpliq*t00a+h00a))
            else if (bndry_flx_local)   then
               if (conserve_dycore)  then
                  condepsf(:ncol,k) = -(cpliq*(state_t(:ncol,k)-t00a  )+zm(:ncol,k)*gravit+state_phis(:ncol))*ntrnprd(:ncol,k) &
                       -(cpice*(state_t(:ncol,k)-t00a  )+zm(:ncol,k)*gravit+state_phis(:ncol))*ntsnprd(:ncol,k)
                  condepsf(:ncol,k) = condepsf(:ncol,k) - &
                       (ntrnprd(:ncol,k)+ntsnprd(:ncol,k))*(cpliq*t00a+h00a)
                  condepsf(:ncol,k) = condepsf(:ncol,k) + &
                       mdqr(:ncol,k)*(cpwv*(state_t(:ncol,k)-t00a)+zm(:ncol,k)*gravit+state_phis(:ncol)+(cpliq*t00a+h00a))
               else if (conserve_physics) then
                  condepsf(:ncol,k) =-condeps_ref(:ncol,k)*(ntrnprd(:ncol,k)+ntsnprd(:ncol,k))
                  condepsf(:ncol,k) = condepsf(:ncol,k)+condeps_ref(:ncol,k)*mdqr(:ncol,k)
               endif
            endif

            ! residual column water change: integrates to surface evaporation
            dcwatr  (:ncol)   = dcwatr(:ncol)  + mdqr(:ncol,k)*state_pdel(:ncol,k)/gravit
         enddo
      else
         mdqr    (:ncol,:)=mdq  (:ncol,:)
         dcwatr  (:ncol)  =dcwat(:ncol)
         condepsf(:ncol,:)=0._r8
         condepss(:ncol,:)=0._r8
         do k=1,pver
            if      (conserve_physics.or..not.l_nolocdcpttend)  then
               condepss(:ncol,k) = condeps_ref(:ncol,k)*mdq(:ncol,k)
            else if (conserve_dycore ) then
               condcp  (:ncol,k) = dvap (:ncol,k)*cpwv +dliq(:ncol,k)*cpliq+dice(:ncol,k)*cpice
               condepss(:ncol,k) = condcp(:ncol,k)*(state_t(:ncol,k)-t00a) &
                    +(zm(:ncol,k)*gravit+state_phis(:ncol))*mdq(:ncol,k)
               condepss(:ncol,k) = condepss(:ncol,k)+(cpliq*t00a+h00a)*mdq(:ncol,k)
            endif
            if      (bndry_flx_surface) then
               condcp  (:ncol,k) = dvap (:ncol,k)*cpwv +dliq(:ncol,k)*cpliq+dice(:ncol,k)*cpice
               condepsf(:ncol,k) = condcp(:ncol,k)*&
                    (tprc(:ncol)-t00a)+state_phis(:ncol)*mdq(:ncol,k)+dvap(:ncol,k)*cpwv*(tevp(:ncol)-tprc(:ncol))
               condepsf(:ncol,k) = condepsf(:ncol,k)+(cpliq*t00a+h00a)*mdq(:ncol,k)
            else if (bndry_flx_local)   then
               condepsf(:ncol,k) = condepss(:ncol,k)
               if (conserve_dycore .and.l_nolocdcpttend) &
                    condepsf(:ncol,k) = condepsf(:ncol,k)+((cpliq-cpair)*t00a+h00a)*mdq(:ncol,k)
            endif
         enddo
      endif

      if (conserve) then ! partition arbitrarily based on sign match
         ! EFLX_OUT here: work array for part of input EFLX not accounted for by NTSN/RNPR
         eflx_out(:ncol) = eflx(:ncol)*dt
         do k = 1, pver
            where(is_invalid(:ncol).eq.0)
               eflx_out(:ncol) = eflx_out(:ncol) - state_pdel(:ncol,k)/gravit*condepsf(:ncol,k)
            elsewhere
               eflx_out(:ncol) = 0._r8
            endwhere
         enddo
         dcqm(:ncol)=0._r8
         do k=1,pver
            where(mdqr(:ncol,k)*dcwatr(:ncol).gt.0._r8)
               dcqm(:ncol)=dcqm(:ncol)+mdqr(:ncol,k)*state_pdel(:ncol,k)/gravit
            endwhere
         enddo
         where(abs(dcwatr(:ncol)).gt.rtiny)
            dcqm(:ncol)=dcwatr(:ncol)/dcqm(:ncol)
         elsewhere
            dcqm(:ncol)=0._r8
         endwhere
         do k=1,pver
            where(mdqr(:ncol,k)*dcwatr(:ncol).gt.0._r8)
               condepsf(:ncol,k) = condepsf(:ncol,k)+eflx_out(:ncol)/dcwatr(:ncol)*mdqr(:ncol,k)*dcqm(:ncol)
            endwhere
            where(is_invalid(:ncol).eq.1)
               condepsf(:ncol,k) = 0._r8
            endwhere
         enddo
      endif

      ! boundary flux of energy due to mass sources (diagnostic)
      mflx_out(:ncol) = 0._r8
      do k = 1, pver
         where(is_invalid(:ncol).eq.0)
            ! boundary-flux diagnostic associated with water exchanged (column water gained/lost)
            mflx_out(:ncol) = mflx_out(:ncol) + state_pdel(:ncol,k)/gravit*mdq(:ncol,k)/dt
         endwhere
      enddo

      ! boundary flux of energy due to mass sources (diagnostic)
      eflx_out(:ncol  ) = 0._r8
      do k = 1, pver
         where(is_invalid(:ncol).eq.0)
            ! boundary-flux diagnostic associated with water exchanged (column water gained/lost)
            eflx_out(:ncol) = eflx_out(:ncol) + state_pdel(:ncol,k)/gravit*condepsf(:ncol,k)/dt
         endwhere
      enddo

      ! make local specific enthalpy incremental
      if (conserve)  then
         do k = 1, pver
            condeps_ref(:ncol,k) = condeps_ref(:ncol,k)*mdq(:ncol,k)
         enddo
      endif

      ! new surface pressure
      state_ps(:ncol) = state_pint(:ncol,1)
      do k = 1, pver
         state_ps(:ncol) = state_ps(:ncol) + state_pdel(:ncol,k)*(1._r8 + mdq(:ncol,k))
      end do

      ! heat exchange with condensates
      htx_cond(:ncol,:) = 0._r8
      do k = 1, pver
         do i=1,ncol
            if(l_nolocdcpttend)then
               ! diff. between destination enthalpy and LOCAL     enthalpy (or zero) is distributed in column below
               if (k.eq.1) then
                  condepsf(i,k)=(condepsf(i,k)-condepss(i,k)) &
                       *state_pdel(i,k)/(state_ps(i)-state_pint(i,k))
               else
                  condepsf(i,k)=(condepsf(i,k)-condepss(i,k)) &
                       *state_pdel(i,k)/(state_ps(i)-state_pint(i,k))   &
                       +condepsf(i,k-1)
               endif
            else
               condepsf(i,k)=(condepsf(i,k)-condepss(i,k))/(1._r8+mdq(i,k))
            endif
            htx_cond(i,k) = condepsf(i,k) &
                 ! diff. between LOCAL  enthalpy and reference enthalpy is applied locally
                 +(condepss(i,k)-condeps_ref(i,k))/(1._r8 + mdq(i,k))
         enddo

         pdel_new(:ncol,k) = state_pdel(:ncol,k)*(1._r8 + mdq(:ncol,k))

         ! compute new total pressure variables
         state_pint(:ncol,k+1) = state_pint(:ncol,k  ) + pdel_new(:ncol,k)

      end do

      ! original pressure
      state_ps  (:ncol)   = ps_old  (:ncol)
      state_pint(:ncol,:) = pint_old(:ncol,:)

    end subroutine dme_bflx

  end subroutine dme_adjust_camnor_run

end module dme_adjust_camnor
