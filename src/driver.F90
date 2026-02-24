! test310_driver.F90
! changed on 09/02/19 for namelist input
!
! this does a test of mam4 (with marine organics and coarse-mode
!    carbonaceous species) that exercisess modal_aero_calcsize,
!    modal_aero_wateruptake, and modal_aero_amicphys (all the amicphys processes),
!    with cloud fraction = 0.
!
! there is no benchmark solution for this test case.
!
! this test driver takes namelist input for initial condition
!
!-------------------------------------------------------------------------------

      module driver

      use precision_mod, only: r8 => f8,  fp, f8
      use mam_utils, only: masterproc, endrun, iulog, pcols, pver,plev, begchunk, endchunk,&
                           mdo_gaschem, mdo_cloudchem, mdo_coldstart, mdo_mambox, &
                           mdo_gasaerexch, mdo_rename, mdo_newnuc, mdo_coag

      use constituents, only: pcnst, cnst_name, cnst_get_ind
      use modal_aero_data, only: ntot_amode
      use physics_buffer, only: physics_buffer_desc
      use physics_types, only : physics_state, physics_ptend
      use mam_opt, only : mam_optics_diagnostics, mamoptdiag
#if(defined USE_NC4)
      use netcdf
#endif
      implicit none

      public

! similaire a  GeosCore mam_driv
!      integer :: mdo_gaschem, mdo_cloudchem
!      integer :: mdo_gasaerexch, mdo_rename, mdo_newnuc, mdo_coag
      integer:: loffset, lchnk
      real(r8) :: deltat

      type(physics_buffer_desc), pointer :: pbuf(:)
      type(physics_state) :: physta
      type(physics_ptend) :: ptend


! propres au driver de mambox
      integer, parameter :: lun_outfld = 90
      integer :: mopt_aero_comp, mopt_aero_load, mopt_ait_size
      integer :: mopt_h2so4_uptake
      integer :: i_cldy_sameas_clear
      integer :: iwrite3x_species_flagaa, iwrite3x_units_flagaa
      integer :: iwrite4x_heading_flagbb
      real(r8) :: xopt_cloudf

      ! in the multiple nbc/npoa code, the following are in modal_aero_data
      integer :: lptr_bca_a_amode(ntot_amode) = -999888777
      integer :: lptr_poma_a_amode(ntot_amode) = -999888777

      integer :: species_class(pcnst) = -1

      real(r8) :: tmin,tmax,rhmin,rhmax ! namelist variables

      contains


!-------------------------------------------------------------------------------
      subroutine cambox_main

      use modal_aero_data,         only: ntot_amode
      use physics_buffer, only: physics_buffer_desc

      use modal_aero_data, only: &
         lmassptrcw_amode, nspec_amode, numptrcw_amode, &
         qqcw_get_field

      use modal_aero_initialize_data, only: MAM_init_basics, MAM_ALLOCATE, MAM_cold_start
      use mam_opt, only : mam_init_opt 
                                                                                            
      integer :: nstop
      real*8  :: deltat

      pcols = 1
      pver = 1

      plev = pver




      masterproc = .true.


      write(*,'(/a)') '*** Hello from MAIN ***'

      write(*,'(/a)') '*** main calling MAM_init_basics'
      call MAM_init_basics(pbuf)
      write(*,'(/a)') '*** main call MAM_allocate'
      call MAM_ALLOCATE (physta,ptend )
      write(*,'(/a)') '*** main call MAM_init_opt'
      call MAM_INIT_OPT()

      call gcmambox_init_run (nstop,deltat)

      write(*,'(/a)') '*** main calling cambox_do_run'

      call gcmambox_do_run(nstop,deltat)

      end subroutine cambox_main


!-------------------------------------------------------------------------------
      subroutine gcmambox_init_run( nstop, deltat )

      use chem_mods, only: adv_mass, gas_pcnst, imozart
      use physconst, only: pi, epsilo, latvap, latice, &
                           rh2o, cpair, tmelt, mwdry, r_universal
      use wv_saturation, only: qsat, gestbl

      use modal_aero_data
      use modal_aero_amicphys, only: &
           gaexch_h2so4_uptake_optaa, newnuc_h2so4_conc_optaa, mosaic

      use  modal_aero_initialize_data, only: MAM_init_basics, MAM_ALLOCATE, MAM_cold_start
      implicit none
      integer,  intent(out  ) :: nstop
      real(r8), intent(out  ) :: deltat

      integer :: i
      integer :: k
      integer :: l, ll, loffset, lun
      integer :: l_nh3g, l_so2g, l_soag, l_hno3g, l_hclg, l_h2so4g
      integer :: l_num_a1, l_num_a2, l_nh4_a1, l_nh4_a2, &
                 l_so4_a1, l_so4_a2, l_soa_a1, l_soa_a2
      integer :: l_numa, l_so4a, l_nh4a, l_soaa, l_poma, l_bcxa, l_ncla, &
                 l_dsta, l_no3a, l_clxa, l_caxa, l_co3a, l_moma
      integer :: mode123_empty
      integer :: mopt_aero_loadaa, mopt_aero_loadbb
      integer :: n, nacc, nait

      logical :: ip

      real(r8) :: ev_sat(pcols,pver)
      real(r8) :: qv_sat(pcols,pver)
      real(r8) :: tmn, tmx, trice
      real(r8) :: tmpa, tmpq




      iwrite3x_species_flagaa = 1
      iwrite3x_units_flagaa   = 10
      iwrite4x_heading_flagbb = 1
      mopt_ait_size       = 2
      xopt_cloudf         = 0.6_r8
      i_cldy_sameas_clear = 0

      write(*,'(/a)') '*** main call MAM_cold_start'
      call MAM_cold_start (physta,nstop=nstop, deltat=deltat, tmin = tmin, tmax =tmax, rhmin =rhmin, rhmax =rhmax)!


! require gestbl to build saturation vapor pressure table.

      call  qsat( physta%t(1:pcols,1:pver), physta%pmid(1:pcols,1:pver), &
                  ev_sat(1:pcols,1:pver), qv_sat(1:pcols,1:pver) )
!test !
      qv_sat = 4.35E-3
      physta%qv(:,:) = physta%relhum(:,:)*qv_sat(:,:)
      physta%q(:,:,1) = physta%qv(:,:)
       
      return
      end subroutine gcmambox_init_run


!-------------------------------------------------------------------------------
      subroutine gcmambox_do_run( nstop, deltat)

      use chem_mods, only: adv_mass, gas_pcnst, imozart
      use physconst, only: mwdry,rga
      use physics_types, only: physics_state, physics_ptend
      use physics_buffer, only: physics_buffer_desc, pbuf_get_chunk
      use mam_utils, only: l_h2so4g, l_nh3g, l_so2g, l_hno3g, l_hclg, l_soag
      use modal_aero_data
      use modal_aero_calcsize, only: modal_aero_calcsize_sub
      use modal_aero_amicphys, only: modal_aero_amicphys_intr, &
          gaexch_h2so4_uptake_optaa, newnuc_h2so4_conc_optaa, mosaic
      use modal_aero_wateruptake, only: modal_aero_wateruptake_dr, load_pbuf, unload_pbuf
      use gaschem_simple, only: gaschem_simple_sub
      use cloudchem_simple, only: cloudchem_simple_sub
      use wv_saturation, only : qsat
      use radconstants, only : nswbands, nlwbands
      use mam_opt , only : mam_aero_sw, mam_aero_lw, mamoptdiag 

      implicit none

      integer,  intent(in   ) :: nstop

      real(r8), intent(in   ) :: deltat

!      type(physics_buffer_desc), pointer :: pbuf2d(:,:)  ! full physics buffer

      integer, parameter :: nqtendbb = 4
      integer, parameter :: iqtend_cond = 1
      integer, parameter :: iqtend_rnam = 2
      integer, parameter :: iqtend_nnuc = 3
      integer, parameter :: iqtend_coag = 4
      integer, parameter :: nqqcwtendbb = 1
      integer, parameter :: iqqcwtend_rnam = 1

      integer :: i, icalcaer_flag, iwaterup_flag
      integer :: istep
      integer :: itmpa, itmpb
      integer :: k
      integer :: l, l2, ll
!      integer :: l_h2so4g, l_nh3g, l_so2g, l_hno3g, l_hclg, l_soag
      integer :: l_num_a1, l_nh4_a1, l_so4_a1
      integer :: l_num_a2, l_nh4_a2, l_so4_a2
      integer :: lmz_h2so4g, lmz_nh3g, lmz_so2g, &
                 lmz_hno3g, lmz_hclg, lmz_soag
      integer :: lmz_num_a1, lmz_nh4_a1, lmz_so4_a1
      integer :: lmz_num_a2, lmz_nh4_a2, lmz_so4_a2
      integer :: lchnk, loffset, lun
      integer :: latndx(pcols), lonndx(pcols)
      integer :: n, nacc, nait, nstep

      logical :: aero_mmr_flag
      logical :: h2o_mmr_flag
      logical :: dotend(pcnst)

      real(r8) ::temp, relh

      real(r8) :: ev_sat(pcols,pver), qv_sat(pcols,pver)

      real(r8) :: cld_ncol(pver)
      real(r8) :: del_h2so4_aeruptk(pcols,pver)
      real(r8) :: del_h2so4_gasprod(pcols,pver)
      real(r8) :: dqdt(pcols,pver,pcnst)        ! Tracer MR tendency array
      real(r8) :: dvmrdt_bb(pcols,pver,gas_pcnst,nqtendbb)   ! mixing ratio changes
      real(r8) :: dvmrcwdt_bb(pcols,pver,gas_pcnst,nqqcwtendbb) ! mixing ratio changes
      real(r8) :: dvmrdt_cond(pcols,pver,gas_pcnst)   ! mixing ratio changes from renaming
      real(r8) :: dvmrcwdt_cond(pcols,pver,gas_pcnst) ! mixing ratio changes from renaming
      real(r8) :: dvmrdt_nnuc(pcols,pver,gas_pcnst)   ! mixing ratio changes from renaming
      real(r8) :: dvmrcwdt_nnuc(pcols,pver,gas_pcnst) ! mixing ratio changes from renaming
      real(r8) :: dvmrdt_coag(pcols,pver,gas_pcnst)   ! mixing ratio changes from renaming
      real(r8) :: dvmrcwdt_coag(pcols,pver,gas_pcnst) ! mixing ratio changes from renaming
      real(r8) :: dvmrdt_rnam(pcols,pver,gas_pcnst)   ! mixing ratio changes from renaming
      real(r8) :: dvmrcwdt_rnam(pcols,pver,gas_pcnst) ! mixing ratio changes from renaming
      real(r8) :: h2so4_pre_gaschem(pcols,pver) ! grid-avg h2so4(g) mix ratio before gas chem (mol/mol)
      real(r8) :: h2so4_aft_gaschem(pcols,pver) ! grid-avg h2so4(g) mix ratio after  gas chem (mol/mol)
      real(r8) :: h2so4_clear_avg(  pcols,pver) ! average clear sub-area h2so4(g) mix ratio (mol/mol)
      real(r8) :: h2so4_clear_fin(  pcols,pver) ! final   clear sub-area h2so4(g) mix ratio (mol/mol)
      real(r8) :: tau_gaschem_simple(pcols,pver)
      real(r8) :: tmpa, tmpb, tmpc
      real(r8) :: tmpveca(999)
      real(r8) :: told, tnew
      real(r8) :: uptkrate_h2so4(   pcols,pver) ! h2so4(g) uptake (by aerosols) rate (1/s)
      real(r8) :: vmr(pcols,pver,gas_pcnst)     ! gas & aerosol volume mixing ratios
      real(r8) :: vmr_svaa(pcols,pver,gas_pcnst)
      real(r8) :: vmr_svbb(pcols,pver,gas_pcnst)
      real(r8) :: vmr_svcc(pcols,pver,gas_pcnst)
      real(r8) :: vmr_svdd(pcols,pver,gas_pcnst)
      real(r8) :: vmr_svee(pcols,pver,gas_pcnst)
      real(r8) :: vmrcw(pcols,pver,gas_pcnst)   ! gas & aerosol volume mixing ratios
      real(r8) :: vmrcw_svaa(pcols,pver,gas_pcnst)
      real(r8) :: vmrcw_svbb(pcols,pver,gas_pcnst)
      real(r8) :: vmrcw_svcc(pcols,pver,gas_pcnst)
      real(r8) :: vmrcw_svdd(pcols,pver,gas_pcnst)
      real(r8) :: vmrcw_svee(pcols,pver,gas_pcnst)

!optical properties 

      real(r8)  :: tauxar(pcols,pver,nswbands)  ! aerosol extinction optical depth
      real(r8)  :: wa(pcols,pver,nswbands)      ! aerosol single scattering albedo * tau
      real(r8)  :: ga(pcols,pver,nswbands)      ! aerosol asymmetry parameter * wa
      real(r8)  :: fa(pcols,pver,nswbands)      ! aerosol forward scattered fraction * ga
      real(r8)  :: taux_lw(pcols,pver,nlwbands) ! aerosol LW absorption optical depth

!
! for netcdf file
!
      character (len = *), parameter :: FILE_NAME = "mam_output.nc"
      integer :: ncid, nstep_dimid, mode_dimid
      integer :: error
      integer :: dimids(2), varid(60)
      character (8) :: date
      real(r8), dimension(nstop,ntot_amode) :: tmp_dgn_a, &
                               tmp_dgn_awet, tmp_num_aer
real(r8) :: tmp_so4_aer(nstop, ntot_amode)
 real(r8) :: tmp_nh4_aer(nstop, ntot_amode)
real(r8) :: tmp_no3_aer(nstop, ntot_amode)
real(r8) :: tmp_soa_aer(nstop, ntot_amode)
real(r8) :: tmp_pom_aer(nstop, ntot_amode)
real(r8) :: tmp_bc_aer(nstop, ntot_amode)
real(r8) :: tmp_dust_aer(nstop, ntot_amode)
real(r8) :: tmp_co3_aer(nstop, ntot_amode)
real(r8) :: tmp_nacl_aer(nstop, ntot_amode)
real(r8) :: tmp_ca_aer(nstop, ntot_amode)
real(r8) :: tmp_cl_aer(nstop, ntot_amode)
real(r8) :: tmp_wat_aer(nstop, ntot_amode)
real(r8) :: tmp_hygro_aer(nstop, ntot_amode)
! Optical diagnostics - species optical depths (AOD)
real(r8), dimension(nstop,ntot_amode) :: tmp_aod_sulfate, tmp_aod_bc, &
                                         tmp_aod_pom, tmp_aod_soa, &
                                         tmp_aod_dust, tmp_aod_seasalt, &
                                         tmp_aod_mode

! Optical diagnostics - species single scattering albedo
real(r8), dimension(nstop,ntot_amode) :: tmp_ssa_sulfate, tmp_ssa_bc, &
                                         tmp_ssa_pom, tmp_ssa_soa, &
                                         tmp_ssa_dust, tmp_ssa_seasalt, &
                                         tmp_ssa_mode

! Optical diagnostics - species asymmetry parameter
real(r8), dimension(nstop,ntot_amode) :: tmp_asm_sulfate, tmp_asm_bc, &
                                         tmp_asm_pom, tmp_asm_soa, &
                                         tmp_asm_dust, tmp_asm_seasalt, &
                                         tmp_asm_mode

! LW optical diagnostics - per-mode total LW absorption optical depth
! and per-species contributions.  Reference band: ilw=9 (atmospheric window,
! ~8-12 um, 1000-1390 cm-1), analogous to isw=10 (550 nm) for the SW.
integer, parameter :: ilw_diag = 9   ! LW diagnostic band index
real(r8), dimension(nstop,ntot_amode) :: tmp_aod_lw_mode, &
                                         tmp_aod_lw_sulfate, &
                                         tmp_aod_lw_bc, &
                                         tmp_aod_lw_pom, &
                                         tmp_aod_lw_soa, &
                                         tmp_aod_lw_dust, &
                                         tmp_aod_lw_seasalt

#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
! All three MOSAIC species
real(r8), dimension(nstop,ntot_amode) :: tmp_aod_mom, tmp_aod_no3, tmp_aod_nh4
real(r8), dimension(nstop,ntot_amode) :: tmp_ssa_mom, tmp_ssa_no3, tmp_ssa_nh4
real(r8), dimension(nstop,ntot_amode) :: tmp_asm_mom, tmp_asm_no3, tmp_asm_nh4
#elif ( defined MODAL_AERO_4MODE_MOM )
! Only MOM (marine organic matter)
real(r8), dimension(nstop,ntot_amode) :: tmp_aod_mom
real(r8), dimension(nstop,ntot_amode) :: tmp_ssa_mom
real(r8), dimension(nstop,ntot_amode) :: tmp_asm_mom
#endif

real(r8), dimension(nstop)            :: tmp_h2so4, tmp_hno3, tmp_nh3, &
                                         tmp_soag, tmp_hcl, tmp_so2 ,tmp_relh, tmp_temp
real(r8), dimension(nstop,ntot_amode) :: qtend_cond_aging_so4, &
                                               qtend_cond_aging_soa, &
                                               qtend_rename_so4, &
                                               qtend_rename_soa, &
                                               qtend_newnuc_so4, &
                                               qtend_newnuc_soa, &
                                               qtend_coag_so4, &
                                               qtend_coag_soa
      real(r8), dimension(nstop)            :: qtend_cond_aging_h2so4, &
                                               qtend_cond_aging_soag,  &
                                               qtend_rename_h2so4,     &
                                               qtend_rename_soag,      &
                                               qtend_newnuc_h2so4,     &
                                               qtend_newnuc_soag,      &
                                               qtend_coag_h2so4,       &
                                               qtend_coag_soag

#if(defined USE_NC4)
!
! output comparison results
!
      ! Create the netCDF file. The nf90_clobber parameter tells
      ! netCDF to overwrite this file, if it already exists.
      call check( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) )

      ! Define the dimensions. NetCDF will hand back an ID for each.
      call check( nf90_def_dim(ncid, "nsteps", nstop, nstep_dimid) )
      call check( nf90_def_dim(ncid, "mode", ntot_amode, mode_dimid) )

      ! The dimids array is used to pass the IDs of the dimensions of
      ! the variables.
      dimids  = (/nstep_dimid, mode_dimid/)

! Original Variables Definition

call check(nf90_def_var(ncid, "num_aer", &
          NF90_DOUBLE, dimids, varid(1)) )
call check(nf90_def_var(ncid, "so4_aer", &
          NF90_DOUBLE, dimids, varid(2)) )
call check(nf90_def_var(ncid, "nh4_aer", &
          NF90_DOUBLE, dimids, varid(3)) )
call check(nf90_def_var(ncid, "no3_aer", &
          NF90_DOUBLE, dimids, varid(4)) )
call check(nf90_def_var(ncid, "soa_aer", &
          NF90_DOUBLE, dimids, varid(5)) )
call check(nf90_def_var(ncid, "pom_aer", &
          NF90_DOUBLE, dimids, varid(6)) )
call check(nf90_def_var(ncid, "bc_aer", &
          NF90_DOUBLE, dimids, varid(7)) )
call check(nf90_def_var(ncid, "dst_aer", &
          NF90_DOUBLE, dimids, varid(8)) )
call check(nf90_def_var(ncid, "co3_aer", &
          NF90_DOUBLE, dimids, varid(9)) )
call check(nf90_def_var(ncid, "ncl_aer", &
          NF90_DOUBLE, dimids, varid(10)) )
call check(nf90_def_var(ncid, "ca_aer", &
          NF90_DOUBLE, dimids, varid(11)) )
call check(nf90_def_var(ncid, "cl_aer", &
          NF90_DOUBLE, dimids, varid(12)) )
call check(nf90_def_var(ncid, "wat_aer", &
          NF90_DOUBLE, dimids, varid(13)) )



! Gas phase variables - added hno3, nh3, soag, hcl, so2
call check(nf90_def_var(ncid, "h2so4_gas", &
          NF90_DOUBLE, dimids(1), varid(14)) )
call check(nf90_def_var(ncid, "hno3_gas", &
          NF90_DOUBLE, dimids(1), varid(15)) )
call check(nf90_def_var(ncid, "nh3_gas", &
          NF90_DOUBLE, dimids(1), varid(16)) )
call check(nf90_def_var(ncid, "soag_gas", &
          NF90_DOUBLE, dimids(1), varid(17)) )
call check(nf90_def_var(ncid, "hcl_gas", &
          NF90_DOUBLE, dimids(1), varid(18)) )
call check(nf90_def_var(ncid, "so2_gas", &
          NF90_DOUBLE, dimids(1), varid(19)) )

call check(nf90_def_var(ncid, "Dgn_mode", &
          NF90_DOUBLE, dimids, varid(20)) )

call check(nf90_def_var(ncid, "hygro_aer", &
          NF90_DOUBLE, dimids, varid(21)) )

call check(nf90_def_var(ncid, "relhum", &
          NF90_DOUBLE, dimids(1), varid(22)) )
call check(nf90_def_var(ncid, "temp", &
          NF90_DOUBLE, dimids(1), varid(23)) )

! Optical depth (AOD) for each species per mode
call check(nf90_def_var(ncid, "aod_sulfate", &
          NF90_DOUBLE, dimids, varid(24)) )
call check(nf90_def_var(ncid, "aod_bc", &
          NF90_DOUBLE, dimids, varid(25)) )
call check(nf90_def_var(ncid, "aod_pom", &
          NF90_DOUBLE, dimids, varid(26)) )
call check(nf90_def_var(ncid, "aod_soa", &
          NF90_DOUBLE, dimids, varid(27)) )
call check(nf90_def_var(ncid, "aod_dust", &
          NF90_DOUBLE, dimids, varid(28)) )
call check(nf90_def_var(ncid, "aod_seasalt", &
          NF90_DOUBLE, dimids, varid(29)) )
call check(nf90_def_var(ncid, "aod_mode", &
          NF90_DOUBLE, dimids, varid(30)) )

! Single scattering albedo for each species per mode
call check(nf90_def_var(ncid, "ssa_sulfate", &
          NF90_DOUBLE, dimids, varid(31)) )
call check(nf90_def_var(ncid, "ssa_bc", &
          NF90_DOUBLE, dimids, varid(32)) )
call check(nf90_def_var(ncid, "ssa_pom", &
          NF90_DOUBLE, dimids, varid(33)) )
call check(nf90_def_var(ncid, "ssa_soa", &
          NF90_DOUBLE, dimids, varid(34)) )
call check(nf90_def_var(ncid, "ssa_dust", &
          NF90_DOUBLE, dimids, varid(35)) )
call check(nf90_def_var(ncid, "ssa_seasalt", &
          NF90_DOUBLE, dimids, varid(36)) )
call check(nf90_def_var(ncid, "ssa_mode", &
          NF90_DOUBLE, dimids, varid(37)) )

! Asymmetry parameter for each species per mode
call check(nf90_def_var(ncid, "asm_sulfate", &
          NF90_DOUBLE, dimids, varid(38)) )
call check(nf90_def_var(ncid, "asm_bc", &
          NF90_DOUBLE, dimids, varid(39)) )
call check(nf90_def_var(ncid, "asm_pom", &
          NF90_DOUBLE, dimids, varid(40)) )
call check(nf90_def_var(ncid, "asm_soa", &
          NF90_DOUBLE, dimids, varid(41)) )
call check(nf90_def_var(ncid, "asm_dust", &
          NF90_DOUBLE, dimids, varid(42)) )
call check(nf90_def_var(ncid, "asm_seasalt", &
          NF90_DOUBLE, dimids, varid(43)) )
call check(nf90_def_var(ncid, "asm_mode", &
          NF90_DOUBLE, dimids, varid(44)) )

#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
! Marine organic matter (MOM)
call check(nf90_def_var(ncid, "aod_mom", &
          NF90_DOUBLE, dimids, varid(45)) )
call check(nf90_def_var(ncid, "ssa_mom", &
          NF90_DOUBLE, dimids, varid(46)) )
call check(nf90_def_var(ncid, "asm_mom", &
          NF90_DOUBLE, dimids, varid(47)) )

! Nitrate (NO3)
call check(nf90_def_var(ncid, "aod_no3", &
          NF90_DOUBLE, dimids, varid(48)) )
call check(nf90_def_var(ncid, "ssa_no3", &
          NF90_DOUBLE, dimids, varid(49)) )
call check(nf90_def_var(ncid, "asm_no3", &
          NF90_DOUBLE, dimids, varid(50)) )

! Ammonium (NH4)
call check(nf90_def_var(ncid, "aod_nh4", &
          NF90_DOUBLE, dimids, varid(51)) )
call check(nf90_def_var(ncid, "ssa_nh4", &
          NF90_DOUBLE, dimids, varid(52)) )
call check(nf90_def_var(ncid, "asm_nh4", &
          NF90_DOUBLE, dimids, varid(53)) )
#endif
  
call check(nf90_put_att(ncid, varid(1), "units", "#/kg-air") )
call check(nf90_put_att(ncid, varid(2), "units", "kg-aer/kg-air") )
call check(nf90_put_att(ncid, varid(3), "units", "kg-aer/kg-air") )
call check(nf90_put_att(ncid, varid(4), "units", "kg-gas/kg-air") )
call check(nf90_put_att(ncid, varid(5), "units", "kg-gas/kg-air") )
call check(nf90_put_att(ncid, varid(6), "units", "kg-aer/kg-air") )
call check(nf90_put_att(ncid, varid(7), "units", "kg-aer/kg-air") )
call check(nf90_put_att(ncid, varid(8), "units", "kg-aer/kg-air") )
call check(nf90_put_att(ncid, varid(9), "units", "kg-aer/kg-air") )
call check(nf90_put_att(ncid, varid(10), "units", "kg-aer/kg-air") )
call check(nf90_put_att(ncid, varid(11), "units", "kg-aer/kg-air") )
call check(nf90_put_att(ncid, varid(12), "units", "kg-aer/kg-air") )
call check(nf90_put_att(ncid, varid(13), "units", "kg-aer/kg-air") )
!gas phase
call check(nf90_put_att(ncid, varid(14), "units", "kg-gas/kg-air") )
call check(nf90_put_att(ncid, varid(15), "units", "kg-gas/kg-air") )
call check(nf90_put_att(ncid, varid(16), "units", "kg-gas/kg-air") )
call check(nf90_put_att(ncid, varid(17), "units", "kg-gas/kg-air") )
call check(nf90_put_att(ncid, varid(18), "units", "kg-gas/kg-air") )
call check(nf90_put_att(ncid, varid(19), "units", "kg-gas/kg-air") )
!properties
call check(nf90_put_att(ncid, varid(20), "units", "micro-m") )
call check(nf90_put_att(ncid, varid(21), "units", "none") )

call check(nf90_put_att(ncid, varid(22), "units", "none") )
call check(nf90_put_att(ncid, varid(23), "units", "K") )

! Add units attributes for optical properties
call check(nf90_put_att(ncid, varid(24), "units", "none") )
call check(nf90_put_att(ncid, varid(24), "long_name", "Sulfate AOD at 550nm") )
call check(nf90_put_att(ncid, varid(25), "units", "none") )
call check(nf90_put_att(ncid, varid(25), "long_name", "BC AOD at 550nm") )
call check(nf90_put_att(ncid, varid(26), "units", "none") )
call check(nf90_put_att(ncid, varid(26), "long_name", "POM AOD at 550nm") )
call check(nf90_put_att(ncid, varid(27), "units", "none") )
call check(nf90_put_att(ncid, varid(27), "long_name", "SOA AOD at 550nm") )
call check(nf90_put_att(ncid, varid(28), "units", "none") )
call check(nf90_put_att(ncid, varid(28), "long_name", "Dust AOD at 550nm") )
call check(nf90_put_att(ncid, varid(29), "units", "none") )
call check(nf90_put_att(ncid, varid(29), "long_name", "Sea salt AOD at 550nm") )
call check(nf90_put_att(ncid, varid(30), "units", "none") )
call check(nf90_put_att(ncid, varid(30), "long_name", "Total mode AOD at 550nm") )

call check(nf90_put_att(ncid, varid(31), "units", "none") )
call check(nf90_put_att(ncid, varid(31), "long_name", "Sulfate SSA at 550nm") )
call check(nf90_put_att(ncid, varid(32), "units", "none") )
call check(nf90_put_att(ncid, varid(32), "long_name", "BC SSA at 550nm") )
call check(nf90_put_att(ncid, varid(33), "units", "none") )
call check(nf90_put_att(ncid, varid(33), "long_name", "POM SSA at 550nm") )
call check(nf90_put_att(ncid, varid(34), "units", "none") )
call check(nf90_put_att(ncid, varid(34), "long_name", "SOA SSA at 550nm") )
call check(nf90_put_att(ncid, varid(35), "units", "none") )
call check(nf90_put_att(ncid, varid(35), "long_name", "Dust SSA at 550nm") )
call check(nf90_put_att(ncid, varid(36), "units", "none") )
call check(nf90_put_att(ncid, varid(36), "long_name", "Sea salt SSA at 550nm") )
call check(nf90_put_att(ncid, varid(37), "units", "none") )
call check(nf90_put_att(ncid, varid(37), "long_name", "Total mode SSA at 550nm") )

call check(nf90_put_att(ncid, varid(38), "units", "none") )
call check(nf90_put_att(ncid, varid(38), "long_name", "Sulfate asymmetry parameter at 550nm") )
call check(nf90_put_att(ncid, varid(39), "units", "none") )
call check(nf90_put_att(ncid, varid(39), "long_name", "BC asymmetry parameter at 550nm") )
call check(nf90_put_att(ncid, varid(40), "units", "none") )
call check(nf90_put_att(ncid, varid(40), "long_name", "POM asymmetry parameter at 550nm") )
call check(nf90_put_att(ncid, varid(41), "units", "none") )
call check(nf90_put_att(ncid, varid(41), "long_name", "SOA asymmetry parameter at 550nm") )
call check(nf90_put_att(ncid, varid(42), "units", "none") )
call check(nf90_put_att(ncid, varid(42), "long_name", "Dust asymmetry parameter at 550nm") )
call check(nf90_put_att(ncid, varid(43), "units", "none") )
call check(nf90_put_att(ncid, varid(43), "long_name", "Sea salt asymmetry parameter at 550nm") )
call check(nf90_put_att(ncid, varid(44), "units", "none") )
call check(nf90_put_att(ncid, varid(44), "long_name", "Total mode asymmetry parameter at 550nm") )

! -----------------------------------------------------------------------
! LW absorption optical depth variables (reference band: ilw_diag=9, ~8-12 um)
! varid(54-60): mode total + 6 species
! -----------------------------------------------------------------------
call check(nf90_def_var(ncid, "aod_lw_mode",    NF90_DOUBLE, dimids, varid(54)) )
call check(nf90_def_var(ncid, "aod_lw_sulfate", NF90_DOUBLE, dimids, varid(55)) )
call check(nf90_def_var(ncid, "aod_lw_bc",      NF90_DOUBLE, dimids, varid(56)) )
call check(nf90_def_var(ncid, "aod_lw_pom",     NF90_DOUBLE, dimids, varid(57)) )
call check(nf90_def_var(ncid, "aod_lw_soa",     NF90_DOUBLE, dimids, varid(58)) )
call check(nf90_def_var(ncid, "aod_lw_dust",    NF90_DOUBLE, dimids, varid(59)) )
call check(nf90_def_var(ncid, "aod_lw_seasalt", NF90_DOUBLE, dimids, varid(60)) )

#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
call check(nf90_put_att(ncid, varid(45), "units", "none") )
call check(nf90_put_att(ncid, varid(45), "long_name", "Marine organic matter AOD at 550nm") )
call check(nf90_put_att(ncid, varid(46), "units", "none") )
call check(nf90_put_att(ncid, varid(46), "long_name", "Marine organic matter SSA at 550nm") )
call check(nf90_put_att(ncid, varid(47), "units", "none") )
call check(nf90_put_att(ncid, varid(47), "long_name", "Marine organic matter asymmetry parameter at 550nm") )

call check(nf90_put_att(ncid, varid(48), "units", "none") )
call check(nf90_put_att(ncid, varid(48), "long_name", "Nitrate AOD at 550nm") )
call check(nf90_put_att(ncid, varid(49), "units", "none") )
call check(nf90_put_att(ncid, varid(49), "long_name", "Nitrate SSA at 550nm") )
call check(nf90_put_att(ncid, varid(50), "units", "none") )
call check(nf90_put_att(ncid, varid(50), "long_name", "Nitrate asymmetry parameter at 550nm") )

call check(nf90_put_att(ncid, varid(51), "units", "none") )
call check(nf90_put_att(ncid, varid(51), "long_name", "Ammonium AOD at 550nm") )
call check(nf90_put_att(ncid, varid(52), "units", "none") )
call check(nf90_put_att(ncid, varid(52), "long_name", "Ammonium SSA at 550nm") )
call check(nf90_put_att(ncid, varid(53), "units", "none") )
call check(nf90_put_att(ncid, varid(53), "long_name", "Ammonium asymmetry parameter at 550nm") )
#endif

! LW optical depth attributes (reference band ilw_diag=9, atmospheric window ~8-12 um)
call check(nf90_put_att(ncid, varid(54), "units", "none") )
call check(nf90_put_att(ncid, varid(54), "long_name", "Total mode LW absorption optical depth (band 9, ~8-12 um)") )
call check(nf90_put_att(ncid, varid(55), "units", "none") )
call check(nf90_put_att(ncid, varid(55), "long_name", "Sulfate LW absorption optical depth (band 9, ~8-12 um)") )
call check(nf90_put_att(ncid, varid(56), "units", "none") )
call check(nf90_put_att(ncid, varid(56), "long_name", "BC LW absorption optical depth (band 9, ~8-12 um)") )
call check(nf90_put_att(ncid, varid(57), "units", "none") )
call check(nf90_put_att(ncid, varid(57), "long_name", "POM LW absorption optical depth (band 9, ~8-12 um)") )
call check(nf90_put_att(ncid, varid(58), "units", "none") )
call check(nf90_put_att(ncid, varid(58), "long_name", "SOA LW absorption optical depth (band 9, ~8-12 um)") )
call check(nf90_put_att(ncid, varid(59), "units", "none") )
call check(nf90_put_att(ncid, varid(59), "long_name", "Dust LW absorption optical depth (band 9, ~8-12 um)") )
call check(nf90_put_att(ncid, varid(60), "units", "none") )
call check(nf90_put_att(ncid, varid(60), "long_name", "Sea salt LW absorption optical depth (band 9, ~8-12 um)") )

! Add global attribute
      call check( nf90_put_att(ncid, NF90_GLOBAL, &
                               "Created_by", "ETHz") )
      call date_and_time(date)
      call check( nf90_put_att(ncid, NF90_GLOBAL, &
                               "Created_date", date) )

      ! End define mode. This tells netCDF we are done defining
      ! metadata.
      call check( nf90_enddef(ncid) )

#endif

      tmp_dgn_a              = 0._r8 ; tmp_dgn_awet          = 0._r8
      tmp_num_aer            = 0._r8 ; tmp_so4_aer           = 0._r8
      tmp_soa_aer            = 0._r8 ; tmp_h2so4             = 0._r8
      tmp_soag               = 0._r8
      tmp_relh =0._r8
      tmp_temp=0._r8
      qtend_cond_aging_so4   = 0._r8 ; qtend_cond_aging_soa  = 0._r8
      qtend_rename_so4       = 0._r8 ; qtend_rename_soa      = 0._r8
      qtend_newnuc_so4       = 0._r8 ; qtend_newnuc_soa      = 0._r8
      qtend_coag_so4         = 0._r8 ; qtend_coag_soa        = 0._r8
      qtend_cond_aging_h2so4 = 0._r8 ; qtend_cond_aging_soag = 0._r8
      qtend_rename_h2so4     = 0._r8 ; qtend_rename_soag     = 0._r8
      qtend_newnuc_h2so4     = 0._r8 ; qtend_newnuc_soag     = 0._r8
      qtend_coag_h2so4       = 0._r8 ; qtend_coag_soag       = 0._r8

      tmp_aod_sulfate  = 0._r8 ; tmp_aod_bc       = 0._r8
      tmp_aod_pom      = 0._r8 ; tmp_aod_soa      = 0._r8
      tmp_aod_dust     = 0._r8 ; tmp_aod_seasalt  = 0._r8
      tmp_aod_mode     = 0._r8

      tmp_ssa_sulfate  = 0._r8 ; tmp_ssa_bc       = 0._r8
      tmp_ssa_pom      = 0._r8 ; tmp_ssa_soa      = 0._r8
      tmp_ssa_dust     = 0._r8 ; tmp_ssa_seasalt  = 0._r8
      tmp_ssa_mode     = 0._r8

      tmp_asm_sulfate  = 0._r8 ; tmp_asm_bc       = 0._r8
      tmp_asm_pom      = 0._r8 ; tmp_asm_soa      = 0._r8
      tmp_asm_dust     = 0._r8 ; tmp_asm_seasalt  = 0._r8
      tmp_asm_mode     = 0._r8

      tmp_aod_lw_mode     = 0._r8
      tmp_aod_lw_sulfate  = 0._r8 ; tmp_aod_lw_bc       = 0._r8
      tmp_aod_lw_pom      = 0._r8 ; tmp_aod_lw_soa      = 0._r8
      tmp_aod_lw_dust     = 0._r8 ; tmp_aod_lw_seasalt  = 0._r8

#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
tmp_aod_mom      = 0._r8 ; tmp_aod_no3      = 0._r8 ; tmp_aod_nh4      = 0._r8
tmp_ssa_mom      = 0._r8 ; tmp_ssa_no3      = 0._r8 ; tmp_ssa_nh4      = 0._r8
tmp_asm_mom      = 0._r8 ; tmp_asm_no3      = 0._r8 ; tmp_asm_nh4      = 0._r8
#endif




      lchnk = begchunk


      latndx = -1
      lonndx = -1


      nacc = modeptr_accum
      l_num_a1 = numptr_amode(nacc)
      l_so4_a1 = lptr_so4_a_amode(nacc)
      l_nh4_a1 = lptr_nh4_a_amode(nacc)

      nait = modeptr_aitken
      l_num_a2 = numptr_amode(nait)
      l_so4_a2 = lptr_so4_a_amode(nait)
      l_nh4_a2 = lptr_nh4_a_amode(nait)

      lmz_h2so4g = l_h2so4g - (imozart-1)
      lmz_so2g   = l_so2g   - (imozart-1)
      lmz_nh3g   = l_nh3g   - (imozart-1)
      lmz_hno3g  = l_hno3g  - (imozart-1)
      lmz_hclg   = l_hclg   - (imozart-1)
      lmz_soag   = l_soag   - (imozart-1)

      lmz_num_a1 = l_num_a1 - (imozart-1)
      lmz_so4_a1 = l_so4_a1 - (imozart-1)
      lmz_nh4_a1 = l_nh4_a1 - (imozart-1)

      lmz_num_a2 = l_num_a2 - (imozart-1)
      lmz_so4_a2 = l_so4_a2 - (imozart-1)
      lmz_nh4_a2 = l_nh4_a2 - (imozart-1)

print*, 'starting time integration' 
main_time_loop:&
do nstep = 1, nstop
      istep = nstep
      if (nstep == 1) tnew = 0.0_r8
      told = tnew
      tnew = told + deltat

      if (nstep == 1) then
        temp = tmin
        relh = rhmin
      end if
!      if (nstep > 1 .and. tmax > tmin ) then
      if (nstep > 1 ) then
        temp = temp + (tmax -tmin)/(nstop -1)
        physta%t(1,1) = temp
      end if
!      if (nstep > 1 .and. rhmax > rhmin ) then
      if (nstep > 1) then
        relh = relh + (rhmax -rhmin)/(nstop -1)
       end if

       physta%t(:,:) = temp
       physta%relhum(:,:) = relh
       call  qsat( physta%t(1:pcols,1:pver), physta%pmid(1:pcols,1:pver), &
                  ev_sat(1:pcols,1:pver), qv_sat(1:pcols,1:pver) )

       physta%qv(:,:) = physta%relhum(:,:)*qv_sat(:,:)
       physta%q(:,:,1) = physta%qv(:,:)

       print*, 'FAB QVSAT DE QSAT', qv_sat(:,:)
!
! calcsize
!
      write(*,'(/a,i8)') 'cambox_do_run doing calcsize, istep=', istep
      print*, lchnk,pcols

! *** new calcsize interface ***
! load state
      physta%lchnk = lchnk
      physta%ncol  = pcols

! load pbuf
      call load_pbuf( pbuf, lchnk, pcols, &
         physta%cld, physta%qqcw, physta%dgncur_a, physta%dgncur_awet,  physta%qaerwat, physta%wetdens, physta%hygro)
      ! load ptend
      ptend%lq    =.false.
      ptend%q     = 0.
      call modal_aero_calcsize_sub( physta, ptend, deltat, pbuf, &
         do_adjust_in=.true., do_aitacc_transfer_in=.true. )

! unload pbuf
      call unload_pbuf( pbuf, lchnk, pcols, &
         physta%cld, physta%qqcw, physta%dgncur_a, physta%dgncur_awet,  physta%qaerwat, physta%wetdens, physta%hygro )

!apply tendencies
      do l = 1, pcnst
         if ( .not. ptend%lq(l) ) cycle
         do k = 1, pver
         do i = 1, pcols
            physta%q(i,k,l) = physta%q(i,k,l) + ptend%q(i,k,l)*deltat
            physta%q(i,k,l) = max( physta%q(i,k,l), 0.0_r8 )
         end do
         end do
      end do

!
! *****wateruptake******
!
     write(*,'(/a,i8)') 'cambox_do_run doing wateruptake, istep=', istep

     call load_pbuf( pbuf, lchnk, pcols, &
        physta%cld, physta%qqcw, physta%dgncur_a, physta%dgncur_awet,  physta%qaerwat, physta%wetdens, physta%hygro )

!
     call modal_aero_wateruptake_dr( physta, pbuf,deltat,nstep )

     call unload_pbuf( pbuf, lchnk, pcols, &
         physta%cld, physta%qqcw, physta%dgncur_a, physta%dgncur_awet,  physta%qaerwat, physta%wetdens, physta%hygro)


!
! switch from q & qqcw to vmr and vmrcw
!
      loffset = imozart - 1
      vmr = 0.0_r8
      vmrcw = 0.0_r8
      do l = imozart, pcnst
         l2 = l - loffset
         vmr(  :,:,l2) = physta%q(  :,:,l)*mwdry/adv_mass(l2)
         vmrcw(:,:,l2) = physta%qqcw(:,:,l)*mwdry/adv_mass(l2)
      end do

IF (nstep > 1) then
!
! gaschem_simple
!
      lun = 6
      write(*,'(/a,i8)') 'cambox_do_run doing gaschem simple, istep=', istep
      vmr_svaa   = vmr
      vmrcw_svaa = vmrcw
      h2so4_pre_gaschem(:,:) = vmr(:,:,lmz_h2so4g)

! global avg ~= 13 d = 1.12e6 s, daytime avg ~= 5.6e5, noontime peak ~= 3.7e5
      tau_gaschem_simple = 3.e4  ! 3.0e5  ! so2 gas-rxn timescale (s)

      if (mdo_gaschem > 0) then
         call gaschem_simple_sub(                       &
            lchnk,    pcols,     nstep,               &
            loffset,  deltat,                        &
            vmr,                tau_gaschem_simple      )
!      else
         ! assumed constant gas chemistry production rate (mol/mol)
!         vmr(:,:,lmz_h2so4g) = vmr(:,:,lmz_h2so4g) + 1.e-16_r8*deltat
      end if

      h2so4_aft_gaschem(:,:) = vmr(:,:,lmz_h2so4g)

!
! cloudchem_simple
!
      write(*,'(/a,i8)') &
         'cambox_do_run doing cloudchem simple, istep=', istep
      vmr_svbb = vmr
      vmrcw_svbb = vmrcw

      if (mdo_cloudchem > 0 ) then

!      call cloudchem_simple_sub(                  &
!         lchnk,    pcols,     nstep,               &
!         loffset,  deltat,                        &
!         vmr,      vmrcw,    cld_ncol             )

      end if ! (mdo_cloudchem > 0 .and. maxval( cld_ncol(:,:) ) > 1.0e-6_r8) then


!
! gasaerexch
!
      write(*,'(/a,i8)') 'cambox_do_run doing gasaerexch, istep=', istep
      vmr_svcc = vmr
      vmrcw_svcc = vmrcw

      dvmrdt_bb = 0.0_r8 ; dvmrcwdt_bb = 0.0_r8

!      if (nstep > 1) then ! want to write  the initial state in output
      call modal_aero_amicphys_intr(              &
         mdo_gasaerexch,     mdo_rename,          &
         mdo_newnuc,         mdo_coag,            &
         lchnk,    pcols,     nstep,               &
         loffset,  deltat,                        &
         latndx,   lonndx,                        &
         physta%t,   physta%pmid, physta%pdel,    &
         physta%zm,  physta%pblh,                 &
         physta%qv,  physta%cld ,                 &
         vmr,                vmrcw,               &   ! after  cloud chem
         vmr_svaa,                                &   ! before gas chem
         vmr_svbb,           vmrcw_svbb,          &   ! before cloud chem
#if ( defined( CAMBOX_ACTIVATE_THIS ) )
         nqtendbb,           nqqcwtendbb,         &
         dvmrdt_bb,          dvmrcwdt_bb,         &
#endif
         physta%dgncur_a,     physta%dgncur_awet, &
         physta%wetdens,      physta%qaerwat      )


      dvmrdt_cond(  :,:,:) = dvmrdt_bb(  :,:,:,iqtend_cond)
      dvmrdt_rnam(  :,:,:) = dvmrdt_bb(  :,:,:,iqtend_rnam)
      dvmrdt_nnuc(  :,:,:) = dvmrdt_bb(  :,:,:,iqtend_nnuc)
      dvmrdt_coag(  :,:,:) = dvmrdt_bb(  :,:,:,iqtend_coag)
      dvmrcwdt_cond(:,:,:) = 0.0_r8
      dvmrcwdt_rnam(:,:,:) = dvmrcwdt_bb(:,:,:,iqqcwtend_rnam)
      dvmrcwdt_nnuc(:,:,:) = 0.0_r8
      dvmrcwdt_coag(:,:,:) = 0.0_r8

END IF
!
! done
!
      write(*,'(/a,i8)') 'cambox_do_run step done, istep=', istep

!
! switch from vmr & vmrcw to q & qqcw
!
      loffset = imozart - 1
      do l = imozart, pcnst
         l2 = l - loffset
         physta%q(    :,:,l)  = vmr(  :,:,l2) * adv_mass(l2)/mwdry
         physta%qqcw( :,:,l)  = vmrcw(:,:,l2) * adv_mass(l2)/mwdry
      end do
!
! calculate optical properties 
!

       call  mam_aero_sw(physta, tauxar, wa, ga, fa, mamoptdiag)

       write(*,'(/a,i8)') 'modal optical properties done, istep=', istep 
       print*,'tauxar', tauxar 
       print*, 'wa/tauxar', wa/tauxar  
       do l=1,4
       
       print*,physta%pdeldry(:,:)*rga, mamoptdiag(l)%vext_sulfate(:,:,10),  mamoptdiag(l)%vext_mode(:,:,10)*physta%pdeldry(:,:)*rga  
       end do

       call  mam_aero_lw(physta, taux_lw, mamoptdiag)

       write(*,'(/a,i8)') 'modal LW optical properties done, istep=', istep
       print*,'taux_lw (band',ilw_diag,')', taux_lw(:,:,ilw_diag)
       ! store the data of each time step for netcdf output
!

      if (l_h2so4g > 0) tmp_h2so4 (nstep) = physta%q(1,1,l_h2so4g)* adv_mass(l_h2so4g-loffset)/mwdry *1E9
      if (l_hno3g > 0) tmp_hno3(nstep) = physta%q(1,1,l_hno3g)* adv_mass(l_hno3g-loffset)/mwdry *1E9
      if (l_nh3g > 0)  tmp_nh3(nstep)  = physta%q(1,1,l_nh3g)* adv_mass(l_nh3g-loffset)/mwdry *1E9
      if (l_soag > 0)  tmp_soag(nstep) = physta%q(1,1,l_soag)* adv_mass(l_soag-loffset)/mwdry *1E9
      if (l_hclg > 0)  tmp_hcl(nstep)  = physta%q(1,1,l_hclg)* adv_mass(l_hclg-loffset)/mwdry *1E9
      if (l_so2g > 0)  tmp_so2(nstep)  = physta%q(1,1,l_so2g) * adv_mass(l_so2g-loffset)/mwdry *1E9

      tmp_dgn_a(nstep,1:ntot_amode)        = physta%dgncur_a(1,1,1:ntot_amode)
      tmp_dgn_awet(nstep,1:ntot_amode)     = physta%dgncur_awet(1,1,1:ntot_amode)
      tmp_dgn_a(nstep,1:ntot_amode)  = tmp_dgn_awet(nstep,1:ntot_amode)  ! use wet in th output

      do i = 1, ntot_amode
       tmp_num_aer(nstep,i) = physta%q(1,1,numptr_amode(i))
       if(lptr_so4_a_amode(i)>0)  tmp_so4_aer(nstep,i) = physta%q(1,1,lptr_so4_a_amode(i))
       if(lptr_nh4_a_amode(i)>0)  tmp_nh4_aer(nstep,i) = physta%q(1,1,lptr_nh4_a_amode(i))
       if(lptr_no3_a_amode(i)>0)  tmp_no3_aer(nstep,i) = physta%q(1,1,lptr_no3_a_amode(i))
       if(lptr_soa_a_amode(i)>0)  tmp_soa_aer(nstep,i) = physta%q(1,1,lptr_soa_a_amode(i))
       if(lptr_pom_a_amode(i)>0)  tmp_pom_aer(nstep,i) = physta%q(1,1,lptr_pom_a_amode(i))
       if(lptr_bc_a_amode(i)>0)   tmp_bc_aer(nstep,i)  = physta%q(1,1,lptr_bc_a_amode(i))
       if(lptr_dust_a_amode(i)>0) tmp_dust_aer(nstep,i) = physta%q(1,1,lptr_dust_a_amode(i))
       if(lptr_co3_a_amode(i)>0)  tmp_co3_aer(nstep,i) = physta%q(1,1,lptr_co3_a_amode(i))
       if(lptr_ca_a_amode(i)>0)   tmp_ca_aer(nstep,i) = physta%q(1,1,lptr_ca_a_amode(i))
       if(lptr_nacl_a_amode(i)>0) tmp_nacl_aer(nstep,i) = physta%q(1,1,lptr_nacl_a_amode(i))
       if(lptr_cl_a_amode(i)>0) tmp_cl_aer(nstep,i) = physta%q(1,1,lptr_cl_a_amode(i))
       tmp_wat_aer(nstep,i) = physta%qaerwat(1,1,i)
       tmp_hygro_aer(nstep,i) = physta%hygro(1,1,i)
       end do

       tmp_relh(nstep) = physta%relhum(1,1)
       tmp_temp(nstep) = physta%t(1,1)


       ! Store optical properties for each mode
! Convert extinction (m2/kg-air) to AOD (dimensionless) using pdeldry*rga
do i = 1, ntot_amode
   ! AOD = extinction * (pdeldry * rga)
   ! pdeldry is in Pa, rga = 1/g = 0.102 kg/m2/Pa
   tmp_aod_sulfate(nstep,i)  = mamoptdiag(i)%vext_sulfate(1,1,10) * physta%pdeldry(1,1) * rga
   tmp_aod_bc(nstep,i)       = mamoptdiag(i)%vext_bc(1,1,10) * physta%pdeldry(1,1) * rga
   tmp_aod_pom(nstep,i)      = mamoptdiag(i)%vext_pom(1,1,10) * physta%pdeldry(1,1) * rga
   tmp_aod_soa(nstep,i)      = mamoptdiag(i)%vext_soa(1,1,10) * physta%pdeldry(1,1) * rga
   tmp_aod_dust(nstep,i)     = mamoptdiag(i)%vext_dust(1,1,10) * physta%pdeldry(1,1) * rga
   tmp_aod_seasalt(nstep,i)  = mamoptdiag(i)%vext_seasalt(1,1,10) * physta%pdeldry(1,1) * rga
!   tmp_aod_mode(nstep,i)     = mamoptdiag(i)%vext_mode(1,1) * physta%pdeldry(1,1) * rga
   tmp_aod_mode(nstep,i)     = mamoptdiag(i)%tauxar(1,1,10)
   
   ! Single scattering albedo (dimensionless)
   tmp_ssa_sulfate(nstep,i)  = mamoptdiag(i)%vssa_sulfate(1,1,10)
   tmp_ssa_bc(nstep,i)       = mamoptdiag(i)%vssa_bc(1,1,10)
   tmp_ssa_pom(nstep,i)      = mamoptdiag(i)%vssa_pom(1,1,10)
   tmp_ssa_soa(nstep,i)      = mamoptdiag(i)%vssa_soa(1,1,10)
   tmp_ssa_dust(nstep,i)     = mamoptdiag(i)%vssa_dust(1,1,10)
   tmp_ssa_seasalt(nstep,i)  = mamoptdiag(i)%vssa_seasalt(1,1,10)
   tmp_ssa_mode(nstep,i)     = mamoptdiag(i)%vssa_mode(1,1,10)
!   tmp_ssa_mode(nstep,i)     = mamoptdiag(i)%ssav(1,1,10)   
   ! Asymmetry parameter (dimensionless)
   tmp_asm_sulfate(nstep,i)  = mamoptdiag(i)%vasm_sulfate(1,1,10)
   tmp_asm_bc(nstep,i)       = mamoptdiag(i)%vasm_bc(1,1,10)
   tmp_asm_pom(nstep,i)      = mamoptdiag(i)%vasm_pom(1,1,10)
   tmp_asm_soa(nstep,i)      = mamoptdiag(i)%vasm_soa(1,1,10)
   tmp_asm_dust(nstep,i)     = mamoptdiag(i)%vasm_dust(1,1,10)
   tmp_asm_seasalt(nstep,i)  = mamoptdiag(i)%vasm_seasalt(1,1,10)
   tmp_asm_mode(nstep,i)     = mamoptdiag(i)%vasm_mode(1,1,10)

   ! LW absorption optical depth (reference band ilw_diag, ~8-12 um)
   ! tauxar_lw already has units of optical depth (dimensionless)
   ! vext_lw_* is m2/kg-air -> multiply by pdeldry*rga to get optical depth
   tmp_aod_lw_mode    (nstep,i) = mamoptdiag(i)%tauxar_lw(1,1,ilw_diag)
   tmp_aod_lw_sulfate (nstep,i) = mamoptdiag(i)%vext_lw_sulfate (1,1,ilw_diag) * physta%pdeldry(1,1) * rga
   tmp_aod_lw_bc      (nstep,i) = mamoptdiag(i)%vext_lw_bc      (1,1,ilw_diag) * physta%pdeldry(1,1) * rga
   tmp_aod_lw_pom     (nstep,i) = mamoptdiag(i)%vext_lw_pom     (1,1,ilw_diag) * physta%pdeldry(1,1) * rga
   tmp_aod_lw_soa     (nstep,i) = mamoptdiag(i)%vext_lw_soa     (1,1,ilw_diag) * physta%pdeldry(1,1) * rga
   tmp_aod_lw_dust    (nstep,i) = mamoptdiag(i)%vext_lw_dust    (1,1,ilw_diag) * physta%pdeldry(1,1) * rga
   tmp_aod_lw_seasalt (nstep,i) = mamoptdiag(i)%vext_lw_seasalt (1,1,ilw_diag) * physta%pdeldry(1,1) * rga
end do
#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
do i = 1, ntot_amode
   ! Marine organic matter (MOM)
   tmp_aod_mom(nstep,i)  = mamoptdiag(i)%vext_mom(1,1,10) * physta%pdeldry(1,1) * rga
   tmp_ssa_mom(nstep,i)  = mamoptdiag(i)%vssa_mom(1,1,10)
   tmp_asm_mom(nstep,i)  = mamoptdiag(i)%vasm_mom(1,1,10)
   
   ! Nitrate (NO3)
   tmp_aod_no3(nstep,i)  = mamoptdiag(i)%vext_no3(1,1,10) * physta%pdeldry(1,1) * rga
   tmp_ssa_no3(nstep,i)  = mamoptdiag(i)%vssa_no3(1,1,10)
   tmp_asm_no3(nstep,i)  = mamoptdiag(i)%vasm_no3(1,1,10)
   
   ! Ammonium (NH4)
   tmp_aod_nh4(nstep,i)  = mamoptdiag(i)%vext_nh4(1,1,10) * physta%pdeldry(1,1) * rga
   tmp_ssa_nh4(nstep,i)  = mamoptdiag(i)%vssa_nh4(1,1,10)
   tmp_asm_nh4(nstep,i)  = mamoptdiag(i)%vasm_nh4(1,1,10)
end do
#endif

end do main_time_loop

#if (defined USE_NC4)
!
! Write the data to the file.
!

call check( nf90_put_var(ncid, varid(1), &
            tmp_num_aer(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(2), &
            tmp_so4_aer(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(3), &
            tmp_nh4_aer(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(4), &
            tmp_no3_aer(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(5), &
            tmp_soa_aer(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(6), &
            tmp_pom_aer(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(7), &
            tmp_bc_aer(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(8), &
            tmp_dust_aer(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(9), &
            tmp_co3_aer(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(10), &
            tmp_nacl_aer(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(11), &
            tmp_ca_aer(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(12), &
            tmp_cl_aer(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(13), &
            tmp_wat_aer(1:nstop,1:ntot_amode)) )

!  Write gas phase data - added hno3, nh3, soag, hcl, so2
call check( nf90_put_var(ncid, varid(14), &
            tmp_h2so4(1:nstop)) )
call check( nf90_put_var(ncid, varid(15), &
            tmp_hno3(1:nstop)) )
call check( nf90_put_var(ncid, varid(16), &
            tmp_nh3(1:nstop)) )
call check( nf90_put_var(ncid, varid(17), &
            tmp_soag(1:nstop)) )
call check( nf90_put_var(ncid, varid(18), &
            tmp_hcl(1:nstop)) )
call check( nf90_put_var(ncid, varid(19), &
            tmp_so2(1:nstop)) )

call check( nf90_put_var(ncid, varid(20), &
            tmp_dgn_a(1:nstop,1:ntot_amode)) )

call check( nf90_put_var(ncid, varid(21), &
            tmp_hygro_aer(1:nstop,1:ntot_amode)) )

call check( nf90_put_var(ncid, varid(22), &
            tmp_relh(1:nstop)) )

call check( nf90_put_var(ncid, varid(23), &
            tmp_temp(1:nstop)) )


! Write optical depth (AOD)
call check( nf90_put_var(ncid, varid(24), &
            tmp_aod_sulfate(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(25), &
            tmp_aod_bc(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(26), &
            tmp_aod_pom(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(27), &
            tmp_aod_soa(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(28), &
            tmp_aod_dust(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(29), &
            tmp_aod_seasalt(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(30), &
            tmp_aod_mode(1:nstop,1:ntot_amode)) )

! Write single scattering albedo (SSA)
call check( nf90_put_var(ncid, varid(31), &
            tmp_ssa_sulfate(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(32), &
            tmp_ssa_bc(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(33), &
            tmp_ssa_pom(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(34), &
            tmp_ssa_soa(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(35), &
            tmp_ssa_dust(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(36), &
            tmp_ssa_seasalt(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(37), &
            tmp_ssa_mode(1:nstop,1:ntot_amode)) )

! Write asymmetry parameter
call check( nf90_put_var(ncid, varid(38), &
            tmp_asm_sulfate(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(39), &
            tmp_asm_bc(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(40), &
            tmp_asm_pom(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(41), &
            tmp_asm_soa(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(42), &
            tmp_asm_dust(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(43), &
            tmp_asm_seasalt(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(44), &
            tmp_asm_mode(1:nstop,1:ntot_amode)) )

! Write LW absorption optical depth (reference band ilw_diag=9, ~8-12 um)
call check( nf90_put_var(ncid, varid(54), &
            tmp_aod_lw_mode(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(55), &
            tmp_aod_lw_sulfate(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(56), &
            tmp_aod_lw_bc(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(57), &
            tmp_aod_lw_pom(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(58), &
            tmp_aod_lw_soa(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(59), &
            tmp_aod_lw_dust(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(60), &
            tmp_aod_lw_seasalt(1:nstop,1:ntot_amode)) )

#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
! Marine organic matter (MOM)
call check( nf90_put_var(ncid, varid(45), &
            tmp_aod_mom(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(46), &
            tmp_ssa_mom(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(47), &
            tmp_asm_mom(1:nstop,1:ntot_amode)) )

! Nitrate (NO3)
call check( nf90_put_var(ncid, varid(48), &
            tmp_aod_no3(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(49), &
            tmp_ssa_no3(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(50), &
            tmp_asm_no3(1:nstop,1:ntot_amode)) )

! Ammonium (NH4)
call check( nf90_put_var(ncid, varid(51), &
            tmp_aod_nh4(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(52), &
            tmp_ssa_nh4(1:nstop,1:ntot_amode)) )
call check( nf90_put_var(ncid, varid(53), &
            tmp_asm_nh4(1:nstop,1:ntot_amode)) )

#endif


    ! Close the file. This frees up any internal netCDF resources
      ! associated with the file, and flushes any buffers.
      call check( nf90_close(ncid) )

#endif

      return
      end subroutine gcmambox_do_run




      subroutine check(status)
      integer, intent(in) :: status
#if(defined USE_NC4)
      if(status /= nf90_noerr) then
         print *, trim(nf90_strerror(status))
         stop "Stopped"
      end if
#endif
      end subroutine check


!-------------------------------------------------------------------------------

      end module driver
