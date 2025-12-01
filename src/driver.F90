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
                           mdo_gaschem, mdo_cloudchem,&
                           mdo_gasaerexch, mdo_rename, mdo_newnuc, mdo_coag

      use constituents, only: pcnst, cnst_name, cnst_get_ind
      use modal_aero_data, only: ntot_amode
      use physics_buffer, only: physics_buffer_desc
      use physics_types, only : physics_state, physics_ptend

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
      
      integer :: nstop
      real*8  :: deltat 
      
      pcols = 1
      pver = 1
 
      plev = pver

      mdo_gaschem=1
      mdo_cloudchem =1

      mdo_gasaerexch=1
      mdo_rename=1
      mdo_newnuc=1
      mdo_coag=1

      masterproc = .true. 


      write(*,'(/a)') '*** Hello from MAIN ***'

      write(*,'(/a)') '*** main calling MAM_init_basics'
      call MAM_init_basics(pbuf)
      write(*,'(/a)') '*** main call MAM_allocate'
      call MAM_ALLOCATE (physta,ptend )
     
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

! iniialise ( q(:,:,1 ) call gestbl to build saturation vapor pressure table.
      tmn   = 173.16_r8
      tmx   = 375.16_r8
      trice =  20.00_r8
      ip    = .true.
      call gestbl(tmn     ,tmx     ,trice   ,ip      ,epsilo  , &
                  latvap  ,latice  ,rh2o    ,cpair   ,tmelt )

      call  qsat( physta%t(1:pcols,1:pver), physta%pmid(1:pcols,1:pver), &
                  ev_sat(1:pcols,1:pver), qv_sat(1:pcols,1:pver) )
      physta%qv(:,:) = physta%relhum(:,:)*qv_sat(:,:)
      physta%q(:,:,1) = physta%qv(:,:)


      return
      end subroutine gcmambox_init_run


!-------------------------------------------------------------------------------
      subroutine gcmambox_do_run( nstop, deltat) 

      use chem_mods, only: adv_mass, gas_pcnst, imozart
      use physconst, only: mwdry
      use physics_types, only: physics_state, physics_ptend
      use physics_buffer, only: physics_buffer_desc, pbuf_get_chunk
      use mam_utils, only: l_h2so4g, l_nh3g, l_so2g, l_hno3g, l_hclg, l_soag
      use modal_aero_data
      use modal_aero_calcsize, only: modal_aero_calcsize_sub
      use modal_aero_amicphys, only: modal_aero_amicphys_intr, &
          gaexch_h2so4_uptake_optaa, newnuc_h2so4_conc_optaa, mosaic
      use modal_aero_wateruptake, only: modal_aero_wateruptake_dr
      use gaschem_simple, only: gaschem_simple_sub
      use cloudchem_simple, only: cloudchem_simple_sub
      use wv_saturation, only : qsat
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

!
! for netcdf file
!
      character (len = *), parameter :: FILE_NAME = "mam_output.nc"
      integer :: ncid, nstep_dimid, mode_dimid
      integer :: error
      integer :: dimids(2), varid(30)
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
      if (nstep > 1 .and. tmax > tmin ) then 
        temp = temp + (tmax -tmin)/(nstop -1)
        physta%t(1,1) = temp
      end if 
      if (nstep > 1 .and. rhmax > rhmin ) then
        relh = relh + (rhmax -rhmin)/(nstop -1)
       end if  

       physta%t(:,:) = temp
       physta%relhum(:,:) = relh
       call  qsat( physta%t(1:pcols,1:pver), physta%pmid(1:pcols,1:pver), &
                  ev_sat(1:pcols,1:pver), qv_sat(1:pcols,1:pver) )
       
       physta%qv(:,:) = physta%relhum(:,:)*qv_sat(:,:)
       physta%q(:,:,1) = physta%qv(:,:)
        

      print*, 'T , H  !! ', physta%t, physta%q(:,:,1)
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

print*, 'FAB qaerwat',physta%qaerwat(1,1,:), physta%hygro

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
      write(lun,'(/a,i8)') 'cambox_do_run doing gaschem simple, istep=', istep
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
      else
         ! assumed constant gas chemistry production rate (mol/mol)
         vmr(:,:,lmz_h2so4g) = vmr(:,:,lmz_h2so4g) + 1.e-16_r8*deltat
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
    
      ! Close the file. This frees up any internal netCDF resources
      ! associated with the file, and flushes any buffers.
      call check( nf90_close(ncid) )

#endif

      return
      end subroutine gcmambox_do_run


!-------------------------------------------------------------------------------
      subroutine load_pbuf( pbuf, lchnk, ncol, &
         cld, qqcw, dgncur_a, dgncur_awet, qaerwat, wetdens, hygro )

      use mam_utils, only: pcols,pver
      use constituents, only : pcnst
      use chem_mods, only: adv_mass, gas_pcnst, imozart
      use physconst, only: mwdry

      use modal_aero_data, only:  &
         lmassptrcw_amode, nspec_amode, numptrcw_amode, &
         qqcw_get_field

      use physics_buffer, only: physics_buffer_desc, &
         pbuf_get_index, pbuf_get_field

      type(physics_buffer_desc), pointer :: pbuf(:)  ! physics buffer for a chunk

      integer,  intent(in   ) :: lchnk, ncol

      real(r8), intent(in   ) :: cld(pcols,pver)    ! stratiform cloud fraction
      real(r8), intent(in   ) :: qqcw(pcols,pver,pcnst)  ! Cloudborne aerosol MR array
      real(r8), intent(in   ) :: dgncur_a(pcols,pver,ntot_amode)
      real(r8), intent(in   ) :: dgncur_awet(pcols,pver,ntot_amode)
      real(r8), intent(in   ) :: qaerwat(pcols,pver,ntot_amode)
      real(r8), intent(in   ) :: wetdens(pcols,pver,ntot_amode)
      real(r8), intent(in   ) :: hygro(pcols,pver,ntot_amode)

      integer :: idx, l, ll, n

      real(r8), pointer :: fldcw(:,:)
      real(r8), pointer :: ycld(:,:)
      real(r8), pointer :: ydgnum(:,:,:)
      real(r8), pointer :: ydgnumwet(:,:,:)
      real(r8), pointer :: yqaerwat(:,:,:)
      real(r8), pointer :: ywetdens(:,:,:)
      real(r8), pointer :: yhygro(:,:,:)

      idx = pbuf_get_index( 'CLD' )
      call pbuf_get_field( pbuf, idx, ycld )
      ycld(:,:) = 0.0_r8
      ycld(1:ncol,:) = cld(1:ncol,:)

      idx = pbuf_get_index( 'DGNUM' )
      call pbuf_get_field( pbuf, idx, ydgnum )
      ydgnum(:,:,:) = 0.0_r8
      ydgnum(1:ncol,:,:) = dgncur_a(1:ncol,:,:)

      idx = pbuf_get_index( 'DGNUMWET' )
      call pbuf_get_field( pbuf, idx, ydgnumwet )
      ydgnumwet(:,:,:) = 0.0_r8
      ydgnumwet(1:ncol,:,:) = dgncur_awet(1:ncol,:,:)

      idx = pbuf_get_index( 'QAERWAT' )
      call pbuf_get_field( pbuf, idx, yqaerwat )
      yqaerwat(:,:,:) = 0.0_r8
      yqaerwat(1:ncol,:,:) = qaerwat(1:ncol,:,:)

      idx = pbuf_get_index( 'WETDENS_AP' )
      call pbuf_get_field( pbuf, idx, ywetdens )
      ywetdens(:,:,:) = 0.0_r8
      ywetdens(1:ncol,:,:) = wetdens(1:ncol,:,:)

      idx = pbuf_get_index( 'HYGROM' )
      call pbuf_get_field( pbuf, idx, yhygro )
      yhygro(:,:,:) = 0.0_r8
      yhygro(1:ncol,:,:) = hygro(1:ncol,:,:)

      do n = 1, ntot_amode
      do ll = 0, nspec_amode(n)
         l = numptrcw_amode(n)
         if (ll > 0) l = lmassptrcw_amode(ll,n)
         fldcw => qqcw_get_field( pbuf, l, lchnk )
         fldcw(:,:) = 0.0_r8
         fldcw(1:ncol,:) = qqcw(1:ncol,:,l)
      end do
      end do
   
      print*,'end pbuf'
      return
      end subroutine load_pbuf


!-------------------------------------------------------------------------------
      subroutine unload_pbuf( pbuf, lchnk, ncol, &
         cld, qqcw, dgncur_a, dgncur_awet, qaerwat, wetdens, hygro )

      use mam_utils, only: pcols,pver
      use constituents, only : pcnst
      use chem_mods, only: adv_mass, gas_pcnst, imozart
      use physconst, only: mwdry

      use modal_aero_data, only:  &
         lmassptrcw_amode, nspec_amode, numptrcw_amode, &
         qqcw_get_field

      use physics_buffer, only: physics_buffer_desc, &
         pbuf_get_index, pbuf_get_field

      type(physics_buffer_desc), pointer :: pbuf(:)  ! physics buffer for a chunk

      integer,  intent(in   ) :: lchnk, ncol

      real(r8), intent(in   ) :: cld(pcols,pver)    ! stratiform cloud fraction

      real(r8), intent(inout) :: qqcw(pcols,pver,pcnst)  ! Cloudborne aerosol MR array
      real(r8), intent(inout) :: dgncur_a(pcols,pver,ntot_amode)
      real(r8), intent(inout) :: dgncur_awet(pcols,pver,ntot_amode)
      real(r8), intent(inout) :: qaerwat(pcols,pver,ntot_amode)
      real(r8), intent(inout) :: wetdens(pcols,pver,ntot_amode)
      real(r8), intent(inout) :: hygro(pcols,pver,ntot_amode)

      integer :: i, idx, k, l, ll, n
      real(r8) :: tmpa

      real(r8), pointer :: fldcw(:,:)
      real(r8), pointer :: ycld(:,:)
      real(r8), pointer :: ydgnum(:,:,:)
      real(r8), pointer :: ydgnumwet(:,:,:)
      real(r8), pointer :: yqaerwat(:,:,:)
      real(r8), pointer :: ywetdens(:,:,:)
      real(r8), pointer :: yhygro(:,:,:)

      idx = pbuf_get_index( 'CLD' )
      call pbuf_get_field( pbuf, idx, ycld )
! cld should not have changed, so check for changes rather than unloading it
!     cld(1:ncol,:) = ycld(1:ncol,:)
      tmpa = maxval( abs( cld(1:ncol,:) - ycld(1:ncol,:) ) )
      if (tmpa /= 0.0_r8) then
         write(*,*) '*** unload_pbuf cld change error - ', tmpa
         stop
      end if

      idx = pbuf_get_index( 'DGNUM' )
      call pbuf_get_field( pbuf, idx, ydgnum )
      dgncur_a(1:ncol,:,:) = ydgnum(1:ncol,:,:)

      idx = pbuf_get_index( 'DGNUMWET' )
      call pbuf_get_field( pbuf, idx, ydgnumwet )
      dgncur_awet(1:ncol,:,:) = ydgnumwet(1:ncol,:,:)

      idx = pbuf_get_index( 'QAERWAT' )
      call pbuf_get_field( pbuf, idx, yqaerwat )
      qaerwat(1:ncol,:,:) = yqaerwat(1:ncol,:,:)

      idx = pbuf_get_index( 'WETDENS_AP' )
      call pbuf_get_field( pbuf, idx, ywetdens )
      wetdens(1:ncol,:,:) = ywetdens(1:ncol,:,:)

      idx = pbuf_get_index( 'HYGROM' )
      call pbuf_get_field( pbuf, idx, yhygro )
      hygro(1:ncol,:,:) = yhygro(1:ncol,:,:)


      do n = 1, ntot_amode
      do ll = 0, nspec_amode(n)
         l = numptrcw_amode(n)
         if (ll > 0) l = lmassptrcw_amode(ll,n)
         fldcw => qqcw_get_field( pbuf, l, lchnk )
         qqcw(1:ncol,:,l) = fldcw(1:ncol,:)
      end do
      end do


      return
      end subroutine unload_pbuf


!-------------------------------------------------------------------------------


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
