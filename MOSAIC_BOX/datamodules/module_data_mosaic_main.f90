      module module_data_mosaic_main

      use module_data_mosaic_kind, only:  r8

      implicit none

      integer, parameter ::   &
                ngas_com = 44,   &
                ngas_urb = 19,   &
                ngas_bio =  7,   &
                ngas_mar = 11

      integer, parameter ::   &
                naer_tot = 28,		   &  ! total num of 3-D variables/bin
                naerbin  = 44 ! 110 ! 41760  	      ! ( 48 size)*(29 wbc)*(30 kappa)
!               naerbin  = 3240  	      ! ( 24 size)*(15 wbc)*( 9 kappa)
!               naerbin  = 90000 	      ! (100 size)*(30 wbc)*(30 kappa)

      integer, parameter ::   &
                ncld_tot = 13,		   &  ! + 8 = total num of 3-D variables/bin
                ncldbin  =  4,		   &  ! num of cloud bins
                ncld     = 22		! num of dynamic cloud species/bin

      integer, parameter :: ngas_max = ngas_com + ngas_urb + ngas_bio + ngas_mar
      integer, parameter :: naer_max = naer_tot*naerbin
      integer, parameter :: ncld_max = ncld_tot*ncldbin

      integer, parameter :: ntot_max = ngas_max + naer_max + ncld_max

      integer, parameter :: kmaxd = 100		! num of vertical levels in a column mode

      integer, save ::   &
      		naerbin_used=0,   &   ! num of aerosol bins being used
      		ncldbin_used=0,   &   ! num of  cloud  bins being used
      		ntot_used=ngas_max    ! portion of cnn array being used

      integer, save ::   &
      		iwrite_gas,   &
      		iwrite_aer_bin,   &
      		iwrite_aer_dist,   &
      		iwrite_aer_species

      integer, save ::   &
      		lun_inp,   &
      		lun_gas,   &
      		lun_aer(naerbin),   &
      		lun_aer_status(naerbin),   &
      		lun_aeroptic,   &
      		lun_drydist,   &
      		lun_wetdist,   &
      		lun_species, &
			lun_totaer

      integer, save ::   &
                lun_sect_170,   & 
                lun_sect_171,   & 
                lun_sect_172,   & 
                lun_sect_180,   & 
                lun_sect_183,   & 
                lun_sect_184,   & 
                lun_sect_185,   & 
                lun_sect_186,   & 
                lun_sect_188,   & 
                lun_sect_190

      integer, save ::   &
      		ipmcmos = 0,   &   ! if > 0, do emissions, dilution, air density,
      		                   ! and relative humidity as in partmc_mosaic 
      		istate_pblh = 0    ! used for dilution from pblh changes

      real(r8), parameter :: press0_pa = 1.01325d5  ! pressure of 1 atm [Pa]
      real(r8), parameter :: mw_air = 28.966d0      ! dry-air mean molecular weight [g/mol]

      character(len=64), save ::   &
      		inputfile,   &
      		gas_output,   &
      		aer_output(naerbin),   &
      		aeroptic_output,   &
      		drydist_output,   &
      		wetdist_output,   &
      		species_output

!------------------------------------------------------------------------
! Global Species Indices
!
      integer, save ::   &
       kh2so4,      khno3,       khcl,        knh3,        kno,       &
       kno2,        kno3,        kn2o5,       khono,       khno4,     &
       ko3,         ko1d,        ko3p,        koh,         kho2,      &
       kh2o2,       kco,         kso2,        kch4,        kc2h6,     &
       kch3o2,      kethp,       khcho,       kch3oh,      kanol,     &
       kch3ooh,     kethooh,     kald2,       khcooh,      krcooh,    &
       kc2o3,       kpan,        karo1,       karo2,       karo3,     &
       karo4,       kalk1,       kole1,       kapi1,       kapi2,     &
       kapi3,       kapi4,       klim1,       klim2,                  &
       kpar,        kaone,       kmgly,       keth,        kolet,     &
       kolei,       ktol,        kxyl,        kcres,       kto2,      &
       kcro,        kopen,       konit,       krooh,       kro2,      &
       kano2,       knap,        kxo2,        kxpar,                  &
       kisop,       kisoprd,     kisopp,      kisopn,      kisopo2,   &
       kapi,        klim,                                             &
       kdms,        kmsa,        kdmso,       kdmso2,      kch3so2h,  &
       kch3sch2oo,  kch3so2,     kch3so3,     kch3so2ch2oo,kch3so2oo, &
       ksulfhox

      integer, save ::   &
       knum_a,      kdpdry_a,    ksigmag_a,  kjhyst_a,    kwater_a,   &
       kso4_a,      kno3_a,      kcl_a,      knh4_a,      koc_a,      &
       kmsa_a,      kco3_a,      kna_a,      kca_a,       kbc_a,      &
       koin_a,      karo1_a,     karo2_a,    karo3_a,     karo4_a,    &
       kalk1_a,     kole1_a,     kapi1_a,    kapi2_a,     kapi3_a,    &
       kapi4_a,     klim1_a,     klim2_a

      integer, save ::   &
       knum_c,      kwater_c,    &
       kso4_c,      kno3_c,      kcl_c,      knh4_c,      koc_c,      &
       kmsa_c,      kco3_c,      kna_c,      kca_c,       kbc_c,      &
       koin_c,      karo1_c,     karo2_c,    karo3_c,     karo4_c,    &
       kalk1_c,     kole1_c,     kapi1_c,    kapi2_c,     kapi3_c,    &
       kapi4_c,     klim1_c,     klim2_c



!-------------------------------------------------------------------------

      character(len=40), save :: species(ntot_max)

      integer, save ::	mmode, mgas, maer, mcld

      integer, save ::	maeroptic, mshellcore

      real(r8), save :: cnn(ntot_max)

      real(r8), save :: tNi, tSi, tCli, tNH4i, DN, DS, DCl, DNH4

      real(r8), save :: emission(ntot_max), emit(ntot_max)

      real(r8), save :: avogad, deg2rad, pi, piover4, piover6, third

      integer, save ::   &
       tbeg_dd,   tbeg_mo,   &
       tbeg_hh,   tbeg_mm,   tbeg_ss,   tmar21_sec,   &
       trun_dd,   trun_hh,   trun_mm,   trun_ss,   &
       tbeg_sec,  trun_sec,   &
       nstep,     it,        iprint

      real(r8), save ::   &
       tcur_sec,  tcur_min,  tcur_hrs,		   &  ! time since beginning of year (UTC)
       time_sec,  time_min,  time_hrs,		   &  ! time since beginning of simulation
       time_UTC,  time_UTC_beg,			   &  ! time of day in hrs (UTC)
       tmid_sec,  tsav_sec,  told_sec, time_sec_old,   &
	 t_since_start,                              &
       dt_sec,    dt_min,    dt_aeroptic_min,      &
       rlon,      rlat,                            &
       zalt_m,    cos_sza

      real(r8), save ::   &
       cair_mlc,  cair_molm3,     h2o,       o2,   &
       h2,        ppb,            speed_molec,     &
       te,        pr_atm,         RH,        pblh, &
       cair_mlc_old, cair_molm3_old,               &
       te_old,    pr_atm_old,     RH_old,    pblh_old

      integer, save ::   &
       idaytime,  msolar, mphoto

	integer, save :: icase_soa, iurbanleg, imodel, iraoultslaw, ireactiveuptake

!      real, save ::   &
!       gas(ngas_max), flux(ngas_max,naerbin)


!------------------------------------------------------------------------

      end module module_data_mosaic_main