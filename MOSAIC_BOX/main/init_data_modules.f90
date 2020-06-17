      subroutine init_data_modules
!
!   place various constant or initial values into common
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero
      use module_data_mosaic_cloud

      implicit none


!--------------------------------------------------
!     define fundamental constants...
       avogad	= 6.02217e+23

!      pi	= 3.141592654
!      deg2rad	= 0.017453293
       pi		= 4.0_r8 * atan( 1.0_r8 )
       piover4	= pi/4.0_r8
       piover6	= pi/6.0_r8
       deg2rad	= pi/180.0_r8
       third	= 1.0_r8/3.0_r8

!---------------------------------
! define species indices
!
! species in inorganic chemistry
       kh2so4		=  1
       khno3		=  2
       khcl			=  3
       knh3			=  4
       kno			=  5
       kno2			=  6
       kno3			=  7
       kn2o5		=  8
       khono		=  9
       khno4		= 10
       ko3			= 11
       ko1d			= 12
       ko3p			= 13
       koh			= 14
       kho2			= 15
       kh2o2		= 16
       kco			= 17
       kso2			= 18
!
! species in methane, ethane, formaldehyde chemistry
       kch4			= 19
       kc2h6		= 20
       kch3o2		= 21
       kethp		= 22
       khcho		= 23
       kch3oh		= 24
       kanol		= 25
       kch3ooh		= 26
       kethooh		= 27
       kald2		= 28
       khcooh		= 29
       krcooh		= 30
       kc2o3		= 31
       kpan			= 32
       karo1		= 33	! soa prec 1
       karo2		= 34	! soa prec 2
       karo3		= 35	! soa prec 3
       karo4		= 36	! soa prec 4
       kalk1		= 37	! soa prec 5
       kole1		= 38	! soa prec 6
       kapi1		= 39	! soa prec 7
       kapi2		= 40	! soa prec 8
       kapi3		= 41	! soa prec 9
       kapi4		= 42	! soa prec 10
       klim1		= 43	! soa prec 11
       klim2		= 44	! soa prec 12
!
! species in hc1 mechanism. initialize indices to zero
       kpar			= 45
       kaone		= 46
       kmgly		= 47
       keth			= 48
       kolet		= 49
       kolei		= 50
       ktol			= 51
       kxyl			= 52
       kcres		= 53
       kto2			= 54
       kcro			= 55
       kopen	 	= 56
       konit		= 57
       krooh		= 58
       kro2			= 59
       kano2		= 60
       knap			= 61
       kxo2			= 62
       kxpar		= 63
!
! species in hc2 mechanism. initialize indices to zero
       kisop		= 64
       kisoprd		= 65
       kisopp		= 66
       kisopn		= 67
       kisopo2		= 68
       kapi			= 69
       klim			= 70

! species in dms mechanism. initialize indices to zero
       kdms			= 71
       kmsa			= 72
       kdmso		= 73
       kdmso2		= 74
       kch3so2h		= 75
       kch3sch2oo		= 76
       kch3so2		= 77
       kch3so3		= 78
       kch3so2oo		= 79
       kch3so2ch2oo	= 80
       ksulfhox		= 81

! aerosol bin: (24 species=bin)
       knum_a		=  1
       kdpdry_a		=  2
       ksigmag_a		=  3
       kjhyst_a		=  4
       kwater_a		=  5
       kso4_a		=  6
       kno3_a		=  7
       kcl_a		=  8
       knh4_a		=  9
       kmsa_a		= 10
       karo1_a		= 11
       karo2_a		= 12
       karo3_a		= 13
       karo4_a		= 14
       kalk1_a		= 15
       kole1_a		= 16
       kapi1_a		= 17
       kapi2_a		= 18
       kapi3_a		= 19
       kapi4_a		= 20
       klim1_a		= 21
       klim2_a		= 22
       kco3_a		= 23
       kna_a		= 24
       kca_a		= 25
       koin_a		= 26
       koc_a		= 27
       kbc_a		= 28



! cloud mode: 1-ncldmode (13=mode)
       knum_c		=  1
       kwater_c		=  2
       kso4_c		=  3
       kno3_c		=  4
       kcl_c		=  5
       knh4_c		=  6
       koc_c		=  7
       kmsa_c		=  8
       kco3_c		=  9
       kna_c		= 10
       kca_c		= 11
       kbc_c		= 12
       koin_c		= 13
!
!
! regime-dependent chemistry definitions
!
       iregime		=  1
!
!     GAS
!
       ih2so4		=  1
       ihno3		=  2
       ihcl			=  3
       inh3			=  4
       ino			=  5
       ino2			=  6
       ino3			=  7
       in2o5		=  8
       ihono		=  9
       ihno4		= 10
       io3			= 11
       io1d			= 12
       io3p			= 13
       ioh			= 14
       iho2			= 15
       ih2o2		= 16
       ico			= 17
       iso2			= 18
!
! species in methane, ethane, formaldehyde chemistry
       ich4			= 19
       ic2h6		= 20
       ich3o2		= 21
       iethp		= 22
       ihcho		= 23
       ich3oh		= 24
       ianol		= 25
       ich3ooh		= 26
       iethooh		= 27
       iald2		= 28
       ihcooh		= 29
       ircooh		= 30
       ic2o3		= 31
       ipan			= 32
       iaro1		= 33	! soa prec 1
       iaro2		= 34	! soa prec 2
       iaro3		= 35	! soa prec 1
       iaro4		= 36	! soa prec 2
       ialk1		= 37	! soa prec 3
       iole1		= 38	! soa prec 4
       iapi1		= 39	! soa prec 5
       iapi2		= 40	! soa prec 6
       iapi3		= 41	! soa prec 5
       iapi4		= 42	! soa prec 6
       ilim1		= 43	! soa prec 7
       ilim2		= 44	! soa prec 8
!
! species in hc1 mechanism. initialize indices to zero
       ipar			= 45
       iaone		= 46
       imgly		= 47
       ieth			= 48
       iolet		= 49
       iolei		= 50
       itol			= 51
       ixyl			= 52
       icres		= 53
       ito2			= 54
       icro			= 55
       iopen		= 56
       ionit		= 57
       irooh		= 58
       iro2			= 59
       iano2		= 60
       inap			= 61
       ixo2			= 62
       ixpar		= 63
!
! species in hc2 mechanism. initialize indices to zero
       iisop		= 64
       iisoprd		= 65
       iisopp		= 66
       iisopn		= 67
       iisopo2		= 68
       iapi			= 69
       ilim			= 70

! species in dms mechanism. initialize indices to zero
       idms			= 71
       imsa			= 72
       idmso		= 73
       idmso2		= 74
       ich3so2h		= 75
       ich3sch2oo		= 76
       ich3so2		= 77
       ich3so3		= 78
       ich3so2oo		= 79
       ich3so2ch2oo	= 80
       isulfhox		= 81

! alkylperoxy radical indices for parameterized permutation reactions
       jch3o2		=  1
       jethp		=  2
       jro2			=  3
       jc2o3		=  4
       jano2		=  5
       jnap			=  6
       jisopp		=  7
       jisopn		=  8
       jisopo2		=  9
       jxo2			= 10

! photolyzing species indices
       jphoto_no2		=  1
       jphoto_no3		=  2
       jphoto_hono	=  3
       jphoto_hno3	=  4
       jphoto_hno4	=  5
       jphoto_n2o5	=  6
       jphoto_o3a		=  7
       jphoto_o3b		=  8
       jphoto_h2o2	=  9
       jphoto_hchoa	= 10
       jphoto_hchob	= 11
       jphoto_ch3ooh	= 12
       jphoto_ethooh	= 13
       jphoto_ald2	= 14
       jphoto_aone	= 15
       jphoto_mgly	= 16
       jphoto_open	= 17
       jphoto_rooh	= 18
       jphoto_onit	= 19
       jphoto_isoprd	= 20


!
!     CLOUD
!
! cloud (local): used for total and undissociated species
       iso4_c		=  1
       ino3_c		=  2
       icl_c		=  3
       inh4_c		=  4
       ioc_c		=  5
       imsa_c		=  6
       ico2_c		=  7
       ina_c		=  8
       ica_c		=  9
       ibc_c		= 10
       ioin_c		= 11
       iso2_c		= 12
       ihono_c		= 13
       ih2o2_c		= 14
       ich3ooh_c		= 15
       ihcooh_c		= 16
       ircooh_c		= 17
       ihcho_c		= 18
       io3_c		= 19
       iho2_c		= 20	! don't actually exist in cloud (surf rxn only)
       ino2_c		= 21	! don't actually exist in cloud (surf rxn only)
       ino3r_c		= 22	! don't actually exist in cloud (surf rxn only)
       in2o5_c		= 23	! don't actually exist in cloud (surf rxn only)

! cloud ionic species
       jh_c			=  1
       jnh4_c		=  2
       jna_c		=  3
       jhso4_c		=  4
       jso4_c		=  5
       jno3_c		=  6
       jcl_c		=  7
       jno2_c		=  8
       jhso3_c		=  9
       jso3_c		= 10
       jhco3_c		= 11
       jco3_c		= 12
       jho2_c		= 13
       jhcoo_c		= 14
       jrcoo_c		= 15
       jmsa_c		= 16
       joh_c		= 17
       jch2oh2_c		= 18


      return
      end
