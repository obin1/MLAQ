*Reactions* 16
{1} ISOP + OH       = ISOPP + .08 XO2       : ARR(2.6e-11, 409.)              ;
{2} ISOP + O3       = .6 HCHO + .65 ISOPRD +
                      .27 OH + .07 HO2 + 
                      .07 CO + .39 RCOOH + 
                      .15 ALD2 + .2 XO2 +
                      .2 C2O3               : ARR(1.2e-14, -2013.)            ;
{3} ISOP + NO3      = ISOPN                 : ARR(3.0e-12, -446.)             ;
{4} ISOPRD + hv     = .97 C2O3 + .33 HO2 + 
                      .33 CO + .7 CH3O2 +
                      .2 HCHO + .07 ALD2 +
                      .03 AONE              : photo_ISOPRD                    ;
{5} ISOPRD + OH     = .5 C2O3 + .5 ISOPO2 +
                      .2 XO2                : 3.3e-11                         ;
{6} ISOPRD + O3     = .27 OH + .1 HO2 + 
                      .11 C2O3 + .07 XO2 +
                      .05 CH3O2 + .16 CO + 
                      .15 HCHO + .02 ALD2 +
                      .09 AONE + .85 MGLY +
                      .46 RCOOH             : 7.0e-18                         ;
{7} ISOPRD + NO3    = .07 C2O3 + .07 HNO3 +
                      .64 CO + .28 HCHO + 
                      .93 ONIT + .28 ALD2 +
                      .93 HO2 + .93 XO2 +
                      1.86 PAR              : 1.0e-15                         ;
{8} ISOPP + NO      = .09 ONIT + .91 NO2 + 
                      .91 HO2 + .63 HCHO +
                      .91 ISOPRD + 0.18 PAR : 4.0e-12                         ;
{9} ISOPN + NO      = 1.2 NO2 + .8 ONIT + 
                      .8 ALD2 + .8 HO2 + 
                      .2 ISOPRD + 1.6 PAR   : 4.0e-12                         ;
{10} ISOPO2 + NO    = NO2 + HO2 + .59 CO +
                      .55 ALD2 + .25 HCHO +
                      .34 MGLY + .63 AONE   : 4.0e-12                         ;
{11} ISOPP + HO2    = ROOH                  : ARR(1.7e-13, 1300.)             ;
{12} ISOPN + HO2    = ONIT + 2 PAR          : ARR(1.7e-13, 1300.)             ;
{13} ISOPO2 + HO2   = ROOH                  : ARR(1.7e-13, 1300.)             ;
{14} ISOPP          = ISOPRD                : rk_param(jisopp)                ;
{15} ISOPN          = ALD2 + ONIT + 2 PAR   : rk_param(jisopn)                ;
{16} ISOPO2         = .5 ALD2 + .5 AONE     : rk_param(jisopo2)               ;
*end*