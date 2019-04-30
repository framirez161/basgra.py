import math

def Harvest(year,doy,DAYS_HARVEST,TILV,TILG2,CLAIV,LAI,DELT,CLV,PHEN,HAGERE,CRES,CST):
    HARV   = 0
    NOHARV = 1  
    for i in range (len(DAYS_HARVEST)):    
        if year==DAYS_HARVEST[i][0] and doy==DAYS_HARVEST[i][1]:
            HARV   = 1
            NOHARV = 0	
  
    FRACTV = TILV/(TILG2 + TILV)
    CLAI   = FRACTV * CLAIV
    if (LAI <= CLAI) :
        HARVFR = 0.0
    else:
        HARVFR = 1.0 - CLAI/LAI

    HARVLA    = (HARV   * LAI * HARVFR) / DELT
    HARVLV    = (HARV   * CLV * HARVFR) / DELT
    HARVPH    = (HARV   * PHEN        ) / DELT
    TV1       = (HARVFR * FRACTV) + (1-FRACTV)*HAGERE
    HARVRE    = (HARV   * TV1 * CRES  ) / DELT
    HARVST    = (HARV   * CST         ) / DELT
    GSTUB     =  HARVST * (1-HAGERE)
    HARVTILG2 = (HARV   * TILG2       ) / DELT    
    return(NOHARV,FRACTV,HARVLA,HARVLV,HARVPH,TV1,HARVRE,HARVST,GSTUB,HARVTILG2)

def Biomass(COCRESMX,CLV,CRES,CST):
    CRESMX = COCRESMX*(CLV + CRES + CST)
    RESNOR = max(0,min(1, CRES/CRESMX ))    
    return(CRESMX, RESNOR)

def Phenology(DAVTMP,DAYLP,DAYL,DAYLB,PHEN,DELT,PHENCR,DLMXGE):
    GPHEN = max(0, (DAVTMP-0.01)*0.000144*24 * (min(DAYLP,DAYL)-0.24) )
    DPHEN = 0
    if DAYL < DAYLB:
      DPHEN = PHEN / DELT
    PHENRF = (1 - PHEN)/(1 - PHENCR)
    if PHENRF > 1.0: 
      PHENRF = 1.0
    if PHENRF < 0.0: 
      PHENRF = 0.0
    DAYLGE = max(0,min(1, (DAYL - DAYLB)/(DLMXGE - DAYLB) ))
    return(GPHEN,DPHEN,PHENRF,DAYLGE)

def Foliage1(TBASE,DAVTMP,DAYLGE,SLAMAX,FSLAMIN,RESNOR):
    EFFTMP = max(TBASE, DAVTMP)
    LERV   =          max(0, (-0.76 + 0.52*EFFTMP)/1000)
    LERG   = DAYLGE * max(0, (-5.46 + 2.80*EFFTMP)/1000)
    SLAMIN = SLAMAX * FSLAMIN
    SLANEW = SLAMAX - RESNOR*(SLAMAX-SLAMIN)
    return(EFFTMP,LERV,LERG,SLAMIN,SLANEW)

def LUECO2TM(DAVTMP,RUBISC,CO2A,PARAV,KLUETILG,FRACTV,K):
#=============================================================================
# Calculate LUEMXQ (mol CO2 mol-1 PAR quanta)
# Inputs : PARAV (micromol PAR quanta m-2 s-)
#=============================================================================
    T      = DAVTMP                                            #(degC)
    RUBISCN = RUBISC * (1.E6/550000.)
    EAVCMX =  68000                                            #(J mol-1)
    EAKMC  =  65800                                            #(J mol-1)
    EAKMO  =   1400                                            #(J mol-1)
    KC25   =     20                                            #(mol CO2 mol-1 Rubisco s-1)
    KMC25  =    460                                            #(ppm CO2)
    KMO25  =     33                                            #(% O2)
    KOKC   =      0.21                                         #(-)
    O2     =     21                                            #(% O2)
    R      =      8.314                                        #(J K-1 mol-1)
    CO2I   = 0.7 * CO2A                                        #(ppm CO2)
    VCMAX  = RUBISCN * KC25 * math.exp((1/298-1/(T+273))*EAVCMX/R) #(micromol CO2 m-2 s-1)
    KMC    =          KMC25 * math.exp((1/298-1/(T+273))*EAKMC /R)  #(ppm CO2)
    KMO    =          KMO25 * math.exp((1/298-1/(T+273))*EAKMO /R)  #(% O2)
    GAMMAX = 0.5 * KOKC * KMC * O2 / KMO                       #(ppm CO2)
    PMAX   = VCMAX * (CO2I-GAMMAX) / (CO2I + KMC * (1+O2/KMO)) #(micromol CO2 m-2 s-1)
    TMPFAC = max( 0, min( 1, (T+4)/5) )                   #(-)
    EFF    = TMPFAC * (1/2.1) * (CO2I-GAMMAX) / (4.5*CO2I+10.5*GAMMAX) #(mol CO2 mol-1 PAR quanta)
    LUEMXQ = EFF*PMAX*(1+KLUETILG*(1-FRACTV)) / (EFF*K*PARAV + PMAX) #(mol CO2 mol-1 PAR quanta)Aug 8    
    return(LUEMXQ)

  
def HardeningSink(reHardRedEnd,reHardRedDay,doy,Tsurf,THARDMX,LT50,LT50MN,Hparam,CLV,KRESPHARD,RESNOR):
    reHardRedStart =  reHardRedEnd - reHardRedDay % 365
    doySinceStart  =  doy-reHardRedStart          % 365
    if doySinceStart < (reHardRedDay+0.5*(365-reHardRedDay)) :
        reHardPeriod = max(0, 1-doySinceStart/reHardRedDay )
    else:
        reHardPeriod = 1
  
    if Tsurf>THARDMX or LT50<LT50MN:
        RATEH = 0
    else:
        RATEH = reHardPeriod * Hparam * (THARDMX-Tsurf) * (LT50-LT50MN)

    RESPHARDSI = RATEH * CLV * KRESPHARD * max(0,min(1, RESNOR*5))
    return(RESPHARDSI)


def Growth(PARINT,TRANRF,LUEMXQ,NOHARV,CRES,TCRES,DAVTMP,RESPHARDSI,CRESMX,DELT,TILG2,CST,CSTAVM,SIMAX1T,PHENRF,NELLVM,LERV,TILV,LFWIDV,LERG,LFWIDG,SHAPE,SLANEW,YG):
    PHOT     = PARINT * TRANRF * 12 * LUEMXQ * NOHARV
    RESMOB   = (CRES * NOHARV / TCRES) * max(0,min( 1,DAVTMP/5))
    SOURCE   = RESMOB + PHOT
    RESPHARD = min(SOURCE,RESPHARDSI)
    ALLOTOT  = SOURCE - RESPHARD
    GRESSI   = 0.5 * (RESMOB + max(0,(CRESMX-CRES)/DELT))
    if TILG2 != 0.0: 
        CSTAV  = CST/TILG2 
    else :
        CSTAV  = 0
  
    SINK1T   = max(0, 1 - (CSTAV/CSTAVM)) * SIMAX1T
    NELLVG   = PHENRF * NELLVM 
    GLAISI   = ((LERV*TILV*NELLVM*LFWIDV) + (LERG*TILG2*NELLVG*LFWIDG)) * SHAPE * TRANRF
    GLVSI    = (GLAISI * NOHARV / SLANEW) / YG
    GSTSI    = (SINK1T * TILG2 * TRANRF * NOHARV) / YG    
    return (PHOT,RESMOB,SOURCE,RESPHARD,ALLOTOT,GRESSI,CSTAV,SINK1T,NELLVG,GLAISI,GLVSI,GSTSI)

def Allocation(GLVSI,GSTSI,DAYLGE,ALLOTOT,GRESSI,YG):
    GSHSI = GLVSI + GSTSI
    if DAYLGE >= 0.1 :
    # Situation 1: Growth has priority over storage (spring and growth period)
    # Calculate amount of assimilates allocated to shoot
        ALLOSH = min(ALLOTOT, GSHSI )
        # Calculate amount of assimilates allocated to reserves    
        GRES   = min(ALLOTOT - ALLOSH, GRESSI)
    else:
        # Situation 2: Storage has priority over shoot (autumn)
        # Calculate amount of assimilates allocated to reserves
        GRES   = min( ALLOTOT, GRESSI )
        # Calculate amount of assimilates allocated to shoot
        ALLOSH = min( ALLOTOT - GRES, GSHSI )
        
    # All surplus carbohydrate goes to roots
    ALLORT  = ALLOTOT - ALLOSH - GRES
    if GSHSI == 0:
        GSHSI = 1
    ALLOLV  = GLVSI * (ALLOSH / GSHSI)
    ALLOST  = GSTSI * (ALLOSH / GSHSI)
    GLV     = ALLOLV * YG
    GST     = ALLOST * YG
    GRT     = ALLORT * YG
    RESPGSH = (ALLOLV + ALLOST) * (1-YG)
    RESPGRT =  ALLORT           * (1-YG)
    return (GRES,ALLOLV,ALLOST,GLV,GST,GRT,RESPGSH,RESPGRT)
    
def PlantRespiration(FO2,FO2MX,RESPGRT,RESPGSH,RESPHARD):
    fAer      = max(0,min(1, FO2/FO2MX ))
    RplantAer = fAer * ( RESPGRT + RESPGSH + RESPHARD )
    return(RplantAer)

def Senescence(LAI,LAICR,RDRSCO,RDRSMX,RDRTMIN,RDRTEM,Tsurf,NOHARV,RDRFROST,RDRTOX,CLV,CSTUB,RDRSTUB,TILV,CRT,RDRROOT):
    if (LAI < LAICR) :
        TV1 = 0.0 
    else :
        TV1 = RDRSCO*(LAI-LAICR)/LAICR

    RDRS   = min(TV1, RDRSMX)
    RDRT   = max(RDRTMIN, RDRTEM * Tsurf)
    TV2    = NOHARV * max(RDRS,RDRT,RDRFROST,RDRTOX)
    TV2TIL = NOHARV * max(RDRS,     RDRFROST,RDRTOX)
    DLAI   = LAI    * TV2
    DLV    = CLV    * TV2
    DSTUB  = CSTUB  * RDRSTUB
    DTILV  = TILV   * TV2TIL
    DRT    = CRT    * RDRROOT
    return(RDRS,RDRT,TV2,TV2TIL,DLAI,DLV,DSTUB,DTILV,DRT)

def AnaerobicDamage(PERMgas,TANAER,DELT,LDT50A,LDT50B,LT50,KRDRANAER):
    if PERMgas==0 :
        dTANAER = 1
    else:
        dTANAER = -TANAER / DELT

    LD50 = LDT50A + LDT50B * LT50
    if (TANAER > 0) :
        RDRTOX = KRDRANAER / (1+math.exp(-KRDRANAER*(TANAER-LD50)))
    else:
        RDRTOX = 0
    return (dTANAER,RDRTOX)

def Hardening(KRSR3H,Tsurf,LT50,Dparam,LT50MX,TsurfDiff,DELT,RATEDMX,RESPHARD,CLV,KRESPHARD):
    RSR3H      = 1 / (1+math.exp(-KRSR3H*(Tsurf-LT50)))
    # RDRFROST should be less than 1 to avoid numerical problems
    # (loss of all biomass but keeping positive reserves). We cap it at 0.5.
    RSRDAY     = RSR3H # In previous versions we had RSRDAY = RSR3H^8 which understimated survival
    RDRFROST   = min( 0.5, 1 - RSRDAY )
    RATED      = min( Dparam*(LT50MX-LT50)*(Tsurf+TsurfDiff), (LT50MX-LT50)/DELT )
    DeHardRate = max(0,min( RATEDMX, RATED ))
    HardRate   = RESPHARD / (CLV * KRESPHARD)
    return (RSRDAY,RDRFROST,RATED,DeHardRate,HardRate)


def Foliage2(SLANEW,GLV,Tsurf,TBASE,PHY,NOHARV,TRANRF,DAYLGE,FRACTV,PHENRF,LAITIL,LAIEFT,LAI,FSMAX,RESNOR,TILV,DAVTMP,TOPTGE,TV2TIL,RGENMX,VERN,DAYL,DAYLG1G2,TILG1,RGRTG1G2):
    GLAI    = SLANEW * GLV
    if Tsurf < TBASE: 
        TV1   = 0 
    else :
        TV1   = Tsurf/PHY

    RLEAF   = TV1 * NOHARV * TRANRF * DAYLGE * ( FRACTV + PHENRF*(1-FRACTV) )
    TV2     = LAITIL - LAIEFT*LAI
    if TV2 > FSMAX:
        TV2 = FSMAX
    if TV2 < 0:
        TV2 = 0
    RGRTV   = max( 0       , TV2 * RESNOR * RLEAF )
    GTILV   = TILV  * RGRTV
    TGE     = max( 0       , 1 - (abs(DAVTMP - TOPTGE))/(TOPTGE-TBASE))
    RGRTVG1 = min( 1-TV2TIL, NOHARV * DAYLGE * TGE * RGENMX ) * VERN
    TILVG1  = TILV  * RGRTVG1
    if (DAYL > DAYLG1G2) :
        TILG1G2 = TILG1 * RGRTG1G2
    else:
        TILG1G2 = 0
    return(GLAI,RLEAF,GTILV,TILVG1,TILG1G2)

