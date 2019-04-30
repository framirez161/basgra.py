import math

def MicroClimate(WAPS):
    if WAPS == 0:
        PERMgas = 1
    else:
        PERMgas = 0
    return (PERMgas)

def RainSnowSurfacePool(Wsupply,LAI):
    RNINTC = min( Wsupply, 0.25*LAI )
    PINFIL = Wsupply - RNINTC
    return (RNINTC,PINFIL)

def precForm(DAVTMP,TrainSnow,RAIN):
    if DAVTMP > TrainSnow:
      Pwater = RAIN
      Psnow  = 0
    else:
      Pwater = 0
      Psnow  = RAIN
    return (Pwater,Psnow)
  
def SnowMeltWmaxStore(Bias,Ampl,Freq,doy,DAVTMP,TmeltFreeze,DRYSTOR,DELT,SWret):

   Melt = Bias + Ampl * math.sin( Freq * (doy-(174-91)) )
   if DAVTMP > TmeltFreeze:
     SnowMelt = max(0, min( DRYSTOR/DELT, Melt*(DAVTMP-TmeltFreeze)))
   else:
     SnowMelt = 0
   WmaxStore = DRYSTOR * SWret
   return (SnowMelt,WmaxStore)

  
def WETSTORdynamics(SWrf,TmeltFreeze,DAVTMP,WETSTOR,DELT):
   reFreezeMax = SWrf * (TmeltFreeze-DAVTMP)
   if WETSTOR > 0 and DAVTMP < TmeltFreeze:
     reFreeze = min(WETSTOR/DELT,reFreezeMax)
   else:
     reFreeze = 0
   StayWet = WETSTOR/DELT - reFreeze
   return (reFreeze,StayWet)
  
def LiquidWaterDistribution(StayWet,SnowMelt,Pwater,WmaxStore):
   Wavail  = StayWet + SnowMelt + Pwater
   Wremain = min(Wavail,WmaxStore)
   Wsupply = Wavail - Wremain 
   return (Wremain,Wsupply)
   
  
def SnowDensity(DRYSTOR,WETSTOR,Sdepth):
    SWE = DRYSTOR + WETSTOR
    if Sdepth > 0:
        DENSITY = min(480, SWE/Sdepth)
    else:
        DENSITY = 0
    return (DENSITY)
  
def SnowDepthDecrease(Sdepth,DELT, RHOpack,SnowMelt,DENSITY):
    if Sdepth > 0:
        PackMelt = max(0,min( Sdepth/DELT, Sdepth*RHOpack - SnowMelt/DENSITY))
    else:
        PackMelt = 0
    return (PackMelt)
  
def INFILrunOn(Fdepth,poolInfilLimit,PINFIL):
    if Fdepth <= poolInfilLimit:
      INFIL = PINFIL
    else:
      INFIL = 0
    runOn = PINFIL - INFIL
    return (INFIL,runOn)
  
def SurfacePool(WpoolMax,WAPL,WAPS,runOn,DELT,Fdepth,poolInfilLimit,Frate,Tsurf,LAMBDAice,RHOwater,LatentHeat):
    poolVolRemain = max(0, WpoolMax - WAPL - WAPS)
    poolInfil     = min(runOn,poolVolRemain)
    #poolRUNOFF    = runOn - poolInfil
    poolWavail    = poolInfil + WAPL/DELT
    if poolWavail == 0:
      poolDrain = 0
    elif Fdepth <= poolInfilLimit:
      poolDrain = poolWavail
    else:
      poolDrain = max(0,min(-Frate*1000, poolWavail ))
    if Tsurf>0 and WAPL==0 and WAPS==0:
      PIrate    = 0
    else:
      eta       = LAMBDAice / ( RHOwater * LatentHeat )      # [m2 C-1 day-1]
      PIrate    = (math.sqrt( max(0,(math.pow(0.001*WAPS),2) - 2*eta*Tsurf*DELT)))/DELT - (0.001*WAPS)/DELT # [m day-1]
    if PIrate < 0:
      FREEZEPL  = 0
      THAWPS    = min( WAPS/DELT , -PIrate*1000)
    else:
      FREEZEPL  = max( 0,min( poolInfil + WAPL/DELT - poolDrain*DELT, PIrate*1000))
      THAWPS    = 0
    return (poolInfil,poolDrain,PIrate,FREEZEPL,THAWPS)
  
def DDAYL(pi,doy,LAT):
#=============================================================================
# Calculate day length (d d-1) from Julian day and latitude (LAT, degN)
# Author - Marcel van Oijen (CEH-Edinburgh)
#=============================================================================
    RAD  = pi / 180                                                    # (radians deg-1)
    DEC  = -math.asin (math.sin (23.45*RAD)*math.cos (2*pi*(doy+10)/365))           # (radians)
    DECC = max(math.atan(-1/math.tan(RAD*LAT)),min( math.atan( 1/math.tan(RAD*LAT)),DEC)) # (radians)
    DAYL = 0.5 * ( 1 + 2 * math.asin(math.tan(RAD*LAT)*math.tan(DECC)) / pi )        # (d d-1) 
    return (DAYL)

#def PEVAPINPUT():
#    PEVAP  =     math.exp(-0.5*LAI)  * PET                      # (mm d-1)
#    PTRAN  = (1-math.exp(-0.5*LAI)) * PET                      # (mm d-1)
#    PTRAN  = max( 0., PTRAN-0.5*RNINTC )                   # (mm d-1)

def PENMAN(DTR,DAVTMP,VP,WN,LAI,RNINTC):
  #=============================================================================
  # Calculate potential rates of evaporation and transpiration (mm d-1)
  # Inputs: LAI (m2 m-2), DTR (MJ GR m-2 d-1), RNINTC (mm d-1)
  # Inputs not in header: VP (kPa), WN (m s-1)
  # Outputs: PEVAP & PTRAN (mm d-1)
  # Author - Marcel van Oijen (CEH-Edinburgh)
  #=============================================================================
    DTRJM2 = DTR * 1E6                                    # (J GR m-2 d-1)
    BOLTZM = 5.668E-8                                      # (J m-2 s-1 K-4)
    LHVAP  = 2.4E6                                         # (J kg-1)
    PSYCH  = 0.067                                         # (kPA degC-1))    
    BBRAD  = BOLTZM * math.pow((DAVTMP+273),4) * 86400.    # (J m-2 d-1)
    SVP    = 0.611 * math.exp(17.4 * DAVTMP / (DAVTMP + 239))  # (kPa)
    SLOPE  = 4158.6 * SVP / math.pow((DAVTMP + 239),2)    # (kPA degC-1)
    RLWN   = BBRAD * max(0,0.55*(1-VP/SVP))              # (J m-2 d-1)
    NRADS  = DTRJM2 * (1-0.15) - RLWN                     # (J m-2 d-1)
    NRADC  = DTRJM2 * (1-0.25) - RLWN                     # (J m-2 d-1)
    PENMRS = NRADS * SLOPE/(SLOPE+PSYCH)                   # (J m-2 d-1)
    PENMRC = NRADC * SLOPE/(SLOPE+PSYCH)                   # (J m-2 d-1)
    WDF    = 2.63 * (1.0 + 0.54 * WN)                      # (kg m-2 d-1 kPa-1)
    PENMD  = LHVAP * WDF * (SVP-VP) * PSYCH/(SLOPE+PSYCH)  # (J m-2 d-1)
    PEVAP  =    math.exp(-0.5*LAI)  * (PENMRS + PENMD) / LHVAP # (mm d-1)
    PTRAN  = (1-math.exp(-0.5*LAI)) * (PENMRC + PENMD) / LHVAP # (mm d-1)
    PTRAN  = max( 0, PTRAN-0.5*RNINTC )                   # (mm d-1)
    return (PEVAP, PTRAN)





