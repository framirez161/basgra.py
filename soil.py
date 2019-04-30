import math

def SoilWaterContent(Fdepth,ROOTD,WAL):
    if Fdepth < ROOTD:
        WCL = WAL*0.001 / (ROOTD-Fdepth)
    else :
        WCL = 0
    return  WCL 

def Physics(DAVTMP,Fdepth,Sdepth,gamma):
    if Fdepth > 0:
        Tsurf = DAVTMP / (1 + 10 * (Sdepth / Fdepth) )
        fPerm = 0
    else :
        Tsurf = DAVTMP * math.exp(-gamma*Sdepth)
        fPerm = 1
    return (Tsurf,fPerm)
    
def FrozenSoil(Fdepth,ROOTD,WAS,WCFC,WCL,Tsurf,LAMBDAsoil,RHOwater,LatentHeat,DELT):
     # Determining the amount of solid water that contributes in transportation of heat to surface 'WCeff'
     if Fdepth > ROOTD:
       WCeff = WCFC
     elif Fdepth > 0:
       WCeff = (0.001*WAS) / Fdepth
     else:
       WCeff = WCL

     # Calculating potential frost rate 'PFrate'
     #  if ((Fdepth == 0.).and.(Tsurf>0.)) : # No soil frost present AND no frost starting
     if Fdepth == 0 and Tsurf>0 or WCeff == 0: # No soil frost present AND no frost starting
       PFrate = 0
     else :
       alpha  = LAMBDAsoil / ( RHOwater * WCeff * LatentHeat)
       PFrate = math.sqrt( math.max(0, math.pow(Fdepth,2) - 2*alpha*Tsurf) ) - Fdepth

     if PFrate >= 0 and Fdepth > 0 and Fdepth < ROOTD:
       Frate = PFrate * (0.001*WAS/Fdepth) / WCFC # Soil frost increasing
     elif (PFrate+Fdepth/DELT) < 0 :
       Frate = -Fdepth / DELT                     # Remaining soil frost thaws away
     else :
       Frate = PFrate
       
     return (WCeff,PFrate,Frate)  

 
def FRDRUNIR(WCFC,ROOTD,Fdepth,WCST,INFIL,poolDrain,WAL,DELT,EVAP,TRAN,Frate,WAS,DRATE,IRRIGF):
    WAFC   = 1000 * WCFC * max(0,(ROOTD-Fdepth))                      # (mm)
    WAST   = 1000 * WCST * max(0,(ROOTD-Fdepth))                      # (mm)
    INFILTOT = INFIL + poolDrain
    if Fdepth < ROOTD :
        FREEZEL = max(0, min( WAL/DELT + (INFILTOT - EVAP - TRAN), (Frate/(ROOTD-Fdepth))*WAL))                 # (mm d-1)
    else:
        REEZEL = 0                                                     # (mm d-1)
    if Fdepth > 0 and Fdepth <= ROOTD :
        THAWS   = max(0,min( WAS/DELT, -Frate*WAS/Fdepth ))              # (mm d-1)
    else :
        THAWS   = 0                                                      # (mm d-1)
        DRAIN  = max(0,min( DRATE, (WAL-WAFC)/DELT + (INFILTOT - EVAP - TRAN - FREEZEL + THAWS) ))                # (mm d-1)
        RUNOFF = max(0,            (WAL-WAST)/DELT + (INFILTOT - EVAP - TRAN - FREEZEL + THAWS - DRAIN) )         # (mm d-1)
        IRRIG  = IRRIGF *  (       (WAFC-WAL)/DELT - (INFILTOT - EVAP - TRAN - FREEZEL + THAWS - DRAIN - RUNOFF)) # (mm d-1)
        
    return (FREEZEL,THAWS,DRAIN,RUNOFF,IRRIG)

def O2status(O2,ROOTD,FGAS):
    FO2 = O2 / (ROOTD * FGAS * 1000/22.4)
    return (FO2)
  
def O2fluxes(RplantAer,KRTOTAER,FO2MX,ROOTD,FGAS,PERMgas,O2,DELT):
    O2OUT = RplantAer * KRTOTAER * 1/12 * 1
    O2MX  = FO2MX * ROOTD * FGAS * 1000/22.4
    O2IN  = PERMgas * ( (O2MX-O2) + O2OUT*DELT )  
    return (O2OUT,O2IN)

