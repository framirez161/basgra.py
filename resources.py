import math

def Light(DAYL,PAR,DTR,K,LAI):
    if DAYL > 0 :
        PARAV = PAR * (1E6/(24*3600)) / DAYL
    else:
        PARAV = 0
    PARINT = PAR * (1 - math.exp(-1.0*K*LAI))
    DTRINT = DTR * (1 - math.exp(-0.75*K*LAI))
    return (PARAV,PARINT,DTRINT)
  
def EVAPTRTRF(Fdepth,ROOTD,WAL,WCAD,PEVAP,WCFC,WCWP,PTRAN,TRANCO,WCST,WCWET,DELT):
    if (Fdepth < ROOTD) :
        WCL = WAL*0.001 / (ROOTD-Fdepth)
    else:
        WCL = 0                                                 
    WAAD = 1000 * WCAD * (ROOTD-Fdepth)                          # (mm)
    EVAP = PEVAP * max(0, min(1, (WCL-WCAD)/(WCFC-WCAD) ))      # (mm d-1)
    WCCR = WCWP + max( 0.01, PTRAN/(PTRAN+TRANCO) * (WCFC-WCWP) ) # (m3 m-3)
    
    if WCL > WCCR and WCST != WCWET:
        FR = max(0, min(1, (WCST-WCL)/(WCST- WCWET)))       
    else:
        FR = max(0, min(1, (WCL-WCWP)/(WCCR-WCWP)))              
                                                     # (mm mm-1)
    TRAN = PTRAN * FR                                # (mm d-1)
    if EVAP+TRAN > 0 :
        AVAILF = min( 1, ((WAL-WAAD)/DELT) / (EVAP+TRAN) )         
    else:
        AVAILF = 0                                                
                                                       # (mm mm-1)
    EVAP = EVAP * AVAILF                               # (mm d-1)
    TRAN = TRAN * AVAILF                               # (mm d-1)
    if (PTRAN > 0.) :
        TRANRF = TRAN / PTRAN                          # (-)
    else:
        TRANRF = 1                                     # (-)

    return(WCL,FR,AVAILF,EVAP,TRAN,TRANRF)

def ROOTDG(Fdepth,ROOTD,WAL,ROOTDM,WCWP,RRDMAX,DELT,WCFC):
    if (Fdepth < ROOTD) :
        WCL = WAL*0.001 / (ROOTD-Fdepth)
    else:
        WCL = 0
                                                     # (m3 m-3)
    if ROOTD < ROOTDM and WCL > WCWP:
        RROOTD = min( RRDMAX, (ROOTDM-ROOTD)/DELT )
    else:
        RROOTD = 0

    EXPLOR = 1000 * RROOTD * WCFC    
    return (WCL,RROOTD,EXPLOR)

