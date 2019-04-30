# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 21:19:01 2019
@author: Fredy
"""

import math
import soil as soil_
import environment as env
import resources as res
import plant as plan


#def mainBasgra(PARAMS,MATRIX_WEATHER,DAYS_HARVEST,NDAYS,NOUT,y):
def mainBasgra(PARAMS,MATRIX_WEATHER,ASING,DAYS_HARVEST):
    #-------------------------------------------------------------------------------
    # This is the BASic GRAss model originally written in MATLAB/Simulink by Marcel
    # van Oijen, Mats Hoglind, Stig Morten Thorsen and Ad Schapendonk.
    # 2011-07-13: Translation to FORTRAN by David Cameron and Marcel van Oijen.
    # 2014-03-17: Extra category of tillers added
    # 2014-04-03: Vernalization added
    # 2014-04-03: Lower limit of temperature-driven leaf senescence no longer zero
    #-------------------------------------------------------------------------------

    # Parameters
    pa = [i[ASING] for i in PARAMS] 
    
    LOG10CLVI  = pa[1]
    LOG10CRESI = pa[2]
    LOG10CRTI  = pa[3]
    CSTI	   = pa[4]
    LOG10LAII  = pa[5]	   
    PHENI	   = pa[6] 
    TILTOTI	   = pa[7] 
    FRTILGI	   = pa[8]
    LT50I      = pa[9]
    
    # Process parameters 
    CLAIV     = pa[10]	   
    COCRESMX  =	pa[11]
    CSTAVM	  = pa[12]
    DAYLB	  =	pa[13]  
    DAYLP	  =	pa[14]
    DLMXGE	  = pa[15]
    FSLAMIN   = pa[16]	 
    FSMAX     = pa[17]  
    HAGERE    =	pa[18]	  
    K         =	pa[19]  
    LAICR	  = pa[20]	 
    LAIEFT    = pa[21]	   
    LAITIL	  =	pa[22] 
    LFWIDG	  =	pa[23]
    LFWIDV	  = pa[24] 
    NELLVM	  = pa[25]	 
    PHENCR    = pa[26]	  
    PHY	      =	pa[27]	  
    RDRSCO	  =	pa[28]  
    RDRSMX	  = pa[29]	 
    RDRTEM    = pa[30] 
    RGENMX	  =	pa[31]  
    ROOTDM	  =	pa[32]	  
    RRDMAX	  = pa[33]
    RUBISC    = pa[34]	   
    SHAPE	  =	pa[35]	  
    SIMAX1T	  =	pa[36]
    SLAMAX    = pa[37]
    TBASE     = pa[38]	   
    TCRES     = pa[39]	  
    TOPTGE	  =	pa[40]	  
    TRANCO	  = pa[41]	 
    YG        = pa[42]
    
    LAT       = pa[43]
    WCI       = pa[44]
    FWCAD     = pa[45]
    FWCWP     = pa[46]
    FWCFC     = pa[47]
    FWCWET    = pa[48]
    WCST      = pa[49]
    WpoolMax  = pa[50]
    
    Dparam	     = pa[51]			
    FGAS	     = pa[52]			
    FO2MX	     = pa[53]			
    gamma	     = pa[54]			
    Hparam	     = pa[55]			
    KRDRANAER    = pa[56]			
    KRESPHARD    = pa[57]			
    KRSR3H	     = pa[58]			
    KRTOTAER     = pa[59]			
    KSNOW	     = pa[60]			
    LAMBDAsoil   = pa[61]			
    LDT50A	     = pa[62]			
    LDT50B	     = pa[63]			
    LT50MN	     = pa[64]			
    LT50MX	     = pa[65]			
    RATEDMX	     = pa[66]			
    reHardRedDay = pa[67]			
    RHOnewSnow	 = pa[68]			
    RHOpack	     = pa[69]			
    SWret	     = pa[70]			
    SWrf	     = pa[71]			
    THARDMX	     = pa[72]			
    TmeltFreeze	 = pa[73]			
    TrainSnow	 = pa[74]			
    TsurfDiff	 = pa[75]
    KLUETILG	 = pa[76]
    FRTILGG1I	 = pa[77]
    DAYLG1G2     = pa[78]
    RGRTG1G2     = pa[79]
    RDRTMIN      = pa[80]
    TVERN        = pa[81]
    
    # Parameter transformations 
    CLVI  = math.pow(10,LOG10CLVI)
    CRESI = math.pow(10,LOG10CRESI)
    CRTI  = math.pow(10,LOG10CRTI)
    LAII  = math.pow(10,LOG10LAII)
    
    WCAD  = FWCAD  * WCST
    WCWP  = FWCWP  * WCST
    WCFC  = FWCFC  * WCST
    WCWET = FWCWET * WCST
         
    # Calendar & weather
    YEARI  = [i[1] for i in MATRIX_WEATHER]    
    DOYI   = [i[2] for i in MATRIX_WEATHER]
    TI     = [i[3] for i in MATRIX_WEATHER]
    GRI    = [i[9] for i in MATRIX_WEATHER]
    TMMNI  = [i[4] for i in MATRIX_WEATHER]
    TMMXI  = [i[5] for i in MATRIX_WEATHER]

    #ifdef weathergen  
    #  RAINI = MATRIX_WEATHER(:,6)
    #  PETI  = MATRIX_WEATHER(:,7)
    #else
    RHI = [i[6] for i in MATRIX_WEATHER]
    
    RAINI = [i[7] for i in MATRIX_WEATHER]
    WNI   = [i[8] for i in MATRIX_WEATHER]
    #endif
    
    #module parameters_plant
    # Initial constants
    CLVDI  = 0
    YIELDI = 0
    CSTUBI = 0
    # Process parameters
    RDRROOT      =  0.001
    RDRSTUB      =  0.2
    reHardRedEnd = 91    
    
    #module parameters_site    
    # Simulation period and time step
    DELT   =   1.0
    # Atmospheric conditions
    CO2A   = 350
    # Soil
    DRATE  =  50
    # Soil - WINTER PARAMETERS
    LAMBDAice      = 1.9354e+005
    LatentHeat     = 335000
    poolInfilLimit =      0.2
    RHOwater       =   1000
    
    # Soil initial constants
    DRYSTORI = 0
    FdepthI  = 0
    SDEPTHI  = 0
    TANAERI  = 0
    WAPLI    = 0
    WAPSI    = 0
    WASI     = 0
    WETSTORI = 0
    
    # Management: harvest dates and irrigation
    IRRIGF = 0 
    
    # Mathematical constants
    pi   = 3.141592653589793
    Freq = 2*pi / 365
    Kmin = 4
    Ampl = 0.625
    Bias = Kmin + Ampl
    
    # Initial constants
    CLV     = CLVI
    CLVD    = CLVDI
    YIELD   = YIELDI
    YIELD_LAST = YIELDI
    CRES    = CRESI
    CRT     = CRTI
    CST     = CSTI
    CSTUB   = CSTUBI
    DRYSTOR = DRYSTORI
    Fdepth  = FdepthI
    LAI     = LAII
    LT50    = LT50I
    O2      = FGAS * ROOTDM * FO2MX * 1000/22.4
    PHEN    = PHENI
    ROOTD   = ROOTDM
    Sdepth  = SDEPTHI
    TANAER  = TANAERI
    TILG1   = TILTOTI *       FRTILGI *    FRTILGG1I
    TILG2   = TILTOTI *       FRTILGI * (1-FRTILGG1I)
    TILV    = TILTOTI * (1 - FRTILGI)
    VERN    = 0
    WAL     = 1000 * ROOTDM * WCI
    WAPL    = WAPLI
    WAPS    = WAPSI
    WAS     = WASI
    WETSTOR = WETSTORI
    
    y=[]
    for day in range(1,len(DOYI)):
        y.append([])
  
    for day in range(1,len(DOYI)):
    
      # Environment
      #set_weather_day()
        year   = float(YEARI[day]) # day of the year (d)        
        doy    = float(DOYI[day])  # day of the year (d)
        RAIN   = float(RAINI[day]) # precipitation (mm d-1)	
        GR     = float(GRI[day])   # irradiation (MJ m-2 d-1)	
        TMMN   = float(TMMNI[day]) # minimum (or average) temperature (degrees Celsius)
        TMMX   = float(TMMXI[day]) # maximum (or average) temperature (degrees Celsius)
        
        # kPa
        T      = float(TI[day])    

        VP     = math.exp(17.27*T/(T+239)) * 0.6108 * float(RHI[day]) / 100   # vapour pressure (kPa)
        
        WN     = float(WNI[day])   # mean wind speed (m s-1)
        DAVTMP = T ############
        DTR    = GR * math.exp(-KSNOW*DRYSTOR)
        PAR    = 0.5*4.56*DTR
      
        WCL =                soil_.SoilWaterContent(Fdepth,ROOTD,WAL)        
        Tsurf,fPerm =        soil_.Physics(DAVTMP,Fdepth,Sdepth,gamma)
        WCeff,PFrate,Frate = soil_.FrozenSoil(Fdepth,ROOTD,WAS,WCFC,WCL,Tsurf,LAMBDAsoil,RHOwater,LatentHeat,DELT)
        
        Pwater,Psnow = env.precForm(DAVTMP,TrainSnow,RAIN)
        
        #WaterSnow 
        SnowMelt,WmaxStore = env.SnowMeltWmaxStore(Bias,Ampl,Freq,doy,DAVTMP,TmeltFreeze,DRYSTOR,DELT,SWret)
        reFreeze,StayWet   = env.WETSTORdynamics(SWrf,TmeltFreeze,DAVTMP,WETSTOR,DELT)
        Wremain,Wsupply    = env.LiquidWaterDistribution(StayWet,SnowMelt,Pwater,WmaxStore)
        DENSITY            = env.SnowDensity(DRYSTOR,WETSTOR,Sdepth)
        PackMelt           = env.SnowDepthDecrease(Sdepth,DELT, RHOpack,SnowMelt,DENSITY)
    
        RNINTC,PINFIL      = env.RainSnowSurfacePool(Wsupply,LAI)
        INFIL,runOn        = env.INFILrunOn(Fdepth,poolInfilLimit,PINFIL)
        
        poolInfil,poolDrain,PIrate,FREEZEPL,THAWPS = env.SurfacePool(WpoolMax,WAPL,WAPS,runOn,DELT,Fdepth,poolInfilLimit,Frate,Tsurf,LAMBDAice,RHOwater,LatentHeat)
          
        PERMgas            = env.MicroClimate(WAPS)
        DAYL               = env.DDAYL(pi,doy,LAT)
        
        #ifdef weathergen  
        #  call PEVAPINPUT ()
        #else        
        PEVAP, PTRAN       = env.PENMAN(DTR,DAVTMP,VP,WN,LAI,RNINTC)
        #endif
        
        # Resources
        PARAV,PARINT,DTRINT = res.Light(DAYL,PAR,DTR,K,LAI)
        WCL,FR,AVAILF,EVAP,TRAN,TRANRF = res.EVAPTRTRF(Fdepth,ROOTD,WAL,WCAD,PEVAP,WCFC,WCWP,PTRAN,TRANCO,WCST,WCWET,DELT)
        WCL,RROOTD,EXPLOR   = res.ROOTDG(Fdepth,ROOTD,WAL,ROOTDM,WCWP,RRDMAX,DELT,WCFC)
        
        # Soil
        FREEZEL,THAWS,DRAIN,RUNOFF,IRRIG = soil_.FRDRUNIR(WCFC,ROOTD,Fdepth,WCST,INFIL,poolDrain,WAL,DELT,EVAP,TRAN,Frate,WAS,DRATE,IRRIGF)
        FO2                 = soil_.O2status(O2,ROOTD,FGAS)
        
        # Plant
        NOHARV,FRACTV,HARVLA,HARVLV,HARVPH,TV1,HARVRE,HARVST,GSTUB,HARVTILG2 = plan.Harvest(year,doy,DAYS_HARVEST,TILV,TILG2,CLAIV,LAI,DELT,CLV,PHEN,HAGERE,CRES,CST)
        CRESMX, RESNOR = plan.Biomass(COCRESMX,CLV,CRES,CST)
        GPHEN,DPHEN,PHENRF,DAYLGE = plan.Phenology(DAVTMP,DAYLP,DAYL,DAYLB,PHEN,DELT,PHENCR,DLMXGE)
        EFFTMP,LERV,LERG,SLAMIN,SLANEW = plan.Foliage1(TBASE,DAVTMP,DAYLGE,SLAMAX,FSLAMIN,RESNOR)
        LUEMXQ = plan.LUECO2TM(DAVTMP,RUBISC,CO2A,PARAV,KLUETILG,FRACTV,K)
        RESPHARDSI = plan.HardeningSink(reHardRedEnd,reHardRedDay,doy,Tsurf,THARDMX,LT50,LT50MN,Hparam,CLV,KRESPHARD,RESNOR)
        PHOT,RESMOB,SOURCE,RESPHARD,ALLOTOT,GRESSI,CSTAV,SINK1T,NELLVG,GLAISI,GLVSI,GSTSI = plan.Growth(PARINT,TRANRF,LUEMXQ,NOHARV,CRES,TCRES,DAVTMP,RESPHARDSI,CRESMX,DELT,TILG2,CST,CSTAVM,SIMAX1T,PHENRF,NELLVM,LERV,TILV,LFWIDV,LERG,LFWIDG,SHAPE,SLANEW,YG)
        GRES,ALLOLV,ALLOST,GLV,GST,GRT,RESPGSH,RESPGRT = plan.Allocation(GLVSI,GSTSI,DAYLGE,ALLOTOT,GRESSI,YG)
        RplantAer = plan.PlantRespiration(FO2,FO2MX,RESPGRT,RESPGSH,RESPHARD)
        dTANAER,RDRTOX = plan.AnaerobicDamage(PERMgas,TANAER,DELT,LDT50A,LDT50B,LT50,KRDRANAER)
        RSRDAY,RDRFROST,RATED,DeHardRate,HardRate = plan.Hardening(KRSR3H,Tsurf,LT50,Dparam,LT50MX,TsurfDiff,DELT,RATEDMX,RESPHARD,CLV,KRESPHARD)
        RDRS,RDRT,TV2,TV2TIL,DLAI,DLV,DSTUB,DTILV,DRT = plan.Senescence(LAI,LAICR,RDRSCO,RDRSMX,RDRTMIN,RDRTEM,Tsurf,NOHARV,RDRFROST,RDRTOX,CLV,CSTUB,RDRSTUB,TILV,CRT,RDRROOT)
        GLAI,RLEAF,GTILV,TILVG1,TILG1G2 = plan.Foliage2(SLANEW,GLV,Tsurf,TBASE,PHY,NOHARV,TRANRF,DAYLGE,FRACTV,PHENRF,LAITIL,LAIEFT,LAI,FSMAX,RESNOR,TILV,DAVTMP,TOPTGE,TV2TIL,RGENMX,VERN,DAYL,DAYLG1G2,TILG1,RGRTG1G2)
        
        # Soil 2
        O2OUT,O2IN = soil_.O2fluxes(RplantAer,KRTOTAER,FO2MX,ROOTD,FGAS,PERMgas,O2,DELT)

        #================
        # Outputs
        #================
        y[1].append(year + (doy-0.5)/366) # "Time" = Decimal year [approximation]
        y[2].append(year)
        y[3].append(doy)
        y[4].append(DAVTMP)
        
        y[5].append(CLV)
        y[6].append(CLVD)
        y[7].append(YIELD_LAST)
        y[8].append(CRES)
        y[9].append(CRT)
        y[10].append(CST)
        y[11].append(CSTUB)
        y[12].append(DRYSTOR)
        y[13].append(Fdepth)
        y[14].append(LAI)
        y[15].append(LT50)
        y[16].append(O2)
        y[17].append(PHEN)
        y[18].append(ROOTD)
        y[19].append(Sdepth)
        y[20].append(TANAER)
        y[21].append(TILG1 + TILG2)
        y[22].append(TILV)
        y[23].append(WAL)
        y[24].append(WAPL)
        y[25].append(WAPS)
        y[26].append(WAS)
        y[27].append(WETSTOR)
        
        # Extra derived variables for calibration
        y[28].append((CLV+CST+CSTUB)/0.45 + CRES/0.40)                                 # "DM"      = Aboveground dry matter in g m-2
        y[29].append((CRES/0.40) / ((CLV+CST+CSTUB)/0.45 + CRES/0.40))                 # "RES"     = Reserves in g g-1 aboveground dry matter
        y[30].append(LERG)                               #
        y[31].append(NELLVG)                             #
        y[32].append(RLEAF)                              #
        y[33].append(LAI / (CLV/0.45))                   # "SLA"     = m2 leaf area g-1 dry matter vegetative tillers
        y[34].append(TILG1 + TILG2 + TILV)               # "TILTOT"  = Total tiller number in # m-2
        y[35].append((TILG1+TILG2) / (TILG1+TILG2+TILV)) # "FRTILG"  = Fraction of tillers that is generative
        y[36].append( TILG1        / (TILG1+TILG2+TILV)) # "FRTILG1" = Fraction of tillers that is in TILG1
        y[37].append(       TILG2  / (TILG1+TILG2+TILV)) # "FRTILG2" = Fraction of tillers that is in TILG2
        y[38].append(RDRT)
        y[39].append(VERN)
        
        CLV     = CLV     + GLV   - DLV    - HARVLV
        CLVD    = CLVD            + DLV
        YIELD   = (HARVLV + HARVST*HAGERE) / 0.45 + HARVRE/0.40
        if (YIELD>0):
            YIELD_LAST = YIELD
        CRES    = CRES    + GRES  - RESMOB - HARVRE
        CRT     = CRT     + GRT   - DRT
        CST     = CST     + GST           - HARVST
        CSTUB   = CSTUB   + GSTUB - DSTUB
        DRYSTOR = DRYSTOR + reFreeze + Psnow - SnowMelt
        Fdepth  = Fdepth  + Frate
        LAI     = LAI     + GLAI - DLAI   - HARVLA
        LT50    = LT50    + DeHardRate - HardRate
        O2      = O2      + O2IN - O2OUT
        PHEN    = min(1, PHEN + GPHEN - DPHEN - HARVPH)
        ROOTD   = ROOTD   + RROOTD
        Sdepth  = Sdepth  + Psnow/RHOnewSnow - PackMelt
        TANAER  = TANAER  + dTANAER
        TILG1   = TILG1           + TILVG1 - TILG1G2
        TILG2   = TILG2                    + TILG1G2 - HARVTILG2
        TILV    = TILV    + GTILV - TILVG1           - DTILV   
        if (DAVTMP<TVERN):
            VERN = 1
        WAL     = WAL  + THAWS  - FREEZEL  + poolDrain + INFIL +EXPLOR+IRRIG-DRAIN-RUNOFF-EVAP-TRAN
        WAPL    = WAPL + THAWPS - FREEZEPL + poolInfil - poolDrain
        WAPS    = WAPS - THAWPS + FREEZEPL
        WAS     = WAS  - THAWS  + FREEZEL
        WETSTOR = WETSTOR + Wremain - WETSTOR
        
    return(y)


  