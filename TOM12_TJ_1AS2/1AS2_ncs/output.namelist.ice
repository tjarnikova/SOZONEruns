&NAMICERUN
 JPL=5          ,
 NLAY_I=2          ,
 NLAY_S=1          ,
 CN_ICERST_IN="restart_ice_in                  ",
 CN_ICERST_INDIR=".                                                                                                                                                                                                                                                               ",
 CN_ICERST_OUT="restart_ice                     ",
 CN_ICERST_OUTDIR=".                                                                                                                                                                                                                                                               ",
 LN_LIMDYN=T,
 RN_AMAX= 0.99900000000000000     ,
 LN_LIMDIAHSB=F,
 LN_LIMDIAOUT=T,
 LN_ICECTL=F,
 IICEPRT=10         ,
 JICEPRT=10         ,
 /
&NAMICEITD
 NN_CATBND=2          ,
 RN_HIMEAN=  2.0000000000000000     ,
 /
&NAMICEHDF
 NN_CONVFRQ=5          ,
 /
&NAMICETHD
 RN_HNEWICE= 0.10000000000000001     ,
 LN_FRAZIL=F,
 RN_MAXFRAZB=  1.0000000000000000     ,
 RN_VFRAZB= 0.41699999999999998     ,
 RN_CFRAZB=  5.0000000000000000     ,
 RN_HIMIN= 0.10000000000000001     ,
 RN_BETAS= 0.66000000000000003     ,
 RN_KAPPA_I=  1.0000000000000000     ,
 NN_CONV_DIF=  50.000000000000000     ,
 RN_TERR_DIF=  1.0000000000000000E-004,
 NN_ICE_THCON=1          ,
 NN_MONOCAT=0          ,
 LN_IT_QNSICE=T,
 /
&NAMICESAL
 NN_ICESAL=2          ,
 RN_ICESAL=  4.0000000000000000     ,
 RN_SAL_GD=  5.0000000000000000     ,
 RN_TIME_GD=  1730000.0000000000     ,
 RN_SAL_FL=  2.0000000000000000     ,
 RN_TIME_FL=  864000.00000000000     ,
 RN_SIMAX=  20.000000000000000     ,
 RN_SIMIN= 0.10000000000000001     ,
 /
&NAMICEITDME
 RN_CS= 0.50000000000000000     ,
 RN_FSNOWRDG= 0.50000000000000000     ,
 RN_FSNOWRFT= 0.50000000000000000     ,
 RN_GSTAR= 0.14999999999999999     ,
 RN_ASTAR=  5.0000000000000003E-002,
 RN_HSTAR=  100.00000000000000     ,
 LN_RAFTING=T,
 RN_HRAFT= 0.75000000000000000     ,
 RN_CRAFT=  5.0000000000000000     ,
 RN_POR_RDG= 0.29999999999999999     ,
 NN_PARTFUN=1          ,
 /
