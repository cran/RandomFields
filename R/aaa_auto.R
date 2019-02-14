# This file has been created automatically by 'rfGenerateConstants'


 ## from  ../../RandomFieldsUtils/RandomFieldsUtils/src/AutoRandomFieldsUtils.h



 MAXUNITS 	<- as.integer(4)
 MAXCHAR 	<- as.integer(18)
 RFOPTIONS 	<- "RFoptions"
 isGLOBAL 	<- as.integer(NA)




 ## from  src/AutoRandomFields.h


 MAXCHAR_RF 	<- as.integer(MAXCHAR)
 METHODMAXCHAR 	<- as.integer(MAXCHAR)
 MAXPARAM 	<- as.integer(20)
 MAXSUB 	<- as.integer(10)

 MAXCEDIM 	<- as.integer(13)
 MAXTBMSPDIM 	<- as.integer(4)
 MAXMPPDIM 	<- as.integer(4)
 MAXMPPVDIM 	<- as.integer(10)
 MAXVDIM 	<- as.integer(MAXMPPVDIM)
 MAXHYPERDIM 	<- as.integer(4)
 MAXVARIODIM 	<- as.integer(20)
 MAXTBMVDIM 	<- as.integer(5)
 MAXGAUSSVDIM 	<- as.integer(10)
 MAXACCEPTED 	<- as.double(1e300)
 MAXVARIANTS 	<- as.integer(6)
 MAXSYSTEMS 	<- as.integer(1)
 MAXDATANAMES 	<- as.integer(5)

 MAXBOXCOXVDIM 	<- as.integer(10)

 MAXCOOORDNAMES 	<- as.integer(4)
 units_none 	<- as.integer(0)
 units_km 	<- as.integer(1)
 units_miles 	<- as.integer(2)
 units_time 	<- as.integer(3)
 units_user 	<- as.integer(4)

 nr_units 	<- as.integer((units_user+1))


 coord_auto 	<- as.integer(0)
 coord_keep 	<- as.integer(1)
 cartesian 	<- as.integer(2)
 earth 	<- as.integer(3)
 sphere 	<- as.integer(4)
 gnomonic 	<- as.integer(5)
 orthographic 	<- as.integer(6)
 coord_mix 	<- as.integer(7)


 nr_coord_sys 	<- as.integer((coord_mix+1))

 reportcoord_always 	<- as.integer(0)
 reportcoord_warnings_orally 	<- as.integer(1)
 reportcoord_warnings 	<- as.integer(2)
 reportcoord_none 	<- as.integer(3)

 nr_reportcoord_modes 	<- as.integer((reportcoord_none+1))

 careless 	<- as.integer(0)
 sloppy 	<- as.integer(1)
 easygoing 	<- as.integer(2)
 normal 	<- as.integer(3)
 precise 	<- as.integer(4)
 pedantic 	<- as.integer(5)
 neurotic 	<- as.integer(6)

 nr_modes 	<- as.integer((neurotic+1))

 output_sp 	<- as.integer(0)
 output_rf 	<- as.integer(1)
 output_geor 	<- as.integer(2)

 nr_output_modes 	<- as.integer((output_geor+1))

 radians 	<- as.integer(0)
 degree 	<- as.integer(1)

 last_angle_mode 	<- as.integer(degree)



 PARAM_DEP 	<- as.integer(-1)
 PREVMODEL_DEP 	<- as.integer(-2)
 SUBMODEL_DEP 	<- as.integer(-3)
 MISMATCH 	<- as.integer(-4)
 UNSET 	<- as.integer(-5)



 XONLY 	<- as.integer(0)
 KERNEL 	<- as.integer(1)
 PREVMODEL_D 	<- as.integer(2)
 SUBMODEL_D 	<- as.integer(3)
 PARAMDEP_D 	<- as.integer(4)
 KEEPCOPY_DOM 	<- as.integer(5)
 DOMAIN_MISMATCH 	<- as.integer(6)


 FIRST_DOMAIN 	<- as.integer(XONLY)
 LAST_DOMAINUSER 	<- as.integer(KERNEL)
 LAST_DOMAIN 	<- as.integer(DOMAIN_MISMATCH)




RC_ISOTROPIC <- ISOTROPIC 	<- as.integer(0)
RC_DOUBLEISOTROPIC <- DOUBLEISOTROPIC 	<- as.integer(1)
 VECTORISOTROPIC 	<- as.integer(2)
 SYMMETRIC 	<- as.integer(3)
RC_CARTESIAN_COORD <- CARTESIAN_COORD 	<- as.integer(4)
RC_GNOMONIC_PROJ <- GNOMONIC_PROJ 	<- as.integer(5)
RC_ORTHOGRAPHIC_PROJ <- ORTHOGRAPHIC_PROJ 	<- as.integer(6)
 SPHERICAL_ISOTROPIC 	<- as.integer(7)
 SPHERICAL_SYMMETRIC 	<- as.integer(8)
 SPHERICAL_COORD 	<- as.integer(9)
RC_EARTH_ISOTROPIC <- EARTH_ISOTROPIC 	<- as.integer(10)
 EARTH_SYMMETRIC 	<- as.integer(11)
 EARTH_COORD 	<- as.integer(12)
 CYLINDER_COORD 	<- as.integer(13)
RC_UNREDUCED <- UNREDUCED 	<- as.integer(14)
 PREVMODEL_I 	<- as.integer(15)
 SUBMODEL_I 	<- as.integer(16)
 PARAMDEP_I 	<- as.integer(17)
 KEEPCOPY_ISO 	<- as.integer(18)
 ISO_MISMATCH 	<- as.integer(19)




 FIRST_CARTESIAN 	<- as.integer(ISOTROPIC)
 FIRST_ISOUSER 	<- as.integer(FIRST_CARTESIAN)
 LAST_REDUCEDXDIM_CART 	<- as.integer(DOUBLEISOTROPIC)
 FIRST_PROJECTION 	<- as.integer(GNOMONIC_PROJ)
 LAST_PROJECTION 	<- as.integer(ORTHOGRAPHIC_PROJ)
 LAST_CARTESIAN 	<- as.integer(LAST_PROJECTION)
 FIRST_SPHERICAL 	<- as.integer(SPHERICAL_ISOTROPIC)
 FIRST_ANYSPHERICAL 	<- as.integer(FIRST_SPHERICAL)
 LAST_TRUESPHERICAL 	<- as.integer(SPHERICAL_COORD)
 FIRST_EARTH 	<- as.integer(EARTH_ISOTROPIC)
 LAST_EARTH 	<- as.integer(EARTH_COORD)
 LAST_ANYSPHERICAL 	<- as.integer(LAST_EARTH)


 LAST_ISOUSER 	<- as.integer(UNREDUCED)
 LAST_ISO 	<- as.integer(ISO_MISMATCH)



 MON_UNSET 	<- as.integer(UNSET)
 MON_MISMATCH 	<- as.integer(MISMATCH)
 MON_SUB_DEP 	<- as.integer(SUBMODEL_DEP)
 MON_PREV_DEP 	<- as.integer(PREVMODEL_DEP)
 MON_PARAMETER 	<- as.integer(PARAM_DEP)
 NOT_MONOTONE 	<- as.integer(0)
 MONOTONE 	<- as.integer(1)
 GNEITING_MON 	<- as.integer(2)
 NORMAL_MIXTURE 	<- as.integer(3)
 COMPLETELY_MON 	<- as.integer(4)
 BERNSTEIN 	<- as.integer(5)


 MONOTONE_TOTAL 	<- as.integer((BERNSTEIN-MON_UNSET+1))



 MAXFIELDS 	<- as.integer(10)
 MODEL_USER 	<- as.integer((MAXFIELDS+0))
 MODEL_COV 	<- as.integer((MAXFIELDS+1))
 MODEL_COVMATRIX 	<- as.integer((MAXFIELDS+2))
 MODEL_VARIOGRAM 	<- as.integer((MAXFIELDS+3))
 MODEL_PSEUDO 	<- as.integer((MAXFIELDS+4))
 MODEL_FCTN 	<- as.integer((MAXFIELDS+5))
 MODEL_DISTR 	<- as.integer((MAXFIELDS+6))
 MODEL_CALC 	<- as.integer((MAXFIELDS+7))


 LAST_MODEL_USER 	<- as.integer((MAXFIELDS+9))
 FIRST_INTERNAL 	<- as.integer((LAST_MODEL_USER+1))
 MODEL_AUX 	<- as.integer((FIRST_INTERNAL+0))
 MODEL_INTERN 	<- as.integer((FIRST_INTERNAL+1))
 MODEL_SPLIT 	<- as.integer((FIRST_INTERNAL+2))
 MODEL_GUI 	<- as.integer((FIRST_INTERNAL+3))
 MODEL_MLE 	<- as.integer((FIRST_INTERNAL+4))
 MODEL_MLESPLIT 	<- as.integer((FIRST_INTERNAL+5))
 MODEL_LSQ 	<- as.integer((FIRST_INTERNAL+6))
 MODEL_BOUNDS 	<- as.integer((FIRST_INTERNAL+7))
 MODEL_KRIGE 	<- as.integer((FIRST_INTERNAL+8))
 MODEL_PREDICT 	<- as.integer((FIRST_INTERNAL+9))
 MODEL_ERR 	<- as.integer((FIRST_INTERNAL+10))
 MODEL_MAX 	<- as.integer(MODEL_ERR)


 TcfType 	<- as.integer(0)
 PosDefType 	<- as.integer(1)
 VariogramType 	<- as.integer(2)
 NegDefType 	<- as.integer(3)
 PointShapeType 	<- as.integer(4)
 ShapeType 	<- as.integer(5)
 TrendType 	<- as.integer(6)
 RandomOrShapeType 	<- as.integer(7)
 ManifoldType 	<- as.integer(8)
 ProcessType 	<- as.integer(9)
 GaussMethodType 	<- as.integer(10)
 NormedProcessType 	<- as.integer(11)
 BrMethodType 	<- as.integer(12)
 SmithType 	<- as.integer(13)
 SchlatherType 	<- as.integer(14)
 PoissonType 	<- as.integer(15)
 PoissonGaussType 	<- as.integer(16)
 RandomType 	<- as.integer(17)
 InterfaceType 	<- as.integer(18)
 MathDefType 	<- as.integer(19)
 OtherType 	<- as.integer(20)
 BadType 	<- as.integer(21)
 SameAsPrevType 	<- as.integer(22)
 LikelihoodType 	<- as.integer(23)
 EvaluationType 	<- as.integer(24)
 MixedInputType 	<- as.integer(25)
 CharInputType 	<- as.integer(26)
 NN1 	<- as.integer(27)
 NN2 	<- as.integer(28)
 NN3 	<- as.integer(29)
 NN4 	<- as.integer(30)


















 LASTTYPE 	<- as.integer(NN4)

 nOptimiser 	<- as.integer(8)
 nNLOPTR 	<- as.integer(15)
 nLikelihood 	<- as.integer(4)
 nNamesOfNames 	<- as.integer(12)
 nVAR2COV_METHODS 	<- as.integer(4)
 nProjections 	<- as.integer(2)

 DetTrendEffect 	<- as.integer(0)
 FixedTrendEffect 	<- as.integer(1)
 FixedEffect 	<- as.integer(2)
 DataEffect 	<- as.integer(3)
 RandomEffect 	<- as.integer(4)
 ErrorEffect 	<- as.integer(5)
 effect_error 	<- as.integer(6)








 VARPARAM 	<- as.integer(0)
 SIGNEDVARPARAM 	<- as.integer(1)
 SDPARAM 	<- as.integer(2)
 SIGNEDSDPARAM 	<- as.integer(3)
 SCALEPARAM 	<- as.integer(4)
 DIAGPARAM 	<- as.integer(5)
 ANISOPARAM 	<- as.integer(6)
 INTEGERPARAM 	<- as.integer(7)
 ANYPARAM 	<- as.integer(8)
 TRENDPARAM 	<- as.integer(9)
 NUGGETVAR 	<- as.integer(10)
 CRITICALPARAM 	<- as.integer(11)
 DONOTVERIFYPARAM 	<- as.integer(12)
 ONLYRETURN 	<- as.integer(13)
 FORBIDDENPARAM 	<- as.integer(14)
 UNKNOWNPARAM 	<- as.integer(15)
 VARONLYMLE 	<- as.integer(16)
 CRITONLYMLE 	<- as.integer(17)
 ONLYMLE 	<- as.integer(18)
 IGNOREPARAM 	<- as.integer(19)
 STANDARD 	<- as.integer(20)
 INCLUDENOTRETURN 	<- as.integer(21)
 INTERNALPARAMETERS 	<- as.integer(22)
 ALLPARAMETERS 	<- as.integer(23)
 NOPARAMETERS 	<- as.integer(24)


 LASTRETURNED 	<- as.integer(FORBIDDENPARAM)
 FIRSTONLYMLE 	<- as.integer(VARONLYMLE)
 LASTONLYMLE 	<- as.integer(ONLYMLE)
 LASTUSERSORTOF 	<- as.integer(FORBIDDENPARAM)
 LASTSORTOF 	<- as.integer(NOPARAMETERS)





 CircEmbed 	<- as.integer(0)
 CircEmbedCutoff 	<- as.integer(1)
 CircEmbedIntrinsic 	<- as.integer(2)
 TBM 	<- as.integer(3)
 SpectralTBM 	<- as.integer(4)
 Direct 	<- as.integer(5)
 Sequential 	<- as.integer(6)
 Trendproc 	<- as.integer(7)
 Average 	<- as.integer(8)
 Nugget 	<- as.integer(9)
 RandomCoin 	<- as.integer(10)
 Hyperplane 	<- as.integer(11)
 Specific 	<- as.integer(12)
 Nothing 	<- as.integer(13)
 Forbidden 	<- as.integer(14)













 original 	<- as.integer(0)
 mle_conform 	<- as.integer(1)
 all_origins 	<- as.integer(2)



 INTERNAL_PARAM 	<- "internal"



 GETMODEL_AS_SAVED 	<- as.integer(0)
 GETMODEL_DEL_NATSC 	<- as.integer(1)
 GETMODEL_SOLVE_NATSC 	<- as.integer(2)
 GETMODEL_DEL_MLE 	<- as.integer(3)
 GETMODEL_SOLVE_MLE 	<- as.integer(4)

 MIXED_X_NAME 	<- "X"
 MIXED_BETA_NAME 	<- "beta"

 COVARIATE_C_NAME 	<- "data"
 COVARIATE_X_NAME 	<- "x"
 COVARIATE_RAW_NAME 	<- "raw"
 COVARIATE_ADDNA_NAME 	<- "addNA"

 CONST_A_NAME 	<- "x"

 COVARIATE_DATA_NAME 	<- "data"

 MINMAX_PMIN 	<- as.integer(1)
 MINMAX_PMAX 	<- as.integer(2)
 MINMAX_TYPE 	<- as.integer(3)
 MINMAX_NAN 	<- as.integer(4)
 MINMAX_MIN 	<- as.integer(5)
 MINMAX_MAX 	<- as.integer(6)
 MINMAX_OMIN 	<- as.integer(7)
 MINMAX_OMAX 	<- as.integer(8)
 MINMAX_ROWS 	<- as.integer(9)
 MINMAX_COLS 	<- as.integer(10)
 MINMAX_BAYES 	<- as.integer(11)
 MINMAX_ENTRIES 	<- as.integer(MINMAX_BAYES)



 XLIST_X 	<- as.integer(0)
 XLIST_Y 	<- as.integer(1)
 XLIST_T 	<- as.integer(2)
 XLIST_GRID 	<- as.integer(3)
 XLIST_SPATIALDIM 	<- as.integer(4)
 XLIST_TIME 	<- as.integer(5)
 XLIST_DIST 	<- as.integer(6)
 XLIST_RESTOT 	<- as.integer(7)
 XLIST_L 	<- as.integer(8)
 XLIST_UNITS 	<- as.integer(9)
 XLIST_NEWUNITS 	<- as.integer(10)
 XLIST_ENTRIES 	<- as.integer((XLIST_NEWUNITS+1))


 PROJ_SPACE 	<- as.integer(-1)
 PROJ_TIME 	<- as.integer(-2)
 PROJECTIONS 	<- as.integer(2)

 VAR2COV_EXTREMAL 	<- as.integer(-1)
 VAR2COV_ORIGIN 	<- as.integer(-2)
 VAR2COV_CENTER 	<- as.integer(-3)
 VAR2COV_METHODS 	<- as.integer(3)


 INTERN_SHOW 	<- as.integer(2)


 METHOD_VARIOGRAM 	<- as.integer(0)
 METHOD_PSEUDO 	<- as.integer(1)
 METHOD_COVARIANCE 	<- as.integer(2)
 METHOD_PSEUDOMADOGRAM 	<- as.integer(3)
 METHOD_MADOGRAM 	<- as.integer(4)

 EV_FFT_EV 	<- as.integer(0)
 EV_FFT_N 	<- as.integer(1)
 EV_EV 	<- as.integer(0)
 EV_SDSQ 	<- as.integer(1)
 EV_N 	<- as.integer(2)

 BRP_UNIF 	<- as.integer(0)
 BRP_COV 	<- as.integer(1)
 BRP_KRIG 	<- as.integer(2)



 VAR2COV_X 	<- as.integer(0)
 VAR2COV_C 	<- as.integer(1)
 VAR2COV_SPECIAL 	<- as.integer(1)


 POISSON_SCATTER_OPTIM 	<- as.integer(0)
 POISSON_SCATTER_STANDARD 	<- as.integer(1)
 POISSON_SCATTER_ANY 	<- as.integer(2)
 nPOISSON_SCATTER 	<- as.integer((POISSON_SCATTER_ANY+1))






 ## from  src/AutoRandomFields.cc

RC_DOMAIN_NAMES <- DOMAIN_NAMES <-
c( "single variable","kernel","framework dependent","submodel dependent","parameter dependent","<keep copy>","mismatch" )


RC_OPTIMISER_NAMES <- OPTIMISER_NAMES <-
c( "optim","optimx","soma","nloptr","GenSA","minqa","pso","DEoptim" )


RC_NLOPTR_NAMES <- NLOPTR_NAMES <-
c( "NLOPT_GN_DIRECT","NLOPT_GN_DIRECT_L","NLOPT_GN_DIRECT_L_RAND","NLOPT_GN_DIRECT_NOSCAL","NLOPT_GN_DIRECT_L_NOSCAL","NLOPT_GN_DIRECT_L_RAND_NOSCAL","NLOPT_GN_ORIG_DIRECT","NLOPT_GN_ORIG_DIRECT_L","NLOPT_LN_PRAXIS","NLOPT_GN_CRS2_LM","NLOPT_LN_COBYLA","NLOPT_LN_NELDERMEAD","NLOPT_LN_SBPLX","NLOPT_LN_BOBYQA","NLOPT_GN_ISRES" )


RC_LIKELIHOOD_NAMES <- LIKELIHOOD_NAMES <-
c( "auto","full","composite","tesselation" )


RC_ISO_NAMES <- ISO_NAMES <-
c( "isotropic","space-isotropic","vector-isotropic","symmetric","cartesian system","gnomonic projection","orthographic projection","spherical isotropic","spherical symmetric","spherical system","earth isotropic","earth symmetric","earth system","cylinder system","non-dimension-reducing","framework dependent","submodel dependent","parameter dependent","<internal keep copy>","<mismatch>" )


RC_TYPE_NAMES <- TYPE_NAMES <-
c( "tail correlation","positive definite","variogram","negative definite","point-shape function","shape function","trend","distribution or shape","of manifold type","process","method for Gauss process","normed process (non-negative values with maximum value being 0 or 1)","method for Brown-Resnick process","Smith","Schlather","Poisson","PoissonGauss","distribution family","interface","mathematical operator","other type","badly defined","<same as previous>","likelihood","evaluation","mixed input","string input","<special I>","<special II>","<special III>","<special IV>" )


RC_MONOTONE_NAMES <- MONOTONE_NAMES <-
c( "not set","mismatch in monotonicity","submodel dependent monotonicity","previous model dependent monotonicity","parameter dependent monotonicity","not monotone","monotone","Gneiting-Schaback class","normal mixture","completely monotone","Bernstein" )


 SORT_ORIGIN_NAMES <-
c( "original","MLE conform","all" )


 MODE_NAMES <-
c( "careless","sloppy","easygoing","normal","precise","pedantic","neurotic" )


 OUTPUTMODE_NAMES <-
c( "sp","RandomFields","geoR" )


 ANGLE_NAMES <-
c( "radians","degree" )


 REPORTCOORD_NAMES <-
c( "always","warn","important","never" )


 UNITS_NAMES <-
c( "","km","miles","<time>","<user defined>" )


 COORD_SYS_NAMES <-
c( "auto","keep","cartesian","earth","sphere","gnomonic","orthographic","coordinate system mix" )


 COORD_NAMES_CART <-
c( "x","y","z","T" )


 COORD_NAMES_GENERAL <-
c( "coords.x","coords.T" )


 COORD_NAMES_EARTH <-
c( "longitude","latitude","height","time" )


 CARTESIAN_SYS_NAMES <-
c( "cartesian","gnomonic","orthographic" )


 TYPEOF_PARAM_NAMES <-
c( "variance","covariance","sd","signed sd","scale","diagonal","aniso","integer","unspecified","trend","nugget","critical to estimate","never verified","internally ignored","forbidden to be estimated","unkown","variance (used internally)","critical to estimate (used internally)","only used by MLE","neither used by MLE nor returned","standard","include never returned","internal","all","none" )


 EQ_NAMES <-
c( "==","!=","<=","<",">=",">" )


 NAMES_OF_NAMES <-
c( "EQ_NAMES","ISO_NAMES","DOMAIN_NAMES","TYPE_NAMES","MONOTONE_NAMES","MODE_NAMES","OUTPUTMODE_NAMES","REPORTCOORD_NAMES","UNITS_NAMES","COORD_SYS_NAMES","CARTESIAN_SYS_NAMES","TYPEOF_PARAM_NAMES" )


 PROJECTION_NAMES <-
c( "space","time" )


 RMCOV_X <-
c( "origin","center","extremals","all" )


 POISSON_SCATTER_NAMES <-
c( "optimized","standard","any" )



list2RMmodel_Names <- c('R.acos', 'R.acosh', 'R.asin', 'R.asinh', 'R.atan', 'R.atan2', 'R.atanh', 'R.c', 'R.cbrt', 'R.ceil', 'R.const', 'R.cos', 'R.cosh', 'R.div', 'R.erf', 'R.erfc', 'R.exp', 'R.exp2', 'R.expm1', 'R.fabs', 'RFboxcox', 'RFcalc', 'RFcov', 'RFcovmatrix', 'RFcrossvalidate', 'RFddistr', 'R.fdim', 'RFdistr', 'RFearth2cartesian', 'RFearth2dist', 'RFempiricalcovariance', 'RFempiricalmadogram', 'RFempiricalvariogram', 'RFfctn', 'RFfit', 'RFformula', 'RFfractaldim', 'RFgetMethodNames', 'RFgetModel', 'RFgetModelInfo', 'RFgetModelInfo_model', 'RFgetModelInfo_register', 'RFgetModelNames', 'RFgridDataFrame', 'RFgui', 'RFhessian', 'RFhurst', 'RFinterpolate', 'RFlikelihood', 'RFlinearpart', 'R.floor', 'RFmadogram', 'R.fmax', 'R.fmin', 'R.fmod', 'RFoldstyle', 'RFoptions', 'RFOPTIONS', 'RFpar', 'RFparameters', 'RFpdistr', 'RFplotEmpVariogram', 'RFplotModel', 'RFplotSimulation', 'RFplotSimulation1D', 'RFpointsDataFrame', 'RFpseudomadogram', 'RFpseudovariogram', 'RFqdistr', 'RFratiotest', 'RFrdistr', 'RFsimulate', 'RFspatialGridDataFrame', 'RFspatialPointsDataFrame', 'RFspDataFrame2conventional', 'RFspDataFrame2dataArray', 'RFvariogram', 'R.gamma', 'R.hypot', 'R.is', 'R.lat', 'R.lgamma', 'R.log', 'R.log1p', 'R.log2', 'R.lon', 'RMangle', 'RMaskey', 'RMave', 'RMball', 'RMbcw', 'RMbernoulli', 'RMbessel', 'RMbicauchy', 'RMbigneiting', 'RMbistable', 'RMbiwm', 'RMblend', 'RMbr2bg', 'RMbr2eg', 'RMbrownresnick', 'RMbubble', 'RMcardinalsine', 'RMcauchy', 'RMcauchytbm', 'RMchoquet', 'RMcircular', 'RMconstant', 'RMcov', 'RMcovariate', 'RM_COVARIATE', 'RMCOV_X', 'RMcoxisham', 'RMcubic', 'RMcurlfree', 'RMcutoff', 'RMdagum', 'RMdampedcos', 'RMdeclare', 'RM_DECLARE', 'RM_DEFAULT', 'RMdelay', 'RMderiv', 'RMdewijsian', 'RM_DISTR', 'RMdivfree', 'RMeaxxa', 'RMepscauchy', 'RMetaxxa', 'RMexp', 'RMexponential', 'RMfbm', 'RMfixcov', 'RMflatpower', 'RMfractdiff', 'RMfractgauss', 'RMgauss', 'RMgencauchy', 'RMgenfbm', 'RMgengneiting', 'RMgennsst', 'RMgneiting', 'RMgneitingdiff', 'RMhandcock', 'RMhyperbolic', 'RMiaco', 'RMid', 'RMidmodel', 'RMintexp', 'RMintrinsic', 'R.minus', 'RMkolmogorov', 'RMlgd', 'RMlsfbm', 'RMm2r', 'RMm3b', 'RMma', 'RMmastein', 'RMmatern', 'RMmatrix', 'RMmppplus', 'RMmps', 'RMmqam', 'RMmult', 'RM_MULT', 'RMmultiquad', 'RMnatsc', 'RMnsst', 'RMnugget', 'RM_NUGGET', 'RMparswm', 'RMparswmX', 'RMpenta', 'RMplus', 'RM_PLUS', 'RMpolygon', 'RMpolynome', 'RMpower', 'RMpoweredexp', 'RMprod', 'RMqam', 'RMqexp', 'RMrational', 'RMrotat', 'RMrotation', 'RMS', 'RMscale', 'RMschlather', 'RMschur', 'RMsign', 'RMsinepower', 'RMspheric', 'RMstable', 'RMstein', 'RMstp', 'RMsum', 'RMtbm', 'RMtent', 'RMtrafo', 'RMtrend', 'RM_TREND', 'RMtruncsupport', 'R.mult', 'RMuser', 'RM_USER', 'RMvector', 'RMwave', 'RMwendland', 'RMwhittle', 'R.p', 'RPaverage', 'RPbernoulli', 'RPbrmixed', 'RPbrorig', 'RPbrownresnick', 'RPbrshifted', 'RPchi2', 'RPcirculant', 'RPcoins', 'RPcutoff', 'RPdirect', 'RPgauss', 'RPhyperplane', 'RPintrinsic', 'RPloggaussnormed', 'R.plus', 'RPnugget', 'RPopitz', 'R.pow', 'RPpoisson', 'RPschlather', 'RPsequential', 'RPsmith', 'RPspecific', 'RPspectral', 'RPt', 'RPtbm', 'RPtrend', 'RRdeterm', 'RRdistr', 'R.remainder', 'RRgauss', 'RRloc', 'RRmcmc', 'R.round', 'RRrectangular', 'RRspheric', 'RRunif', 'R.sin', 'R.sinh', 'R.sqrt', 'R.tan', 'R.tanh', 'R.trunc')

list2RMmodel_oldNames <- c('#', '>', '>', '+', '*', '$', '$power', 'ave', 'shape.ave', 'bcw', 'locstatfbm', 'bessel', 'bigneiting', 'bernoulli', 'biWM', 'bistable', 'blend', 'brownresnick', 'br2bg', 'br2eg', 'bubble', 'cauchy', 'circular', 'CDeWijsian', 'constant', 'fixcov', 'coxisham', 'cubic', 'curlfree', 'cutoff', 'dagum', 'dampedcosine', 'deriv', 'DeWijsian', 'divfree', 'epsC', 'exponential', 'Exp', 'extremalgauss', 'FD', 'flatpower', 'fractalB', 'fractgauss', 'gauss', 'genB', 'gencauchy', 'bicauchy', 'gengneiting', 'gengneiting', 'gengneit_intern', 'gneiting', 'gennsst', 'gennsst_intern', 'hyperbolic', 'iacocesare', 'identity', 'kolmogorov', 'lgd1', 'mastein', 'ma1', 'ma2', 'M', 'matern', 'mqam', 'multiquad', 'natsc', 'nsst', 'parsWM', 'penta', 'power', 'Pow', 'prod', 'qam', 'qexponential', 'scale', 'schur', 'shift', 'sinepower', 'spherical', 'stable', 'Stein', 'steinst1', 'stp', 'shape.stp', 'tbm', 'sum', 'U', 'cov', 'vector', 'wave', 'whittle', 'nugget', 'missing', 'null', 'trend', 'select', 'angle', 'ball', 'covariate', 'declare', 'EAxxA', 'EtAxxA', 'id', 'trafo', 'mult_inverse', 'polygon', 'rational', 'rotat', 'Rotat', 'scatter', 'sign', 'setparam', 'm2r', 'm3b', 'r3binner', 'mps', 'truncsupport', 'arcsqrt', 'determ', 'distr', 'normal', 'setDistr', 'loc', 'mcmc', 'rectangular', 'spheric', 'unif', 'MCMC_PGS', 'zhou', 'ballani', 'standardShape', '++', 'statShape', 'Cov', 'CovMatrix', 'Dummy', 'get', 'Fctn', 'Distr', 'loglikelihood', 'linearpart', 'predict', 'Pseudovariogram', 'Variogram', 'Simulate', '$proc', 'plusproc', 'prodproc', 'trafoproc', 'mppplusproc', 'multproc', 'matrixproc', 'covproc', 'trend', 'average', 'coins', 'averageIntern', 'circulant', 'cutoff', 'cutoffIntern', 'intrinsic', 'intrinsIntern', 'direct', 'hyperplane', 'hyperIntern', 'nugget', 'nuggetIntern', 'sequential', 'spectral', 'spectralIntern', 'specific', 'tbm', 'tbmIntern', 'loggaussnormed', 'brorig', 'brorigIntern', 'brmixed', 'brmixedIntern', 'brshifted', 'brshiftIntern', 'brownresnick', 'binaryprocess', 'gauss.process', 'poisson', 'extremalgauss', 'extremalt', 'smith', 'chi2', 't', 'minus', 'plus', 'div', 'mult', 'const', 'p', 'c', 'is', '.asin', '.atan', '.atan2', '.cos', '.sin', '.tan', '.asinh', '.atanh', '.cosh', '.sinh', '.tanh', '.log', '.expm1', '.log1p', '.exp2', '.log2', '.hypot', '.cbrt', '.ceil', '.floor', '.fmod', '.round', '.trunc', '.erfc', '.lgamma', '.remainder', '.fdim', '.fmax', '.fmin', '.gamma', '.exp', '.erf', '.fabs', '.acos', '.acosh', '.pow', '.sqrt')

rfgui_Names1 <- c('RMaskey', 'RMbcw', 'RMbessel', 'RMcardinalsine', 'RMcauchy', 'RMcauchytbm', 'RMcircular', 'RMconstant', 'RMcubic', 'RMdagum', 'RMdampedcos', 'RMepscauchy', 'RMexp', 'RMfractdiff', 'RMfractgauss', 'RMgauss', 'RMgencauchy', 'RMgengneiting', 'RMgneiting', 'RMgneitingdiff', 'RMhyperbolic', 'RMlgd', 'RMlsfbm', 'RMparswm', 'RMparswmX', 'RMpenta', 'RMpoweredexp', 'RMqexp', 'RMspheric', 'RMstable', 'RMwave', 'RMwendland')

rfgui_Names2 <- c('RMaskey', 'RMbcw', 'RMbessel', 'RMcardinalsine', 'RMcauchy', 'RMcircular', 'RMconstant', 'RMcubic', 'RMdagum', 'RMdampedcos', 'RMepscauchy', 'RMexp', 'RMgauss', 'RMgencauchy', 'RMgengneiting', 'RMgneiting', 'RMhyperbolic', 'RMlgd', 'RMlsfbm', 'RMparswm', 'RMparswmX', 'RMpenta', 'RMpoweredexp', 'RMqexp', 'RMspheric', 'RMstable', 'RMwave', 'RMwendland')
