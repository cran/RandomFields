#include "AutoRandomFields.h"

const char
*dummyconstantnotused = "just to get the formatting in emacs right",
  *DOMAIN_NAMES[LAST_DOMAIN + 1] = { //RC
  "single variable", "kernel", "framework dependent", "submodel dependent",
  "parameter dependent", "<keep copy>", "mismatch"},


  *OPTIMISER_NAMES[nOptimiser] = { // RC
    "optim", "optimx", "soma", "nloptr", "GenSA", "minqa", "pso", "DEoptim"},
  
  *NLOPTR_NAMES[nNLOPTR] = { // RC
    "NLOPT_GN_DIRECT", "NLOPT_GN_DIRECT_L", 
    "NLOPT_GN_DIRECT_L_RAND", "NLOPT_GN_DIRECT_NOSCAL", 
    "NLOPT_GN_DIRECT_L_NOSCAL", "NLOPT_GN_DIRECT_L_RAND_NOSCAL", 
    "NLOPT_GN_ORIG_DIRECT", "NLOPT_GN_ORIG_DIRECT_L",
    "NLOPT_LN_PRAXIS", "NLOPT_GN_CRS2_LM",
    "NLOPT_LN_COBYLA", "NLOPT_LN_NELDERMEAD", 
    "NLOPT_LN_SBPLX", "NLOPT_LN_BOBYQA", 
    "NLOPT_GN_ISRES"},

  *LIKELIHOOD_NAMES[nLikelihood] = { // RC 
    "auto", "full", "composite", "tesselation"},

  *ISO_NAMES[LAST_ISO + 1] =  { // RC
    "isotropic", "space-isotropic", "vector-isotropic", "symmetric",
    "cartesian system", "gnomonic projection", // 5
    "orthographic projection",
    "spherical isotropic", "spherical symmetric", "spherical system", 
    "earth isotropic", // 10
    "earth symmetric",  "earth system",
    // "symmetric log-cartesian", "log-cartesian",
    "cylinder system",
    "non-dimension-reducing", "framework dependent", // 15
    "submodel dependent",
    "parameter dependent", "<internal keep copy>", "<mismatch>"},
  
  *TYPE_NAMES[LASTTYPE + 1] = { // RC
    "tail correlation", "positive definite", "variogram", "negative definite",
    "point-shape function", "shape function", "trend", "distribution or shape",
    "of manifold type", "process", "method for Gauss process",  // 10
    "normed process (non-negative values with maximum value being 0 or 1)",
    "method for Brown-Resnick process",
    "Smith", "Schlather", "Poisson", "PoissonGauss", // 16
    "distribution family",
    "interface", 
    "mathematical operator", 
    "other type", // 20
    //
    "badly defined", "<same as previous>", "likelihood", "evaluation",// 24
    
    "mixed input", "string input", "<special I>", "<special II>",
    "<special III>", "<special IV>"//  31
  }, 

  
  *MONOTONE_NAMES[MONOTONE_TOTAL] = { // RC
    "not set", "mismatch in monotonicity", "submodel dependent monotonicity",
    "previous model dependent monotonicity",
    "parameter dependent monotonicity",
    "not monotone", "monotone", "Gneiting-Schaback class", 
    "normal mixture", 
    "completely monotone",  
    "Bernstein"},

  *SORT_ORIGIN_NAMES[1 + (int) all_origins] = {
    "original", "MLE conform", "all"
  },
 
  *MODE_NAMES[nr_modes] = {
    "careless", "sloppy", "easygoing", "normal", 
    "precise", "pedantic", "neurotic"},

  *OUTPUTMODE_NAMES[nr_output_modes] = {
    "sp", "RandomFields", "geoR"},

  *ANGLE_NAMES[last_angle_mode + 1] = {
    "radians", "degree"},

  *REPORTCOORD_NAMES[nr_reportcoord_modes] = {
    "always", "warn", "important", "never"},

  *UNITS_NAMES[nr_units] = {
    "", "km", "miles", "<time>", "<user defined>"},

  *COORD_SYS_NAMES[nr_coord_sys] = {
    "auto", "keep", "cartesian", "earth",
    "sphere", "gnomonic", "orthographic", "coordinate system mix"},

  *COORD_NAMES_CART[4] = {
    "x", "y", "z", "T"},
  
  *COORD_NAMES_GENERAL[2] = {
    "coords.x", "coords.T"},
  
  *COORD_NAMES_EARTH[4] = {
    "longitude", "latitude", "height", "time"},
  
  *CARTESIAN_SYS_NAMES[3] = {
    "cartesian", "gnomonic", "orthographic"},

  
  *TYPEOF_PARAM_NAMES[LASTSORTOF + 1] = {
    "variance", "covariance", "sd", "signed sd", "scale", 
    "diagonal", "aniso", "integer", "unspecified",  "trend",
    "nugget", "critical to estimate", "never verified",    
    "internally ignored", "forbidden to be estimated", "unkown",
    "variance (used internally)",
    "critical to estimate (used internally)", "only used by MLE", 
    "neither used by MLE nor returned",
    "standard", "include never returned", "internal", "all",
    "none"},

  *EQ_NAMES[6] = {"==", "!=", "<=", "<", ">=", ">"},
  
  
  *NAMES_OF_NAMES[nNamesOfNames] = {// see also LIST_OF_NAMES in init.general.cc
    "EQ_NAMES", "ISO_NAMES", // never change ordering   
    "DOMAIN_NAMES", // from here on changings are ok
    "TYPE_NAMES", "MONOTONE_NAMES",
    "MODE_NAMES", "OUTPUTMODE_NAMES", "REPORTCOORD_NAMES",
    "UNITS_NAMES", "COORD_SYS_NAMES", "CARTESIAN_SYS_NAMES",
    "TYPEOF_PARAM_NAMES"},
  
  *PROJECTION_NAMES[nProjections] = {
    "space", "time"},

  *RMCOV_X[nVAR2COV_METHODS] = { // c ==  COVARIATE_C_NAME 
  "origin", "center", "extremals", "all"},  

  *POISSON_SCATTER_NAMES[nPOISSON_SCATTER] = {// must be in the ordering of
  "optimized", "standard", "any"}; // addPGS/local in extremes.cc and RMS.cc
 

  
