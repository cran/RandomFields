#include "AutoRandomFields.h"

const char
*dummyconstantnotused = "just to get the formatting in emacs right",
  *DOMAIN_NAMES[LAST_DOMAIN + 1] = { //RC
  "single variable", "kernel", "framework dependent", "mismatch"},

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

  *ISONAMES[LAST_ISO + 1] =  { // RC
    "isotropic", "space-isotropic", "zero-space-isotropic", 
    "vector-isotropic", "symmetric", "cartesian system",
    "gnomonic projection", "orthographic projection",
    "spherical isotropic", "spherical symmetric", "spherical system", 
    "earth isotropic", "earth symmetric",  "earth system", 
    "cylinder system",
    "non-dimension-reducing", "parameter dependent", "<mismatch>"},
  
  *TYPENAMES[OtherType + 1] = { // RC
    "tail correlation", "positive definite", "variogram", 
    "negative definite", "process", 
    "method for Gauss process", "method for Brown-Resnick process",
    "point-shape function", "distribution family", "shape function",
    "trend", "interface",  "distribution or shape", "undefined", 
    "<math definition>", "other type"},
  
  *MONOTONE_NAMES[BERNSTEIN + 1 - (int) MISMATCH] = { // RC
    "mismatch in monotonicity", "submodel dependent monotonicity",
    "previous model dependent monotonicity",
    "parameter dependent monotonicity",
    "not monotone", "monotone", "Gneiting-Schaback class", 
    "normal mixture", 
    "completely monotone",  
    "Bernstein"},
  
  *MODENAMES[nr_modes] = {
    "careless", "sloppy", "easygoing", "normal", 
    "precise", "pedantic", "neurotic"},

  *OUTPUTMODENAMES[nr_output_modes] = {
    "sp", "RandomFields", "geoR"},

  *REPORTCOORDNAMES[nr_reportcoord_modes] = {
    "always", "warn", "important", "never"},

  *UNITS_NAMES[nr_units] = {
    "", "km", "miles", "<time>", "<user defined>"},

  *COORD_SYS_NAMES[nr_coord_sys] = {
    "auto", "keep", "cartesian", "earth",
    "sphere", "gnomonic", "orthographic", "coordinate system mix"},
  
  *CARTESIAN_SYSTEMS[3] = {
    "cartesian", "gnomonic", "orthographic"},

  
  *TYPEOF_PARAM_NAMES[FORBIDDENPARAM + 1] = {
    "variance", "covariance", "sd", "signed sd", "scale", 
    "diagonal", "aniso", "integer", "unspecified",  "trend",
    "nugget", "mixed variance", "critical to estimate", 
    "internally ignored", "never varified", 
    "never returned", "forbidden"},

  *EQNAMES[6] = {"==", "!=", "<=", "<", ">=", ">"},
  
  
  *NAMES_OF_NAMES[12] = {"EQNAMES", "ISONAMES",  // never change ordering   
		       "DOMAIN_NAMES", // from here on changings are ok
		       "TYPENAMES", "MONOTONE_NAMES",
		       "MODENAMES", "OUTPUTMODENAMES", "REPORTCOORDNAMES",
		       "UNITS_NAMES", "COORD_SYS_NAMES", "CARTESIAN_SYSTEMS",
		       "TYPEOF_PARAM_NAMES"},

  *PROJECTION_NAMES[PROJECTIONS] = {
    "space", "time"};

 
  

  
