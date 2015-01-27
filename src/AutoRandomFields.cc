#include "AutoRandomFields.h"

const char
*dummyconstantnotused = "just to get the formatting in emacs right",
  *DOMAIN_NAMES[LAST_DOMAIN + 1] = { //RC
  "single variable", "kernel", "framework dependent", "mismatch"},
  
  *ISONAMES[LAST_ISO + 1] =  { // RC
    "isotropic", "space-isotropic", "zero-space-isotropic", 
    "vector-isotropic", "symmetric", "cartesian system",
    "gnomonic projection", "orthographic projection",
    "spherical isotropic", "spherical system", 
    "earth isotropic", "earth system", 
    "cylinder system",
    "non-dimension-reducing", "parameter dependent", "<mismatch>"},
  
  *TYPENAMES[OtherType + 1] = { // RC
    "tail correlation", "positive definite", "variogram", 
    "negative definite", "process", 
    "method for Gauss process", "method for Brown-Resnick process",
    "point-shape function", "distribution family", "shape function",
    "trend", "interface",  "undefined", 
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

  *TYPEOF_PARAM_NAMES[14] = {
    "var", "signed var", "sd", "signed sd", "scale",
    "diag", "aniso", "integer", "unspecfd", "trend",  
    "nugget", "mixed var", "regress", "any"},



  *NAMES_OF_NAMES[11] = {"DOMAIN_NAMES", "ISONAMES",  // never change ordering
		     "TYPENAMES", "MONOTONE_NAMES",
		     "MODENAMES", "OUTPUTMODENAMES", "REPORTCOORDNAMES",
		     "UNITS_NAMES", "COORD_SYS_NAMES", "CARTESIAN_SYSTEMS",
		     "TYPEOF_PARAM_NAMES"};

  
