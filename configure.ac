AC_INIT([RandomFields], 1.3) dnl package name, version

         

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
   echo "could not determine R_HOME"
   exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`

GMRF_EXT_lib_path="/usr/lib"
GMRF_lib_path="/usr/local/lib"
# 
HAVE_GMRF="" ## HAVE_GMRF undefined 
AC_ARG_ENABLE(GMRF, [  --enable-GMRF           include the GMRF method by Havard Rue],
    [AC_ARG_WITH([GMRF-lib],
	AC_HELP_STRING([--with-GMRF-lib=LIB_PATH],
		[the location of GMRF libraries]),
		[GMRF_lib_path=$withval])

    AC_CHECK_HEADERS(GMRFLib/GMRFLib.h)
    if test "${ac_cv_header_GMRFLib_h}" = no; then
         AC_MSG_ERROR(["GMRF headers GMRFLib.h not found"])
    fi

    AC_ARG_WITH([GMRF-EXT-lib],
	AC_HELP_STRING([--with-GMRF-EXT-lib=LIB_PATH],
		[the location of external GMRF libraries]),
		[GMRF_EXT_lib_path=$withval])

    AC_SUBST([LIBS],["-L${GMRF_lib_path}"])
    AC_CHECK_LIB(GMRFLib, GMRFLib_set_error_handler_off, [],
                 [AC_MSG_ERROR("GMRFLib library not found on ${GMRF_lib_path} ")])
    AC_CHECK_LIB(GMRFLib-geo, main, [],
                 [AC_MSG_ERROR("GMRFLib-geo library not found on ${GMRF_lib_path} ")])

    AC_SUBST([LIBS],[" -L${GMRF_EXT_lib_path} ${LIBS}"])
#   AC_CHECK_LIB(lapack, main, ,
#     [AC_MSG_ERROR("lapack library not found on ${GMRF_EXT_lib_path} -- include from GMRFLib/EXTLIBS")])
    AC_SUBST([LIBS], [" -llapack ${LIBS}"])
    AC_CHECK_LIB(metis, main, ,
      [AC_MSG_ERROR("metis library not found on ${GMRF_EXT_lib_path} --include from GMRFLib/EXTLIBS")])
    AC_CHECK_LIB(taucs, main, , 
      [AC_MSG_ERROR("taucs library not found on ${GMRF_EXT_lib_path} -- include from GMRFLib/EXTLIBS")])

#    AC_CHECK_LIB(blas, main,  ,
#       [AC_MSG_ERROR("blas library not found on ${GMRF_EXT_lib_path} -- inlcude from GMRFLib/EXTLIBS")])
    AC_SUBST([LIBS], [" -lblas ${LIBS}"])
#    AC_CHECK_LIB(z, main, ,
#        [AC_MSG_ERROR("z library not found on ${GMRF_EXT_lib_path} -- include from GMRFLib/EXTLIBS")])
    AC_SUBST([LIBS], [" -lz ${LIBS}"])
#    AC_CHECK_LIB(gslcblas, main, ,
#        [AC_MSG_ERROR("gslcblas library not found on ${GMRF_EXT_lib_path} -- include from GMRFLib/EXTLIBS")])
    AC_SUBST([LIBS], [" -lgslcblas ${LIBS}"])
# AC_CHECK_LIB(gsl, gsl_version_, ,
#     [AC_MSG_ERROR("gsl library not found on ${GMRF_EXT_lib_path} -- include from GMRFLib/EXTLIBS")])
    AC_SUBST([LIBS], [" -lgsl ${LIBS}"])

    HAVE_GMRF="#define INCLUDEMARKOV 1" 
    
    AC_SUBST([LIBS],[" -L${GMRF_lib_path} -lGMRFLib -lGMRFLib-geo -L/usr/lib  -L${GMRF_EXT_lib_path} -lgsl -lgslcblas -lz -lblas  -L${GMRF_EXT_lib_path} -ltaucs -lmetis -llapack"])
])


AC_SUBST(HAVE_GMRF)
AC_CONFIG_FILES([src/includeMarkov.h src/Makevars])
AC_OUTPUT

cat src/includeMarkov.h
cat src/Makevars

dnl     AC_CHECK_LIB(foo, fun, [HAVE_FOO=TRUE], [HAVE_FOO=FALSE])
dnl     AC_SUBST(HAVE_FOO)
dnl     ......
dnl     AC_CONFIG_FILES([foo.R])
dnl     AC_OUTPUT
dnl
dnl in configure.ac (assuming Autoconf 2.50 or better).
dnl
dnl The definition of the respective R function in foo.R.in could be
dnl
dnl     foo <- function(x) {
dnl         if(!@HAVE_FOO@)
dnl           stop("Sorry, library 'foo' is not available"))
dnl         ...

