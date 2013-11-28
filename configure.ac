# autoconf configure.ac > configure && chmod +x configure

# Urversion vo 2.0.66
# ab 2.0.67 geaendert.

AC_INIT([RandomFields], 1.3) dnl package name, version

        
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
   echo "could not determine R_HOME"
   exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`


AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT