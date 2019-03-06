AC_INIT([RandomFields], 3.3.4)

if test "${CXX}" = gcc || "${CXX}" = clang; then
dnl  RANDOMFIELD_CXXFLAGS="-mavx"
  RANDOMFIELDS_CXXFLAGS=""
dnl  AC_MSG_RESULT([${OBJC} found, so that the package is compiled with 'avx'])
else  
  RANDOMFIELDS_CXXFLAGS=""
  AC_MSG_RESULT([Neither 'gcc' nor 'clang' found, so that AVX will not be used])
fi


AC_SUBST(RANDOMFIELDS_CXXFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
