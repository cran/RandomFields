AC_INIT([RandomFields], 3.3.4)


if "${R_HOME}/bin/R" CMD config CXX11 > /dev/null ; then
  RANDOMFIELDS_CXX=`"${R_HOME}/bin/R" CMD config CXX11`
elif "${R_HOME}/bin/R" CMD config CXX1X > /dev/null ; then
  RANDOMFIELDS_CXX=`"${R_HOME}/bin/R" CMD config CXX1X`
else
  RANDOMFIELDS_CXX=`"${R_HOME}/bin/R" CMD config CXX`
fi

if test "x$TARGET" = "x" ; then
  case `basename $RANDOMFIELDS_CXX` in
  g++ | clang++ )
    if uname -m | grep -iq sparc ; then
      echo Since this is a sparc processor we set TARGET=default
      TARGET=default
    else
      TARGET=native
    fi
    ;;
  *)
    echo unknown compiler: $RANDOMFIELDS_CXX therefore we set TARGET=default
    TARGET=default
    ### Now we also should be sure that there is something like -std=c++0x??
    ;;
  esac
fi

if test "x$RANDOMFIELDS_TARGET" != "x" ; then
  TARGET=$RANDOMFIELDS_TARGET
fi

if test "x$1" != "x" ; then
case $1 in
native | generic | default | debug | empty )
  TARGET=$1
  shift 1
  ;;
esac
fi


echo Using target: $TARGET 
echo Using further args: $@


if test $TARGET = native ; then
  CXX_FLAGS="-O2 -march=native $CXX_FLAGS"
elif test $TARGET = generic ; then
  CXX_FLAGS="-mtune=generic -msse2 -O3 $CXX_FLAGS"
elif test $TARGET = default ; then
  CXX_FLAGS="$CXX_FLAGS"
elif test $TARGET = debug ; then
  CXX_FLAGS="-mtune=generic -msse2 -g $CXX_FLAGS"
elif test $TARGET = empty ; then
  : ## don't do anything
fi


echo Config in `pwd`
echo Config: using CXX_FLAG=$CXX_FLAGS
echo Config: using ADDITIONAL_LIBS=$ADDITIONAL_LIBS
echo Config: using ADDITIONAL_CPP_FLAGS=$ADDITIONAL_CPP_FLAGS
echo Config: using CUDA_PATH=$CUDA_PATH
echo Config: using CUDA_OBJECTS=$CUDA_OBJECTS



AC_SUBST(RANDOMFIELDS_CXXFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
