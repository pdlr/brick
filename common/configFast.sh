export INSTALLDIR="/var/tmp"
export COMMONFLAGS="-Wall -Wextra -Wunsafe-loop-optimizations -Wshadow -O2"

export CPPFLAGS="-I$INSTALLDIR"
export CFLAGS="$COMMONFLAGS"
export CXXFLAGS="$COMMONFLAGS"

./configure --prefix=$INSTALLDIR
