export INSTALLDIR="/var/tmp"
export COMMONFLAGS="-Wall -Werror -Wextra -Wunsafe-loop-optimizations -Wshadow -O2"

export CPPFLAGS="-I$INSTALLDIR/include"
export LDFLAGS="-L$INSTALLDIR/lib"
export CFLAGS="$COMMONFLAGS"
export CXXFLAGS="$COMMONFLAGS"

./configure --prefix=$INSTALLDIR
