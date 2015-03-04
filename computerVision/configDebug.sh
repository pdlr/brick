export INSTALLDIR="/var/tmp"
# export COMMONFLAGS="-Wall -Werror -Wextra -Wunsafe-loop-optimizations -Wshadow -g"
export COMMONFLAGS="-Wall -Werror -Wextra -Wshadow -g -DBRICK_NUMERIC_CHECKBOUNDS=1"

export CPPFLAGS="-I$INSTALLDIR/include"
export LDFLAGS="-L$INSTALLDIR/lib"
export CFLAGS="$COMMONFLAGS"
export CXXFLAGS="$COMMONFLAGS"

./configure --prefix=$INSTALLDIR
