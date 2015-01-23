export INSTALLDIR="/var/tmp"
export COMMONFLAGS="-Wall -Wextra -Werror -Wshadow -g -DBRICK_NUMERIC_CHECKBOUNDS=1"

export CPPFLAGS="-I$INSTALLDIR/include"
export LDFLAGS="-L$INSTALLDIR/lib"
export CFLAGS="$COMMONFLAGS"
export CXXFLAGS="$COMMONFLAGS -std=c++11"

./configure --prefix=$INSTALLDIR
