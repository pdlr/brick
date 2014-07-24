export INSTALLDIR="/var/tmp"
export COMMONFLAGS="-Wall -Wextra -Werror -Wshadow -O2"

export CPPFLAGS="-I$INSTALLDIR/include"
export LDFLAGS="-L$INSTALLDIR/lib"
export CFLAGS="$COMMONFLAGS"
export CXXFLAGS="$COMMONFLAGS"

./configure --prefix=$INSTALLDIR
