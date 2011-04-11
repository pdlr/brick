export COMMONFLAGS="-Wall -Wextra -Wunsafe-loop-optimizations -Wshadow -g"

export CFLAGS="$COMMONFLAGS"

export CXXFLAGS="$COMMONFLAGS"

./configure --prefix=/var/tmp
