# g++ -std=gnu++23 -Ofast -march=native -funroll-loops test.cpp /usr/lib/x86_64-linux-gnu/libprofiler.so /usr/lib/x86_64-linux-gnu/libtcmalloc.so -o test
g++ -std=gnu++23 -Ofast -march=native -funroll-loops test.cpp -o test
./test