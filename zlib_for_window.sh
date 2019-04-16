/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc -O3 -Wall -c -o adler32.o adler32.c
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc -O3 -Wall -c -o compress.o compress.c
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc -O3 -Wall -c -o crc32.o crc32.c
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc -O3 -Wall -c -o deflate.o deflate.c
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc -O3 -Wall -c -o gzclose.o gzclose.c
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc -O3 -Wall -c -o gzlib.o gzlib.c
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc -O3 -Wall -c -o gzread.o gzread.c
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc -O3 -Wall -c -o gzwrite.o gzwrite.c
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc -O3 -Wall -c -o infback.o infback.c
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc -O3 -Wall -c -o inffast.o inffast.c
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc -O3 -Wall -c -o inflate.o inflate.c
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc -O3 -Wall -c -o inftrees.o inftrees.c
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc -O3 -Wall -c -o trees.o trees.c
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc -O3 -Wall -c -o uncompr.o uncompr.c
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc -O3 -Wall -c -o zutil.o zutil.c

/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-ar rcs libz.a adler32.o compress.o crc32.o deflate.o gzclose.o gzlib.o gzread.o gzwrite.o infback.o inffast.o inflate.o inftrees.o trees.o uncompr.o zutil.o

/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-windres --define GCC_WINDRES -o zlibrc.o win32/zlib1.rc

/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc -shared -Wl,--out-implib,libz.dll.a -o zlib1.dll win32/zlib.def adler32.o compress.o crc32.o deflate.o gzclose.o gzlib.o gzread.o gzwrite.o infback.o inffast.o inflate.o inftrees.o trees.o uncompr.o zutil.o  zlibrc.o
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-strip zlib1.dll
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc  -O3 -Wall -I. -c -o example.o test/example.c
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc  -o example.exe example.o libz.a
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-strip example.exe
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc  -O3 -Wall -I. -c -o minigzip.o test/minigzip.c
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc  -o minigzip.exe minigzip.o libz.a
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-strip minigzip.exe
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc  -o example_d.exe example.o libz.dll.a
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-strip example_d.exe
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-gcc  -o minigzip_d.exe minigzip.o libz.dll.a
/usr/local/Cellar/mingw-w64/5.0.4_1/toolchain-x86_64/bin/x86_64-w64-mingw32-strip minigzip_d.exe



