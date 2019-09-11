CXXFLAGS=-Wall -O3 -std=c++11 -DNDEBUG -march=native
ZLIB=/mnt/lustre/groups/ukbiobank/Edinburgh_Data/Software/PRSice-cpp_development/PRSice.code/lib/zlib-1.2.11/build/libz.a 
CXX=/opt/apps/compilers/gcc/6.2.0/bin/g++
INCLUDES := -I inc/ -isystem lib/ -isystem lib/zlib-1.2.11/
THREAD := -Wl,--whole-archive -lpthread
SERVER := -L /usr/lib/x86_64-redhat-linux5E/lib64
GCC := -Wl,--no-whole-archive  -static-libstdc++ -static-libgcc -static
CSRC := src/*.c
CPPSRC := src/*.cpp
OBJ := gzstream.o bgen_lib.o binaryplink.o genotype.o misc.o dcdflib.o regression.o snp.o binarygen.o commander.o main.o plink_common.o prsice.o region.o reporter.o fastlm.o

%.o: src/%.c
		$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

%.o: src/%.cpp
		$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

PRSice: $(OBJ)
		$(CXX) $(CXXFLAGS) $(INCLUDES) $(SERVER)  $^ $(ZLIB) $(THREAD) $(GCC) -o $@
