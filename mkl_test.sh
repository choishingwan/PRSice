source /opt/intel/mkl/bin/mklvars.sh intel64


g++ -DEIGEN_USE_MKL_ALL -DNDEBUG  -DMKL_LP64 -m64 -I${MKLROOT}/include -march=native -O3 -std=c++11 -isystem lib/ -I inc/ src/*.cpp ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -lpthread -lm -ldl -lz -o PRSice_mac

