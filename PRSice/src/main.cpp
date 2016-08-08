#include <iostream>
#include <cstdlib>
#include <string>
#include <boost/ptr_container/ptr_vector.hpp>
#include "snp.h"
#include "commander.h"
#include "prsice.h"



int main(int argc, char *argv[])
{
    Commander commander;
    try{
        commander(argc, argv);
    }
    catch(...){
    }
    boost::ptr_vector<SNP> snp_list;
    SNP::read_snp(commander, snp_list);
    PRSice prsice = PRSice();
    prsice.run();
    return 0;
}
