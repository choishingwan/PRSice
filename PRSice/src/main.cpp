#include <iostream>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <boost/ptr_container/ptr_vector.hpp>
#include "commander.h"


int main(int argc, char *argv[])
{
    Commander commander = Commander();
    try{
        if(!commander.initialize(argc, argv)) return 0; //only require the usage information
    }
    catch (const std::runtime_error& error){
        std::cerr << error.what() << std::endl;
    }
//    boost::ptr_vector<SNP> snp_list;
//    SNP::read_snp(commander, snp_list);
//    PRSice prsice = PRSice();
//    prsice.run();
    return 0;
}
