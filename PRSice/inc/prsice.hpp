#ifndef PRSICE_H
#define PRSICE_H

#include <string>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <map>
#include <stdio.h>
#include "commander.hpp"
#include "misc.hpp"
#include "plink.hpp"
#include "snp.hpp"
#include "region.hpp"

//This should be the class to handle all the procedures
class PRSice
{
    public:
        PRSice();
        virtual ~PRSice();
        void run(const Commander &c_commander, Region &region);
        void process(const std::string &c_input, const Commander &c_commander, Region &region);
    protected:
    private:
        void get_inclusion(std::map<std::string, bool> &inclusion, const std::string &target_bim_name,
                       std::vector<SNP> &snp_list, const std::map<std::string, size_t> &snp_index);
        void score(const std::map<std::string, bool> &inclusion, const std::string target_name,
                   const std::map<std::string, size_t> &snp_index, std::vector<SNP> &snp_list);
};

#endif // PRSICE_H
