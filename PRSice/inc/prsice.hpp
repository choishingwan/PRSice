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
};

#endif // PRSICE_H
