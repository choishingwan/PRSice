// This file is part of PRSice2.0, copyright (C) 2016-2017
// Shing Wan Choi, Jack Euesden, Cathryn M. Lewis, Paul F. O’Reilly
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef SRC_GENOTYPEFACTORY_HPP_
#define SRC_GENOTYPEFACTORY_HPP_
#include "binarygen.hpp"
#include "binaryplink.hpp"
#include "commander.hpp"
#include "genotype.hpp"
#include "misc.hpp"
#include "reporter.hpp"
#include <fstream>
#include <string>
#include <vector>

class GenomeFactory
{
private:
    const std::unordered_map<std::string, int> file_type{{"bed", 0},
                                                         {"ped", 1},
                                                         {"bgen", 2}};

public:
    Genotype* createGenotype(const Commander& commander, Reporter& reporter,
                             const bool is_ref = false)
    {
        const std::string type =
            (is_ref) ? commander.ref_type() : commander.target_type();
        int code = 0;
        if (file_type.find(type) != file_type.end()) {
            code = file_type.at(type);
        }
        else
        {
            throw std::invalid_argument("ERROR: Only support bgen and bed");
        }
        switch (code)
        {
        case 0: return new BinaryPlink(commander, is_ref, reporter);
        case 2: return new BinaryGen(commander, is_ref,reporter);
        default:
            throw std::invalid_argument("ERROR: Only support bgen and bed");
        }
    }
};

#endif
