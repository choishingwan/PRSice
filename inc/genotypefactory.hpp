// This file is part of PRSice-2, copyright (C) 2016-2019
// Shing Wan Choi, Paul F. O’Reilly
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
    const std::unordered_map<std::string, int> file_type {{"bed", 0},
                                                          {"ped", 1},
                                                          {"bgen", 2}};

public:
    Genotype* createGenotype(const Commander& commander, Reporter& reporter,
                             const bool is_ref = false)
    {
        const std::string type =
            (is_ref) ? commander.ref_type() : commander.target_type();
        int code = 0;
        if (file_type.find(type) != file_type.end())
        { code = file_type.at(type); }
        else
        {
            throw std::invalid_argument("ERROR: Only support bgen and bed");
        }
        switch (code)
        {
        case 0:
        {
            if (is_ref)
            {
                return new BinaryPlink(
                    commander.ref_list(), commander.ref_name(),
                    commander.thread(), commander.ignore_fid(),
                    commander.nonfounders(), commander.keep_ambig(), is_ref,
                    &reporter);
            }
            else
            {
                return new BinaryPlink(
                    commander.target_list(), commander.target_name(),
                    commander.thread(), commander.ignore_fid(),
                    commander.nonfounders(), commander.keep_ambig(), is_ref,
                    &reporter);
            }
        }
        case 2:
        {
            if (is_ref)
            {
                // for reference bgen, we ALWAYS hard code it. so we don't
                // really use the hard_coded flag
                return new BinaryGen(
                    commander.ref_list(), commander.ref_name(),
                    commander.pheno_file(), commander.out(),
                    commander.get_ref_hard_threshold(),
                    commander.ref_dose_thres(), commander.thread(),
                    commander.use_inter(), commander.hard_coded(),
                    commander.no_regress(), commander.ignore_fid(),
                    commander.nonfounders(), commander.keep_ambig(), is_ref,
                    &reporter);
            }
            else
            {
                return new BinaryGen(
                    commander.target_list(), commander.target_name(),
                    commander.pheno_file(), commander.out(),
                    commander.get_target_hard_threshold(),
                    commander.target_hard_threshold(), commander.thread(),
                    commander.use_inter(), commander.hard_coded(),
                    commander.no_regress(), commander.ignore_fid(),
                    commander.nonfounders(), commander.keep_ambig(), is_ref,
                    &reporter);
            }
        }
        default:
            throw std::invalid_argument("ERROR: Only support bgen and bed");
        }
    }
};

#endif
