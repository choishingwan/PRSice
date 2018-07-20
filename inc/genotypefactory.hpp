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
    std::unordered_map<std::string, int> file_type{{"bed", 0},
                                                   {"ped", 1},
                                                   {"bgen", 2}};

public:
    Genotype* createGenotype(const std::string& prefix, const std::string& type,
                             const std::string& multi_input, const int thread,
                             const bool ignore_fid, const bool keep_nonfounder,
                             const bool keep_ambig, Reporter& reporter,
                             const Commander& commander,
                             const bool is_ref = false)
    {
        std::vector<std::string> external_sample = misc::split(prefix, ",");
        std::string sample_file = "";
        std::string binary_file = prefix;
        const bool intermediate = commander.intermediate();
        if (external_sample.size() > 1) {
            sample_file = external_sample[1];
            binary_file = external_sample[0];
        }
        int code =
            (file_type.find(type) != file_type.end()) ? file_type[type] : 0;
        switch (code)
        {
        case 0:
        {
            std::string message =
                "Loading Genotype file: " + binary_file + " (bed)\n";
            if (!multi_input.empty())
                message = "Loading Genotype info from file (bed) \n";
            if (!sample_file.empty()) {
                message.append("With external fam file: " + sample_file + "\n");
            }
            reporter.report(message);

            return new BinaryPlink(binary_file, sample_file, multi_input,
                                   thread, ignore_fid, keep_nonfounder,
                                   keep_ambig);
        }
        case 2:
        {
            std::string message =
                "Loading Genotype file: " + binary_file + " (bgen)\n";
            if (!multi_input.empty())
                message = "Loading Genotype info from file (bgen) \n";
            if (!sample_file.empty()) {
                message.append("With sample file: " + sample_file + "\n");
            }
            reporter.report(message);
            if (commander.pheno_file().empty() && !commander.no_regress()
                && sample_file.empty())
            {
                throw std::runtime_error("ERROR: You must provide a phenotype "
                                         "file for bgen format!\n");
            }
            if (sample_file.empty()) sample_file = commander.pheno_file();
            return new BinaryGen(binary_file, sample_file, multi_input, thread,
                                 ignore_fid, keep_nonfounder, keep_ambig,
                                 is_ref, intermediate);
        }
        default:
            throw std::invalid_argument("ERROR: Only support bgen and bed");
        }
    }
};

#endif
