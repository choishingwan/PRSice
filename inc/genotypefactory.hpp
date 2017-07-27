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
#include "genotype.hpp"
#include "binaryplink.hpp"
#include "binarygen.hpp"

class GenomeFactory
{
    private:
        std::unordered_map<std::string, int> file_type { { "bed", 0 }, { "ped", 1 }, { "bgen", 2 } };
    public:
        Genotype* createGenotype(Commander &commander,
                const std::string &prefix, const std::string &type,
                bool verbose)
        {
            fprintf(stderr, "Loading Genotype file: %s ", prefix.c_str());
            int code = (file_type.find(type) != file_type.end()) ? file_type[type] : 0;
            switch (code)
            {
                case 1:
                    /*
                     return std::unique_ptr<Genotype>(new Plink(prefix, commander.num_auto(),
                     commander.no_x(), commander.no_y(), commander.no_xy(), commander.no_mt(),
                     commander.thread(), verbose));
                     */
                case 2:
                    fprintf(stderr, "(bgen)\n");
                    if(commander.pheno_file().empty() && !commander.no_regress())
                    {
                        throw std::runtime_error("ERROR: You must provide a phenotype file for bgen format!\n");
                    }
                    if(!commander.filter_hard_threshold())
                    {
                        fprintf(stderr, "WARNING: Hard threshold was not given, will use default %f\n", commander.hard_threshold());
                    }
                    return new BinaryGen(prefix,
                            commander.pheno_file(),
                            commander.has_pheno_col(),
                            commander.remove_sample_file(),
                            commander.keep_sample_file(),
                            commander.extract_snp_file(),
                            commander.exclude_snp_file(),
                            commander.ignore_fid(), commander.num_auto(),
                            commander.no_x(), commander.no_y(),
                            commander.no_xy(), commander.no_mt(),
							commander.keep_ambig(),
                            commander.thread(), verbose);
                default:
                case 0:
                    fprintf(stderr, "(bed)\n");
                    return new BinaryPlink(prefix,
                            commander.remove_sample_file(),
                            commander.keep_sample_file(),
                            commander.extract_snp_file(),
                            commander.exclude_snp_file(),
                            commander.ignore_fid(), commander.num_auto(),
                            commander.no_x(), commander.no_y(),
                            commander.no_xy(), commander.no_mt(),
							commander.keep_ambig(),
                            commander.thread(), verbose);

            }
        }

};

#endif
