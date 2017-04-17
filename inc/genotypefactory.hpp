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

                    return new BinaryGen(prefix,
                            commander.pheno_file(),
                            commander.has_pheno_col(),
                            commander.remove_sample_file(),
                            commander.keep_sample_file(),
                            commander.ignore_fid(), commander.num_auto(),
                            commander.no_x(), commander.no_y(),
                            commander.no_xy(), commander.no_mt(),
                            commander.thread(), verbose);
                default:
                case 0:
                    fprintf(stderr, "(bed)\n");
                    return new BinaryPlink(prefix,
                            commander.remove_sample_file(),
                            commander.keep_sample_file(),
                            commander.ignore_fid(), commander.num_auto(),
                            commander.no_x(), commander.no_y(),
                            commander.no_xy(), commander.no_mt(),
                            commander.thread(), verbose);

            }
        }

};

#endif
