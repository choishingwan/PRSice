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


#include "plink_common.hpp"
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

#include "IITree.h"
#include "cgranges.h"
#include "commander.hpp"
#include "genotype.hpp"
#include "genotypefactory.hpp"
#include "prsice.hpp"
#include "region.hpp"
#include "reporter.hpp"
int main(int argc, char* argv[])
{
    // initialize reporter, use to generate log
    Reporter reporter;
    try
    {
        // initialize commander, use to parse command line arguments
        // initialization should help setting the default values
        Commander commander;
        try
        {
            if (!commander.init(argc, argv, reporter))
            {
                return 0; // only require the usage information
            }
        }
        catch (...)
        {
            return -1; // all error messages should have printed
        }
        bool verbose = true;
        // parse the exclusion range and put it into the exclusion object
        // Generate the exclusion region
        std::vector<IITree<size_t, size_t>> exclusion_regions;
        Region::generate_exclusion(exclusion_regions,
                                   commander.exclusion_range());
        // this allow us to generate the appropriate object (i.e. binaryplink /
        // binarygen)
        double maf, geno, info;
        bool init_ref = false;
        // load the filtering parameters for the target file
        bool maf_filter = commander.target_maf(maf);
        bool geno_filter = commander.target_geno(geno);
        bool info_filter = commander.target_info(info);
        bool hard_coded = commander.hard_coded();
        GenomeFactory factory;
        Genotype *target_file = nullptr, *reference_file = nullptr;
        try
        {
            // initialize the target object using the factory
            target_file = factory.createGenotype(commander, reporter);
            // load base file into memory
            const std::string base_name = misc::remove_extension<std::string>(
                misc::base_name<std::string>(commander.base_name()));
            std::string message = "Start processing " + base_name + "\n";
            message.append(
                "==================================================");
            reporter.report(message);
            target_file->snp_extraction(commander.extract_file(),
                                        commander.exclude_file());
            target_file->read_base(
                commander.base_name(), commander.index(), commander.has_col(),
                commander.bar_levels(), commander.lower(), commander.inter(),
                commander.upper(), exclusion_regions,
                commander.maf_base_control(), commander.maf_base_case(),
                commander.base_info_score(), commander.fastscore(),
                commander.no_full(), commander.beta(), commander.is_index(),
                commander.keep_ambig());
            // no longer need the exclusion region object
            // then we will read in the sample information
            message = "Loading Genotype info from target\n";
            message.append(
                "==================================================");
            reporter.report(message);
            target_file->load_samples(commander.keep_sample_file(),
                                      commander.remove_sample_file(),
                                      commander.delim(), verbose);
            // For bgen or any format that require intermediate generation, we
            // need to know if a reference file is going to be used, therefore
            // decide if we are going to generate the target intermediate or the
            // reference intermediate
            if (commander.use_ref()) target_file->expect_reference();
            // read in SNPs included in the GWAS summary statistic
            // don't store or do anything, but just check which SNPs are
            // included so that we can ignore SNPs not found in GWAS
            // when we do geno and maf
            // Finally, we can read in the SNP information
            target_file->load_snps(commander.out(), exclusion_regions, verbose);
            // now load the reference file
            // initialize the memory map file
            target_file->init_mmap();
            if (commander.use_ref()
                && (!commander.no_clump() || commander.use_ref_maf()))
            {
                message = "Start processing reference\n";
                reporter.report(message);
                reference_file =
                    factory.createGenotype(commander, reporter, true);
                init_ref = true;
                message = "Loading Genotype info from reference\n";
                message.append(
                    "==================================================");
                reporter.report(message);
                reference_file->load_samples(commander.ref_keep_file(),
                                             commander.ref_remove_file(),
                                             commander.delim(), verbose);
                // load the reference file
                reference_file->load_snps(commander.out(), exclusion_regions,
                                          verbose, target_file);
            }
            exclusion_regions.clear();
            // with the reference file read, we can start doing filtering and
            // calculate relevent metric

            if (maf_filter || geno_filter || info_filter
                || (commander.use_inter()
                    && (hard_coded || !commander.use_ref())))
            {
                message =
                    "Calculate MAF and perform filtering on target SNPs\n";
                message.append(
                    "==================================================");
                reporter.report(message);
                // only calculate the MAF if we need to
                // We want to only invoke the MAF calculation if we need to
                // i.e after clumping, to speed up the process
                target_file->calc_freqs_and_intermediate(
                    maf, geno, info, maf_filter, geno_filter, info_filter,
                    hard_coded, true);
            }
            maf_filter = commander.ref_maf(maf);
            geno_filter = commander.ref_geno(geno);
            info_filter = commander.ref_info(info);
            hard_coded = commander.ref_hard_threshold();
            bool has_ref_maf = false;
            if (init_ref)
            {
                reference_file->init_mmap();
                if (maf_filter || geno_filter || info_filter
                    || commander.use_inter())
                {
                    // we only go through the reference file if we are
                    // 1. Need the reference MAF
                    // 2. Need to filter the reference file (need hard code
                    // info)
                    // 3. Need to generate an intermediate file for clumping
                    message = "Calculate MAF and perform filtering on "
                              "reference SNPs\n";
                    message.append("======================================="
                                   "===========");
                    reporter.report(message);
                    reference_file->calc_freqs_and_intermediate(
                        maf, geno, info, maf_filter, geno_filter, info_filter,
                        hard_coded, true, target_file);
                    has_ref_maf = true;
                }
            }
            // now should get the correct MAF and should have filtered the
            // SNPs accordingly Generate Region flag information
            Region region(commander, &reporter);
            std::unordered_map<std::string, std::vector<size_t>> snp_in_sets;
            std::vector<IITree<size_t, size_t>> gene_sets;
            size_t num_regions = region.generate_regions(
                gene_sets, snp_in_sets, target_file->max_chr());
            std::vector<std::string> region_names = region.get_names();
            target_file->add_flags(gene_sets, snp_in_sets, num_regions,
                                   commander.genome_wide_background());

            gene_sets.clear();
            // start processing other files before doing clumping
            PRSice prsice(commander, num_regions > 2, &reporter);
            prsice.pheno_check(commander.pheno_file(), commander.pheno_col(),
                               commander.is_binary());
            // Store relevant parameters to the target object
            target_file->set_info(commander);
            if (!commander.no_clump())
            {
                // now go through the snp vector an define the
                // windows so that we can jump directly to the
                // relevant SNPs immediately when doing clumping
                target_file->build_clump_windows();
                // get the sort by p index vector for target
                // so that we can still find out the relative coordinates of
                // each SNPs This is only required for clumping
                if (!target_file->sort_by_p())
                {
                    std::string error_message =
                        "No SNPs left for PRSice processing";
                    reporter.report(error_message);
                    return -1;
                }
                // now perform clumping
                target_file->efficient_clumping(
                    commander.use_ref() ? *reference_file : *target_file);
                // immediately free the memory
            }
            if (init_ref)
            {
                if (commander.use_ref_maf() && !has_ref_maf)
                {
                    message = "Calculate MAF based on Reference\n";
                    message.append("======================================="
                                   "===========");
                    reporter.report(message);
                    reference_file->calc_freqs_and_intermediate(
                        maf, geno, info, maf_filter, geno_filter, info_filter,
                        hard_coded, true, target_file);
                }
                delete reference_file;
            }
            // We now have problem:
            // MAF has not been calculated, and can therefore cause a
            // problem in PRS calculation

            // can do the update structure here
            // Use sparse matrix for space and speed
            // Column = set, row = SNPs (because EIGEN is column major)
            // need to also know the number of threshold included
            std::vector<size_t> region_membership;
            std::vector<size_t> region_start_idx;
            std::vector<size_t>::const_iterator background_start_idx,
                background_end_idx;
            target_file->prepare_prsice();
            target_file->build_membership_matrix(
                region_membership, region_start_idx, num_regions,
                commander.out(), region_names, commander.print_snp());
            background_start_idx = region_membership.cbegin();
            std::advance(background_start_idx,
                         static_cast<long>(region_start_idx[1]));
            background_end_idx = region_membership.cbegin();
            if (num_regions > 2)
            {
                std::advance(background_end_idx,
                             static_cast<long>(region_start_idx[2]));
            }
            // we can now quickly check if any of the region are empty
            bool has_empty_region = false;
            std::ofstream empty_region;
            std::string empty_region_name = commander.out() + ".xregion";
            // region_start_idx size always = num_regions
            // check regions to see if there are any empty regions
            for (size_t i = 2; i < region_start_idx.size(); ++i)
            {
                size_t cur_idx = region_start_idx[i];
                if (i + 1 >= region_start_idx.size())
                {
                    // this is the last set
                    if (cur_idx == region_membership.size())
                    {
                        if (!has_empty_region)
                        {
                            empty_region.open(empty_region_name.c_str());
                            if (!empty_region.is_open())
                            {
                                std::string error_message =
                                    "Error: Cannot open file: "
                                    + empty_region_name + " to write!";
                                reporter.report(error_message);
                                return -1;
                            }
                            has_empty_region = true;
                        }
                        empty_region << region_names[i] << std::endl;
                    }
                }
                else if (cur_idx == region_start_idx[i + 1])
                {
                    if (!has_empty_region)
                    {
                        empty_region.open(empty_region_name.c_str());
                        if (!empty_region.is_open())
                        {
                            std::string error_message =
                                "Error: Cannot open file: " + empty_region_name
                                + " to write!";
                            reporter.report(error_message);
                            return -1;
                        }
                        has_empty_region = true;
                    }
                    empty_region << region_names[i] << std::endl;
                }
            }
            if (has_empty_region) empty_region.close();
            const size_t num_pheno = prsice.num_phenotype();

            // Initialize the progress bar
            prsice.init_process_count(commander, num_regions,
                                      target_file->num_threshold());
            for (size_t i_pheno = 0; i_pheno < num_pheno; ++i_pheno)
            {
                fprintf(stderr, "Processing the %zu th phenotype\n",
                        i_pheno + 1);
                prsice.init_matrix(commander, i_pheno, commander.delim(),
                                   *target_file);
                fprintf(stderr, "Preparing Output Files\n");
                prsice.prep_output(commander.out(), commander.all_scores(),
                                   commander.has_prevalence(), *target_file,
                                   region_names, i_pheno,
                                   commander.no_regress());
                // go through each region
                fprintf(stderr, "\nStart Processing\n");
                for (size_t i_region = 0; i_region < num_regions; ++i_region)
                {
                    // always skip background region
                    if (i_region == 1) continue;
                    if (!prsice.run_prsice(commander, i_pheno, i_region,
                                           region_membership, region_start_idx,
                                           *target_file))
                    {
                        // did not run
                        continue;
                    }
                    if (!commander.no_regress())
                    {
                        // if we performed regression, we'd like to generate
                        // the output file (.prsice)
                        prsice.output(commander, region_names, i_pheno,
                                      i_region);
                    }
                    else
                    {
                        prsice.no_regress_out(commander, region_names, i_pheno,
                                              i_region);
                    }
                }
                if (!commander.no_regress() && commander.perform_set_perm())
                {
                    // only perform permutation if regression is performed
                    // and user request it
                    prsice.run_competitive(*target_file, background_start_idx,
                                           background_end_idx, commander,
                                           i_pheno);
                }
            }
            prsice.print_progress(true);
            fprintf(stderr, "\n");
            if (!commander.no_regress())
                // now generate the summary file
                prsice.summarize(commander);
        }
        catch (const std::invalid_argument& ia)
        {
            reporter.report(ia.what());
            return -1;
        }
        catch (const std::runtime_error& error)
        {
            reporter.report(error.what());
            return -1;
        }
        catch (const std::out_of_range& error)
        {
            reporter.report(error.what());
            return -1;
        }
        catch (const std::exception& ex)
        {
            reporter.report(ex.what());
        }
        catch (...)
        {
            std::string error_message =
                "Error: Bad Allocation exception detected. "
                "This is likely due to insufficient memory "
                "for PRSice. You can try re-running PRSice "
                "with more memory.";
            reporter.report(error_message);
        }
        delete target_file;
    }
    catch (const std::exception& ex)
    {
        reporter.report(ex.what());
    }
    catch (...)
    {
        std::string error_message = "Error: Bad Allocation exception detected. "
                                    "This is likely due to insufficient memory "
                                    "for PRSice. You can try re-running PRSice "
                                    "with more memory.";
        reporter.report(error_message);
    }
    return 0;
}
