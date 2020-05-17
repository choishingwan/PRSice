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


#include "IITree.h"
#include "cgranges.h"
#include "commander.hpp"
#include "genotype.hpp"
#include "genotypefactory.hpp"
#include "pipeline_functions.hpp"
#include "plink_common.hpp"
#include "prsice.hpp"
#include "region.hpp"
#include "reporter.hpp"
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>


int main(int argc, char* argv[])
{
    const std::string separator =
        "==================================================";
    // initialize reporter, use to generate log
    Reporter reporter;
    try
    {
        // initialize commander, use to parse command line arguments
        // initialization should help setting the default values
        Commander commander;
        try
        {
            if (!commander.process_command(argc, argv, reporter))
            {
                return 0; // only require the usage information
            }
        }
        catch (...)
        {
            return -1; // all error messages should have printed
        }
        // parse the exclusion range and put it into the exclusion object
        // Generate the exclusion region
        std::vector<IITree<size_t, size_t>> exclusion_regions;
        Region::generate_exclusion(exclusion_regions,
                                   commander.exclusion_range());
        GenomeFactory factory;
        Genotype *target_file = nullptr, *reference_file = nullptr;
        try
        {
            // initialize the target object using the factory
            target_file = factory.createGenotype(commander.get_target(),
                                                 commander.get_pheno(),
                                                 commander.delim(), reporter);
            initialize_target(exclusion_regions, commander, target_file,
                              reporter);
            if (commander.use_ref() && commander.need_ref())
            {
                reference_file = factory.createGenotype(
                    commander.get_reference(), commander.get_pheno(),
                    commander.delim(), reporter);
                initialize_reference(exclusion_regions, commander, target_file,
                                     reference_file, reporter);
            }
            exclusion_regions.clear();
            target_file->calc_freqs_and_intermediate(commander.get_target_qc(),
                                                     commander.out(), true);
            if (reference_file != nullptr)
            {
                reference_file->set_thresholds(commander.get_ref_qc());
                reference_file->calc_freqs_and_intermediate(
                    commander.get_ref_qc(), commander.out(), true, target_file);
            }
            if (target_file->num_snps() == 0)
            {
                reporter.report("No SNPs left for PRSice processing");
                return -1;
            }
            const auto [region_names, num_regions] =
                add_gene_set_info(commander, target_file, reporter);
            auto prs_instruction = commander.get_prs_instruction();
            auto pheno_info = commander.get_pheno();
            const bool no_regress = prs_instruction.no_regress;
            PRSice::pheno_check(no_regress, pheno_info, reporter);

            if (!commander.get_clump_info().no_clump)
            {
                target_file->build_clump_windows(
                    commander.get_clump_info().distance);
                target_file->sort_by_p();
                // now perform clumping
                target_file->clumping(commander.get_clump_info(),
                                      commander.use_ref() ? *reference_file
                                                          : *target_file,
                                      commander.get_prs_instruction().thread);
            }
            // immediately free the memory
            if (reference_file != nullptr) { delete reference_file; }
            if (commander.ultra_aggressive())
            { target_file->load_genotype_to_memory(); }
            target_file->prepare_prsice();
            // from now on, we are not allow to sort the m_existed_snps
            auto snp_file = commander.print_snp()
                                ? misc::load_ostream(commander.out() + ".snp")
                                : nullptr;
            // vector containing the index for each SNP in each set
            // structure is [vec of Set][vec of SNP]
            auto region_membership = target_file->build_membership_matrix(
                num_regions, region_names, commander.print_snp(),
                *snp_file.get());
            // we can now quickly check if any of the region are empty
            try
            {
                print_empty_region(commander.out(), region_membership,
                                   region_names);
            }
            catch (const std::runtime_error& er)
            {
                reporter.report(er.what());
                return -1;
            }
            // Initialize the progress bar
            // one progress bar per phenotype
            // one extra progress bar for competitive permutation
            assert(target_file->get_set_thresholds().size() == num_regions);

            const auto [max_fid, max_iid] = target_file->get_max_id_length();
            const size_t num_pheno = pheno_info.pheno_col_idx.size();
            // prsice and summary file will be per run
            // all score and best file will be per phenotype
            // this is mainly because of the size of the file and the way we
            // need to generate them
            // with prsice and summary file, we can do row wise output, so it is
            // ok for us to keep using it
            // but for all and best, we are doing column-wise output, which need
            // expensive seek operations
            std::unique_ptr<std::ostream> summary_file = nullptr;
            auto prsice_out = misc::load_ostream(commander.out() + ".prsice");
            bool has_prevalence = !pheno_info.prevalence.empty();
            print_prsice_header(has_prevalence, no_regress, prsice_out);
            auto perm_info = commander.get_perm();
            if (!no_regress)
            {
                summary_file = misc::load_ostream(commander.out() + ".summary");
                print_summary_header(has_prevalence, perm_info.run_set_perm,
                                     perm_info.run_perm, summary_file);
            }
            for (size_t i_pheno = 0; i_pheno < num_pheno; ++i_pheno)
            {
                if (pheno_info.skip_pheno[i_pheno])
                {
                    reporter.simple_report("Skipping the "
                                           + std::to_string(i_pheno + 1)
                                           + " th phenotype");
                    continue;
                }
                reporter.simple_report("Processing the "
                                       + std::to_string(i_pheno + 1)
                                       + " th phenotype");
                PRSice prsice(commander.get_prs_instruction(),
                              commander.get_p_threshold(), perm_info,
                              commander.out(), pheno_info.binary[i_pheno],
                              &reporter);
                prsice.init_progress_count(target_file->get_set_thresholds());
                std::unique_ptr<std::ostream> best_file = nullptr,
                                              all_score_file = nullptr;
                if (!no_regress)
                {
                    prsice.init_matrix(
                        pheno_info.cov_colname, pheno_info.col_index_of_cov,
                        pheno_info.col_index_of_factor_cov,
                        pheno_info.pheno_file, pheno_info.cov_file,
                        pheno_info.pheno_col[i_pheno], commander.delim(),
                        pheno_info.pheno_col_idx[i_pheno],
                        pheno_info.ignore_fid, *target_file);
                    best_file = misc::load_ostream(
                        commander.out()
                        + ((num_pheno > 1) ? "." + pheno_info.pheno_col[i_pheno]
                                           : "")
                        + ".best");
                    prsice.prep_best_output(*target_file, region_membership,
                                            region_names, max_fid, max_iid,
                                            best_file);
                }
                if (commander.all_scores())
                {
                    all_score_file = misc::load_ostream(
                        commander.out()
                        + ((num_pheno > 1) ? "." + pheno_info.pheno_col[i_pheno]
                                           : "")
                        + ".all_score");
                    prsice.prep_all_score_output(
                        *target_file, region_membership, region_names, max_fid,
                        max_iid, all_score_file);
                }
                // go through each region
                fprintf(stderr, "\nStart Processing\n");
                for (size_t i_region = 0; i_region < num_regions; ++i_region)
                {
                    // always skip background region and empty regions
                    if (i_region == 1 || region_membership[i_region].empty())
                        continue;
                    prsice.run_prsice(i_pheno, i_region, region_membership,
                                      commander.all_scores(), *target_file);
                    prsice.output(region_names, i_pheno, i_region);
                }
                if (!no_regress)
                {
                    prsice.print_best(*target_file, region_names, i_pheno);
                    if (perm_info.run_set_perm && region_names.size() > 2)
                    {
                        // only perform permutation if regression is performed
                        // and user request it
                        assert(region_membership.size() >= 2);
                        prsice.run_competitive(
                            *target_file, region_membership[1].begin(),
                            region_membership[1].end(), i_pheno);
                    }
                }
            }
            // prsice.print_progress(true);
            // fprintf(stderr, "\n");
            /*
            if (!commander.get_prs_instruction().no_regress)
                // now generate the summary file
                prsice.summarize();
                */
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
