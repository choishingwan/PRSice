#ifndef PIPELINE_FUNCTIONS_HPP
#define PIPELINE_FUNCTIONS_HPP
#include "IITree.h"
#include "genotype.hpp"
#include "genotypefactory.hpp"
#include "misc.hpp"
#include "region.hpp"
#include "reporter.hpp"
#include <exception>
#include <ostream>
#include <string>
#include <vector>

inline void
print_empty_region(const std::string& out,
                   const std::vector<std::vector<size_t>>& region_membership,
                   const std::vector<std::string>& region_names)
{
    bool has_empty_region = false;
    std::ostringstream empty_regions;
    std::string empty_region_name = out + ".xregion";
    // region_start_idx size always = num_regions
    // check regions to see if there are any empty regions
    for (size_t region_idx = 2; region_idx < region_membership.size();
         ++region_idx)
    {
        if (region_membership[region_idx].empty())
        {
            has_empty_region = true;
            empty_regions << region_names[region_idx] << std::endl;
        }
    }
    if (has_empty_region)
    {
        std::ofstream empty_region_file;
        empty_region_file.open(empty_region_name.c_str());
        if (!empty_region_file.is_open())
        {
            throw std::runtime_error(
                "Error: Cannot open file: " + empty_region_name + " to write!");
        }
        empty_region_file << empty_regions.str();
        empty_region_file.close();
    }
}

inline void initialize_genotype(
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    const Commander& commander, Genotype* current_file, Reporter& reporter,
    Genotype* target_file)
{
    bool is_ref = true;
    std::string type = "reference";
    const std::string separator =
        "==================================================";
    current_file = &current_file->keep_nonfounder(commander.nonfounders())
                        .keep_ambig(commander.keep_ambig())
                        .intermediate(commander.use_inter())
                        .set_prs_instruction(commander.get_prs_instruction())
                        .set_weight();

    if (target_file == nullptr)
    {
        is_ref = false;
        type = "target";
        const std::string base_name = commander.get_base_name();
        reporter.report("Start processing " + base_name + "\n" + separator);
        current_file->snp_extraction(commander.extract_file(),
                                     commander.exclude_file());
        auto [filter_count, dup_rs_id] = current_file->read_base(
            commander.get_base(), commander.get_base_qc(),
            commander.get_p_threshold(), exclusion_regions);
        current_file->print_base_stat(filter_count, dup_rs_id, commander.out(),
                                      commander.get_base_qc().info_score);
        current_file->set_thresholds(commander.get_target_qc());
    }
    if (is_ref)
    {
        current_file = &current_file->reference();
        current_file->set_thresholds(commander.get_ref_qc());
    }
    reporter.report("Loading Genotype info from " + type + "\n" + separator);
    if (!is_ref && commander.use_ref()) current_file->expect_reference();
    current_file->load_samples();
    const bool verbose = true;
    current_file->load_snps(commander.out(), exclusion_regions, verbose);
}

inline void
initialize_target(const std::vector<IITree<size_t, size_t>>& exclusion_regions,
                  const Commander& commander, Genotype* target_file,
                  Reporter& reporter)
{
    initialize_genotype(exclusion_regions, commander, target_file, reporter,
                        nullptr);
}
inline void initialize_reference(
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    const Commander& commander, Genotype* target_file, Genotype* reference_file,
    Reporter& reporter)
{
    initialize_genotype(exclusion_regions, commander, reference_file, reporter,
                        target_file);
}

inline std::tuple<std::vector<std::string>, size_t>
add_gene_set_info(const Commander& commander, Genotype* target_file,
                  Reporter& reporter)
{
    Region region(commander.get_set(), &reporter);
    const size_t num_regions = region.generate_regions(
        target_file->included_snps_idx(), target_file->included_snps(),
        target_file->max_chr());
    target_file->add_flags(region.get_gene_sets(), num_regions,
                           commander.get_set().full_as_background);
    return {region.get_names(), num_regions};
}

#endif // PIPELINE_FUNCTIONS_HPP
