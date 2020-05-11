#ifndef MOCK_BINARYGEN_HPP
#define MOCK_BINARYGEN_HPP
#include "binarygen.hpp"
#include "reporter.hpp"

class mock_binarygen : public ::BinaryGen
{
public:
    mock_binarygen() {}
    mock_binarygen(GenoFile& geno, Phenotype& pheno, const std::string& delim,
                   Reporter* reporter)
        : BinaryGen(geno, pheno, delim, reporter)
    {
    }


    size_t test_get_sex_col(const std::string& header,
                            const std::string& format_line)
    {
        return get_sex_col(header, format_line);
    }
    void test_handle_pheno_header(std::unique_ptr<std::istream>& sample)
    {
        handle_pheno_header(sample);
    }
    void set_reporter(Reporter* reporter) { m_reporter = reporter; }
    void set_sample_size(uintptr_t sample_size)
    {
        m_unfiltered_sample_ct = sample_size;
    }
    void update_sample(uintptr_t sample_size)
    {
        set_sample_size(sample_size);
        m_sample_ct = m_unfiltered_sample_ct;
        m_founder_ct = m_unfiltered_sample_ct;
        init_sample_vectors();
        for (size_t i = 0; i < m_unfiltered_sample_ct; ++i)
        {
            SET_BIT(i, m_sample_for_ld.data());
            SET_BIT(i, m_calculate_prs.data());
        }
        post_sample_read_init();
    }

    void set_founder_vector(const std::vector<bool>& founder)
    {
        m_founder_ct = 0;
        for (size_t i = 0; i < founder.size(); ++i)
        {
            if (founder[i])
            {
                SET_BIT(i, m_sample_for_ld.data());
                ++m_founder_ct;
            }
        }
    }
    void test_read_genotype(uintptr_t* genotype, SNP& snp)
    {
        read_genotype(genotype, snp);
    }
    void add_select_sample(const std::string& in)
    {
        m_sample_selection_list.insert(in);
    }
    void change_sample_selection(bool remove) { m_remove_sample = remove; }

    std::vector<uintptr_t> sample_for_ld() const { return m_sample_for_ld; }
    std::vector<uintptr_t> calculate_prs() const { return m_calculate_prs; }
    std::vector<Sample_ID> sample_id() const { return m_sample_id; }
    void gen_bgen_header(const std::string& file_name,
                         uint32_t number_of_snp_blocks,
                         uint32_t number_of_samples, std::string free_data,
                         uint32_t flags)
    {
        std::ofstream dummy(file_name, std::ofstream::binary);
        genfile::bgen::Context context;
        context.number_of_variants = number_of_snp_blocks;
        context.number_of_samples = number_of_samples;
        context.free_data = free_data;
        context.flags = flags;
        genfile::bgen::write_offset(dummy, context.free_data.length() + 20);
        genfile::bgen::write_header_block(dummy, context);
        dummy.close();
    }
    size_t test_transverse_bgen_for_snp(
        const std::vector<IITree<size_t, size_t>>& exclusion_regions,
        const std::string mismatch_snp_record_name, const size_t file_idx,
        std::unique_ptr<std::istream> bgen_file,
        std::unordered_set<std::string>& duplicated_snps,
        std::unordered_set<std::string>& processed_snps,
        std::vector<bool>& retain_snp, bool& chr_error, bool& sex_error,
        Genotype* genotype)
    {
        return transverse_bgen_for_snp(
            exclusion_regions, mismatch_snp_record_name, file_idx,
            std::move(bgen_file), duplicated_snps, processed_snps, retain_snp,
            chr_error, sex_error, genotype);
    }
    void set_context(genfile::bgen::Context context, size_t idx)
    {
        if (m_context_map.size() < idx + 1) { m_context_map.resize(idx); }
        m_context_map[idx] = context;
    }
    void load_context(std::istream& in_file, size_t idx = 0)
    {
        if (m_context_map.size() <= idx + 1) { m_context_map.resize(idx + 1); }
        genfile::bgen::Context context;
        uint32_t offset;
        genfile::bgen::read_offset(in_file, &offset);
        genfile::bgen::read_header_block(in_file, &context);
        context.offset = offset;
        m_context_map[idx] = context;
    }
    void test_init_chr(int num_auto = 22, bool no_x = false, bool no_y = false,
                       bool no_xy = false, bool no_mt = false)
    {
        init_chr(num_auto, no_x, no_y, no_xy, no_mt);
    }
    std::tuple<AlleleCounts, double, double>
    test_plink_generator(std::unique_ptr<std::istream> bgen_file,
                         const std::streampos& bytepos)
    {
        auto cur_idx = 0ul;
        PLINK_generator setter(&m_calculate_prs, m_tmp_genotype.data(),
                               m_hard_threshold, m_dose_threshold);
        // we use tellg to get the location of the variant info, so don't need
        // to do offset jump
        bgen_file->seekg(bytepos);
        // bgen_file->seekg(offset + 4);
        genfile::bgen::read_and_parse_genotype_data_block<PLINK_generator>(
            *bgen_file, m_context_map[cur_idx], setter, &m_buffer1, &m_buffer2);
        AlleleCounts ct;
        setter.get_count(ct.homcom, ct.het, ct.homrar, ct.missing);
        double impute = setter.info_score(INFO::IMPUTE2);
        double mach = setter.info_score(INFO::MACH);
        return {ct, impute, mach};
    }

    std::string gen_mock_snp(const std::vector<std::vector<double>>& geno_prob,
                             std::vector<SNP>& input, uint32_t n_sample,
                             genfile::OrderType& phased,
                             genfile::bgen::Layout& layout,
                             genfile::bgen::Compression& compress)
    {
        if (geno_prob.size() != input.size())
        {
            throw std::runtime_error("Expect each SNP to have their own prob");
        }
        genfile::bgen::Context context;
        context.flags |= compress;
        context.flags |= layout;
        context.number_of_samples = n_sample;
        context.magic = "bgen";
        context.number_of_variants = static_cast<uint32_t>(input.size());
        context.free_data = "Test with automatically generated mock data";
        std::ostringstream mock_file, mock_w_offset;
        genfile::bgen::write_header_block(mock_file, context);
        genfile::bgen::write_little_endian_integer(mock_w_offset,
                                                   context.header_size());
        auto header = mock_w_offset.tellp();
        std::streampos byte_pos;
        for (size_t i = 0; i < input.size(); ++i)
        {
            auto&& snp = input[i];
            std::string chr = std::to_string(snp.chr());
            gen_ostringstream_for_snp(geno_prob[i], snp.rs(), snp.rs(), chr,
                                      snp.loc(), snp.ref(), snp.alt(), phased,
                                      mock_file, byte_pos, context);
            snp.update_file(snp.get_file_idx(), byte_pos + header, false);
        }
        mock_w_offset.write(mock_file.str().data(), mock_file.str().size());
        return (mock_w_offset.str());
    }
    bool test_calc_freq_gen_inter(const QCFiltering& filter_info,
                                  const std::string& prefix)
    {
        return calc_freq_gen_inter(filter_info, prefix, this);
    }
    std::string gen_mock_snp(const std::vector<SNP>& input,
                             uint32_t number_individual,
                             genfile::OrderType& phased,
                             genfile::bgen::Layout& layout,
                             genfile::bgen::Compression& compress)
    {
        genfile::bgen::Context context;
        context.flags |= compress;
        context.flags |= layout;
        context.number_of_samples = number_individual;
        context.magic = "bgen";
        context.number_of_variants = static_cast<uint32_t>(input.size());
        context.free_data = "Test with automatically generated mock data";
        std::ostringstream mock_file, mock_w_offset;
        genfile::bgen::write_header_block(mock_file, context);
        genfile::bgen::write_little_endian_integer(mock_w_offset,
                                                   context.header_size());
        auto multi = (phased == genfile::ePerUnorderedGenotype) ? 3ul : 4ul;
        auto temp = std::vector<double>(context.number_of_samples * multi, 0.0);
        std::streampos byte_pos;
        for (auto&& snp : input)
        {
            std::string chr =
                snp.chr() > m_autosome_ct ? "chrX" : std::to_string(snp.chr());
            gen_ostringstream_for_snp(temp, snp.rs(), snp.rs(), chr, snp.loc(),
                                      snp.ref(), snp.alt(), phased, mock_file,
                                      byte_pos, context);
        }
        mock_w_offset.write(mock_file.str().data(), mock_file.str().size());
        return (mock_w_offset.str());
    }
    void add_file(const std::string& name)
    {
        m_genotype_file_names.clear();
        m_genotype_file_names.push_back(name);
    }
    void gen_ostringstream_for_snp(
        const std::vector<double>& geno_prob, const std::string& snpid,
        const std::string& rsid, const std::string chr, const size_t loc,
        const std::string a1, const std::string a2,
        const genfile::OrderType& phased, std::ostringstream& mock_file,
        std::streampos& byte_pos, genfile::bgen::Context& context)
    {
        std::vector<genfile::byte_t> buffer, buffer2;
        genfile::byte_t* end = genfile::bgen::write_snp_identifying_data(
            &buffer, context, snpid, rsid, chr, static_cast<uint32_t>(loc), 2,
            [&a1, &a2](std::size_t i) {
                if (i == 0) { return a1; }
                return a2;
            });
        mock_file.write(reinterpret_cast<char*>(&buffer[0]), end - &buffer[0]);
        byte_pos = mock_file.tellp();
        int bits_per_probability = 16;
        genfile::bgen::GenotypeDataBlockWriter writer(
            &buffer, &buffer2, context, bits_per_probability);
        writer.initialise(context.number_of_samples, 2);
        gen_prob_vector(geno_prob, context.number_of_samples, phased, writer);
        mock_file.write(reinterpret_cast<char const*>(writer.repr().first),
                        writer.repr().second - writer.repr().first);
    }

    void gen_prob_vector(const std::vector<double>& geno_prob,
                         const size_t number_of_individuals,
                         const genfile::OrderType& type,
                         genfile::bgen::GenotypeDataBlockWriter& writer)
    {
        for (size_t i = 0; i < number_of_individuals; ++i)
        {

            writer.set_sample(i);
            if (type == genfile::ePerUnorderedGenotype)
            {
                writer.set_number_of_entries(2, 3, type, genfile::eProbability);
                for (uint32_t j = 0; j < 3; ++j)
                { writer.set_value(j, geno_prob[i * 3 + j]); }
            }
            else
            {
                writer.set_number_of_entries(
                    2, 4, genfile::ePerPhasedHaplotypePerAllele,
                    genfile::eProbability);
                for (uint32_t j = 0; j < 4; ++j)
                { writer.set_value(j, geno_prob[i * 4 + j]); }
            }
        }
        writer.finalise();
    }
    size_t num_info_filter() const { return m_num_info_filter; }
    size_t num_geno_filter() const { return m_num_geno_filter; }
    size_t num_maf_filter() const { return m_num_maf_filter; }
    size_t num_miss_filter() const { return m_num_miss_filter; }
    void manual_load_snp(SNP cur)
    {
        m_existed_snps_index[cur.rs()] = m_existed_snps.size();
        m_existed_snps.emplace_back(cur);
    }
    std::vector<SNP> existed_snps() const { return m_existed_snps; }
    std::vector<std::string> genotype_file_names() const
    {
        return m_genotype_file_names;
    }
    static void round_probs_to_simplex(double* p, std::size_t n,
                                       std::size_t number_of_bits)
    {
        double const scale =
            uint64_t(0xFFFFFFFFFFFFFFFF) >> (64 - number_of_bits);
        std::multimap<double, double*, std::greater<double>> fractional_parts;
        double total_fractional_part = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            *(p + i) *= scale;
            double const fractional_part = *(p + i) - std::floor(*(p + i));
            fractional_parts.insert(std::make_pair(fractional_part, (p + i)));
            total_fractional_part += fractional_part;
        }
        std::size_t const upper = std::floor(0.5 + total_fractional_part);
        std::multimap<double, double*, std::greater<double>>::const_iterator
            i = fractional_parts.begin(),
            end_i = fractional_parts.end();
        for (std::size_t count = 0; i != end_i; ++i, ++count)
        {
            if (count < upper)
            { *(i->second) = std::ceil(*(i->second)) / scale; }
            else
            {
                *(i->second) = std::floor(*(i->second)) / scale;
            }
        }
    }

    static void generate_samples(
        const size_t n_entries, const size_t n_sample, const QCFiltering qc,
        std::vector<uintptr_t>& plink_genotype, std::vector<double>& data_prob,
        std::vector<bool>& founder, genfile::bgen::Layout version,
        genfile::OrderType phasing, double& exp_mach, double& exp_impute,
        uint32_t& ref_ct, uint32_t& het_ct, uint32_t& alt_ct, uint32_t& miss_ct)
    {
        ref_ct = 0;
        het_ct = 0;
        alt_ct = 0;
        miss_ct = 0;
        exp_mach = 0.0;
        exp_impute = 0.0;
        // just test 16 as bgen library should have test run other options
        const size_t bits_per_probability = 16;
        const uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(n_sample);
        const uintptr_t unfiltered_sample_ctv2 = 2 * unfiltered_sample_ctl;
        plink_genotype.resize(unfiltered_sample_ctv2, 0);
        std::vector<size_t> expected_geno(n_sample, 0);
        std::vector<double> cur_geno_prob(3, 0.0), hap_prob(4, 0.0);
        std::random_device rnd_device;
        std::mt19937 mersenne_engine {
            rnd_device()}; // Generates random integers
        std::uniform_real_distribution<double> dist {0, 1};
        std::uniform_int_distribution<size_t> rand_miss {1, 10};
        // round to 4 decimal places so that BGEN's imprecision would not cause
        // our test to break
        auto prob_generation = [&dist, &mersenne_engine, &version, &phasing]() {
            // return dist(mersenne_engine);
            switch (version)
            {
            case genfile::bgen::e_Layout1:
            {
                double tmp = dist(mersenne_engine);
                std::vector<double> prob = {tmp * tmp, tmp * (1.0 - tmp) * 2.0,
                                            (1.0 - tmp) * (1.0 - tmp)};
                for (auto&& p : prob)
                { p = std::floor(0.5 + p * 32768.0) / 32768.0; }
                return prob;
            }
            case genfile::bgen::e_Layout2:
                switch (phasing)
                {
                case genfile::OrderType::ePerUnorderedGenotype:
                {
                    double tmp = dist(mersenne_engine);
                    std::vector<double> prob = {tmp * tmp,
                                                tmp * (1.0 - tmp) * 2.0,
                                                (1.0 - tmp) * (1.0 - tmp)};
                    round_probs_to_simplex(prob.data(), 3,
                                           bits_per_probability);
                    return prob;
                }
                case genfile::OrderType::ePerPhasedHaplotypePerAllele:
                {
                    double tmp = dist(mersenne_engine);
                    double tmp2 = dist(mersenne_engine);
                    std::vector<double> prob = {tmp, 1.0 - tmp, tmp2,
                                                1.0 - tmp2};
                    round_probs_to_simplex(prob.data(), 2,
                                           bits_per_probability);
                    round_probs_to_simplex(prob.data() + 2, 2,
                                           bits_per_probability);
                    return prob;
                }
                default: return std::vector<double> {}; break;
                }
            default: return std::vector<double> {}; break;
            }
        };
        auto missing_prob = [&rand_miss, &mersenne_engine]() {
            return rand_miss(mersenne_engine) <= 2;
        };
        double sum_prob = 0.0, impute_sum = 0.0;
        misc::RunningStat statistic;
        size_t geno_idx = 0;
        for (size_t i = 0; i < n_sample; ++i)
        {
            double exp = 0.0, tmp = 0.0;
            expected_geno[i] = 3;
            if (missing_prob())
            {
                if (founder[i])
                {
                    ++miss_ct;
                    SET_BIT(geno_idx, plink_genotype.data());
                    geno_idx += 2;
                }
            }
            else
            {
                if (n_entries == 3)
                {
                    cur_geno_prob = prob_generation();
                    data_prob[i * n_entries] = cur_geno_prob[0];
                    data_prob[i * n_entries + 1] = cur_geno_prob[1];
                    data_prob[i * n_entries + 2] = cur_geno_prob[2];
                }
                else
                {
                    hap_prob = prob_generation();
                    data_prob[i * n_entries] = hap_prob[0];
                    data_prob[i * n_entries + 1] = hap_prob[1];
                    data_prob[i * n_entries + 2] = hap_prob[2];
                    data_prob[i * n_entries + 3] = hap_prob[3];
                    cur_geno_prob[0] = hap_prob[0] * hap_prob[2];
                    cur_geno_prob[1] =
                        hap_prob[0] * hap_prob[3] + hap_prob[1] * hap_prob[2];
                    cur_geno_prob[2] = hap_prob[1] * hap_prob[3];
                }
                exp = cur_geno_prob[1] + cur_geno_prob[2] * 2.0;
                tmp = cur_geno_prob[1] + cur_geno_prob[2] * 4.0;
                for (auto&& g : cur_geno_prob) { sum_prob += g; }
                const double prob1 = cur_geno_prob[0] * 2 + cur_geno_prob[1];
                const double prob2 = cur_geno_prob[2] * 2 + cur_geno_prob[1];
                const double hard_score =
                    (std::fabs(prob1 - std::round(prob1))
                     + std::fabs(prob2 - std::round(prob2)))
                    * 0.5;
                if (hard_score <= qc.hard_threshold)
                {
                    double hard_prob = 0;
                    for (size_t geno = 0; geno < 3; ++geno)
                    {
                        if (cur_geno_prob[geno] > hard_prob
                            && cur_geno_prob[geno] >= qc.dose_threshold)
                        {
                            // +1 because geno ==1 represents missing
                            expected_geno[i] = geno;
                            hard_prob = cur_geno_prob[geno];
                        }
                    }
                }

                statistic.push(exp);
                impute_sum -= tmp - exp * exp;
                if (founder[i])
                {
                    switch (expected_geno[i])
                    {
                    case 0: ++ref_ct; break;
                    case 1:
                        SET_BIT(geno_idx + 1, plink_genotype.data());
                        ++het_ct;
                        break;
                    case 2:
                        ++alt_ct;
                        SET_BIT(geno_idx, plink_genotype.data());
                        SET_BIT(geno_idx + 1, plink_genotype.data());
                        break;
                    case 3:
                        ++miss_ct;
                        SET_BIT(geno_idx, plink_genotype.data());
                        break;
                    }
                    geno_idx += 2;
                }
            }
        }
        double p = statistic.mean() / 2.0;
        double p_all = 2.0 * p * (1.0 - p);
        exp_mach = statistic.var() / p_all;
        exp_impute = 1.0 + ((impute_sum / p_all) / sum_prob);
    }
    static void generate_samples(const size_t n_entries, const size_t n_sample,
                                 const QCFiltering qc,
                                 std::vector<uintptr_t>& plink_genotype,
                                 std::vector<double>& data_prob,
                                 std::vector<bool>& founder,
                                 genfile::bgen::Layout version,
                                 genfile::OrderType phasing)
    {
        double exp_mach, exp_impute;
        uint32_t ref_ct, het_ct, alt_ct, miss_ct;
        generate_samples(n_entries, n_sample, qc, plink_genotype, data_prob,
                         founder, version, phasing, exp_mach, exp_impute,
                         ref_ct, het_ct, alt_ct, miss_ct);
    }
};

#endif // MOCK_BINARYGEN_HPP
