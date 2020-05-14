#ifndef BINARYGEN_SETTERS_HPP
#define BINARYGEN_SETTERS_HPP

#include "bgen_lib.hpp"
#include "misc.hpp"
#include "plink_common.hpp"
#include "storage.hpp"
#include <stdexcept>
#include <zlib.h>

// TODO: Use ref MAf for dosage score too
class PRS_Interpreter
{
public:
    virtual ~PRS_Interpreter() {}
    PRS_Interpreter(std::vector<PRS>* sample_prs,
                    std::vector<uintptr_t>* sample_inclusion,
                    MISSING_SCORE missing)
        : m_sample_prs(sample_prs), m_sample_inclusion(sample_inclusion)
    {
        m_miss_count = m_ploidy * (missing != MISSING_SCORE::SET_ZERO);
        m_setzero = (missing == MISSING_SCORE::SET_ZERO);
        m_centre = (missing == MISSING_SCORE::CENTER);
        m_geno_probs.resize(3, 0.0);
    }
    void set_stat(const double& stat, const double& homcom_weight,
                  const double& het_weight, const double& homrar_weight,
                  const bool flipped)
    {
        // don't use external expected as that doesn't take into account of
        // the weighting
        m_missing.clear();
        m_stat = stat;
        m_cal_expected = 0;
        dose_statistic.clear();
        // to match the encoding in PLINK format, we "unflip" SNPs here
        // otherwise our polygenic score will be going to an opposite
        // direction
        m_weights[0] = homcom_weight;
        m_weights[1] = het_weight;
        m_weights[2] = homrar_weight;
        if (!flipped)
        {
            m_weights[0] = homrar_weight;
            m_weights[2] = homcom_weight;
        }
    }

    void initialise(std::size_t, std::size_t)
    {
        // we will reset the IDs
        // m_prs_sample_i represent the sample index of the result PRS
        // vector which does not contain samples that are removed
        m_prs_sample_i = 0;
    }

    void set_min_max_ploidy(uint32_t, uint32_t, uint32_t, uint32_t) {}

    bool set_sample(std::size_t i)
    {
        m_is_missing = false;
        return IS_SET(m_sample_inclusion->data(), i);
    }
    void set_number_of_entries(std::size_t ploidy, std::size_t,
                               genfile::OrderType phased, genfile::ValueType)
    {
        m_ploidy = ploidy;
        m_phased = (phased == genfile::OrderType::ePerPhasedHaplotypePerAllele);
        m_probs.resize(3 + m_phased, 0.0);
    }

    void set_value(uint32_t idx, double value) { m_probs[idx] = value; }

    void set_value(uint32_t, genfile::MissingValue) { m_is_missing = true; }

    void sample_completed()
    {
        if (!m_is_missing)
        {
            double v11_miss_check = 0.0;
            for (auto&& p : m_probs) { v11_miss_check += p; }
            if (misc::logically_equal(v11_miss_check, 0.0))
            { m_is_missing = true; }
        }
        if (m_phased)
        {
            m_geno_probs[0] = m_probs[0] * m_probs[2];
            m_geno_probs[1] = m_probs[0] * m_probs[3] + m_probs[1] * m_probs[2];
            m_geno_probs[2] = m_probs[1] * m_probs[3];
        }
        else
        {
            m_geno_probs = m_probs;
        }
        m_sum = 0.0;
        for (size_t i = 0; i < 3; ++i)
        { m_sum += m_geno_probs[i] * m_weights[i]; }
        add_prs_score(m_prs_sample_i);
        // go to next sample that we need (not the bgen index)
        ++m_prs_sample_i;
    }
    virtual void add_prs_score(size_t) {}
    void finalise()
    {
        m_adj_score = 0;
        m_miss_score = 0;
        m_miss_count = 0;
        if (m_centre) { m_adj_score = m_stat * dose_statistic.mean(); }
        if (!m_setzero)
        {
            m_miss_count = m_ploidy;
            m_miss_score = m_stat * dose_statistic.mean();
        }
        if (!m_centre) { process_missing(); }
        else
        {
            process_centre_missing();
        }
    }
    virtual void process_missing() {}
    virtual void process_centre_missing() {}

protected:
    std::vector<PRS>* m_sample_prs;
    std::vector<uintptr_t>* m_sample_inclusion;
    std::vector<size_t> m_missing;
    std::vector<double> m_probs;
    std::vector<double> m_geno_probs;
    std::vector<double> m_weights = {0, 0.5, 1};
    misc::RunningStat dose_statistic;
    double m_stat = 0.0;
    double m_sum = 0.0;
    double m_miss_score = 0.0;
    double m_adj_score = 0.0;
    double m_cal_expected = 0.0;
    uint32_t m_prs_sample_i = 0;
    size_t m_ploidy = 2;
    size_t m_miss_count = 0;
    bool m_is_missing = false;
    bool m_phased = false;
    bool m_setzero = false;
    bool m_centre = false;
};
class First_PRS : public PRS_Interpreter
{
public:
    First_PRS(std::vector<PRS>* sample_prs,
              std::vector<uintptr_t>* sample_inclusion, MISSING_SCORE missing)
        : PRS_Interpreter(sample_prs, sample_inclusion, missing)
    {
    }
    virtual ~First_PRS() {}
    void add_prs_score(size_t idx)
    {
        if (m_is_missing) { m_missing.push_back(m_prs_sample_i); }
        // this is not a missing sample and we can either add the prs or
        // assign the PRS
        else
        {

            (*m_sample_prs)[idx].num_snp = m_ploidy;
            (*m_sample_prs)[idx].prs = m_sum * m_stat;
            dose_statistic.push(m_sum);
        }
    }
    void process_centre_missing()
    {
        size_t cur_idx = 0;
        for (size_t i = 0; i < m_sample_prs->size(); ++i)
        {
            if (cur_idx < m_missing.size() && i == m_missing[cur_idx])
            {
                (*m_sample_prs)[i].prs = m_miss_score;
                (*m_sample_prs)[i].num_snp = m_miss_count;
                ++cur_idx;
            }
            else if (m_centre)
            {
                // if it is not missing and we want the centre the
                // score we will need to minus the adjusted score
                // which was 0 before this run
                (*m_sample_prs)[i].prs -= m_adj_score;
            }
        }
    }
    void process_missing()
    {
        // only need to do this if we don't have the expected
        // information
        for (auto&& idx : m_missing)
        {
            (*m_sample_prs)[idx].prs = m_miss_score;
            (*m_sample_prs)[idx].num_snp = m_miss_count;
        }
    }
};
class Add_PRS : public PRS_Interpreter
{
public:
    Add_PRS(std::vector<PRS>* sample_prs,
            std::vector<uintptr_t>* sample_inclusion, MISSING_SCORE missing)
        : PRS_Interpreter(sample_prs, sample_inclusion, missing)
    {
    }
    virtual ~Add_PRS() {}
    void add_prs_score(size_t idx)
    {
        if (m_is_missing) { m_missing.push_back(m_prs_sample_i); }
        // this is not a missing sample and we can either add the prs or
        // assign the PRS
        else
        {

            (*m_sample_prs)[idx].num_snp += m_ploidy;
            (*m_sample_prs)[idx].prs += m_sum * m_stat;
            dose_statistic.push(m_sum);
        }
    }
    void process_centre_missing()
    {
        size_t cur_idx = 0;
        for (size_t i = 0; i < m_sample_prs->size(); ++i)
        {
            if (cur_idx < m_missing.size() && i == m_missing[cur_idx])
            {
                (*m_sample_prs)[i].prs += m_miss_score;
                (*m_sample_prs)[i].num_snp += m_miss_count;
                ++cur_idx;
            }
            else if (m_centre)
            {
                // if it is not missing and we want the centre the
                // score we will need to minus the adjusted score
                // which was 0 before this run
                (*m_sample_prs)[i].prs -= m_adj_score;
            }
        }
    }
    void process_missing()
    {
        // only need to do this if we don't have the expected
        // information
        for (auto&& idx : m_missing)
        {
            (*m_sample_prs)[idx].prs += m_miss_score;
            (*m_sample_prs)[idx].num_snp += m_miss_count;
        }
    }
};


struct PLINK_generator
{
    PLINK_generator(uintptr_t* sample, uintptr_t* genotype,
                    double hard_threshold, double dose_threshold)
        : m_sample(sample)
        , m_genotype(genotype)
        , m_hard_threshold(hard_threshold)
        , m_dose_threshold(dose_threshold)
    {
        m_prob.resize(3);
        m_geno_prob.resize(3);
    }
    void initialise(std::size_t, std::size_t)
    {
        // initialize is called when we meet a new SNP
        // we therefore want to reset the index and shift, use for
        // indicate where we push the byte onto genotype
        m_index = 0;
        m_shift = 0;
        m_homcom_ct = 0;
        m_homrar_ct = 0;
        m_het_ct = 0;
        m_missing_ct = 0;
        m_impute2 = 0;
        m_total_probability = 0.0;
        // we also clean the running stat, so that we can go
        // through another round of calculation of mean and sd
        statistic.clear();
    }
    void set_min_max_ploidy(uint32_t, uint32_t, uint32_t, uint32_t) {}

    bool set_sample(std::size_t i)
    {
        // we set the sample index to i
        m_sample_i = i;
        // we also reset the hard_prob to 0
        m_hard_prob = 0.0;
        // and expected value to 0
        m_missing = false;
        // then we determine if we want to include sample using by
        // consulting the flag on m_sample
        return IS_SET(m_sample, m_sample_i);
    }
    void set_number_of_entries(std::size_t, std::size_t,
                               genfile::OrderType phased, genfile::ValueType)
    {
        m_phased = (phased == genfile::OrderType::ePerPhasedHaplotypePerAllele);
        m_prob.resize(3 + m_phased, 0.0);
    }

    void set_value(uint32_t idx, double value)
    {
        // if the current probability is the highest and the value is
        // higher than the required threshold, we will assign it
        //            if (value > m_hard_prob && value >=
        //            m_hard_threshold)
        //            {
        //                /*
        //                 * Representation of each geno to their binary
        //                 code:
        //                 *   geno    desired binary  decimal
        //                 representation
        //                 *   0           00              0
        //                 *   1           01              1
        //                 *   2           11              3
        //                 *   the binary code 10 is reserved for
        //                 missing sample
        //                 */
        //                // plink do 3 2 0 1
        //                m_geno = (geno == 2) ? 3 : geno;
        //                m_hard_prob = value;
        //            }

        m_prob[idx] = value;
    }

    void set_value(uint32_t, genfile::MissingValue) { m_missing = true; }
    void finalise() {}
    /*!
     * \brief sample_completed is called when each sample is completed.
     * This will set the binary vector accordingly
     */
    void sample_completed()
    {
        // if m_shift is zero, it is when we meet the index for the
        // first time, therefore we want to initialize it
        if (m_shift == 0) m_genotype[m_index] = 0;
        // check missing not required for phased data, as phased data
        // only supported in 1.2+, which represent missing value
        // differently
        if (!m_phased && !m_missing)
        {
            m_missing =
                misc::logically_equal(m_prob[0] + m_prob[1] + m_prob[2], 0.0);
            m_geno_prob = m_prob;
        }
        else
        {
            // convert haplotype prob into genotype prob
            m_geno_prob[0] = m_prob[0] * m_prob[2];
            m_geno_prob[1] = m_prob[0] * m_prob[3] + m_prob[1] * m_prob[2];
            m_geno_prob[2] = m_prob[1] * m_prob[3];
        }

        const double prob1 = m_geno_prob[0] * 2 + m_geno_prob[1];
        const double prob2 = m_geno_prob[2] * 2 + m_geno_prob[1];
        uintptr_t obs_genotype = 1;
        const double hard_score = (std::fabs(prob1 - std::round(prob1))
                                   + std::fabs(prob2 - std::round(prob2)))
                                  * 0.5;
        if (!m_missing && (hard_score <= m_hard_threshold))
        {
            m_hard_prob = 0;
            for (size_t geno = 0; geno < 3; ++geno)
            {
                if (m_geno_prob[geno] > m_hard_prob
                    && m_geno_prob[geno] >= m_dose_threshold)
                {
                    // +1 because geno ==1 represents missing
                    obs_genotype = (geno == 0) ? geno : geno + 1;
                    m_hard_prob = m_geno_prob[geno];
                }
            }
        }
        // we now add the genotype to the vector
        m_genotype[m_index] |= obs_genotype << m_shift;
        switch (obs_genotype)
        {
        case 3: ++m_homrar_ct; break;
        case 2: ++m_het_ct; break;
        case 1: ++m_missing_ct; break;
        case 0: ++m_homcom_ct; break;
        }
        // as the genotype is represented by two bit, we will +=2
        m_shift += 2;
        // if we reach the boundary, we will now add the index and reset
        // the shift
        if (m_shift == BITCT)
        {
            ++m_index;
            m_shift = 0;
        }

        double exp_value = m_geno_prob[1] + m_geno_prob[2] * 2.0;
        double impute2_tmp = m_geno_prob[1] + m_geno_prob[2] * 4.0;
        // we can now push in the expected value for this sample. This
        // can then use for the calculation of the info score
        // only calculate for included samples
        if (IS_SET(m_sample, m_sample_i) && !m_missing)
        {
            statistic.push(exp_value);
            m_impute2 -= impute2_tmp - exp_value * exp_value;
            for (auto&& g : m_geno_prob) m_total_probability += g;
        }
    }
    double info_score(const INFO type) const
    {
        double p = statistic.mean() / 2.0;
        double p_all = 2.0 * p * (1.0 - p);
        if (type == INFO::MACH) { return (statistic.var() / p_all); }
        else if (type == INFO::IMPUTE2)
        {
            return 1.0 + ((m_impute2 / p_all) / m_total_probability);
        }
        throw std::runtime_error("Error: Undefined info score type");
    }
    double expected() const { return statistic.mean(); }
    void get_count(size_t& homcom_ct, size_t& het_ct, size_t& homrar_ct,
                   size_t& missing_ct) const
    {
        homcom_ct = m_homcom_ct;
        het_ct = m_het_ct;
        homrar_ct = m_homrar_ct;
        missing_ct = m_missing_ct;
    }
    void get_count(uint32_t& homcom_ct, uint32_t& het_ct, uint32_t& homrar_ct,
                   uint32_t& missing_ct) const
    {
        homcom_ct = m_homcom_ct;
        het_ct = m_het_ct;
        homrar_ct = m_homrar_ct;
        missing_ct = m_missing_ct;
    }
    ~PLINK_generator() {}

private:
    // the sample inclusion vector, if bit is set, sample is required
    std::vector<double> m_prob;
    std::vector<double> m_geno_prob;
    // is the genotype vector
    uintptr_t* m_sample;
    uintptr_t* m_genotype;
    misc::RunningStat statistic;
    double m_hard_threshold = 0.0;
    double m_dose_threshold = 0.0;
    double m_hard_prob = 0;
    double m_impute2 = 0.0;
    double m_total_probability = 0.0;
    uint32_t m_shift = 0;
    uint32_t m_index = 0;
    size_t m_sample_i = 0;
    uint32_t m_homcom_ct = 0;
    uint32_t m_homrar_ct = 0;
    uint32_t m_het_ct = 0;
    uint32_t m_missing_ct = 0;
    bool m_phased = false;
    bool m_missing = false;
};
#endif // BINARYGEN_SETTERS_HPP
