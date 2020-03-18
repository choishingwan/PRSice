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

#include "prsice.hpp"

void PRSice::produce_null_prs(
    Thread_Queue<std::pair<std::vector<double>, size_t>>& q, Genotype& target,
    std::vector<size_t> background, size_t num_consumer,
    std::map<size_t, std::vector<size_t>>& set_index)
{
    // we need to know the size of the biggest set
    const size_t max_size = set_index.rbegin()->first;
    const size_t num_sample = m_matrix_index.size();
    const size_t num_regress_sample =
        static_cast<size_t>(m_independent_variables.rows());
    size_t processed = 0;
    size_t prev_size = 0;
    // we seed the random number generator
    std::mt19937 g(m_perm_info.seed);
    bool first_run = true;
    while (processed < m_perm_info.num_permutation)
    {
        // sample without replacement
        fisher_yates(background, g, max_size);
        first_run = true;
        prev_size = 0;
        for (auto&& set_size : set_index)
        {
            // for each gene sets size, we calculate the PRS
            target.get_null_score(set_size.first, prev_size, background,
                                  first_run);
            first_run = false;
            // we need to know how many SNPs we have already read, such that
            // we can skip reading this number of SNPs for the next set
            prev_size = set_size.first;
            // we store the PRS in a new vector to avoid crazy error with
            // move semetics and stuff which I have not fully understand
            std::vector<double> prs(num_regress_sample, 0);
            for (size_t sample_id = 0; sample_id < num_sample; ++sample_id)
            {
                // propagate the prs vector
                prs[sample_id] =
                    target.calculate_score(m_matrix_index[sample_id]);
            }
            // then we push the result prs to the queue, which can then
            // picked up by the consumers
            q.emplace(std::make_pair(prs, set_size.first), num_consumer);
            ++m_analysis_done;
            print_progress();
        }
        ++processed;
    }
    // send termination signal to the consumers
    q.completed();
}


void PRSice::consume_prs(
    Thread_Queue<std::pair<std::vector<double>, size_t>>& q,
    const Regress& decomposed, std::map<size_t, std::vector<size_t>>& set_index,
    const std::vector<double>& obs_t_value,
    std::vector<std::atomic<size_t>>& set_perm_res, const bool is_binary)
{
    const Eigen::Index num_regress_sample =
        static_cast<Eigen::Index>(m_matrix_index.size());
    Eigen::MatrixXd independent;
    Eigen::VectorXd prs, beta, effects;
    if (m_perm_info.logit_perm && is_binary)
        independent = m_independent_variables;
    // to avoid false sharing and frequent lock, we wil first store all
    // permutation results within a temporary vector
    double coefficient, standard_error, r2;
    double obs_p = 2.0; // for safety reason, make sure it is out bound
    // results from queue will be stored in the prs_info
    std::pair<std::vector<double>, size_t> prs_info;
    // now listen for producer
    while (!q.pop(prs_info))
    {
        // update the independent variable matrix with the new PRS

        if (is_binary && m_perm_info.logit_perm)
        {
            independent.col(1) = Eigen::Map<Eigen::VectorXd>(
                std::get<0>(prs_info).data(), num_regress_sample);
            Regression::glm(m_phenotype, independent, obs_p, r2, coefficient,
                            standard_error, 1);
        }
        else
        {
            prs = Eigen::Map<Eigen::VectorXd>(std::get<0>(prs_info).data(),
                                              num_regress_sample);
            get_t_value(decomposed, prs, beta, effects, coefficient,
                        standard_error);
        }
        double t_value = std::fabs(coefficient / standard_error);
        auto&& index = set_index[std::get<1>(prs_info)];
        // we register the number of time a more significant / bigger
        // t-value is obtained when compared to the observed t-value
        for (auto&& ref : index)
        {
            // in theory because set_perm_res is now atomic, it should be ok
            if (obs_t_value[ref] < t_value) ++set_perm_res[ref];
        }
    }
}

void PRSice::observe_set_perm(Thread_Queue<size_t>& progress_observer,
                              size_t num_thread)
{
    size_t progress = 0, total_run = 0;
    while (!progress_observer.pop(progress, num_thread))
    {
        m_analysis_done += progress;
        total_run += progress;
        print_progress();
    }
}

// Shuffle the idx vector
// By selecting the first n element from idx, we've got the random selection
// without replacement
void PRSice::fisher_yates(std::vector<size_t>& idx, std::mt19937& g, size_t n)
{
    size_t begin = 0;
    // we will shuffle n where n is the set with the largest size
    // this is the Fisher-Yates shuffle algorithm for random selection
    // without replacement
    size_t num_idx = idx.size() - 1;
    size_t advance_index;
    while (n--)
    {
        std::uniform_int_distribution<size_t> dist(begin, num_idx);
        advance_index = dist(g);
        std::swap<size_t>(idx[begin], idx[advance_index]);
        ++begin;
    }
}

template <typename T>
void PRSice::subject_set_perm(T& progress_observer, Genotype& target,
                              std::vector<size_t> background,
                              std::map<size_t, std::vector<size_t>>& set_index,
                              std::vector<std::atomic<size_t>>& set_perm_res,
                              const std::vector<double>& obs_t_value,
                              const std::random_device::result_type seed,
                              const Regress& decomposed, const size_t num_perm,
                              const bool is_binary)
{
    assert(set_index.size() != 0);
    const size_t max_size = set_index.rbegin()->first;
    const Eigen::Index num_sample =
        static_cast<Eigen::Index>(m_matrix_index.size());
    double coefficient, standard_error, r2, obs_p, t_value;
    Eigen::VectorXd prs = Eigen::VectorXd::Zero(num_sample);
    Eigen::VectorXd beta, effects;
    Eigen::MatrixXd independent;
    if (m_perm_info.logit_perm && is_binary)
    { independent = m_independent_variables; }
    // each thread should have their own cur_prs to avoid overhead
    std::vector<PRS> cur_prs(target.num_sample());
    bool first_run = true;
    std::mt19937 g(seed);
    size_t processed = 0;

    while (processed < num_perm)
    {
        fisher_yates(background, g, max_size);
        //  we have now selected N SNPs from the background. We can then
        //  construct the PRS based on these index
        first_run = true;
        size_t prev_size = 0;
        for (auto&& set_size : set_index)
        {
            target.get_null_score(cur_prs, set_size.first, prev_size,
                                  background, first_run);
            first_run = false;
            prev_size = set_size.first;
            for (Eigen::Index sample_id = 0; sample_id < num_sample;
                 ++sample_id)
            {
                size_t idx = m_matrix_index[static_cast<size_t>(sample_id)];
                if (m_perm_info.logit_perm && is_binary)
                {
                    independent(sample_id, 1) =
                        target.calculate_score(cur_prs, idx);
                }
                else
                {
                    prs(sample_id) = target.calculate_score(cur_prs, idx);
                }
            }
            progress_observer.emplace(1);
            //  we can now perform the glm or linear regression analysis
            if (is_binary && m_perm_info.logit_perm)
            {
                Regression::glm(m_phenotype, m_independent_variables, obs_p, r2,
                                coefficient, standard_error, 1);
            }
            else
            {
                get_t_value(decomposed, prs, beta, effects, coefficient,
                            standard_error);
            }
            t_value = std::fabs(coefficient / standard_error);
            // set_size second contain the indexs to each set with this size
            for (auto&& set_index : set_size.second)
            { set_perm_res[set_index] += (obs_t_value[set_index] < t_value); }
        }
        ++processed;
    }
    progress_observer.completed();
}

void PRSice::run_competitive(
    Genotype& target, const std::vector<size_t>::const_iterator& bk_start_idx,
    const std::vector<size_t>::const_iterator& bk_end_idx,
    const size_t pheno_index)
{
    if (!m_perm_info.run_set_perm) { return; }
    fprintf(stderr, "\n");
    m_reporter->report("\n\nStart competitive permutation\n");
    const bool is_binary = m_pheno_info.binary[pheno_index];
    const size_t num_prs_res = m_prs_summary.size();
    const size_t num_bk_snps =
        static_cast<size_t>(std::distance(bk_start_idx, bk_end_idx));
    // the number of items to skip from the front of prs_summary
    size_t pheno_start_idx = 0;
    // obs_t_value stores the observed t-value
    std::vector<double> obs_t_value;
    // set_index stores the index of sets with "key" size
    std::map<size_t, std::vector<size_t>> set_index;
    bool started = false;
    // start at 1 to avoid the base set
    size_t cur_set_index = 0;
    size_t max_set_size = 0;
    // pre-decompose the Y+ cov matrix to speed up the permutation
    if (is_binary && !m_printed_warning)
    {
        if (!m_perm_info.logit_perm)
        {
            m_reporter->report(
                "Warning: To speed up the permutation, "
                "linear regression instead of logistic "
                "regression were performed within the permutation, and "
                "constructs the null distribution using the absolute "
                "z-scores. This is based on the assumption that "
                "linear Regression & logistic "
                "regression should produce similar absolute z-scores. In "
                "addition, the regression equation changed from "
                "Phenotype~PRS+Covariates to PRS~Phenotype+Covariate. This "
                "two "
                "equations should generate similar z-score for the "
                "independent "
                "variable and will allow us to perform some optimizations "
                "to "
                "speed up the permutation\n\n");
        }
        else
        {
            m_reporter->report("Warning: Using --logit-perm will be "
                               "ridiculously slow\n");
        }
    }
    const Eigen::Index p = m_independent_variables.cols();
    Regress decomposed;
    if (!m_perm_info.logit_perm)
    {
        decomposed.YCov = m_independent_variables;
        decomposed.YCov.col(1) = m_phenotype;
        // TODO The phenotype vector here doesn't seems to represent the
        // phenotype in the reference file?
        decomposed.PQR.compute(decomposed.YCov);
        decomposed.Pmat = decomposed.PQR.colsPermutation();
        decomposed.rank = decomposed.PQR.rank();
        if (decomposed.rank != p)
        {
            decomposed.Rinv =
                decomposed.PQR.matrixQR()
                    .topLeftCorner(decomposed.rank, decomposed.rank)
                    .triangularView<Eigen::Upper>()
                    .solve(Eigen::MatrixXd::Identity(decomposed.rank,
                                                     decomposed.rank));
        }
        get_se_matrix(decomposed.PQR, decomposed.Pmat, decomposed.Rinv, p,
                      decomposed.rank, decomposed.se);
    }
    m_printed_warning = true;
    for (size_t i = 0; i < num_prs_res; ++i)
    {
        if (m_prs_summary[i].has_competitive || m_prs_summary[i].set == "Base")
            continue;
        if (!started)
        {
            pheno_start_idx = i;
            started = true;
        }
        auto&& res = m_prs_summary[i].result;
        set_index[res.num_snp].push_back(cur_set_index);
        ++cur_set_index;
        if (res.num_snp > max_set_size) max_set_size = res.num_snp;
        // ori_t_value will contain the obesrved t-value
        obs_t_value.push_back(std::fabs(res.coefficient / res.se));
    }
    // set_perm_res stores number of perm where a more sig result is
    // obtained initialize here as atomic doesn't have copy and move
    // constructor, making any operator that need to invoke constructor
    // invalid
    std::vector<std::atomic<size_t>> set_perm_res(obs_t_value.size());
    for (auto& set : set_perm_res) { set = 0; }
    if (max_set_size > num_bk_snps)
    {
        for (size_t i = pheno_start_idx; i < num_prs_res; ++i)
        {
            // set them to true so that we will skip them for the next
            // phenotype (though in reality, they will all encounter the
            // same error. Need a better structure here)
            m_prs_summary[i].has_competitive = true;
        }
        m_reporter->report("Error: Insufficient background SNPs for "
                           "competitive analysis. Please ensure you have "
                           "use the correct background. Will now generate skip "
                           "the competitive analysis\n");
        return;
    }
    // now we can run the competitive testing
    // know how many thread we are allowed to use
    int num_thread = m_prs_info.thread;
    const uintptr_t num_regress_sample =
        static_cast<uintptr_t>(m_independent_variables.rows());

    // This is a rough estimate, we might be using more memory than
    // indicated here
    const uintptr_t basic_memory_required_per_thread =
        m_perm_info.logit_perm ? (
            4 * num_regress_sample + 2ULL * static_cast<unsigned long long>(p)
            + 1ULL + num_regress_sample * static_cast<unsigned long long>(p))
                               : num_regress_sample;

    for (; num_thread > 0; --num_thread)
    {
        double* bigstack = nullptr;
        try
        {
            bigstack = new double[basic_memory_required_per_thread
                                  * static_cast<size_t>(num_thread)];
            delete bigstack;
            break;
        }
        catch (const std::bad_alloc&)
        {
        }
    }
    if (num_thread == 0)
    {
        fprintf(stderr, "\n");
        throw std::runtime_error(
            "(DEBUG) Error: Not enough memory left for permutation. "
            "Minimum require memory = "
            + std::to_string(basic_memory_required_per_thread / 1048576)
            + " Mb");
    }
    m_reporter->report("Running permutation with " + misc::to_string(num_thread)
                       + " threads");
    size_t ran_perm = 0;
    if (num_thread > 1)
    {
        //  similar to permutation for empirical p-value calculation, we
        //  employ the producer consumer pattern where one thread is
        //  responsible for reading in the PRS and construct the required
        //  independent variable and other threads are responsible for the
        //  calculation
        if (!target.genotyped_stored())
        {
            ran_perm = m_perm_info.num_permutation;
            Thread_Queue<std::pair<std::vector<double>, size_t>> set_perm_queue;
            std::thread producer(&PRSice::produce_null_prs, this,
                                 std::ref(set_perm_queue), std::ref(target),
                                 std::vector<size_t>(bk_start_idx, bk_end_idx),
                                 num_thread - 1, std::ref(set_index));
            std::vector<std::thread> consumer_store;
            for (int i_thread = 0; i_thread < num_thread - 1; ++i_thread)
            {
                consumer_store.push_back(std::thread(
                    &PRSice::consume_prs, this, std::ref(set_perm_queue),
                    std::cref(decomposed), std::ref(set_index),
                    std::cref(obs_t_value), std::ref(set_perm_res), is_binary));
            }

            producer.join();
            for (auto&& thread : consumer_store) thread.join();
        }
        else
        {
            // we don't need gatherer function, we can just let all the threads
            // run subset of the permutation. This should be much faster
            Thread_Queue<size_t> progress_observer;
            std::thread observer(&PRSice::observe_set_perm, this,
                                 std::ref(progress_observer), num_thread);
            std::vector<std::thread> subjects;
            std::mt19937 rand_gen {m_perm_info.seed};
            std::uniform_int_distribution<unsigned int> dis(
                std::numeric_limits<unsigned int>::min(),
                std::numeric_limits<unsigned int>::max());
            size_t job_per_thread =
                m_perm_info.num_permutation / static_cast<size_t>(num_thread);
            int remain = static_cast<int>(
                static_cast<size_t>(m_perm_info.num_permutation)
                % static_cast<size_t>(num_thread));
            for (int i_thread = 0; i_thread < num_thread; ++i_thread)
            {
                auto seed = dis(rand_gen);
                subjects.push_back(std::thread(
                    &PRSice::subject_set_perm<Thread_Queue<size_t>>, this,
                    std::ref(progress_observer), std::ref(target),
                    std::vector<size_t>(bk_start_idx, bk_end_idx),
                    std::ref(set_index), std::ref(set_perm_res),
                    std::cref(obs_t_value), seed, std::cref(decomposed),
                    job_per_thread + (remain > 0), is_binary));
                ran_perm += job_per_thread + (remain > 0);
                remain--;
            }
            observer.join();
            for (auto&& thread : subjects) thread.join();
        }
    }
    else
    {
        // alternatively, if we only got one thread, we will use the no
        // thread function to reduce threading overhead
        dummy_reporter<size_t> dummy(*this);
        subject_set_perm(dummy, target,
                         std::vector<size_t>(bk_start_idx, bk_end_idx),
                         set_index, set_perm_res, obs_t_value, m_perm_info.seed,
                         decomposed, m_perm_info.num_permutation, is_binary);
        ran_perm = m_perm_info.num_permutation;
    }
    // start_index is the index of m_prs_summary[i], not the actual index
    // on set_perm_res.
    // this will iterate all sets from beginning of current phenotype
    // to the last set within the current phenotype
    // because of the sequence of how set_perm_res is contructed,
    // the results for each set should be sequentially presented in
    // set_perm_res. Index for set_perm_res results are therefore
    // i - start_index
    for (size_t i = pheno_start_idx; i < num_prs_res; ++i)
    {
        auto&& res = m_prs_summary[i].result;
        // we need to minus out the start index from i such that our index
        // start at 0, which is the assumption of set_perm_res
        res.competitive_p =
            (static_cast<double>(set_perm_res[(i - pheno_start_idx)]) + 1.0)
            / (static_cast<double>(ran_perm) + 1.0);
        m_prs_summary[i].has_competitive = true;
    }
}
