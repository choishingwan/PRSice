#include "bgen.hpp"

BGEN::BGEN(std::string prefix,  std::string remove_sample, std::string keep_sample,
        bool ignore_fid, int num_auto=22, bool no_x=false, bool no_y=false, bool no_xy=false,
        bool no_mt=false, const size_t thread=1, bool verbose=false)
{
    if(remove_sample.empty()) m_remove_sample = false;
    else m_remove_sample_list = load_ref(remove_sample, ignore_fid);
    if(keep_sample.empty()) m_keep_sample = false;
    else m_keep_sample_list = load_ref(keep_sample, ignore_fid);

    m_xymt_codes.resize(XYMT_OFFSET_CT);
    init_chr(num_auto, no_x, no_y, no_xy, no_mt);
    m_thread = thread;
    set_genotype_files(prefix);
    m_sample_names = load_samples(ignore_fid);
    m_existed_snps = load_snps();
    if(verbose)
    {
        fprintf(stderr, "%zu people (%zu males, %zu females) included\n", m_unfiltered_sample_ct, m_num_male, m_num_female);
        if(m_num_ambig!=0) fprintf(stderr, "%u ambiguous variants excluded\n", m_num_ambig);
        fprintf(stderr, "%zu variants included\n", m_marker_ct);
    }

    m_founder_ctl = BITCT_TO_WORDCT(m_founder_ct);
    m_founder_ctv3 = BITCT_TO_ALIGNED_WORDCT(m_founder_ct);
    m_founder_ctsplit = 3 * m_founder_ctv3;
}
