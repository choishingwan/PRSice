/*
 * dosagetemplate.cpp
 *
 *  Created on: 16 May 2017
 *      Author: shingwan
 */

#include "dosagetemplate.h"

Dosage_template::Dosage_template(std::string prefix, std::string pheno_file, bool header,
        std::string remove_sample, std::string keep_sample,
        bool ignore_fid, int num_auto = 22, bool no_x = false,
        bool no_y = false, bool no_xy = false, bool no_mt = false,
        const size_t thread = 1, bool verbose = false)
{
    // load the remove and keep sample list here. The result will be stored in the
    // corresponding unordered_set
    if(remove_sample.empty()) m_remove_sample = false;
    else m_remove_sample_list = load_ref(remove_sample, ignore_fid);
    if(keep_sample.empty()) m_keep_sample = false;
    else m_keep_sample_list = load_ref(keep_sample, ignore_fid);

    // This was intended for compatibility with species other than human
    // but currently this does not really work as we have not
    // implemented the required plink functions to handle the haplo
    // problems
    m_xymt_codes.resize(XYMT_OFFSET_CT);
    init_chr(num_auto, no_x, no_y, no_xy, no_mt);
    // Set the number of thread use (mainly for clumping)
    m_thread = thread;
    // For multiple chromosome, this will generate a vector containing file names
    // with # replaced by the chromosome numbers. At the moment, that'd be 1-22
    set_genotype_files(prefix);
    // You should implement the load_samples function, which will setup the
    // m_sample_names vector which contains the Sample struct:
    /*
      struct Sample{
        std::string FID;
        std::string IID;
        std::string pheno;
        double prs; // set as 0 when initialize
        int num_snp; // set as 0 when initialize
        bool included; // if this is true, it will be use for the regression
                       // otherwise, this will be excluded
    };
     */
    m_sample_names = load_samples(ignore_fid);
    // This should set up the SNP vector (by returning it)
    // Construct the SNP using
    // SNP(rs_id, chr_code, loc, A1, A2, prefix, num_line));
    // rs_id is the SNP id for matching with reference and base file
    // chr_code is an int corresponding to the chromosome id
    //          if chr is in chr# format or CHR# format, you can
    //          use get_chrom_code_raw to get the int
    // loc is the location of the SNP
    // A1 is the effective allele
    // A2 is the alternative
    // prefix is the genotype file name (mainly for multi-chromosome input)
    // num_line is the SNP id from the file, e.g. 0 for the first SNP
    //          this will be use for jumping to the specific SNP
    m_existed_snps = load_snps();
    if(verbose)
    {
        fprintf(stderr, "%zu people (%zu males, %zu females) included\n", m_unfiltered_sample_ct, m_num_male, m_num_female);
        if(m_num_ambig!=0) fprintf(stderr, "%u ambiguous variants excluded\n", m_num_ambig);
        fprintf(stderr, "%zu variants included\n", m_marker_ct);
    }
    // You should always include this to initialize the variables correctly
    initialize();
    // You should also initialize this here
    // this is a temporary object required for Plink processing. By pre-initializing
    // it here, this should help to speed things up
    m_tmp_genotype = new uintptr_t[m_unfiltered_sample_ctl*2];

}

Dosage_template::~Dosage_template()
{
    // Don't need to delete m_tmp_genotype here as the base class should do it
}

