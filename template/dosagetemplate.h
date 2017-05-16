/*
 * dosagetemplate.h
 *
 *  Created on: 16 May 2017
 *      Author: shingwan
 */

#ifndef DOSAGETEMPLATE_H_
#define DOSAGETEMPLATE_H_

#include <genotype.hpp>

class Dosage_template: public Genotype
{
    public:
        Dosage_template(std::string prefix, std::string pheno_file, bool header,
                std::string remove_sample, std::string keep_sample,
                bool ignore_fid, int num_auto = 22, bool no_x = false,
                bool no_y = false, bool no_xy = false, bool no_mt = false,
                const size_t thread = 1, bool verbose = false);
        virtual ~Dosage_template();
    private:
        void cleanup(){};
        std::vector<Sample> load_samples(bool ignore_fid)
        {
            // You need to initialize the followings
            // Most importantly, you need to setup the m_unfiltered_sample_ct
            // that is required for proper initialization of the binary arrays

            /*
            m_unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
            m_unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;

            m_sex_male = new uintptr_t[m_unfiltered_sample_ctl];
            std::memset(m_sex_male, 0x0, m_unfiltered_sample_ctl*sizeof(uintptr_t));

            m_founder_info = new uintptr_t[m_unfiltered_sample_ctl];
            std::memset(m_founder_info, 0x0, m_unfiltered_sample_ctl*sizeof(uintptr_t));

            m_sample_exclude = new uintptr_t[m_unfiltered_sample_ctl];
            std::memset(m_sample_exclude, 0x0, m_unfiltered_sample_ctl*sizeof(uintptr_t));
            */
            return std::vector < Sample > (0);
        };
        std::vector<SNP> load_snps()
        {
            // You need to initialize the following:
            // m_chr_order
            // And also properly setup the SNP vectors
            return std::vector < SNP > (0);
        };

        inline void load_raw(uintptr_t* genotype, const uint32_t snp_index,
                       const std::string &file_name)
        {
            // this is the most difficult function for you to write
            // basically, this is converting the dosage file into
            // the plink binary format
            // A semi pseudo code is provided below

            // First the inputs:
            // genotype  -> Where the resulting binary array is store at
            //           -> This is initialized with the correct size
            // snp_index -> The SNP id on the file. this is important
            //           -> for quickly jump to the corresponding SNPs
            //           -> either using seekg from ifstream or fseeko
            //           -> from FILE
            // file_name -> The file where the SNP is located on
            // NOTE: You can define a member variable to store the
            //       file stream so that you don't need to reopen it
            //       for every single SNP. You just need to also retain
            //       information of the current file e.g. m_cur_file
            //       that stores which file is open and compare that
            //       with file_name. Only when it is different do you
            //       need to re-open the file.
            //       I recommend using the ifstream as that's what it
            //       was used here. You can also refer to the binarygen
            //       class for detail example though the code might
            //       be a bit messy

            // std::memset(m_tmp_genotype, 0x0, m_unfiltered_sample_ctl * 2 * sizeof(uintptr_t));
            // here we assume the file is ifstream format
            // dosage_file.seekg(snp_index, std::ios_base::beg); // for jumping to the SNP
            // code for obtaining the dosage information
            // we assume the dosage information is stored in probability vector
            // where the probability vector is a vector of double of double
            // with each index correspond to each individual
            // the below code will properly setup the genotype array
            // you should adjust the following code according to your dosage
            // structure
            /*
            for(size_t i_sample=0; i_sample < probability.size(); ++i_sample)
            {
                auto &&prob = probability[i_sample];
                if(prob.size() != 3) // assume it is hom ref, het, hom alt
                {
                    // this is likely phased
                    fprintf(stderr, "ERROR: Currently don't support phased data\n");
                    fprintf(stderr, "       (It is because the lack of development time)\n");
                    throw std::runtime_error("");
                }
                uintptr_t cur_geno = 1;
                for(size_t g = 0; g < prob.size(); ++g)
                {
                    if(prob[g] >= filter.info_score)
                    {
                        cur_geno = (g==2)? 0: 3-g; // binary code for plink
                        break;
                    }
                }
                // now genotype contain the genotype of this sample after filtering
                // need to bit shift here
                int shift = (i_sample%BITCT*2);
                int index = (i_sample*2)/BITCT;
                genotype[index] |= cur_geno << shift;
            }
            */
        }

        uint32_t load_and_collapse_incl(const uint32_t snp_index,
                const std::string &file_name, uint32_t unfiltered_sample_ct, uint32_t sample_ct,
                const uintptr_t* __restrict sample_include, uintptr_t final_mask,
                uint32_t do_reverse, uintptr_t* __restrict rawbuf, uintptr_t* __restrict mainbuf) {
            assert(unfiltered_sample_ct);
            uint32_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
            if (unfiltered_sample_ct == sample_ct) {
                rawbuf = mainbuf;
            }
            load_raw(rawbuf, snp_index, file_name);

            if (unfiltered_sample_ct != sample_ct) {
                copy_quaterarr_nonempty_subset(rawbuf, sample_include, unfiltered_sample_ct, sample_ct, mainbuf);
            } else {
                mainbuf[(unfiltered_sample_ct - 1) / BITCT2] &= final_mask;
            }
            if (do_reverse) {
                reverse_loadbuf(sample_ct, (unsigned char*)mainbuf);
            }

            // mainbuf should contains the information
            return 0;
        }

        inline void read_genotype(uintptr_t* genotype, const uint32_t snp_index,
                const std::string &file_name)
        {
            // You can simply retain the structure here as long as you use
            // ifstream for your core file
            std::memset(m_tmp_genotype, 0x0, m_unfiltered_sample_ctl * 2 * sizeof(uintptr_t));
            if (load_and_collapse_incl(snp_index, file_name, m_unfiltered_sample_ct, m_founder_ct,
                    m_founder_info, m_final_mask, false,
                    m_tmp_genotype, genotype))
            {
                throw std::runtime_error("ERROR: Cannot read the bed file!");
            }
        };

        void read_score(
                std::vector<std::vector<Sample_lite> > &current_prs_score,
                size_t start_index, size_t end_bound)
        {
            // simplest way will be to reuse the load_and_collapse_incl function though
            // that might not be the most efficient
            // if that's the case, then you can simply use the following codes
            /*
            m_cur_file = "";
            size_t prev =0;
            uint32_t uii;
            uint32_t ujj;
            uint32_t ukk;
            uintptr_t ulii = 0;
            int32_t delta1 = 1;
            int32_t delta2 = 2;
            int32_t deltam = 0;
            uint32_t max_reverse = BITCT_TO_WORDCT(m_unfiltered_marker_ct);
            if(current_prs_score.size() != m_region_size)
            {
                throw std::runtime_error("Size of Matrix doesn't match number of region!!");
            }
            for(size_t i = 0; i < m_region_size; ++i)
            {
                if(current_prs_score[i].size()!=m_unfiltered_sample_ct)
                {
                    throw std::runtime_error("Size of Matrix doesn't match number of samples!!");
                }
            }
            std::vector<bool> in_region(m_region_size);
            // index is w.r.t. partition, which contain all the information
            uintptr_t* genotype = new uintptr_t[m_unfiltered_sample_ctl*2];
            for(size_t i_snp = start_index; i_snp < end_bound; ++i_snp)
            { // for each SNP
                for(size_t i_region = 0; i_region < m_region_size; ++i_region)
                {
                    in_region[i_region]=m_existed_snps[i_snp].in(i_region);
                }
                size_t cur_line = m_existed_snps[i_snp].snp_id();
                std::memset(genotype, 0x0, m_unfiltered_sample_ctl*2*sizeof(uintptr_t));
                std::memset(m_tmp_genotype, 0x0, m_unfiltered_sample_ctl*2*sizeof(uintptr_t));
                if(load_and_collapse_incl(m_existed_snps[i_snp].snp_id(), m_existed_snps[i_snp].file_name(),
                        m_unfiltered_sample_ct, m_founder_ct, m_sample_exclude, m_final_mask,
                        false, m_tmp_genotype, genotype))
                {
                    throw std::runtime_error("ERROR: Cannot read the bed file!");
                }
                uintptr_t* lbptr = genotype;
                uii = 0;
                std::vector<size_t> missing_samples;
                double stat = m_existed_snps[i_snp].stat();
                bool flipped = m_existed_snps[i_snp].is_flipped();
                std::vector<double> genotypes(m_unfiltered_sample_ct);
                int total_num = 0;
                uint32_t sample_idx=0;
                int nmiss = 0;
                // This whole thing was wrong...
                do
                {
                    ulii = ~(*lbptr++);
                    if (uii + BITCT2 > m_unfiltered_sample_ct)
                    {
                        ulii &= (ONELU << ((m_unfiltered_sample_ct & (BITCT2 - 1)) * 2)) - ONELU;
                    }
                    while (ulii)
                    {
                        ujj = CTZLU(ulii) & (BITCT - 2);
                        ukk = (ulii >> ujj) & 3;
                        sample_idx = uii + (ujj / 2);
                        if(ukk==1 || ukk==3) // Because 01 is coded as missing
                        {
                            // 3 is homo alternative
                            //int flipped_geno = snp_list[snp_index].geno(ukk);
                            total_num+=(ukk==3)? 2: ukk;
                            genotypes[sample_idx] = (ukk==3)? 2: ukk;
                        }
                        else // this should be 2
                        {
                            missing_samples.push_back(sample_idx);
                            nmiss++;
                        }
                        ulii &= ~((3 * ONELU) << ujj);
                    }
                    uii += BITCT2;
                } while (uii < m_unfiltered_sample_ct);


                size_t i_missing = 0;
                double maf = ((double)total_num/((double)((int)m_unfiltered_sample_ct-nmiss)*2.0)); // MAF does not count missing
                if(flipped) maf = 1.0-maf;
                double center_score = stat*maf;
                size_t num_miss = missing_samples.size();
                for(size_t i_sample=0; i_sample < m_unfiltered_sample_ct; ++i_sample)
                {
                    if(i_missing < num_miss && i_sample == missing_samples[i_missing])
                    {
                        for(size_t i_region=0; i_region < m_region_size; ++i_region)
                        {
                            if(in_region[i_region])
                            {
                                if(m_scoring == SCORING::MEAN_IMPUTE) current_prs_score[i_region][i_sample].prs += center_score;
                                if(m_scoring != SCORING::SET_ZERO) current_prs_score[i_region][i_sample].num_snp++;
                            }
                        }
                        i_missing++;
                    }
                    else
                    { // not missing sample
                        for(size_t i_region=0; i_region < m_region_size; ++i_region)
                        {
                            if(in_region[i_region])
                            {
                                if(m_scoring == SCORING::CENTER)
                                {
                                    // if centering, we want to keep missing at 0
                                    current_prs_score[i_region][i_sample].prs -= center_score;
                                }
                                current_prs_score[i_region][i_sample].prs += ((flipped)?fabs(genotypes[i_sample]-2):genotypes[i_sample])*stat*0.5;
                                current_prs_score[i_region][i_sample].num_snp++;
                            }
                        }
                    }
                }
            }
            delete [] genotype;
            */
        };
};

#endif /* DOSAGETEMPLATE_H_ */
