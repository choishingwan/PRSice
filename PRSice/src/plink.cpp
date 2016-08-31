//
//  plink.cpp
//  plink
//
//  Created by Shing Wan Choi on 18/08/2016.
//  Copyright © 2016 Shing Wan Choi. All rights reserved.
//

#include "plink.hpp"

#define MULTIPLEX_LD 1920
#define BITCT 64
#define BITCT2 (BITCT / 2)

std::mutex PLINK::clump_mtx;

PLINK::~PLINK(){
	for(size_t i = 0; i < m_genotype.size(); ++i){
		delete [] m_genotype[i];
		delete [] m_missing[i];
	}
}

#ifdef __LP64__
//This is obtained from plink
void PLINK::ld_dot_prod_batch(__m128i* vec1, __m128i* vec2, __m128i* mask1, __m128i* mask2, int32_t* return_vals, uint32_t iters) {
	const __m128i m1 = {FIVEMASK,FIVEMASK};
	const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
	const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
	__m128i loader1;
	__m128i loader2;
	__m128i sum1;
	__m128i sum2;
	__m128i sum11;
	__m128i sum22;
	__m128i sum12;
	__m128i tmp_sum1;
	__m128i tmp_sum2;
	__m128i tmp_sum12;
	__univec acc;
	__univec acc1;
	__univec acc2;
	__univec acc11;
	__univec acc22;
	acc.vi = _mm_setzero_si128();
	acc1.vi = _mm_setzero_si128();
	acc2.vi = _mm_setzero_si128();
	acc11.vi = _mm_setzero_si128();
	acc22.vi = _mm_setzero_si128();
	do {
		loader1 = *vec1++;
	    	loader2 = *vec2++;
	    	sum1 = *mask2++;
	    	sum2 = *mask1++;
	    	sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
	    	sum1 = _mm_and_si128(sum1, loader1);
	    	sum2 = _mm_and_si128(sum2, loader2);
	    	sum11 = _mm_and_si128(sum1, m1);
	    	sum22 = _mm_and_si128(sum2, m1);
	    	loader1 = _mm_andnot_si128(_mm_add_epi64(m1, sum12), _mm_xor_si128(loader1, loader2));
	    	sum12 = _mm_or_si128(sum12, loader1);

	    	sum1 = _mm_add_epi64(_mm_and_si128(sum1, m2), _mm_and_si128(_mm_srli_epi64(sum1, 2), m2));
	    	sum2 = _mm_add_epi64(_mm_and_si128(sum2, m2), _mm_and_si128(_mm_srli_epi64(sum2, 2), m2));
	    	sum12 = _mm_add_epi64(_mm_and_si128(sum12, m2), _mm_and_si128(_mm_srli_epi64(sum12, 2), m2));

	    	loader1 = *vec1++;
	    	loader2 = *vec2++;
	    	tmp_sum1 = *mask2++;
	    	tmp_sum2 = *mask1++;

	    	tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
	    	tmp_sum1 = _mm_and_si128(tmp_sum1, loader1);
	    	tmp_sum2 = _mm_and_si128(tmp_sum2, loader2);
	    	sum11 = _mm_add_epi64(sum11, _mm_and_si128(tmp_sum1, m1));
	    	sum22 = _mm_add_epi64(sum22, _mm_and_si128(tmp_sum2, m1));
	    	loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12), _mm_xor_si128(loader1, loader2));
	    	tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);

	    	sum1 = _mm_add_epi64(sum1, _mm_add_epi64(_mm_and_si128(tmp_sum1, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum1, 2), m2)));
	    	sum2 = _mm_add_epi64(sum2, _mm_add_epi64(_mm_and_si128(tmp_sum2, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum2, 2), m2)));
	    	sum12 = _mm_add_epi64(sum12, _mm_add_epi64(_mm_and_si128(tmp_sum12, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

	    	loader1 = *vec1++;
	    	loader2 = *vec2++;

	   	tmp_sum1 = *mask2++;
	    	tmp_sum2 = *mask1++;

	    	tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
	    	tmp_sum1 = _mm_and_si128(tmp_sum1, loader1);
	    	tmp_sum2 = _mm_and_si128(tmp_sum2, loader2);
	    	sum11 = _mm_add_epi64(sum11, _mm_and_si128(tmp_sum1, m1));
	    	sum22 = _mm_add_epi64(sum22, _mm_and_si128(tmp_sum2, m1));
	    	loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12), _mm_xor_si128(loader1, loader2));
	    	tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);
	    	sum1 = _mm_add_epi64(sum1, _mm_add_epi64(_mm_and_si128(tmp_sum1, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum1, 2), m2)));
	    	sum2 = _mm_add_epi64(sum2, _mm_add_epi64(_mm_and_si128(tmp_sum2, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum2, 2), m2)));
	    	sum11 = _mm_add_epi64(_mm_and_si128(sum11, m2), _mm_and_si128(_mm_srli_epi64(sum11, 2), m2));
	    	sum22 = _mm_add_epi64(_mm_and_si128(sum22, m2), _mm_and_si128(_mm_srli_epi64(sum22, 2), m2));
	    	sum12 = _mm_add_epi64(sum12, _mm_add_epi64(_mm_and_si128(tmp_sum12, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

	    	acc1.vi = _mm_add_epi64(acc1.vi, _mm_add_epi64(_mm_and_si128(sum1, m4), _mm_and_si128(_mm_srli_epi64(sum1, 4), m4)));
	    	acc2.vi = _mm_add_epi64(acc2.vi, _mm_add_epi64(_mm_and_si128(sum2, m4), _mm_and_si128(_mm_srli_epi64(sum2, 4), m4)));
	    	acc11.vi = _mm_add_epi64(acc11.vi, _mm_add_epi64(_mm_and_si128(sum11, m4), _mm_and_si128(_mm_srli_epi64(sum11, 4), m4)));
	    	acc22.vi = _mm_add_epi64(acc22.vi, _mm_add_epi64(_mm_and_si128(sum22, m4), _mm_and_si128(_mm_srli_epi64(sum22, 4), m4)));
	    	acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(sum12, m4), _mm_and_si128(_mm_srli_epi64(sum12, 4), m4)));

	}while (--iters);
	// moved down because we've almost certainly run out of xmm registers
	const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
	#if MULTIPLEX_LD > 960
		acc1.vi = _mm_add_epi64(_mm_and_si128(acc1.vi, m8), _mm_and_si128(_mm_srli_epi64(acc1.vi, 8), m8));
		acc2.vi = _mm_add_epi64(_mm_and_si128(acc2.vi, m8), _mm_and_si128(_mm_srli_epi64(acc2.vi, 8), m8));
		acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
	#else
		acc1.vi = _mm_and_si128(_mm_add_epi64(acc1.vi, _mm_srli_epi64(acc1.vi, 8)), m8);
		acc2.vi = _mm_and_si128(_mm_add_epi64(acc2.vi, _mm_srli_epi64(acc2.vi, 8)), m8);
		acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
	#endif
	acc11.vi = _mm_and_si128(_mm_add_epi64(acc11.vi, _mm_srli_epi64(acc11.vi, 8)), m8);
	acc22.vi = _mm_and_si128(_mm_add_epi64(acc22.vi, _mm_srli_epi64(acc22.vi, 8)), m8);
	return_vals[0] -= ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
	return_vals[1] += ((acc1.u8[0] + acc1.u8[1]) * 0x1000100010001LLU) >> 48;
	return_vals[2] += ((acc2.u8[0] + acc2.u8[1]) * 0x1000100010001LLU) >> 48;
	return_vals[3] += ((acc11.u8[0] + acc11.u8[1]) * 0x1000100010001LLU) >> 48;
	return_vals[4] += ((acc22.u8[0] + acc22.u8[1]) * 0x1000100010001LLU) >> 48;
}

double PLINK::get_r2(const size_t i, const size_t j){
	uintptr_t founder_ct_mld = (m_num_sample + MULTIPLEX_LD - 1) / MULTIPLEX_LD;
	uint32_t founder_ct_mld_m1 = ((uint32_t)founder_ct_mld) - 1;
	uint32_t founder_ct_mld_rem = (MULTIPLEX_LD / 192) - (founder_ct_mld * MULTIPLEX_LD - m_num_sample) / 192;
	uint32_t fixed_missing_ct;
	uint32_t fixed_non_missing_ct;
	uint32_t non_missing_ct;
	uintptr_t founder_ctwd = m_num_sample / BITCT2;
	uintptr_t founder_ctwd12 = founder_ctwd / 12;
	uintptr_t founder_ctwd12_rem = founder_ctwd - (12 * founder_ctwd12);
	uintptr_t lshift_last = 2 * ((0x7fffffc0 - m_num_sample) % BITCT2);
	long_type* vec1 = m_genotype[j];
	long_type* vec2 = m_genotype[i];
	long_type* mask1 = m_missing[j];
	long_type* mask2 = m_missing[i];
	int32_t dp_result[5];
    fixed_missing_ct = m_num_missing[i];
    fixed_non_missing_ct = m_num_sample - fixed_missing_ct;
	non_missing_ct = fixed_non_missing_ct - m_num_missing[j];
	if (fixed_missing_ct && m_num_missing[j]) {
		non_missing_ct += ld_missing_ct_intersect(mask1, mask2, founder_ctwd12, founder_ctwd12_rem, lshift_last);
	}
    dp_result[0] = m_num_sample;
	dp_result[1] = -fixed_non_missing_ct;
	dp_result[2] = (int)m_num_missing[j] - (int)m_num_sample;
	dp_result[3] = dp_result[1];
	dp_result[4] = dp_result[2];
	while (founder_ct_mld_m1--) {
	    ld_dot_prod_batch((__m128i*)vec1, (__m128i*)vec2, (__m128i*)mask1, (__m128i*)mask2, dp_result, MULTIPLEX_LD / 192);
	    vec1 = &(vec1[MULTIPLEX_LD / BITCT2]);
	    vec2 = &(vec2[MULTIPLEX_LD / BITCT2]);
	    mask1 = &(mask1[MULTIPLEX_LD / BITCT2]);
	    mask2 = &(mask2[MULTIPLEX_LD / BITCT2]);
	}
	ld_dot_prod_batch((__m128i*)vec1, (__m128i*)vec2, (__m128i*)mask1, (__m128i*)mask2, dp_result, founder_ct_mld_rem);
	double non_missing_ctd = (double)((int32_t)non_missing_ct);
	double dxx = dp_result[1];
	double dyy = dp_result[2];
	double cov12 = dp_result[0] * non_missing_ctd - dxx * dyy;
	dxx = (dp_result[3] * non_missing_ctd + dxx * dxx) * (dp_result[4] * non_missing_ctd + dyy * dyy);
	dxx = (cov12 * cov12) / dxx;
	return dxx;
}


#else
// This should work for uint32_t but not uint64_t because of the special popcount he used
double PLINK::get_r2(const size_t i, const size_t j){
	double r2 =0.0;
	if(i >= m_genotype.size() || j >=m_genotype.size()) throw std::runtime_error("Out of bound error! In R2 calculation.");
	// Crazy stuff of plink here
	size_t range = (m_required_bit /(m_bit_size))+1;
	long_type loader1, loader2, sum1, sum2, sum11, sum12, sum22;
	long_type final_sum1 = 0;
	long_type final_sum2 = 0;
	long_type final_sum11 = 0;
	long_type final_sum22 = 0;
	long_type final_sum12 = 0;
	double return_vals[5];
	return_vals[0] = (double) m_num_sample;
	return_vals[1] = -(double) m_num_missing[j];
	return_vals[2] = -(double) m_num_missing[i];
	return_vals[3] = return_vals[1];
	return_vals[4] = return_vals[2];
	size_t N =0;
	for(size_t i_geno = 0; i_geno < range;){
	        loader1 = m_genotype[i][i_geno];
	        loader2 = m_genotype[j][i_geno];
	        sum1 = m_missing[j][i_geno];
	        sum2 = m_missing[i][i_geno];
			i_geno++;
			N+= __builtin_popcountll(sum1&sum2)/2;
			sum12 = (loader1 | loader2) & FIVEMASK;
			sum1 = sum1 & loader1;
			sum2 = sum2 & loader2;
			loader1 = (loader1 ^ loader2) & (AAAAMASK - sum12);
			sum12 = sum12 | loader1;
			sum11 = sum1 & FIVEMASK;
			sum22 = sum2 & FIVEMASK;
			sum1 = (sum1 & THREEMASK) + ((sum1 >> 2) & THREEMASK);
			sum2 = (sum2 & THREEMASK) + ((sum2 >> 2) & THREEMASK);
			sum12 = (sum12 & THREEMASK) + ((sum12 >> 2) & THREEMASK);
			long_type tmp_sum1=0 , tmp_sum2=0;
			if(i_geno < range){
				loader1 = m_genotype[i][i_geno];
				loader2 = m_genotype[j][i_geno];
				tmp_sum1 = m_missing[j][i_geno];
				tmp_sum2 = m_missing[i][i_geno];
				N+= __builtin_popcountll(tmp_sum1&tmp_sum2)/2;
			}
			else{
				loader1 = 0;
				loader2 = 0;
			}
		    i_geno++;
		    long_type tmp_sum12 = (loader1 | loader2) & FIVEMASK;
			tmp_sum1 = tmp_sum1 & loader1;
			tmp_sum2 = tmp_sum2 & loader2;
			loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
			tmp_sum12 = tmp_sum12 | loader1;
			sum11 += tmp_sum1 & FIVEMASK;
			sum22 += tmp_sum2 & FIVEMASK;
			sum1 += (tmp_sum1 & THREEMASK) + ((tmp_sum1 >> 2) & THREEMASK);
			sum2 += (tmp_sum2 & THREEMASK) + ((tmp_sum2 >> 2) & THREEMASK);
			sum12 += (tmp_sum12 & THREEMASK) + ((tmp_sum12 >> 2) & THREEMASK);
		    if(i_geno < range){
				loader1 = m_genotype[i][i_geno];
				loader2 = m_genotype[j][i_geno];
				tmp_sum1 = m_missing[j][i_geno];
				tmp_sum2 = m_missing[i][i_geno];
				N+= __builtin_popcountll(tmp_sum1&tmp_sum2)/2;
			}
			else{
				loader1=0;
				loader2=0;
				tmp_sum1=0;
				tmp_sum2=0;
			}
			i_geno++;
			tmp_sum12 = (loader1 | loader2) & FIVEMASK;
			tmp_sum1 = tmp_sum1 & loader1;
			tmp_sum2 = tmp_sum2 & loader2;
			loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
			tmp_sum12 = tmp_sum12 | loader1;
			sum11 += tmp_sum1 & FIVEMASK;
			sum22 += tmp_sum2 & FIVEMASK;
			sum1 += (tmp_sum1 & THREEMASK) + ((tmp_sum1 >> 2) & THREEMASK);
			sum2 += (tmp_sum2 & THREEMASK) + ((tmp_sum2 >> 2) & THREEMASK);
			sum11 = (sum11 & THREEMASK) + ((sum11 >> 2) & THREEMASK);
			sum22 = (sum22 & THREEMASK) + ((sum22 >> 2) & THREEMASK);
			sum12 += (tmp_sum12 & THREEMASK) + ((tmp_sum12 >> 2) & THREEMASK);
			sum1 = (sum1 & OFMASK) + ((sum1 >> 4) & OFMASK);
			sum2 = (sum2 & OFMASK) + ((sum2 >> 4) & OFMASK);
			sum11 = (sum11 & OFMASK) + ((sum11 >> 4) & OFMASK);
			sum22 = (sum22 & OFMASK) + ((sum22 >> 4) & OFMASK);
			sum12 = (sum12 & OFMASK) + ((sum12 >> 4) & OFMASK);
			final_sum1 += (sum1 * ONEZEROMASK) >> 24;
			final_sum2 += (sum2 * ONEZEROMASK) >> 24;
			final_sum11 += (sum11 * ONEZEROMASK) >> 24;
			final_sum22 += (sum22 * ONEZEROMASK) >> 24;
			final_sum12 += (sum12 * ONEZEROMASK) >> 24;
		}

		return_vals[0] -= final_sum12;
		return_vals[1] += final_sum1;
		return_vals[2] += final_sum2;
		return_vals[3] += final_sum11;
		return_vals[4] += final_sum22;

	    double dxx = return_vals[1];
	    double dyy = return_vals[2];
	    double n = N;
	    double cov12 = return_vals[0] * n - dxx * dyy;
	    dxx = (return_vals[3] * n + dxx * dxx) * (return_vals[4] * n + dyy * dyy);
	    if(dxx !=0.0) r2 =(cov12 * cov12) / dxx;
	return r2;
}
#endif

void PLINK::compute_clump( size_t index, size_t i_start, size_t i_end, std::vector<SNP> &snp_list, const std::deque<size_t> &index_check, const double r2_threshold){
	std::cerr << "Computing the clump" << std::endl;
	size_t ref_index = index_check[index];
	std::vector<size_t> self_index; // index we want to push into the current index
	for(size_t i = i_start; i < i_end && i < index_check.size(); ++i){
		size_t target_index = index_check[i];
		if(i != index && snp_list[target_index].get_p_value() > snp_list[ref_index].get_p_value()){
			// only calculate r2 if more significant
			double r2 = get_r2(i, index);
			if(r2 >= r2_threshold) self_index.push_back(target_index);
		}
	}
	PLINK::clump_mtx.lock();
		snp_list[index].add_clump(self_index);
	PLINK::clump_mtx.unlock();
}

void PLINK::clump_thread(const size_t index, const std::deque<size_t> &index_check, std::vector<SNP> &snp_list, const double r2_threshold){
	std::vector<std::thread> thread_store;
	std::cerr << "In clump thread" << std::endl;
	if((index_check.size()-1) < m_thread){
		for(size_t i = 0; i < index_check.size(); ++i){
			if(index_check[i]!=index) thread_store.push_back(std::thread(&PLINK::compute_clump, this, index,i, i+1, std::ref(snp_list), std::cref(index_check), r2_threshold));
		}
	}
	else{
		int num_snp_per_thread =(int)(index_check.size()-1) / (int)m_thread;  //round down
		int remain = (int)(index_check.size()-1) % (int)m_thread;
		int cur_start = 0;
		int cur_end = num_snp_per_thread;
		for(size_t i = 0; i < m_thread; ++i){
			thread_store.push_back(std::thread(&PLINK::compute_clump, this, index, cur_start, cur_end+(remain>0), std::ref(snp_list), std::cref(index_check),r2_threshold ));
			cur_start = cur_end+(remain>0);
			cur_end+=num_snp_per_thread+(remain>0);
			if(cur_end>index_check.size()) cur_end =index_check.size();
			remain--;
		}
	}
	for (size_t i = 0; i < thread_store.size(); ++i) thread_store[i].join();
	thread_store.clear();
}

void PLINK::clumping(std::map<std::string, bool> inclusion, std::vector<SNP> &snp_list, const std::map<std::string, size_t> &snp_index, double p_threshold, double r2_threshold, size_t kb_threshold){
	// Go through all SNPs
	std::deque<size_t> bp_check;
	std::deque<size_t> index_check;
	std::string prev_chr = "";
	size_t require_bp=0, require_index=0; //this is the index on index_check
	bool requiring=false;
	std::cerr << "Start clumping" << std::endl;
	for(size_t i = 0; i < m_snp_list.size(); ++i){
		std::string rs = m_snp_list[i];
		if(inclusion.find(rs)!=inclusion.end() && snp_index.find(rs)!=snp_index.end()){
			// required SNP
			std::cerr << "Require snp:" << rs << std::endl;
			std::string cur_chr = snp_list[snp_index.at(rs)].get_chr();
			if(prev_chr.empty()) prev_chr=cur_chr; // First SNP
			size_t cur_loc = snp_list[snp_index.at(rs)].get_loc();
			if(!requiring && bp_check.size() != 0){ // no point doing this if this is the first snp
				// none of the previous snps are required
				if(prev_chr.compare(cur_chr)!=0){
					// Change chromosome without anything to do, just clean up
					bp_check.clear();
					index_check.clear();
					lerase(m_genotype.size());
				}
				else{
					// Still the same chromosome, just need to remove anything that
					// are out of range
					size_t remove = 0;
					while(cur_loc-bp_check.front() > kb_threshold) remove++;
					if(remove!=0){
						bp_check.erase(bp_check.begin()+remove);
						index_check.erase(index_check.begin()+remove);
						lerase(remove);
					}
				}
			}
			while(bp_check.size() != 0 && requiring && (cur_loc-require_bp > kb_threshold || cur_chr.compare(prev_chr)!=0)){
				// All genotypes are here, now calculate the P-value
				clump_thread(require_index, index_check, snp_list, r2_threshold);
				// now go through the whole remaining index to check if there are any required
				requiring = false;
				for(size_t check = require_index+1; check< index_check.size(); ++check){
					if(snp_list[index_check[check]].get_p_value() < p_threshold){
						requiring = true;
						require_index =check;
						prev_chr = snp_list[index_check[check]].get_chr();
						require_bp = snp_list[index_check[check]].get_loc();
					}
				}
				if(requiring){
					// remove SNPs in the front that are too far away
					size_t num_remove=0;
					for(size_t check; check < require_index; ++check){
						if(require_bp-snp_list[index_check[check]].get_loc() > kb_threshold)num_remove++;
						else break;
					}
					if(num_remove!=0) lerase(num_remove);
				}
			}
			if(prev_chr.compare(cur_chr)!=0){
				// When we get here, we just remove everything
				bp_check.clear();
				index_check.clear();
				lerase(m_genotype.size());
				prev_chr = cur_chr; //update the chromosome
			}
			std::cerr << "Start reading snp" << std::endl;
			//Then read the SNP
			read_snp(1, true);
			std::cerr << "Got the snp" << std::endl;
			bp_check.push_back(cur_loc);
			index_check.push_back(snp_index.at(rs));
			if(!requiring && snp_list[index_check.back()].get_p_value() < p_threshold){
				//only add in if we haven't found ourself an index snp
				require_bp = cur_loc;
				prev_chr = cur_chr;
				require_index = index_check.size()-1;
			}
		}
		else m_bed.seekg(m_num_bytes, m_bed.cur); //This should skip 1 SNP
	}
	if(m_genotype.size() != 0){
		//Still got stuff to work on
		// keep clumping
		while(requiring){
			clump_thread(require_index, index_check, snp_list, r2_threshold);
			requiring = false;
			for(size_t check = require_index+1; check< index_check.size(); ++check){
				if(snp_list[index_check[check]].get_p_value() < p_threshold){
					requiring = true;
					require_index =check;
					require_bp = snp_list[index_check[check]].get_loc();
				}
			}
			if(requiring){
				size_t num_remove=0;
				for(size_t check; check < require_index; ++check){
					if(require_bp-snp_list[index_check[check]].get_loc() > kb_threshold)num_remove++;
					else break;
				}
				if(num_remove!=0) lerase(num_remove);
			}
		}
	}
	//completed
	// Now get the list of SNPs that we want to retain (in index)
	std::map<std::string, bool> include_ref = inclusion; // now update the inclusion such that it only contain the index snps
	inclusion.clear();
	std::vector<size_t> p_sort_order = SNP::sort_by_p(snp_list);
	for(size_t i = 0; i < p_sort_order.size(); ++i){
		if(include_ref.find(snp_list[p_sort_order[i]].get_rs_id()) != include_ref.end() && !snp_list[p_sort_order[i]].clumped() && snp_list[p_sort_order[i]].get_p_value() < p_threshold){
			snp_list[p_sort_order[i]].clump_all(snp_list);
			inclusion[snp_list[p_sort_order[i]].get_rs_id()]=true;
		}
		else if(snp_list[p_sort_order[i]].get_p_value() >= p_threshold) break;
	}
	//Anything remaining should be the required SNPs
	//This is for testing
	for(std::map<std::string, bool>::iterator iter = inclusion.begin();
			iter != inclusion.end(); ++iter){
		std::cout << iter->first << std::endl;
	}
}

std::vector<int> PLINK::get_genotype(int geno) const{
	if(geno >= m_genotype.size()){
		std::string error_message = "Asked for "+std::to_string(geno)+" genotype but contain only "+std::to_string(m_genotype.size());
		throw std::runtime_error(error_message);
	}
	std::vector<int> res(m_num_sample);
	for(size_t i = 0; i < m_num_sample; ++i){
		int index =(i*2)/m_bit_size;
		int info = m_genotype[geno][index] >> (m_bit_size-(i+1)*2)& THREEMASK;
		switch(info){
			case 0:
				res[i] = 0;
				break;
			case 1:
				res[i] = -1;
				break;
			case 2:
				res[i] = 1;
				break;
			case 3:
				res[i] = 2;
				break;
			default:
				throw std::runtime_error("Undefined genotype");
		}
	}
	return res;
}

void PLINK::lerase(int num){
    if(num <0){
        std::string error_message = "Number of removed SNPs cannot be less than 1: "+std::to_string(num);
        throw std::runtime_error(error_message);
    }
    if(num >= m_genotype.size()){
        std::string error_message = "Number of removed SNPs exceed number of SNPs available "+std::to_string(num)+" "+std::to_string(m_genotype.size());
        throw std::runtime_error(error_message);
    }
    for(size_t i = 0; i < num; ++i){
        delete [] m_genotype[i];
        delete [] m_missing[i];
    }
    m_genotype.erase(m_genotype.begin(), m_genotype.begin()+num);
    m_missing.erase(m_missing.begin(), m_missing.begin()+num);
    m_chr_list.erase(m_chr_list.begin(), m_chr_list.begin()+num);
    //m_snp_list.erase(m_snp_list.begin(), m_snp_list.begin()+num);
    m_cm_list.erase(m_cm_list.begin(), m_cm_list.begin()+num);
    m_bp_list.erase(m_bp_list.begin(), m_bp_list.begin()+num);
    m_ref_allele.erase(m_ref_allele.begin(), m_ref_allele.begin()+num);
    m_alt_allele.erase(m_alt_allele.begin(), m_alt_allele.begin()+num);
    m_maf.erase(m_maf.begin(), m_maf.begin()+num);
    m_num_missing.erase(m_num_missing.begin(), m_num_missing.begin()+num);
}

//The return value should be the number of remaining SNPs
int PLINK::read_snp(int num_snp, bool ld){
    if(!m_init) throw std::runtime_error("Class uninitialize! Must initialize before use!");
    if(num_snp <= 0){
        std::string error_message = "Number of required SNPs cannot be less than 1: "+std::to_string(num_snp);
        throw std::runtime_error(error_message);
    }
    std::string line;
    //First get the information of the SNPs
    size_t cur_iter = 0;
    for(;m_snp_iter < m_num_snp & cur_iter<num_snp; ++m_snp_iter){
    		std::getline(m_bim, line);
    		misc::trim(line);
    	    	if(!line.empty()){
    	    		std::vector<std::string> token = misc::split(line);
    	    		if(token.size() >= 6){
    	    			m_chr_list.push_back(token[0]);
    	    			//m_snp_list.push_back(token[1]);
    	    			int temp = misc::convert<int>(token[2]);
    	    			if(temp < 0){
    	    				std::string error_message = "Negative CM: "+line;
    	    				throw std::runtime_error(error_message);
    	    			}
    	    			m_cm_list.push_back(temp);
    	    			temp = misc::convert<int>(token[3]);
    	    			if(temp < 0){
    	    				std::string error_message = "Negative BP: "+line;
    	    				throw std::runtime_error(error_message);
    	    			}
    	    			m_bp_list.push_back(temp);
    	    			m_ref_allele.push_back(token[4]);
    	    			m_alt_allele.push_back(token[5]);
    	    		}
    	    		else throw std::runtime_error("Malformed bim file");
    	    	}else throw std::runtime_error("Malformed bim file");
        cur_iter++;
        char genotype_list[m_num_bytes];
        m_bed.read(genotype_list, m_num_bytes);
        size_t i_genotype = 0;
        size_t total_allele = 0;
        size_t num_missing = 0;
        long_type *genotype = new long_type[(m_required_bit /(m_bit_size))+1];
        long_type *missing = new long_type[(m_required_bit /(m_bit_size))+1];
        for(size_t byte_runner= 0; byte_runner < m_num_bytes;){
#ifdef __LP64__
            long_type current_genotypes = 0ULL;
#else
            long_type current_genotypes=0UL;
#endif
            for(int byte_set = 0; byte_set < sizeof(long_type)/sizeof(char) && byte_runner < m_num_bytes; ++byte_set){
            		//long_type current_byte = static_cast<long_type>(genotype_list[byte_runner]) << ((sizeof(long_type)-1)*CHAR_BIT) >> CHAR_BIT*byte_set;
            		long_type current_byte = static_cast<long_type>(genotype_list[byte_runner]) << ((sizeof(long_type)-1)*CHAR_BIT) >> (((sizeof(long_type)-1)-byte_set)*CHAR_BIT);
            		current_genotypes |= current_byte;
                byte_runner++;
            }
            long_type five_masked_geno = current_genotypes & FIVEMASK;
            long_type inter = (five_masked_geno & (current_genotypes>>1)) ^ five_masked_geno;
            long_type current_missing = inter | (inter << 1);
            genotype[i_genotype] = current_genotypes;
            if(!ld) genotype[i_genotype] = current_genotypes;
            else {
            		genotype[i_genotype] =(current_genotypes &(five_masked_geno <<1));
            		genotype[i_genotype] |= (five_masked_geno^((current_genotypes &(FIVEMASK*2))>>1));
            }
            missing[i_genotype] = ~current_missing;
            total_allele += __builtin_popcountll(current_genotypes & (~current_missing));
            num_missing +=__builtin_popcountll(inter); // because inter only contain one bit for each missing
            i_genotype++;
        }
        m_genotype.push_back(genotype);
        m_missing.push_back(missing);
        double maf = (double)total_allele/((double)m_required_bit-(double)num_missing);
        maf = (maf > 0.5)? 1.0-maf: maf;
        m_maf.push_back(maf);
        m_num_missing.push_back(num_missing);

    }
    return m_num_snp-m_snp_iter;
}

void PLINK::initialize(){
    std::string fam_name = m_prefix+".fam";
    std::string bim_name = m_prefix+".bim";
    std::string bed_name = m_prefix+".bed";
    // Start processing the fam file
    std::ifstream fam;
    fam.open(fam_name.c_str());
    if(!fam.is_open()){
        std::string error_message = "Cannot open fam file: "+fam_name;
        throw std::runtime_error(error_message);
    }
    std::string line;
    while(std::getline(fam, line))
        if(!misc::trimmed(line).empty()) m_num_sample++;
    fam.close();
    // Check whether if the bed file is correct
    bool snp_major = openPlinkBinaryFile(bed_name, m_bed);
    if(!snp_major) throw std::runtime_error("Currently does not support sample major format");
    // Check whether if the bim file is correct
    m_bim.open(bim_name.c_str());
    if(!m_bim.is_open()){
        std::string error_message = "Cannot open bim file: "+bim_name;
        throw std::runtime_error(error_message);
    }
    
    while(std::getline(m_bim, line)){
        if(!misc::trimmed(line).empty()) m_num_snp++;
        m_snp_list.push_back(misc::split(line)[1]); // This is dangerous as we don't check bim file format
    }
    m_bim.clear();
    m_bim.seekg(0);
    m_num_bytes=ceil((double)m_num_sample/4.0);
    m_required_bit = m_num_sample*2;
    m_snp_iter=0;
    m_init = true;
}

bool PLINK::openPlinkBinaryFile(const std::string s, std::ifstream & BIT){
    BIT.open(s.c_str(), std::ios::in | std::ios::binary);
    if(!BIT.is_open()){
        std::string error_message= "Cannot open the bed file: "+s;
        throw std::runtime_error(error_message);
    }
    // 1) Check for magic number
    // 2) else check for 0.99 SNP/Ind coding
    // 3) else print warning that file is too old
    char ch[3];
    BIT.read(ch,3);
    bool bfile_SNP_major = false;
    bool v1_bfile = true;
    // If v1.00 file format
    // Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file
    // check magic number
    if(static_cast<uint32_t>(ch[0])== 108 && static_cast<uint32_t>(ch[1])==27)
        bfile_SNP_major=(static_cast<uint32_t>(ch[2])==1);
    else
        v1_bfile = false;
    // Reset file if < v1
    if ( ! v1_bfile ) {
        std::cerr << "Warning, old BED file <v1.00 : will try to recover..." << std::endl;
        std::cerr << "  but you should --make-bed from PED )" << std::endl;
        BIT.close();
        BIT.clear();
        BIT.open(s.c_str(), std::ios::in | std::ios::binary);
        BIT.read(ch,1);
        uint32_t file_info = static_cast<uint32_t>(ch[0]);
        if(file_info != 1 && file_info!=0){
            std::cerr << std::endl <<   " *** Possible problem: guessing that BED is < v0.99       *** " << std::endl;
            std::cerr <<                " *** High chance of data corruption, spurious results     *** " << std::endl;
            std::cerr <<                " *** Unless you are _sure_ this really is an old BED file *** " << std::endl;
            std::cerr <<                " *** you should recreate PED -> BED                       *** " << std::endl << std::endl;
            bfile_SNP_major = false;
            BIT.close();
            BIT.clear();
            BIT.open(s.c_str(), std::ios::in | std::ios::binary);
        }
        else{
            std::cerr << "Binary PED file is v0.99" << std::endl;
            bfile_SNP_major = (file_info==1);
        }
    }
    return bfile_SNP_major;
}




uint32_t PLINK::ld_missing_ct_intersect(long_type* lptr1, long_type* lptr2, uintptr_t word12_ct, uintptr_t word12_rem, uintptr_t lshift_last) {
  // variant of popcount_longs_intersect()
	uintptr_t tot = 0;
	long_type* lptr1_end2;
#ifdef __LP64__
	const __m128i m1 = {FIVEMASK, FIVEMASK};
	const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
	const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
	const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
	__m128i* vptr1 = (__m128i*)lptr1;
	__m128i* vptr2 = (__m128i*)lptr2;
	__m128i* vend1;
	__m128i loader1;
	__m128i loader2;
	__univec acc;

	while (word12_ct >= 10) {
		word12_ct -= 10;
		vend1 = &(vptr1[60]);
		ld_missing_ct_intersect_main_loop:
		acc.vi = _mm_setzero_si128();
		do {
			loader1 = _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1);
			loader2 = _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1);
			loader1 = _mm_add_epi64(loader1, _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
			loader2 = _mm_add_epi64(loader2, _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
			loader1 = _mm_add_epi64(loader1, _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
			loader2 = _mm_add_epi64(loader2, _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
			loader1 = _mm_add_epi64(_mm_and_si128(loader1, m2), _mm_and_si128(_mm_srli_epi64(loader1, 2), m2));
			loader1 = _mm_add_epi64(loader1, _mm_add_epi64(_mm_and_si128(loader2, m2), _mm_and_si128(_mm_srli_epi64(loader2, 2), m2)));
			acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(loader1, m4), _mm_and_si128(_mm_srli_epi64(loader1, 4), m4)));
		} while (vptr1 < vend1);
		acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
		tot += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
	}
	if (word12_ct) {
		vend1 = &(vptr1[word12_ct * 6]);
		word12_ct = 0;
		goto ld_missing_ct_intersect_main_loop;
	}
	lptr1 = (long_type*)vptr1;
	lptr2 = (long_type*)vptr2;
#else
	uintptr_t* lptr1_end = &(lptr1[word12_ct * 12]);
	uintptr_t tmp_stor;
	uintptr_t loader1;
	uintptr_t loader2;
	while (lptr1 < lptr1_end) {
		loader1 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
		loader2 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
		loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
		loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
		loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
		loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
		loader1 = (loader1 & 0x33333333) + ((loader1 >> 2) & 0x33333333);
		loader1 += (loader2 & 0x33333333) + ((loader2 >> 2) & 0x33333333);
		tmp_stor = (loader1 & 0x0f0f0f0f) + ((loader1 >> 4) & 0x0f0f0f0f);

		loader1 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
		loader2 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
		loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
		loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
		loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
		loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
		loader1 = (loader1 & 0x33333333) + ((loader1 >> 2) & 0x33333333);
		loader1 += (loader2 & 0x33333333) + ((loader2 >> 2) & 0x33333333);
    		tmp_stor += (loader1 & 0x0f0f0f0f) + ((loader1 >> 4) & 0x0f0f0f0f);
    		tot += (tmp_stor * 0x01010101) >> 24;
	}
#endif
	lptr1_end2 = &(lptr1[word12_rem]);
	while (lptr1 < lptr1_end2) {
		tot += popcount2_long((~((*lptr1++) | (*lptr2++))) & FIVEMASK);
	}
	if (lshift_last) {
		tot += popcount2_long(((~((*lptr1) | (*lptr2))) & FIVEMASK) << lshift_last);
	}
	return tot;
}
