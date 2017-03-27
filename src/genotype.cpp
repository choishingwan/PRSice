/*
 * genotype.cpp
 *
 *  Created on: 27 Mar 2017
 *      Author: shingwanchoi
 */

#include "genotype.hpp"

Genotype::Genotype(std::string prefix, const size_t thread, uint32_t species_code, bool verbose) {
	// TODO Auto-generated constructor stub
	if (init_chrom_info(&m_chrom_info)) {
		throw std::runtime_error("Cannot initialize the chromosome information");
	}
}

Genotype::~Genotype() {
	// TODO Auto-generated destructor stub
}


int32_t init_delim_and_species(uint32_t flag_ct, char* flag_buf, uint32_t* flag_map, int32_t argc, char** argv, char* range_delim_ptr, Chrom_info* chrom_info_ptr) {
	uint32_t species_code = SPECIES_DEFAULT;
	uint32_t flag_idx = 0;
	uint32_t retval = 0;
	int32_t cur_arg;
	uint32_t param_ct;
	int32_t ii;
	uint32_t param_idx;
	if (flag_match("autosome-num", &flag_idx, flag_ct, flag_buf)) {
		species_code = SPECIES_UNKNOWN;
		cur_arg = flag_map[flag_idx - 1];
		param_ct = param_count(argc, argv, cur_arg);
		if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
			goto init_delim_and_species_ret_INVALID_CMDLINE_2A;
		}
		if (scan_posint_capped(argv[cur_arg + 1], MAX_CHROM_TEXTNUM, (uint32_t*)(&ii))) {
			sprintf(g_logbuf, "Error: Invalid --autosome-num parameter '%s'.\n", argv[cur_arg + 1]);
			goto init_delim_and_species_ret_INVALID_CMDLINE_WWA;
		}
		chrom_info_ptr->xymt_codes[X_OFFSET] = ii + 1;
		chrom_info_ptr->xymt_codes[Y_OFFSET] = -1;
		chrom_info_ptr->xymt_codes[XY_OFFSET] = -1;
		chrom_info_ptr->xymt_codes[MT_OFFSET] = -1;
		chrom_info_ptr->max_code = ii + 1;
		chrom_info_ptr->autosome_ct = ii;
		set_bit(ii + 1, chrom_info_ptr->haploid_mask);
	}
	if (flag_match("chr-set", &flag_idx, flag_ct, flag_buf)) {
		if (species_flag(&species_code, SPECIES_UNKNOWN)) {
			goto init_delim_and_species_ret_INVALID_CMDLINE;
		}
		cur_arg = flag_map[flag_idx - 1];
		param_ct = param_count(argc, argv, cur_arg);
		if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 5)) {
			goto init_delim_and_species_ret_INVALID_CMDLINE_2A;
		}
		if (scan_int_abs_bounded(argv[cur_arg + 1], MAX_CHROM_TEXTNUM, &ii) || (!ii)) {
			sprintf(g_logbuf, "Error: Invalid --chr-set parameter '%s'.\n", argv[cur_arg + 1]);
			goto init_delim_and_species_ret_INVALID_CMDLINE_WWA;
		}
		if (ii < 0) {
			if (param_ct > 1) {
				logerrprint("Error: --chr-set does not accept multiple parameters in haploid mode.\n");
			}
			ii = -ii;
			chrom_info_ptr->autosome_ct = ii;
			chrom_info_ptr->xymt_codes[X_OFFSET] = -1;
			chrom_info_ptr->xymt_codes[Y_OFFSET] = -1;
			chrom_info_ptr->xymt_codes[XY_OFFSET] = -1;
			chrom_info_ptr->xymt_codes[MT_OFFSET] = -1;
			chrom_info_ptr->max_code = ii;
			fill_all_bits(((uint32_t)ii) + 1, chrom_info_ptr->haploid_mask);
		} else {
			chrom_info_ptr->autosome_ct = ii;
			chrom_info_ptr->xymt_codes[X_OFFSET] = ii + 1;
			chrom_info_ptr->xymt_codes[Y_OFFSET] = ii + 2;
			chrom_info_ptr->xymt_codes[XY_OFFSET] = ii + 3;
			chrom_info_ptr->xymt_codes[MT_OFFSET] = ii + 4;
			set_bit(ii + 1, chrom_info_ptr->haploid_mask);
			set_bit(ii + 2, chrom_info_ptr->haploid_mask);
			for (param_idx = 2; param_idx <= param_ct; param_idx++) {
				if (!strcmp(argv[cur_arg + param_idx], "no-x")) {
					chrom_info_ptr->xymt_codes[X_OFFSET] = -1;
					clear_bit(ii + 1, chrom_info_ptr->haploid_mask);
				} else if (!strcmp(argv[cur_arg + param_idx], "no-y")) {
					chrom_info_ptr->xymt_codes[Y_OFFSET] = -1;
					clear_bit(ii + 2, chrom_info_ptr->haploid_mask);
				} else if (!strcmp(argv[cur_arg + param_idx], "no-xy")) {
					chrom_info_ptr->xymt_codes[XY_OFFSET] = -1;
				} else if (!strcmp(argv[cur_arg + param_idx], "no-mt")) {
					chrom_info_ptr->xymt_codes[MT_OFFSET] = -1;
				} else {
					sprintf(g_logbuf, "Error: Invalid --chr-set parameter '%s'.\n", argv[cur_arg + param_idx]);
				}
			}
			if (chrom_info_ptr->xymt_codes[MT_OFFSET] != -1) {
				chrom_info_ptr->max_code = ii + 4;
			} else if (chrom_info_ptr->xymt_codes[XY_OFFSET] != -1) {
				chrom_info_ptr->max_code = ii + 3;
			} else if (chrom_info_ptr->xymt_codes[Y_OFFSET] != -1) {
				chrom_info_ptr->max_code = ii + 2;
			} else if (chrom_info_ptr->xymt_codes[X_OFFSET] != -1) {
				chrom_info_ptr->max_code = ii + 1;
			} else {
				chrom_info_ptr->max_code = ii;
			}
		}
	}
// set species flag accordingly
	return retval;
}
