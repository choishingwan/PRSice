#include "snp.hpp"


SNP::SNP(const std::string rs_id, const std::string chr, const size_t loc,
		const std::string ref_allele, const std::string alt_allele,
		const double statistic, const double se, const double p_value, long_type* flag,
		const size_t size_of_flag):m_ref_allele(ref_allele), m_alt_allele(alt_allele),
		m_rs_id(rs_id), m_chr(chr), m_loc(loc), m_stat(statistic), m_standard_error(se),
		m_p_value(p_value), m_flags(flag), m_size_of_flag(size_of_flag)  {
    m_bit_size = sizeof(long_type)*CHAR_BIT;
	m_region_clumped = new long_type[((size_of_flag+1)/m_bit_size)+1];
	memset(m_region_clumped, 0x0,(((size_of_flag+1)/m_bit_size)+1)*sizeof(long_type));
}

std::vector<size_t> SNP::sort_by_p(const boost::ptr_vector<SNP> &input){
	std::vector<size_t> idx(input.size());
	std::iota(idx.begin(), idx.end(),0);
	std::sort(idx.begin(), idx.end(), [&input](size_t i1, size_t i2){
		if(input[i1].m_p_value==input[i2].m_p_value){
			if(fabs(input[i1].m_stat)==fabs(input[i2].m_stat)){
				if(input[i1].m_chr.compare(input[i2].m_chr)==0){
					return input[i1].m_loc < input[i2].m_loc;
				}else return input[i1].m_chr.compare(input[i2].m_chr)<0;
			}
			else return fabs(input[i1].m_stat) > fabs(input[2].m_stat);
		}
		else return input[i1].m_p_value < input[i2].m_p_value;
	});
	return idx;
}

SNP::SNP(){
	m_loc=0;
	m_stat=0.0;
	m_standard_error=0.0;
	m_p_value=0.0;
	m_flags=nullptr;
	m_region_clumped=nullptr;
	m_size_of_flag =0;
    m_bit_size = sizeof(long_type)*CHAR_BIT;
}

SNP::~SNP(){
	// This is not exactly safe because if we have also assign this flag
	// to some other SNP, then the flag will be deleted twice
	if(m_flags != nullptr) delete [] m_flags;
	if(m_region_clumped != nullptr) delete [] m_region_clumped;
}

size_t SNP::index_check(const std::string &c_in){
    int temp = atoi(c_in.c_str());
    if(temp < 0) throw std::runtime_error("Index of column cannot be less than 0");
    return temp;
}

size_t SNP::index_check(const std::string &c_in, const std::vector<std::string> &c_header, const std::string &typeOfError){
    for(size_t i = 0; i < c_header.size(); ++i){
        if(c_in.compare(c_header[i])==0){
            return i;
        }
    }
    std::string error_message = typeOfError+": No "+c_in+" colume in input data";
    if(typeOfError.compare("ERROR")==0) throw std::runtime_error(error_message);
    else fprintf(stderr, "%s\n", error_message.c_str());
    return -1; //Cannot find the index
}

std::vector<int> SNP::get_index(const Commander &c_commander, const std::string &c_input){
    // This function should return the index in the following order
    // CHR, A1, A2, STAT, SNP, BP, SE, P
    // Absent field = -1
    // Here we will also check if the field are actually within the file
    std::vector<int> result(9,-1);
    if(c_commander.index()){
        //Index was provided, check if the index is correct, then return the vector
        result[0] = index_check(c_commander.chr());
        result[1] = index_check(c_commander.ref());
        result[2] = index_check(c_commander.alt());
        result[3] = index_check(c_commander.statistic());
        result[4] = index_check(c_commander.snp());
        result[5] = index_check(c_commander.bp());
        result[6] = index_check(c_commander.se());
        result[7] = index_check(c_commander.p());
    }
    else{
        std::ifstream in;
        in.open(c_input.c_str());
        if(!in.is_open()){
            std::string error_message = "Cannot open file: "+ c_input;
            throw std::runtime_error(error_message);
        }
        std::string header_line;
        std::getline(in, header_line);
        in.close();
        if(header_line.empty()){
            std::string error_message = "Empty header line for "+c_input;
            throw std::runtime_error(error_message);
        }
        std::vector<std::string> header = misc::split(header_line);
        result[0] = index_check(c_commander.chr(), header, "WARNING");
        result[1] = index_check(c_commander.ref(), header, "ERROR");
        result[2] = index_check(c_commander.alt(), header, "WARNING");
        result[3] = index_check(c_commander.statistic(), header, "WARNING");
        result[4] = index_check(c_commander.snp(), header, "ERROR");
        result[5] = index_check(c_commander.bp(), header, "WARNING");
        result[6] = index_check(c_commander.se(), header, "WARNING");
        result[7] = index_check(c_commander.p(), header, "ERROR");
        sort( header.begin(), header.end() );
        size_t before = header.size();
        header.erase( unique( header.begin(), header.end() ), header.end() );
        size_t after= header.size();
        if(before!=after){
    			fprintf(stderr, "WARNING: Header contain duplicated elements\n");
    			fprintf(stderr, "         Only the first occurrence is used\n");
    			fprintf(stderr, "         Please do check your input file\n");
        }
    }
    int max_index = -1;
    for(size_t i = 0; i < 9; ++i) max_index = (result[i]> max_index)? result[i]:max_index; //get the maximum index
    result[8] = max_index;
    if(c_commander.index()){
    		std::ifstream in;
    		in.open(c_input.c_str());
    		std::string line;
    		std::getline(in, line);
    		std::vector<std::string> col = misc::split(misc::trimmed(line));
    		if(col.size() < max_index){
    			throw std::runtime_error("ERROR: Number of column in file less than the specified index!");
    		}
    }

    fprintf(stderr,"Column Headers\n");
    fprintf(stderr,"==============================\n");
    if(!c_commander.chr().empty() && result[0] != -1) fprintf(stderr,"Chr            : %s\n", c_commander.chr().c_str());
    fprintf(stderr,"SNP            : %s\n", c_commander.snp().c_str());
    if(!c_commander.bp().empty() && result[5] != -1) fprintf(stderr,"BP             : %s\n", c_commander.bp().c_str());
    fprintf(stderr,"Ref Allele     : %s\n", c_commander.ref().c_str());
    if(!c_commander.alt().empty() && result[2] != -1) fprintf(stderr,"Alt Allele     : %s\n", c_commander.alt().c_str());
    if(!c_commander.statistic().empty() && result[3] != -1) fprintf(stderr,"Statistic      : %s\n", c_commander.statistic().c_str());
    if(!c_commander.se().empty() && result[6] != -1) fprintf(stderr,"Standard Error : %s\n", c_commander.se().c_str());
    fprintf(stderr,"P-value        : %s\n", c_commander.p().c_str());

    return result;
}
