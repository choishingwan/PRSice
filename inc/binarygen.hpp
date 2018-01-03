// This file is part of PRSice2.0, copyright (C) 2016-2017
// Shing Wan Choi, Jack Euesden, Cathryn M. Lewis, Paul F. O’Reilly
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
#ifndef BinaryGEN_H
#define BinaryGEN_H

#include "genotype.hpp"
#include <stdexcept>
#include <zlib.h>

/**
 * Potential problem:
 * Where multi-allelic variants exist in these data, they have been
 * split into a series of bi-allelic variants. This implies that
 * several variants may share the same genomic position but with
 * different alternative alleles.
 */
class BinaryGen : public Genotype
{
public:
    BinaryGen(const std::string& prefix, const std::string& sample_file,
              const size_t thread = 1, const bool ignore_fid = false,
              const bool keep_nonfounder = false,
              const bool keep_ambig = false);
    ~BinaryGen();

private:
    enum FlagMask
    {
        e_NoFlags = 0,
        e_CompressedSNPBlocks = 0x3,
        e_Layout = 0x3C
    };
    enum Layout
    {
        e_Layout0 = 0x0,
        e_Layout1 = 0x4,
        e_Layout2 = 0x8
    };
    enum Structure
    {
        e_SampleIdentifiers = 0x80000000
    };
    enum Compression
    {
        e_NoCompression = 0,
        e_ZlibCompression = 1,
        e_ZstdCompression = 2
    };
    enum OrderType
    {
        eUnknownOrderType = 0,
        eUnorderedList = 1, // a list, treated as unordered
        eOrderedList = 2,   // a list, treated as ordered
        ePerUnorderedGenotype =
            3, // E.g. genotype probabilities; GEN, BGEN v1.1, etc.
        ePerOrderedHaplotype =
            4, // E.g. Phased GT, SHAPEIT or IMPUTE haplotypes
        ePerUnorderedHaplotype = 5, // E.g. Unphased GT, binary PED file.
        ePerPhasedHaplotypePerAllele =
            6,          // E.g. BGEN v1.2-style haplotype probabilities.
        ePerAllele = 7, // E.g. assay intensities.
        ePerSample = 8  // 'B' allele dosage, one value per sample.
    };
    enum ValueType
    {
        eUnknownValueType = 0,
        eProbability = 1,
        eAlleleIndex = 2,
        eDosage = 3
    };
    struct Context
    {
        Context();
        Context(Context const& other);
        Context& operator=(Context const& other);
        uint32_t header_size() const;

    public:
        uint32_t number_of_samples;
        uint32_t number_of_variants;
        std::string magic;
        std::string free_data;
        uint32_t flags;
        uint32_t offset;
    };
    typedef uint8_t byte_t;

    std::unordered_map<std::string, Context> m_context_map;
    std::vector<Sample> gen_sample_vector();
    // check if the sample file is of the sample format specified by bgen
    // or just a simple text file
    bool check_is_sample_format();
    std::vector<SNP> gen_snp_vector(const double geno, const double maf,
                                    const double info_score,
                                    const double hard_threshold,
                                    const bool hard_coded,
                                    const std::string& out_prefix);
    Context get_context(std::string& bgen_name);
    bool check_sample_consistent(const std::string& bgen_name,
                                 const Context& context);


    typedef std::vector<std::vector<double>> Data;


    std::string m_cur_file;
    void prob_to_plink_v11(uintptr_t* genotype, byte_t const* buffer,
                           byte_t const* const end, Context context,
                           const double hard_threshold);
    void prob_to_plink_v12(uintptr_t* genotype, byte_t const* buffer,
                           byte_t const* const end, Context context);
    inline void prob_to_plink(uintptr_t* genotype, std::vector<byte_t> buffer,
                              Context context)
    {
        if ((context.flags & e_Layout) == e_Layout0
            || (context.flags & e_Layout) == e_Layout1)
        {
            prob_to_plink_v11(genotype, &(*buffer)[0],
                              &(*buffer)[0] + buffer->size(), context,
                              m_hard_threshold);
        }
        else
        {
            prob_to_plink_v12(genotype, &(*buffer)[0],
                              &(*buffer)[0] + buffer->size(), context,
                              m_hard_threshold);
        }
    }
    inline void load_raw(uintptr_t* genotype, const std::streampos byte_pos,
                         const std::string& file_name)
    {

        if (m_cur_file.empty() || file_name.compare(m_cur_file) != 0
            || !m_bgen_file.is_open())
        {
            if (m_bgen_file.is_open()) m_bgen_file.close();
            std::string bgen_name = file_name + ".bgen";
            m_bgen_file.open(bgen_name.c_str(), std::ifstream::binary);
            if (!m_bgen_file.is_open()) {
                std::string error_message =
                    "ERROR: Cannot open bgen file: " + file_name;
                throw std::runtime_error(error_message);
            }
            m_cur_file = file_name;
        }
        auto&& context = m_context_map[file_name];
        std::vector<byte_t> buffer1, buffer2;
        m_bgen_file.seekg(byte_pos, std::ios_base::beg);
        read_genotype_data_block(m_bgen_file, context, buffer1);
        uncompress_probability_data(context, buffer1, buffer2);
        prob_to_plink(genotype, buffer2, context);
    };

    inline void read_genotype(uintptr_t* genotype, const SNP& snp,
                              const std::string& file_name)
    {
        // the bgen library seems over complicated
        // try to use PLINK one. The difficulty is ifstream vs FILE
        // this is for the LD calculation
        uintptr_t final_mask = get_final_mask(m_founder_ct);
        // std::fill(m_tmp_genotype.begin(), m_tmp_genotype.end(), 0);
        // std::memset(m_tmp_genotype, 0x0, m_unfiltered_sample_ctl * 2 *
        // sizeof(uintptr_t));
        if (load_and_collapse_incl(snp.byte_pos(), file_name,
                                   m_unfiltered_sample_ct, m_founder_ct,
                                   m_founder_info.data(), final_mask, false,
                                   m_tmp_genotype.data(), genotype))
        {
            throw std::runtime_error("ERROR: Cannot read the bgen file!");
        }
    };

    // borrowed from plink
    uint32_t load_and_collapse_incl(const std::streampos byte_pos,
                                    const std::string& file_name,
                                    uint32_t unfiltered_sample_ct,
                                    uint32_t sample_ct,
                                    const uintptr_t* __restrict sample_include,
                                    uintptr_t final_mask, uint32_t do_reverse,
                                    uintptr_t* __restrict rawbuf,
                                    uintptr_t* __restrict mainbuf)
    {
        assert(unfiltered_sample_ct);
        if (unfiltered_sample_ct == sample_ct) {
            rawbuf = mainbuf;
        }
        load_raw(rawbuf, byte_pos, file_name);

        if (unfiltered_sample_ct != sample_ct) {
            copy_quaterarr_nonempty_subset(rawbuf, sample_include,
                                           unfiltered_sample_ct, sample_ct,
                                           mainbuf);
        }
        else
        {
            mainbuf[(unfiltered_sample_ct - 1) / BITCT2] &= final_mask;
        }
        if (do_reverse) {
            reverse_loadbuf(sample_ct, (unsigned char*) mainbuf);
        }

        // mainbuf should contains the information
        return 0;
    }


    void read_score(std::vector<Sample_lite>& current_prs_score,
                    size_t start_index, size_t end_bound,
                    const size_t region_index);
    void hard_code_score(std::vector<Sample_lite>& current_prs_score,
                         size_t start_index, size_t end_bound,
                         const size_t region_index);
    void dosage_score(std::vector<Sample_lite>& current_prs_score,
                      size_t start_index, size_t end_bound,
                      const size_t region_index);


    std::ifstream m_bgen_file;


    /** DON'T TOUCH      */
    // For our use case, we might be able to directly get the PRS
    // or the expected value without getting the whole vector
    // which I imagine can speed up the bgen read rather quickly
    // however, that'd = rewriting my own parsing

    struct GenotypeDataBlock
    {
        uint32_t numberOfSamples;
        uint16_t numberOfAlleles;
        byte_t ploidyExtent[2];
        byte_t const* ploidy; // Must contain at least N bytes.
        bool phased;
        byte_t bits;
        byte_t const* buffer;
        byte_t const* end;
    };
    GenotypeDataBlock init_genoData(Context const& context,
                                    byte_t const* buffer,
                                    byte_t const* const end);

    uint32_t n_choose_k(uint32_t n, uint32_t k)
    {
        if (k == 0) {
            return 1;
        }
        else if (k == 1)
        {
            return n;
        }
        return (n * n_choose_k(n - 1, k - 1)) / k;
    }
    byte_t const* read_bits_from_buffer(byte_t const* buffer,
                                        byte_t const* const end, uint64_t* data,
                                        int* size, uint8_t const bits)
    {
        assert(bits <= 64 - 8);
        while ((*size) < bits && buffer < end) {
            (*data) |= uint64_t(*(reinterpret_cast<byte_t const*>(buffer++)))
                       << (*size);
            (*size) += 8;
        }
        if ((*size) < bits) {
            throw std::runtime_error(
                "ERROR: BGEN format error! Invalid block size");
        }
        return buffer;
    }

    double parse_bit_representation(uint64_t* data, int* size, int const bits)
    {
        assert(bits <= 32);
        uint64_t bitMask = (0xFFFFFFFFFFFFFFFF >> (64 - bits));
        double const result = (*data & bitMask) / double(bitMask);
        (*size) -= bits;
        (*data) >>= bits;
        return result;
    }
    bool filter_snp(std::vector<byte_t> buffer, Context context,
                    const double geno, const double maf,
                    const double info_score, const double hard_threshold,
                    const bool hard_coded);
    bool filter_snp_v11(byte_t const* buffer, byte_t const* const end,
                        Context context, const double geno, const double maf,
                        const double info_score,
                        const double hard_code_threshold,
                        const bool hard_coded);
    bool filter_snp_v12(byte_t const* buffer, byte_t const* const end,
                        Context context, const double geno, const double maf,
                        const double info_score,
                        const double hard_code_threshold,
                        const bool hard_coded);
    double get_probability_conversion_factor(uint32_t flags)
    {
        uint32_t layout = flags & e_Layout;
        if (layout == e_Layout0) {
            // v1.0-style blocks, deprecated
            return 10000.0;
        }
        else if (layout == e_Layout1)
        {
            // v1.1-style blocks
            return 32768.0;
        }
        else
        {
            // v1.2 style (or other) blocks, these are treated
            // differently and this function does not apply.
            assert(0);
        }
        return -1;
    }


    template <typename IntegerType>
    void read_little_endian_integer(std::istream& in_stream,
                                    IntegerType* integer_ptr)
    {
        byte_t buffer[sizeof(IntegerType)];
        in_stream.read(reinterpret_cast<char*>(buffer), sizeof(IntegerType));
        if (!in_stream) {
            throw std::runtime_error("ERROR: Unable to read bgen file!");
        }
        read_little_endian_integer(buffer, buffer + sizeof(IntegerType),
                                   integer_ptr);
    }

    template <typename IntegerType>
    byte_t const* read_little_endian_integer(byte_t const* buffer,
                                             byte_t const* const end,
                                             IntegerType* integer_ptr)
    {
        assert(end >= buffer + sizeof(IntegerType));
        *integer_ptr = 0;
        for (std::size_t byte_i = 0; byte_i < sizeof(IntegerType); ++byte_i) {
            (*integer_ptr) |=
                IntegerType(*reinterpret_cast<byte_t const*>(buffer++))
                << (8 * byte_i);
        }
        return buffer;
    }

    template <typename T>
    void zlib_uncompress(byte_t const* begin, byte_t const* const end,
                         std::vector<T>* dest)
    {
        uLongf const source_size = (end - begin);
        uLongf dest_size = dest->size() * sizeof(T);
        uncompress(reinterpret_cast<Bytef*>(&dest->operator[](0)), &dest_size,
                   reinterpret_cast<Bytef const*>(begin), source_size);
        assert(result == Z_OK);
        assert(dest_size % sizeof(T) == 0);
        dest->resize(dest_size / sizeof(T));
    }

    // Uncompress the given data, symmetric with zlib_compress.
    // The destination must be large enough to fit the uncompressed data,
    // and it will be resized to exactly fit the uncompressed data.
    template <typename T>
    void zlib_uncompress(std::vector<byte_t> const& source,
                         std::vector<T>* dest)
    {
        byte_t const* begin = &source[0];
        byte_t const* const end = &source[0] + source.size();
        zlib_uncompress(begin, end, dest);
    }

    void uncompress_probability_data(Context const& context,
                                     std::vector<byte_t> const& compressed_data,
                                     std::vector<byte_t>* buffer)
    {
        // compressed_data contains the (compressed or uncompressed) probability
        // data.
        uint32_t const compressionType =
            (context.flags & e_CompressedSNPBlocks);
        if (compressionType != e_NoCompression) {
            byte_t const* begin = &compressed_data[0];
            byte_t const* const end =
                &compressed_data[0] + compressed_data.size();
            uint32_t uncompressed_data_size = 0;
            if ((context.flags & e_Layout) == e_Layout1) {
                uncompressed_data_size = 6 * context.number_of_samples;
            }
            else
            {
                begin = read_little_endian_integer(begin, end,
                                                   &uncompressed_data_size);
            }
            buffer->resize(uncompressed_data_size);
            if (compressionType == e_ZlibCompression) {
                zlib_uncompress(begin, end, buffer);
            }
            else if (compressionType == e_ZstdCompression)
            {
                throw std::runtime_error(
                    "ERROR: zstd compression currently not supported");
                // zstd_uncompress( begin, end, buffer ) ;
            }
            assert(buffer->size() == uncompressed_data_size);
        }
        else
        {
            // copy the data between buffers.
            buffer->assign(compressed_data.begin(), compressed_data.end());
        }
    }

    void read_genotype_data_block(std::istream& aStream, Context const& context,
                                  std::vector<byte_t>* buffer)
    {
        uint32_t payload_size = 0;
        if ((context.flags & e_Layout) == e_Layout2
            || ((context.flags & e_CompressedSNPBlocks) != e_NoCompression))
        {
            read_little_endian_integer(aStream, &payload_size);
        }
        else
        {
            payload_size = 6 * context.number_of_samples;
        }
        buffer->resize(payload_size);
        aStream.read(reinterpret_cast<char*>(&(*buffer)[0]), payload_size);
    }

    bool read_snp_identifying_data(std::istream& aStream,
                                   Context const& context, std::string* SNPID,
                                   std::string* RSID, std::string* chromosome,
                                   uint32_t* SNP_position,
                                   std::vector<std::string>& alleles)
    {
        uint16_t SNPID_size = 0;
        uint16_t RSID_size = 0;
        uint16_t numberOfAlleles = 0;
        uint16_t chromosome_size = 0;
        uint32_t allele_size = 0;
        std::string allele;
        uint32_t const layout = context.flags & e_Layout;

        // If we can't read a valid first field we return false; this will
        // indicate EOF. Any other fail to read is an error and an exception
        // will be thrown.
        if (layout == e_Layout1 || layout == e_Layout0) {
            uint32_t number_of_samples;
            read_little_endian_integer(aStream, &number_of_samples);
            if (number_of_samples != context.number_of_samples) {
                throw std::runtime_error(
                    "ERROR: Mismatched number of samples!");
            }
            read_length_followed_by_data(aStream, &SNPID_size, SNPID);
        }
        else if (layout == e_Layout2)
        {
            read_length_followed_by_data(aStream, &SNPID_size, SNPID);
        }
        else
        {
            assert(0);
        }

        read_length_followed_by_data(aStream, &RSID_size, RSID);
        read_length_followed_by_data(aStream, &chromosome_size, chromosome);
        read_little_endian_integer(aStream, SNP_position);
        if (layout == e_Layout2) {
            read_little_endian_integer(aStream, &numberOfAlleles);
        }
        else
        {
            numberOfAlleles = 2;
        }
        alleles.reserve(numberOfAlleles);
        for (uint16_t i = 0; i < numberOfAlleles; ++i) {
            read_length_followed_by_data(aStream, &allele_size, &allele);
            std::transform(allele.begin(), allele.end(), allele.begin(),
                           ::toupper);
            alleles.push_back(allele);
        }
        if (!aStream) {
            throw std::runtime_error("ERROR: Unable to read bgen file!");
        }
        return true;
    }


    template <typename IntegerType>
    void read_length_followed_by_data(std::istream& in_stream,
                                      IntegerType* length_ptr,
                                      std::string* string_ptr)
    {
        IntegerType& length = *length_ptr;
        read_little_endian_integer(in_stream, length_ptr);
        std::vector<char> buffer(length);
        in_stream.read(&buffer[0], length);
        if (!in_stream) {
            throw std::runtime_error("ERROR: Unable to read bgen file!");
        }
        string_ptr->assign(buffer.begin(), buffer.end());
    }
    //////////////////////
    void read_and_parse_genotype_data_block(std::istream& aStream,
                                            Context context,
                                            std::vector<byte_t>* buffer1,
                                            std::vector<byte_t>* buffer2,
                                            bool quick)
    {
        read_genotype_data_block(aStream, context, buffer1);
        if (quick) return;
        uncompress_probability_data(context, *buffer1, buffer2);
        parse_probability_data(&(*buffer2)[0], &(*buffer2)[0] + buffer2->size(),
                               context, setter);
    }
};

#endif
