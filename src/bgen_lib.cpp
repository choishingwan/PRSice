#include "bgen_lib.hpp"


//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifdef CHAR_BIT
#if (CHAR_BIT != 8)
#error CHAR_BIT "Sorry, this implementation assumes 8-bit bytes. It won't work on your platform"
#endif
#endif

namespace genfile
{
namespace bgen
{
#if DEBUG_BGEN_FORMAT
    namespace impl
    {
        std::string to_hex(std::string const& str)
        {
            std::ostringstream o;
            for (std::size_t i = 0; i < str.size(); ++i)
            {
                if (i % 4 == 0) o << "|";
                o << std::hex << std::setw(2) << std::setfill('0')
                  << static_cast<int>(static_cast<unsigned char>(str[i]));
            }
            return o.str();
        }
    }
#endif

    Context::Context()
        : number_of_samples(0)
        , number_of_variants(0)
        , magic("bgen")
        , free_data("")
        , flags(0)
    {
    }

    Context::Context(Context const& other)
        : number_of_samples(other.number_of_samples)
        , number_of_variants(other.number_of_variants)
        , magic(other.magic)
        , free_data(other.free_data)
        , flags(other.flags)
    {
    }

    Context& Context::operator=(Context const& other)
    {
        number_of_samples = other.number_of_samples;
        number_of_variants = other.number_of_variants;
        magic = other.magic;
        free_data = other.free_data;
        flags = other.flags;
        return *this;
    }

    uint32_t Context::header_size() const { return free_data.size() + 20; }

    void read_offset(std::istream& iStream, uint32_t* offset)
    {
        read_little_endian_integer(iStream, offset);
    }

    void write_offset(std::ostream& oStream, uint32_t const offset)
    {
        write_little_endian_integer(oStream, offset);
    }

    std::size_t read_header_block(std::istream& aStream, Context* context)
    {
        assert(context != 0);
        uint32_t header_size = 0, number_of_snp_blocks = 0,
                 number_of_samples = 0, flags = 0;

        char magic[4];
        std::size_t fixed_data_size = 20;
        std::vector<char> free_data;

        read_little_endian_integer(aStream, &header_size);
        assert(header_size >= fixed_data_size);
        read_little_endian_integer(aStream, &number_of_snp_blocks);
        read_little_endian_integer(aStream, &number_of_samples);
        aStream.read(&magic[0], 4);
        free_data.resize(header_size - fixed_data_size);
        aStream.read(&free_data[0], free_data.size());
        read_little_endian_integer(aStream, &flags);

        if ((magic[0] != 'b' || magic[1] != 'g' || magic[2] != 'e'
             || magic[3] != 'n')
            && (magic[0] != 0 || magic[1] != 0 || magic[2] != 0
                || magic[3] != 0))
        { throw BGenError(); } if (aStream)
        {
            context->number_of_samples = number_of_samples;
            context->number_of_variants = number_of_snp_blocks;
            context->magic.assign(&magic[0], &magic[0] + 4);
            context->free_data.assign(free_data.begin(), free_data.end());
            context->flags = flags;

            return (header_size);
        }
        else
        {
            throw BGenError();
        }
    }

    void write_header_block(std::ostream& aStream, Context const& context)
    {
        uint32_t header_size = context.header_size();
        write_little_endian_integer(aStream, header_size);
        write_little_endian_integer(aStream, context.number_of_variants);
        write_little_endian_integer(aStream, context.number_of_samples);
        aStream.write(context.magic.data(), 4);
        aStream.write(context.free_data.data(), context.free_data.size());
        write_little_endian_integer(aStream, context.flags);
    }

    std::size_t
    write_sample_identifier_block(std::ostream& aStream, Context const& context,
                                  std::vector<std::string> const& sample_ids)
    {
        assert(sample_ids.size() == context.number_of_samples);
        uint32_t block_size = 8;
        for (uint32_t i = 0; i < sample_ids.size(); ++i)
        { block_size += 2 + sample_ids[i].size(); }
        write_little_endian_integer(aStream, block_size);
        write_little_endian_integer(aStream, context.number_of_samples);
#if DEBUG_BGEN_FORMAT
        std::cerr << "genfile::bgen::write_sample_identifier_block(): writing "
                  << sample_ids.size() << " samples...\n";
#endif
        for (uint32_t i = 0; i < sample_ids.size(); ++i)
        {
#if DEBUG_BGEN_FORMAT
            std::cerr << "genfile::bgen::write_sample_identifier_block(): "
                         "writing sample "
                      << sample_ids[i] << ".\n";
#endif
            std::string const& identifier = sample_ids[i];
            assert(identifier.size()
                   <= std::size_t(std::numeric_limits<uint16_t>::max()));
            uint16_t const id_size = uint16_t(identifier.size());
            write_length_followed_by_data(aStream, id_size, identifier);
        }
        return block_size;
    }


    namespace impl
    {
        void check_for_two_alleles(uint16_t numberOfAlleles)
        {
            if (numberOfAlleles != 2)
            {
                std::cerr << "genfile::bgen::impl::check_for_two_alleles: only "
                             "biallelic variants are currently supported.\n";
                assert(0);
            }
        }

        struct TwoAlleleSetter
        {
            TwoAlleleSetter(std::string* allele1, std::string* allele2)
                : m_allele1(allele1), m_allele2(allele2)
            {
                assert(allele1 != 0);
                assert(allele2 != 0);
            }

            void operator()(uint16_t i, std::string const& value)
            {
                if (i == 0) { *m_allele1 = value; }
                else if (i == 1)
                {
                    *m_allele2 = value;
                }
                else
                {
                    assert(0);
                }
            }
            std::string* m_allele1;
            std::string* m_allele2;
        };
    }

    namespace v10
    {
        bool read_snp_identifying_data(
            std::istream& aStream, Context const& context, std::string* SNPID,
            std::string* RSID, std::string* chromosome, uint32_t* SNP_position,
            std::string* first_allele, std::string* second_allele)
        {
            // v1.0-style layout, deprecated
            uint32_t number_of_samples = 0;
            unsigned char max_id_size = 0;
            unsigned char SNPID_size = 0;
            unsigned char RSID_size = 0;
            if (aStream)
            {
                read_little_endian_integer(aStream, &number_of_samples);
                if (!aStream) { return false; }
                if (number_of_samples != context.number_of_samples)
                { throw BGenError(); } }
            if (aStream) { read_little_endian_integer(aStream, &max_id_size); }
            if (aStream)
            {
                read_length_followed_by_data(aStream, &SNPID_size, SNPID);
                assert(SNPID_size <= max_id_size);
                aStream.ignore(max_id_size - SNPID_size);
            }
            if (aStream)
            {
                read_length_followed_by_data(aStream, &RSID_size, RSID);
                assert(RSID_size <= max_id_size);
                aStream.ignore(max_id_size - RSID_size);
            }
            if (aStream)
            {
                unsigned char chromosome_char = 0;
                read_little_endian_integer(aStream, &chromosome_char);
                read_little_endian_integer(aStream, SNP_position);

                switch (chromosome_char)
                {
                case 1: *chromosome = "01"; break;
                case 2: *chromosome = "02"; break;
                case 3: *chromosome = "03"; break;
                case 4: *chromosome = "04"; break;
                case 5: *chromosome = "05"; break;
                case 6: *chromosome = "06"; break;
                case 7: *chromosome = "07"; break;
                case 8: *chromosome = "08"; break;
                case 9: *chromosome = "09"; break;
                case 10: *chromosome = "10"; break;
                case 11: *chromosome = "11"; break;
                case 12: *chromosome = "12"; break;
                case 13: *chromosome = "13"; break;
                case 14: *chromosome = "14"; break;
                case 15: *chromosome = "15"; break;
                case 16: *chromosome = "16"; break;
                case 17: *chromosome = "17"; break;
                case 18: *chromosome = "18"; break;
                case 19: *chromosome = "19"; break;
                case 20: *chromosome = "20"; break;
                case 21: *chromosome = "21"; break;
                case 22: *chromosome = "22"; break;
                case 23: *chromosome = "0X"; break;
                case 24: *chromosome = "0Y"; break;
                case 253: *chromosome = "XY"; break;
                case 254: *chromosome = "MT"; break;
                default: *chromosome = "NA";
                }
                unsigned char allele1_size = 0, allele2_size = 0;
                read_length_followed_by_data(aStream, &allele1_size,
                                             first_allele);
                read_length_followed_by_data(aStream, &allele2_size,
                                             second_allele);
            }
            return true;
        }
    }

    bool read_snp_identifying_data(std::istream& aStream,
                                   Context const& context, std::string* SNPID,
                                   std::string* RSID, std::string* chromosome,
                                   uint32_t* SNP_position,
                                   std::string* first_allele,
                                   std::string* second_allele)
    {
#if DEBUG_BGEN_FORMAT
        std::cerr
            << "genfile::bgen::impl::read_snp_identifying_data(): flags = 0x"
            << std::hex << context.flags << ".\n";
#endif
        uint32_t const layout = context.flags & e_Layout;
        if (layout == e_Layout1 || layout == e_Layout2)
        {
            // forward to v12 version which handles multiple alleles.
            impl::TwoAlleleSetter allele_setter(first_allele, second_allele);
            return bgen::read_snp_identifying_data(
                aStream, context, SNPID, RSID, chromosome, SNP_position,
                &impl::check_for_two_alleles, allele_setter);
        }
        else if (layout == e_Layout0)
        {
            return v10::read_snp_identifying_data(aStream, context, SNPID, RSID,
                                                  chromosome, SNP_position,
                                                  first_allele, second_allele);
        }
        else
        {
            assert(0);
        }
        return true;
    }

    namespace v11
    {
        namespace impl
        {
            double get_probability_conversion_factor(uint32_t flags)
            {
                uint32_t layout = flags & e_Layout;
                if (layout == e_Layout0)
                {
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
        }
    }

    void ignore_genotype_data_block(std::istream& aStream,
                                    Context const& context)
    {
        if ((context.flags & bgen::e_CompressedSNPBlocks) != e_NoCompression)
        {
            uint32_t compressed_data_size = 0;
            read_little_endian_integer(aStream, &compressed_data_size);
            if (compressed_data_size > 0)
            {
                // gcc std::istream::ignore() has a bug / feature in which
                // it peeks at the next char and sets eof() if you ignore all
                // the bytes in the file. This breaks our expected invariant and
                // means subsequent calls to peekg() fail with -1, which stops
                // us getting file size. We deal with this by simply ignoring
                // one less character here.
                aStream.ignore(compressed_data_size - 1);
                aStream.get();
            }
        }
        else
        {
            aStream.ignore(6 * context.number_of_samples);
        }
    }

    void read_genotype_data_block(std::istream& aStream, Context const& context,
                                  std::vector<byte_t>* buffer)
    {
        uint32_t payload_size = 0;
        if ((context.flags & e_Layout) == e_Layout2
            || ((context.flags & e_CompressedSNPBlocks) != e_NoCompression))
        { read_little_endian_integer(aStream, &payload_size); } else
        {
            payload_size = 6 * context.number_of_samples;
        }
        buffer->resize(payload_size);
        aStream.read(reinterpret_cast<char*>(&(*buffer)[0]), payload_size);
    }

    void uncompress_probability_data(Context const& context,
                                     std::vector<byte_t> const& compressed_data,
                                     std::vector<byte_t>* buffer)
    {
        // compressed_data contains the (compressed or uncompressed) probability
        // data.
        uint32_t const compressionType =
            (context.flags & bgen::e_CompressedSNPBlocks);
        if (compressionType != e_NoCompression)
        {
            byte_t const* begin = &compressed_data[0];
            byte_t const* const end =
                &compressed_data[0] + compressed_data.size();
            uint32_t uncompressed_data_size = 0;
            if ((context.flags & e_Layout) == e_Layout1)
            { uncompressed_data_size = 6 * context.number_of_samples; } else
            {
                begin = read_little_endian_integer(begin, end,
                                                   &uncompressed_data_size);
            }
            buffer->resize(uncompressed_data_size);
            if (compressionType == e_ZlibCompression)
            { zlib_uncompress(begin, end, buffer); }
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

    namespace v12
    {
        namespace impl
        {
            // Fill a data field, encoded as a 64-bit integer, with bytes
            // from a buffer, until the data contains at least the given number
            // of bits. (For this to work we require in general that bits <= 56
            // (on an 8-bit byte machine).
            byte_t const* read_bits_from_buffer(byte_t const* buffer,
                                                byte_t const* const end,
                                                uint64_t* data, int* size,
                                                uint8_t const bits)
            {
                assert(bits <= 64 - 8);
                while ((*size) < bits && buffer < end)
                {
                    (*data) |=
                        uint64_t(*(reinterpret_cast<byte_t const*>(buffer++)))
                        << (*size);
                    (*size) += 8;
                }
                if ((*size) < bits) { throw BGenError(); }
                return buffer;
            }

            // Consume the given number of bits from the the least significant
            // end of a data field (encoded as a uint64_t) and interpret them as
            // a floating-point number in the range 0...1
            double parse_bit_representation(uint64_t* data, int* size,
                                            int const bits)
            {
                assert(bits <= 32);
                uint64_t bitMask = (0xFFFFFFFFFFFFFFFF >> (64 - bits));
                double const result = (*data & bitMask) / double(bitMask);
                (*size) -= bits;
                (*data) >>= bits;
                return result;
            }

            namespace
            {
                struct CompareFractionalPart
                {
                    CompareFractionalPart(double* v, std::size_t n)
                        : m_v(v), m_n(n)
                    {
                    }
                    CompareFractionalPart(CompareFractionalPart const& other)
                        : m_v(other.m_v), m_n(other.m_n)
                    {
                    }
                    CompareFractionalPart&
                    operator=(CompareFractionalPart const& other)
                    {
                        m_v = other.m_v;
                        m_n = other.m_n;
                        return *this;
                    }

                    bool operator()(std::size_t a, std::size_t b) const
                    {
                        return (m_v[a] - std::floor(m_v[a]))
                               > (m_v[b] - std::floor(m_v[b]));
                    }

                private:
                    double* m_v;
                    std::size_t m_n;
                };
            }

            void compute_approximate_probabilities(double* p,
                                                   std::size_t* index,
                                                   std::size_t const n,
                                                   int const number_of_bits)
            {
                double const scale =
                    (0xFFFFFFFFFFFFFFFF >> (64 - number_of_bits));
                double total_fractional_part = 0.0;
                double sum = 0.0;
                for (std::size_t i = 0; i < n; ++i)
                {
                    p[i] *= scale;
                    sum += p[i];
                    index[i] = i;
                    total_fractional_part += p[i] - std::floor(p[i]);
                }
                // Suppose the n input numbers sum to 1
                // Each has rounding error of at most machine epsilon
                // the sum thus has maximum error of n * epsilon * scale.
                assert(sum
                       < (scale
                          * (1 + n * std::numeric_limits<double>::epsilon())));

                // We have n numbers which sum to scale ± delta,
                // where delta is the above rounding error < n*epsilon.
                // sum_i floor(p_i) is an integer by definition.
                // Total fractional part is therefore of the form r ± delta
                // where r is an integer. Since scale = sum_i floor(p_i) + r,
                // rounding up r of the p_i's yields a set of integers summing
                // to scale.
                std::size_t const r = std::floor(total_fractional_part + 0.5);
                std::sort(index, index + n, CompareFractionalPart(p, n));

                for (std::size_t i = 0; i < r; ++i)
                { p[index[i]] = std::ceil(p[index[i]]); }
                for (std::size_t i = r; i < n; ++i)
                { p[index[i]] = std::floor(p[index[i]]); } }

            byte_t* write_scaled_probs(uint64_t* data, std::size_t* offset,
                                       double const* probs, std::size_t const n,
                                       int const number_of_bits,
                                       byte_t* destination, byte_t* const end)
            {
                for (std::size_t i = 0; i < (n - 1); ++i)
                {
                    uint64_t const storedValue = uint64_t(probs[i]);
                    *data |= storedValue << (*offset);
                    (*offset) += number_of_bits;
                    if ((*offset) >= 32)
                    {
                        assert((destination + 4) <= end);
                        destination =
                            std::copy(reinterpret_cast<byte_t const*>(data),
                                      reinterpret_cast<byte_t const*>(data) + 4,
                                      destination);
                        (*offset) -= 32;
                        (*data) >>= 32;
                    }
                }
                return destination;
            }
        }
    }
}
}
