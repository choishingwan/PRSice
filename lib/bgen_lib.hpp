
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#ifndef BGEN_REFERENCE_IMPLEMENTATION_HPP
#define BGEN_REFERENCE_IMPLEMENTATION_HPP

#include "memoryread.hpp"
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdint.h>
#include <vector>
#include <zlib.h>
//#include "genfile/snp_data_utils.hpp"
//#include "genfile/get_set.hpp"


/*
 * This file contains a reference implementation of the BGEN file format
 * specification described at:
 * http://www.well.ox.ac.uk/~gav/bgen_format/bgen_format.html
 *
 */

// #define DEBUG_BGEN_FORMAT 1

#if DEBUG_BGEN_FORMAT
#include <iomanip>
#include <iostream>
#include <sstream>
#endif

///////////////////////////////////////////////////////////////////////////////////////////
// INTERFACE
///////////////////////////////////////////////////////////////////////////////////////////

namespace genfile
{
struct MissingValue
{
    bool operator<(MissingValue const& other) const;
    bool operator<=(MissingValue const& other) const;
    bool operator==(MissingValue const& other) const;
};

enum OrderType
{
    eUnknownOrderType = 0,
    eUnorderedList = 1, // a list, treated as unordered
    eOrderedList = 2,   // a list, treated as ordered
    ePerUnorderedGenotype =
        3, // E.g. genotype probabilities; GEN, BGEN v1.1, etc.
    ePerOrderedHaplotype = 4,   // E.g. Phased GT, SHAPEIT or IMPUTE haplotypes
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

typedef uint8_t byte_t;


// Compress the given data into the given destination buffer.  The destination
// will be resized to fit the compressed data.  (Since the capacity of dest may
// be larger than its size, to save memory you may need to copy the contents of
// dest elsewhere after calling this function).
//
// If offset is nonzero, compressed data will be written starting at position
// [offset]. The first [offset] bytes will be untouched.
inline void zlib_compress(byte_t const* buffer, byte_t const* const end,
                          std::vector<byte_t>* dest,
                          std::size_t const offset = 0,
                          int const compressionLevel = Z_BEST_COMPRESSION)
{
    assert(dest != 0);
    assert(compressionLevel >= 0 && compressionLevel <= Z_BEST_COMPRESSION);
    uLongf const source_size = (end - buffer);
    uLongf compressed_size = compressBound(source_size);
    dest->resize(compressed_size + offset);
    int result =
        compress2(reinterpret_cast<Bytef*>(
                      const_cast<uint8_t*>(&(dest->operator[](0)) + offset)),
                  &compressed_size, reinterpret_cast<Bytef const*>(buffer),
                  source_size, compressionLevel);
    assert(result == Z_OK);
    dest->resize(compressed_size + offset);
}

/*
    void zstd_compress(
        byte_t const* buffer,
        byte_t const* const end,
        std::vector< byte_t >* dest,
        std::size_t const offset = 0,
        int const compressionLevel = 22
    ) ;
*/
// Compress the given data into the given destination buffer.  The destination
// will be resized to fit the compressed data.  (Since the capacity of dest may
// be larger than its size, to save memory you may need to copy the contents of
// dest elsewhere after calling this function).
template <typename T>
void zlib_compress(std::vector<T> const& source, std::vector<byte_t>* dest,
                   int const compressionLevel = Z_BEST_COMPRESSION)
{
    byte_t const* begin = reinterpret_cast<byte_t const*>(&source[0]);
    byte_t const* const end =
        reinterpret_cast<byte_t const*>(&source[0] + source.size());
    return zlib_compress(begin, end, dest, 0, compressionLevel);
}

template <typename T>
void zlib_uncompress(byte_t const* begin, byte_t const* const end,
                     std::vector<T>* dest)
{
    uLongf const source_size = (end - begin);
    uLongf dest_size = dest->size() * sizeof(T);
    int result =
        uncompress(reinterpret_cast<Bytef*>(&dest->operator[](0)), &dest_size,
                   reinterpret_cast<Bytef const*>(begin), source_size);
    assert(result == Z_OK);
    assert(dest_size % sizeof(T) == 0);
    dest->resize(dest_size / sizeof(T));
}
/*
    template< typename T >
    void zstd_uncompress( byte_t const* begin, byte_t const* const end,
   std::vector< T >* dest ) { std::size_t const source_size = ( end - begin ) ;
        std::size_t const dest_size = dest->size() * sizeof( T ) ;
        std::size_t const uncompressed_size = ZSTD_getDecompressedSize(
   reinterpret_cast< void const* >( begin ), source_size ) ; std::size_t const
   result = ZSTD_decompress( reinterpret_cast< void* >( &dest->operator[]( 0 )
   ), dest_size, reinterpret_cast< void const* >( begin ), source_size ) ;
        assert( result == uncompressed_size ) ;
        dest->resize( dest_size / sizeof( T )) ;
    }
*/
// Uncompress the given data, symmetric with zlib_compress.
// The destination must be large enough to fit the uncompressed data,
// and it will be resized to exactly fit the uncompressed data.
template <typename T>
void zlib_uncompress(std::vector<byte_t> const& source, std::vector<T>* dest)
{
    byte_t const* begin = &source[0];
    byte_t const* const end = &source[0] + source.size();
    zlib_uncompress(begin, end, dest);
}
std::ostream& operator<<(std::ostream& o, MissingValue const& v);

namespace bgen
{
#if DEBUG_BGEN_FORMAT
    namespace impl
    {
        std::string to_hex(std::string const& str);

        template <typename I, typename I2>
        std::string to_hex(I i, I2 end_i)
        {
            std::ostringstream o;
            for (std::size_t count = 0; i < end_i; ++i, ++count)
            {
                if (count % 4 == 0) o << "|";
                o << std::hex << std::setw(2) << std::setfill('0')
                  << static_cast<int>(static_cast<unsigned char>(*i));
            }
            return o.str();
        }
    }
#endif

    namespace impl
    {
        // n choose k implementation
        // a faster implementation is of course possible, (e.g. table lookup)
        // but this is not a bottleneck.
        template <typename Integer>
        Integer n_choose_k(Integer n, Integer k)
        {
            if (k == 0) { return 1; }
            else if (k == 1)
            {
                return n;
            }
            return (n * n_choose_k(n - 1, k - 1)) / k;
        }
    }

    // class thrown when errors are detected
    struct BGenError : public virtual std::exception
    {
        ~BGenError() throw() {}
        char const* what() const throw() { return "BGenError"; }
    };

    // integer types
    typedef ::uint32_t uint32_t;
    typedef ::uint16_t uint16_t;

    // Header flag definitions
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

    // Structure containing information from the header block.
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

    // Read the offset from the start of the stream.
    void read_offset(std::istream& iStream, uint32_t* offset);
    // Write an offset value to the stream.
    void write_offset(std::ostream& oStream, uint32_t const offset);

    // Read a header block from the supplied stream,
    // filling the fields of the supplied context object.
    // Return the number of bytes read.
    std::size_t read_header_block(std::istream& aStream, Context* context);

    // Write a bgen header block to the supplied stream,
    // taking data from the fields of the supplied context object.
    void write_header_block(std::ostream& aStream, Context const& context);

    // Read a sample identifier block from the given stream.
    // The setter object passed in must be a unary function or
    // function object that takes a string.  It must be callable as
    // setter.set_value( value ) ;
    // where value is of type std::string.  It will be called once for
    // each sample identifier in the block, in the order they occur.
    template <typename SampleSetter>
    std::size_t read_sample_identifier_block(std::istream& aStream,
                                             Context const& context,
                                             SampleSetter setter);

    // Write the sample identifiers contained in the
    // given vector to the stream.
    std::size_t
    write_sample_identifier_block(std::ostream& aStream, Context const& context,
                                  std::vector<std::string> const& sample_ids);

    // Attempt to read identifying information for the next variant in the file.
    // This function will return true if SNP data was successfully read or false
    // if the initial reads met an EOF.  It will throw an BGenError if only a
    // subset of fields can be read before EOF. The two setter objects must be
    // callable as set_number_of_alleles( n ) set_allele( i, a ) where n is an
    // unsigned integer representing the number of alleles, i is an unsigned
    // integer in the range 0..n-1, and a is of type std::string representing
    // the ith allele.
    template <typename NumberOfAllelesSetter, typename AlleleSetter>
    bool read_snp_identifying_data(std::istream& aStream,
                                   Context const& context, std::string* SNPID,
                                   std::string* RSID, std::string* chromosome,
                                   uint32_t* SNP_position,
                                   NumberOfAllelesSetter set_number_of_alleles,
                                   AlleleSetter set_allele);

    // Read identifying data fields for the next variant in the file, assuming 2
    // alleles. This function forwards to the generic multi-allele version,
    // above, and will throw a BGenError if the number of alleles is different
    // than 2.
    bool read_snp_identifying_data(std::istream& aStream,
                                   Context const& context, std::string* SNPID,
                                   std::string* RSID, std::string* chromosome,
                                   uint32_t* SNP_position,
                                   std::string* first_allele,
                                   std::string* second_allele);

    // Write identifying data fields for the given variant.
    // The buffer will be resized to fit.
    // Return a pointer to past-the-end of the data written.
    template <typename AlleleGetter>
    byte_t* write_snp_identifying_data(
        std::vector<byte_t>* buffer, Context const& context, std::string SNPID,
        std::string RSID, std::string chromosome, uint32_t position,
        uint16_t const number_of_alleles, AlleleGetter get_allele);

    // Ignore (and seek forward past) the genotype data block contained in the
    // given stream. (The main purpose of this function is to encapsulate the
    // slightly complicated rules for reading or computing the compressed data
    // size needed to skip the right number of bytes).
    void ignore_genotype_data_block(std::istream& aStream,
                                    Context const& context);

    // Low-level function which reads raw probability data from a genotype data
    // block contained in the input stream into a supplied buffer. The buffer
    // will be resized to fit the data (incurring an allocation if the buffer is
    // not large enough.) (The main purpose of this function is to encapsulate
    // the slightly complicated rules for reading or computing the compressed
    // data size.  Where applicable this function first reads the four bytes
    // indicating the compressed data size; these are discarded and do not
    // appear in the buffer).
    void read_genotype_data_block(std::istream& aStream, Context const& context,
                                  std::vector<byte_t>* buffer1);
    void read_genotype_data_block(FileRead& aStream,
                                  const std::string& file_name,
                                  Context const& context,
                                  std::vector<byte_t>* buffer1,
                                  const std::streampos idx);

    // Low-level function which uncompresses probability data stored in the
    // genotype data block contained in the first buffer into a second buffer
    // (or just copies it over if the probability data is not compressed.) The
    // second buffer will be resized to fit the result (incurring an allocation
    // if the buffer is not large enough). (The main purpose of this function is
    // to encapsulate the rules for reading or computing the uncompressed data
    // size.  Where applicable this function first consumes four bytes of the
    // input buffer and interprets them as the uncompressed data size, before
    // uncompressing the rest.) Usually bgen files are stored compressed.  If
    // the data is not compressed, this function simply copies the source buffer
    // to the target buffer.
    void uncompress_probability_data(Context const& context,
                                     std::vector<byte_t> const& buffer1,
                                     std::vector<byte_t>* buffer2);

    // template< typename Setter >
    // parse uncompressed genotype probability data stored in the given buffer.
    // Values are returned as doubles or as missing values using the
    // setter object provided.  The setter must support the following
    // expressions:
    //
    // - setter.initialise( N, K ) ;
    // where N and K are convertible from std::size_t and represent the number
    // of samples and number of alleles present for the variant.
    //
    // - setter.set_min_max_ploidy( a, b, c, d ) (OPTIONAL)
    // if present this is called with a, b = minimum and maximum ploidy among
    // samples and c,d = minimum and maximum number of probabilities per sample.
    //
    // - setter.set_sample( i )
    // where i is convertible from std::size_t and represents the index of a
    // sample between 0 and N-1. This function should return a value convertible
    // to bool, representing a hint as to whether the setter wants to use the
    // values for this sample or not. If the return value is false, the
    // implementation may choose not to report any data for this sample, whence
    // the next call, if any, will be set_sample(i+1).
    //
    // - setter.set_number_of_entries( P, Z, order_type, value_type ) ;
    // - or setter.set_number_of_entries( Z, order_type, value_type ) ;
    // where
    // - P is an unsigned integer reflecting the ploidy of this sample at this
    // variant
    // - Z is an unsigned integer reflecting the number of probability values
    // comprising the data for this sample
    // order_type is of type genfile::OrderType and equal to either
    // ePerUnorderedGenotype if unphased data is stored or
    // ePerPhasedHaplotypePerAllele if phased data is stored. and value_type is
    // of type genfile::ValueType and always equal to eProbability.
    //
    // This method is called once per sample (for which set_sample is true)
    // and informs the setter the number of probability values present for this
    // sample (depending on the ploidy, the number of alleles, and whether the
    // data is phased.)  E.g. In the common case of a diploid sample at a
    // biallelic variant, this will be called with Z=3.
    //
    // Users can implement either the first or second versions of this function
    // (or both, in which case the first version will be called).
    //
    // For bgen <= 1.2, all data is interpreted as probabilities so the
    // value_type is always eProbability. For bgen <= 1.1, all data is unphased
    // so the order_type is ePerUnorderedGenotype.
    //
    // - setter.set_value( i, value ) ;
    // where value is of type double or genfile::MissingValue.
    // This is called Z times and reflects the probability data stored for this
    // sample. Samples with missing data have Z missing values; those without
    // missing data have Z double values.
    //
    // - setter.finalise()
    // If this method is present it is called once at the end of the call.
    //
    template <typename Setter>
    void parse_probability_data(byte_t const* buffer, byte_t const* const end,
                                Context const& context, Setter& setter);

    // Utility function which wraps the above steps for reading probability data
    // into a single function. Concretely this function: 1: calls
    // read_genotype_data_block(), reading probability data from the given
    // stream. 2: calls uncompress_probability_data() to uncompress the data
    // where necessary. 3: calls parse_probability_data to parse it, returning
    // values using the setter object provided. The buffers are used as
    // intermediate storage and will be resized to fit data as needed.
    template <typename Setter>
    void read_and_parse_genotype_data_block(std::istream& aStream,
                                            Context const& context,
                                            Setter& setter,
                                            std::vector<byte_t>* buffer1,
                                            std::vector<byte_t>* buffer2);
    template <typename Setter>
    void read_and_parse_genotype_data_block(
        FileRead& aStream, const std::string& file_name, Context const& context,
        Setter& setter, std::vector<byte_t>* buffer1,
        std::vector<byte_t>* buffer2, std::streampos idx);
}
}

///////////////////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION
///////////////////////////////////////////////////////////////////////////////////////////

#if (defined(__BYTE_ORDER) && __BYTE_ORDER == __LITTLE_ENDIAN)              \
    || defined(__LITTLE_ENDIAN) || defined(__ARMEL) || defined(__THUMBEL__) \
    || defined(__AARCH64EL__) || defined(_MIPSEL) || defined(__MIPSEL)      \
    || defined(__MIPSEL__)                                                  \
    || (defined(__BYTE_ORDER__)                                             \
        && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__))                     \
    || (defined(__FLOAT_WORD_ORDER__)                                       \
        && (__FLOAT_WORD_ORDER__ == __ORDER_LITTLE_ENDIAN__))
#define BGEN_BIG_ENDIAN 0
#define BGEN_LITTLE_ENDIAN 1
#elif (defined(__BYTE_ORDER) && __BYTE_ORDER == __BIG_ENDIAN)                \
    || defined(__BIG_ENDIAN) || defined(__ARMEB) || defined(__THUMBEB__)     \
    || defined(__AARCH64EB__) || defined(_MIPSEB) || defined(__MIPSEB)       \
    || defined(__MIPSEB__)                                                   \
    || (defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)) \
    || (defined(__FLOAT_WORD_ORDER__)                                        \
        && (__FLOAT_WORD_ORDER__ == __ORDER_BIG_ENDIAN__))
#define BGEN_BIG_ENDIAN 1
#define BGEN_LITTLE_ENDIAN 0
#else
#error "Unable to determine architecture endian-ness"
#endif

#if !BGEN_LITTLE_ENDIAN
#error \
    "BGEN support on big endian machines is currently untested.  Please remove this line if you want to try it."
#endif

namespace genfile
{
namespace bgen
{
    // Read an integer stored in little-endian format into an integer stored in
    // memory.
    template <typename IntegerType>
    inline byte_t const* read_little_endian_integer(byte_t const* buffer,
                                                    byte_t const* const end,
                                                    IntegerType* integer_ptr)
    {
        assert(end >= buffer + sizeof(IntegerType));
#if BGEN_LITTLE_ENDIAN
        *integer_ptr =
            IntegerType(*reinterpret_cast<IntegerType const*>(buffer));
        buffer += sizeof(IntegerType);
#elif BGEN_BIG_ENDIAN
        *integer_ptr = 0;
        for (std::size_t byte_i = 0; byte_i < sizeof(IntegerType); ++byte_i)
        {
            (*integer_ptr) |=
                IntegerType(*reinterpret_cast<byte_t const*>(buffer++))
                << (8 * byte_i);
        }
#else
#error "unknown endianness"
#endif
        return buffer;
    }

    // Read an integer stored in little-endian format into an integer stored in
    // memory. The stream is assumed to have sizeof( Integertype ) readable
    // bytes.
    template <typename IntegerType>
    void read_little_endian_integer(std::istream& in_stream,
                                    IntegerType* integer_ptr)
    {
        byte_t buffer[sizeof(IntegerType)];
        in_stream.read(reinterpret_cast<char*>(buffer), sizeof(IntegerType));
        if (!in_stream) { throw BGenError(); }
        read_little_endian_integer(buffer, buffer + sizeof(IntegerType),
                                   integer_ptr);
    }

    template <typename IntegerType>
    void read_little_endian_integer(FileRead& in_stream,
                                    const std::string& file_name,
                                    IntegerType* integer_ptr,
                                    const std::streampos idx)
    {
        byte_t buffer[sizeof(IntegerType)];
        in_stream.read(file_name, idx, sizeof(IntegerType),
                       reinterpret_cast<char*>(buffer));
        // if (!in_stream) { throw BGenError(); }
        read_little_endian_integer(buffer, buffer + sizeof(IntegerType),
                                   integer_ptr);
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
        if (!in_stream) { throw BGenError(); }
        string_ptr->assign(buffer.begin(), buffer.end());
    }

    // Write an integer to the buffer in little-endian format.
    template <typename IntegerType>
    byte_t* write_little_endian_integer(byte_t* buffer, byte_t* const end,
                                        IntegerType const integer)
    {
        assert(end >= buffer + sizeof(IntegerType));
        for (std::size_t byte_i = 0; byte_i < sizeof(IntegerType); ++byte_i)
        { *buffer++ = (integer >> (8 * byte_i)) & 0xff; }
        return buffer;
    }

    // Write data contained in a std::string to the buffer, preceded
    // by a length of the given integral type in little-endian format.
    // Return past-the-end of what was written.
    template <typename IntegerType>
    byte_t* write_length_followed_by_data(byte_t* buffer, byte_t* const end,
                                          IntegerType length,
                                          std::string const data_string)
    {
        assert(end >= buffer + length + sizeof(IntegerType));
        assert(length <= data_string.size());
        buffer = write_little_endian_integer(buffer, end, length);
        buffer = std::copy(data_string.begin(), data_string.begin() + length,
                           buffer);
        return buffer;
    }

    // Write an integer to the stream in little-endian format.
    // The stream is assumed to have sizeof( Integertype ) bytes writeable.
    template <typename IntegerType>
    void write_little_endian_integer(std::ostream& out_stream,
                                     IntegerType const integer)
    {
        byte_t buffer[sizeof(IntegerType)];
        write_little_endian_integer(buffer, buffer + sizeof(IntegerType),
                                    integer);
        out_stream.write(reinterpret_cast<char const*>(buffer),
                         sizeof(IntegerType));
    }

    // Write data containd in a std::string to the stream,
    // Preceded by a length of the given integral type in little-endian format.
    template <typename IntegerType>
    void write_length_followed_by_data(std::ostream& out_stream,
                                       IntegerType length,
                                       std::string const data_string)
    {
        assert(length <= data_string.size());
        write_little_endian_integer(out_stream, length);
        out_stream.write(data_string.data(), length);
    }

    template <typename SampleSetter>
    std::size_t read_sample_identifier_block(std::istream& aStream,
                                             Context const& context,
                                             SampleSetter setter)
    {
        uint32_t block_size = 0;
        uint32_t number_of_samples = 0;
        uint16_t identifier_size;
        std::string identifier;
        std::size_t bytes_read = 0;

        read_little_endian_integer(aStream, &block_size);
        read_little_endian_integer(aStream, &number_of_samples);
        bytes_read += 8;
        assert(number_of_samples == context.number_of_samples);

        for (uint32_t i = 0; i < number_of_samples; ++i)
        {
            read_length_followed_by_data(aStream, &identifier_size,
                                         &identifier);
            if (aStream)
            {
                bytes_read += sizeof(identifier_size) + identifier_size;
                setter(identifier);
            }
            else
            {
                throw BGenError();
            }
        }
        assert(bytes_read == block_size);
        return bytes_read;
    }

    template <typename NumberOfAllelesSetter, typename AlleleSetter>
    bool read_snp_identifying_data(std::istream& aStream,
                                   Context const& context, std::string* SNPID,
                                   std::string* RSID, std::string* chromosome,
                                   uint32_t* SNP_position,
                                   NumberOfAllelesSetter set_number_of_alleles,
                                   AlleleSetter set_allele)
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
        if (layout == e_Layout1 || layout == e_Layout0)
        {
            uint32_t number_of_samples;
            try
            {
                read_little_endian_integer(aStream, &number_of_samples);
            }
            catch (BGenError const&)
            {
                return false;
            }
            if (number_of_samples != context.number_of_samples)
            { throw BGenError(); }
            read_length_followed_by_data(aStream, &SNPID_size, SNPID);
        }
        else if (layout == e_Layout2)
        {
            try
            {
                read_length_followed_by_data(aStream, &SNPID_size, SNPID);
            }
            catch (BGenError const&)
            {
                return false;
            }
        }
        else
        {
            assert(0);
        }

        read_length_followed_by_data(aStream, &RSID_size, RSID);
        read_length_followed_by_data(aStream, &chromosome_size, chromosome);
        read_little_endian_integer(aStream, SNP_position);
        if (layout == e_Layout2)
        { read_little_endian_integer(aStream, &numberOfAlleles); }
        else
        {
            numberOfAlleles = 2;
        }
        set_number_of_alleles(numberOfAlleles);
        for (uint16_t i = 0; i < numberOfAlleles; ++i)
        {
            read_length_followed_by_data(aStream, &allele_size, &allele);
            set_allele(i, allele);
        }
        if (!aStream)
        {
#if DEBUG_BGEN_FORMAT
            std::cerr << "bgen: layout = " << layout
                      << ", alleles = " << numberOfAlleles << ".\n"
                      << std::flush;
            std::cerr << *SNPID << ", " << *RSID << ", " << *chromosome << ", "
                      << *SNP_position << ".\n"
                      << std::flush;
#endif
            throw BGenError();
        }
        return true;
    }

    namespace
    {
        // TODO: make this C++-03 compatible.
        template <typename Setter>
        struct has_set_min_max_ploidy
        {
            template <typename U,
                      void (U::*)(uint32_t, uint32_t, uint32_t, uint32_t)>
            struct SFINAE
            {
            };
            template <typename U>
            static uint8_t Test(SFINAE<U, &U::set_min_max_ploidy>*);
            template <typename U>
            static uint32_t Test(...);
            static const bool Yes = sizeof(Test<Setter>(0)) == sizeof(uint8_t);
            static const bool No = sizeof(Test<Setter>(0)) == sizeof(uint16_t);
        };

        template <bool C>
        struct tag
        {
        };

        template <typename Setter>
        void call_set_min_max_ploidy(Setter& setter, uint32_t min_ploidy,
                                     uint32_t max_ploidy,
                                     uint32_t numberOfAlleles, bool phased,
                                     tag<true> const&)
        {
            uint32_t min_count =
                phased ? (min_ploidy * numberOfAlleles)
                       : impl::n_choose_k(min_ploidy + numberOfAlleles - 1,
                                          numberOfAlleles - 1);
            uint32_t max_count =
                phased ? (max_ploidy * numberOfAlleles)
                       : impl::n_choose_k(max_ploidy + numberOfAlleles - 1,
                                          numberOfAlleles - 1);
            setter.set_min_max_ploidy(min_ploidy, max_ploidy, min_count,
                                      max_count);
        }

        template <typename Setter>
        void call_set_min_max_ploidy(Setter& setter, uint32_t min_ploidy,
                                     uint32_t max_ploidy,
                                     uint32_t numberOfAlleles, bool phased,
                                     tag<false> const&)
        {
            // do nothing
        }

        template <typename Setter>
        void call_set_min_max_ploidy(Setter& setter, uint32_t min_ploidy,
                                     uint32_t max_ploidy,
                                     uint32_t numberOfAlleles, bool phased)
        {
            call_set_min_max_ploidy(setter, min_ploidy, max_ploidy,
                                    numberOfAlleles, phased,
                                    tag<has_set_min_max_ploidy<Setter>::Yes>());
        }

        template <typename Setter>
        struct has_finalise
        {
            template <typename U, void (U::*)()>
            struct SFINAE
            {
            };
            template <typename U>
            static uint8_t Test(SFINAE<U, &U::finalise>*);
            template <typename U>
            static uint32_t Test(...);
            static const bool Yes = sizeof(Test<Setter>(0)) == sizeof(uint8_t);
            static const bool No = sizeof(Test<Setter>(0)) == sizeof(uint16_t);
        };

        template <typename Setter>
        void call_finalise(Setter& setter, tag<true> const&)
        {
            setter.finalise();
        }

        template <typename Setter>
        void call_finalise(Setter& setter, tag<false> const&)
        {
            // do nothing
        }

        template <typename Setter>
        void call_finalise(Setter& setter)
        {
            call_finalise(setter, tag<has_finalise<Setter>::Yes>());
        }
    }

    namespace impl
    {
        struct ProbabilityDataWriterBase
        {
            virtual ~ProbabilityDataWriterBase() {};

            virtual void initialise(uint32_t nSamples, uint16_t nAlleles,
                                    byte_t* buffer, byte_t* const end) = 0;
            virtual bool set_sample(std::size_t i) = 0;
            virtual void set_number_of_entries(uint32_t ploidy,
                                               uint32_t number_of_entries,
                                               OrderType const order_type,
                                               ValueType const value_type) = 0;
            virtual void set_value(uint32_t entry_i,
                                   genfile::MissingValue const value) = 0;
            virtual void set_value(uint32_t entry_i, double const value) = 0;
            virtual void finalise() = 0;

            virtual std::pair<byte_t const*, byte_t const*> repr() const = 0;
        };
    }

    namespace v11
    {
        namespace impl
        {
            template <typename FloatType>
            FloatType convert_from_integer_representation(uint16_t number,
                                                          FloatType factor)
            {
                FloatType result = number;
                result /= factor;
                return result;
            }

            template <typename FloatType>
            uint16_t convert_to_integer_representation(FloatType number,
                                                       FloatType factor)
            {
                number *= factor;
                number = std::min(std::max(number, 0.0), 65535.0);
                return static_cast<uint16_t>(std::floor(number + 0.5));
            }

            double get_probability_conversion_factor(uint32_t flags);
        }

        struct ProbabilityDataWriter
            : public genfile::bgen::impl::ProbabilityDataWriterBase
        {
            enum Missing
            {
                eNotSet = 0,
                eMissing = 1,
                eNotMissing = 2
            };
            enum State
            {
                eUninitialised = 0,
                eInitialised = 1,
                eSampleSet = 2,
                eNumberOfEntriesSet = 3,
                eValueSet = 4,
                eBaked = 5,
                eFinalised = 6
            };

            ~ProbabilityDataWriter() {}

            ProbabilityDataWriter()
                : m_state(eUninitialised), m_sample_i(0), m_missing(eNotSet)
            {
            }

            void initialise(uint32_t nSamples, uint16_t nAlleles,
                            byte_t* buffer, byte_t* const end)
            {
                assert(nAlleles == 2);
                m_p = m_buffer = buffer;
                m_end = end;
                m_number_of_samples = nSamples;
                m_state = eInitialised;
            }

            bool set_sample(std::size_t i)
            {
                assert(m_state == eInitialised || m_state == eBaked
                       || m_state == eSampleSet);
                assert((m_sample_i == 0 && i == 0) || (i == m_sample_i + 1));
                if (m_state == eSampleSet)
                {
                    // last sample had no data, write zeroes
                    m_values[0] = m_values[1] = m_values[2] = 0.0;
                    bake(&m_values[0]);
                }
                m_sample_i = i;
                m_state = eSampleSet;
                return true;
            }

            void set_number_of_entries(uint32_t ploidy,
                                       uint32_t number_of_entries,
                                       OrderType const order_type,
                                       ValueType const value_type)
            {
                assert(m_state == eSampleSet);
                if (ploidy != uint32_t(2)) { throw BGenError(); }
                assert(number_of_entries == uint32_t(3));
                assert(order_type == ePerUnorderedGenotype);
                m_entry_i = 0;
                m_missing = eNotSet;
                m_state = eNumberOfEntriesSet;
            }

            void set_value(uint32_t entry_i, genfile::MissingValue const value)
            {
                assert(m_state == eNumberOfEntriesSet || m_state == eValueSet);
                assert(m_entry_i < 3);
                assert(m_entry_i == 0 || m_missing == eMissing);
                m_values[m_entry_i++] = 0.0;
                m_missing = eMissing;
                if (m_entry_i == 3)
                {
                    bake(&m_values[0]);
                    m_state = eBaked;
                }
                else
                {
                    m_state = eValueSet;
                }
            }

            void set_value(uint32_t entry_i, double const value)
            {
                assert(m_state == eNumberOfEntriesSet || m_state == eValueSet);
                assert(m_missing == eNotSet || m_missing == eNotMissing);
                assert(m_entry_i < 3);
                m_values[m_entry_i++] = value;
                if (value != 0.0) { m_missing = eNotMissing; }
                if (m_entry_i == 3)
                {
                    bake(&m_values[0]);
                    m_state = eBaked;
                }
                else
                {
                    m_state = eValueSet;
                }
            }

            void finalise()
            {
                assert((m_number_of_samples == 0 && m_state == eInitialised)
                       || m_state == eBaked);
                m_state = eFinalised;
            }

            std::pair<byte_t const*, byte_t const*> repr() const
            {
                return std::make_pair(m_buffer, m_p);
            }

        private:
            byte_t* m_buffer;
            byte_t* m_p;
            byte_t* m_end;
            State m_state;
            uint32_t m_number_of_samples;
            std::size_t m_sample_i;
            std::size_t m_entry_i;
            Missing m_missing;
            double m_values[3];

        private:
            void bake(double* values)
            {
                assert(m_p + 6 <= m_end);
                m_p = write_little_endian_integer(
                    m_p, m_end,
                    impl::convert_to_integer_representation(*values++,
                                                            32768.0));
                m_p = write_little_endian_integer(
                    m_p, m_end,
                    impl::convert_to_integer_representation(*values++,
                                                            32768.0));
                m_p = write_little_endian_integer(
                    m_p, m_end,
                    impl::convert_to_integer_representation(*values++,
                                                            32768.0));
            }
        };

        template <typename Setter>
        void parse_probability_data(byte_t const* buffer,
                                    byte_t const* const end,
                                    Context const& context, Setter& setter)
        {
            if (end != buffer + 6 * context.number_of_samples)
            { throw BGenError(); }
            setter.initialise(context.number_of_samples, 2);
            uint32_t const ploidy = 2;
            call_set_min_max_ploidy(setter, 2ul, 2ul, 2ul, false);
            double const probability_conversion_factor =
                impl::get_probability_conversion_factor(context.flags);
            for (uint32_t i = 0; i < context.number_of_samples; ++i)
            {
                if (setter.set_sample(i))
                {
                    setter.set_number_of_entries(
                        ploidy, 3, ePerUnorderedGenotype, eProbability);
                    assert(end >= buffer + 6);
                    for (std::size_t g = 0; g < 3; ++g)
                    {
                        uint16_t prob;
                        buffer = read_little_endian_integer(buffer, end, &prob);
                        setter.set_value(
                            g, impl::convert_from_integer_representation(
                                   prob, probability_conversion_factor));
                    }
                    setter.sample_completed();
                }
                else
                {
                    buffer += 6;
                }
            }
            call_finalise(setter);
        }
    }

    namespace v12
    {
        namespace impl
        {
            struct BitParser
            {
                BitParser(byte_t const* buffer, byte_t const* const end,
                          int const bits)
                    : m_buffer(buffer)
                    , m_end(end)
                    , m_bits(bits)
                    , m_bitMask(((uint64_t)(0xFFFFFFFFFFFFFFFF) >> (64 - bits)))
                    , m_denominator(m_bitMask)
                    , m_shift(0)
                {
                    assert(bits > 0 && bits <= 32);
#if BGEN_BIG_ENDIAN
                    read_little_endian_integer(buffer, end, &m_data);
#endif
                }

                // check we can consume n more values
                bool check(std::size_t n) const
                {
                    // We need enough bytes to deal with the current shift value
                    // plus enough to deal with the requested bits.
                    std::size_t const bitsNeeded = n * m_bits + m_shift;
                    std::size_t const bytesNeededFromBuffer =
                        (bitsNeeded + 7) / 8;
                    return (m_buffer + bytesNeededFromBuffer) <= m_end;
                }

                // consume and return next value
                double next()
                {
#if BGEN_LITTLE_ENDIAN
                    // Machine endianness matches stored endianness so no
                    // reordering of bytes is needed.
                    // We travel through the data 32 bits at a time
                    // Each time we use a shift to get the appropriate bits from
                    // the current 32-bit word.
                    double value =
                        ((*reinterpret_cast<uint64_t const*>(m_buffer)
                          >> m_shift)
                         & m_bitMask)
                        / m_denominator;
#else // BGEN_BIG_ENDIAN
                    double value =
                        ((m_data >> m_shift) & m_bitMask) / m_denominator;
#endif
                    m_shift += m_bits;
                    if (m_shift > 31)
                    {
                        // m_bits is at most 32
                        // so m_shift can now be a maximum of 63.
                        // std::cerr << "m_shift = " << m_shift << ", moving
                        // buffer.\n" ;
                        m_buffer += 4;
                        m_shift -= 32;
#if BGEN_BIG_ENDIAN
                        // marshal data through a uint64.
                        read_little_endian_integer(buffer, end, &m_data);
#endif
                    }
                    return value;
                }

            private:
                byte_t const* m_buffer;
                byte_t const* const m_end;
                int const m_bits;
                uint64_t const m_bitMask;
                double const m_denominator;
                int m_shift;
#if BGEN_BIG_ENDIAN
                uint64_t m_data;
#endif
            };

            // Optimised bit parser implementations
            template <int bits>
            struct SpecialisedBitParser;

            // Specialisation for 8 bit data
            template <>
            struct SpecialisedBitParser<8>
            {
                SpecialisedBitParser(byte_t const* buffer,
                                     byte_t const* const end)
                    : m_buffer(buffer), m_end(end)
                {
                }

                // check we can consume n more values
                bool check(std::size_t n) const
                {
                    return (m_buffer + n) <= m_end;
                }

                double next()
                {
                    return double(
                               (*reinterpret_cast<uint8_t const*>(m_buffer++)))
                           / 255.0;
                }

            private:
                byte_t const* m_buffer;
                byte_t const* const m_end;
            };

            // Specialisation for 16 bit data
            template <>
            struct SpecialisedBitParser<16>
            {
                SpecialisedBitParser(byte_t const* buffer,
                                     byte_t const* const end)
                    : m_buffer(buffer), m_end(end)
                {
                }

                // check we can consume n more values
                bool check(std::size_t n) const
                {
                    return (m_buffer + 2 * n) <= m_end;
                }

                double next()
                {
#if BGEN_LITTLE_ENDIAN
                    double const value =
                        double(*reinterpret_cast<uint16_t const*>(m_buffer))
                        / 65535.0;
#else // BGEN_BIG_ENDIAN
      // machine is big-endian, get bytes in right order.
                    uint16_t data =
                        uint16_t(*m_buffer) | uint16_t(*(m_buffer + 1)) << 8;
                    double const value = double(data) / 65535.0;
#endif
                    m_buffer += 2;
                    return value;
                }

            private:
                byte_t const* m_buffer;
                byte_t const* const m_end;
            };

            // Round a point on the unit simplex (expressed as n floating-point
            // probabilities) to a point representable with the given number of
            // bits. precondition: p points to n doubles between 0 and 1 that
            // sum to 1 (up to floating-point precision). postcondition: the
            // values p points to are n integer values that sum to
            // 2^(number_of_bits)-1.
            void compute_approximate_probabilities(double* p,
                                                   std::size_t* index,
                                                   std::size_t const n,
                                                   int const number_of_bits);

            // Write data encoding n probabilities, given in probs, that sum to
            // 1, starting at the given offset in data, to the given buffer.
            byte_t* write_scaled_probs(uint64_t* data, std::size_t* offset,
                                       double const* probs, std::size_t const n,
                                       int const number_of_bits,
                                       byte_t* destination, byte_t* const end);
        }

        struct GenotypeDataBlock
        {
        public:
            GenotypeDataBlock();

            GenotypeDataBlock(Context const& context_, byte_t const* buffer,
                              byte_t const* const end);

            void initialise(Context const& context, byte_t const* buffer,
                            byte_t const* const end);

        public:
            Context const* context;
            uint32_t numberOfSamples;
            uint16_t numberOfAlleles;
            byte_t ploidyExtent[2];
            byte_t const* ploidy; // Must contain at least N bytes.
            bool phased;
            byte_t bits;
            byte_t const* buffer;
            byte_t const* end;

        private:
            // forbid copying.
            GenotypeDataBlock(GenotypeDataBlock const& other);
            GenotypeDataBlock& operator=(GenotypeDataBlock const& other);
        };

        inline GenotypeDataBlock::GenotypeDataBlock()
            : context(0)
            , numberOfSamples(0)
            , numberOfAlleles(0)
            , ploidy(0)
            , phased(false)
            , bits(0)
            , buffer(0)
            , end(0)
        {
        }

        inline GenotypeDataBlock::GenotypeDataBlock(Context const& context,
                                                    byte_t const* buffer,
                                                    byte_t const* const end)
        {
            initialise(context, buffer, end);
        }

        inline void GenotypeDataBlock::initialise(Context const& context_,
                                                  byte_t const* buffer,
                                                  byte_t const* const end)
        {
            if (end < buffer + 8) { throw BGenError(); }

            uint32_t N = 0;
            buffer = read_little_endian_integer(buffer, end, &N);
            if (N != context_.number_of_samples) { throw BGenError(); }
            if (end < buffer + N + 2) { throw BGenError(); }

            buffer = read_little_endian_integer(buffer, end, &numberOfAlleles);
            buffer = read_little_endian_integer(buffer, end, &ploidyExtent[0]);
            buffer = read_little_endian_integer(buffer, end, &ploidyExtent[1]);

            // Keep a pointer to the ploidy and move buffer past the ploidy
            // information
            this->context = &context_;
            this->numberOfSamples = N;
            this->ploidy = buffer;
            buffer += N;
            // Get the phased flag and number of bits
            this->phased = ((*buffer++) & 0x1);
            this->bits = *reinterpret_cast<byte_t const*>(buffer++);
            this->buffer = buffer;
            this->end = end;
        }


        template <typename Setter>
        void parse_probability_data(byte_t const* buffer,
                                    byte_t const* const end,
                                    Context const& context, Setter& setter)
        {
            parse_probability_data(GenotypeDataBlock(context, buffer, end),
                                   setter);
        }

        template <typename Setter>
        void parse_probability_data(GenotypeDataBlock const& pack,
                                    Setter& setter)
        {
            Context const& context = *(pack.context);
            // We optimise the most common and simplest-to- parse cases.
            // These are the case where all samples are diploid, and/or where
            // the number of bits is a multiple of 8.
            // This if statement chooses an appropriate implementation.
            if (pack.ploidyExtent[0] == 2 && pack.ploidyExtent[1] == 2
                && pack.numberOfAlleles == 2)
            {
                switch (pack.bits)
                {
                case 8:
                    parse_probability_data_diploid_biallelic(
                        pack,
                        impl::SpecialisedBitParser<8>(pack.buffer, pack.end),
                        context, setter);
                    break;
                case 16:
                    parse_probability_data_diploid_biallelic(
                        pack,
                        impl::SpecialisedBitParser<16>(pack.buffer, pack.end),
                        context, setter);
                    break;
                default:
                    parse_probability_data_diploid_biallelic(
                        pack, impl::BitParser(pack.buffer, pack.end, pack.bits),
                        context, setter);
                    break;
                }
            }
            else
            {
                switch (pack.bits)
                {
                case 8:
                    parse_probability_data_general(
                        pack,
                        impl::SpecialisedBitParser<8>(pack.buffer, pack.end),
                        context, setter);
                    break;
                case 16:
                    parse_probability_data_general(
                        pack,
                        impl::SpecialisedBitParser<16>(pack.buffer, pack.end),
                        context, setter);
                    break;
                default:
                    parse_probability_data_general(
                        pack, impl::BitParser(pack.buffer, pack.end, pack.bits),
                        context, setter);
                    break;
                }
            }
        }

        template <typename Setter, typename BitParser>
        void parse_probability_data_diploid_biallelic(
            GenotypeDataBlock const& pack, BitParser valueConsumer,
            Context const& context, Setter& setter)
        {
            assert(pack.numberOfAlleles == 2);
            assert(pack.ploidyExtent[0] == 2);
            assert(pack.ploidyExtent[1] == 2);

            // These values are specific bit combinations and should not be
            // changed.
            enum SampleStatus
            {
                eIgnore = 0,
                eSetThisSample = 1,
                eSetAsMissing = 3
            };
            byte_t const* ploidy_p = pack.ploidy;
#if DEBUG_BGEN_FORMAT
            std::cerr << "parse_probability_data_v12(): numberOfSamples = "
                      << numberOfSamples << ", phased = " << phased << ".\n";
#endif

            setter.initialise(pack.numberOfSamples, uint32_t(2));
            call_set_min_max_ploidy(setter, uint32_t(2), uint32_t(2), 2,
                                    pack.phased);

            {
                if (pack.phased)
                {
                    for (uint32_t i = 0; i < pack.numberOfSamples;
                         ++i, ++ploidy_p)
                    {
                        bool const missing = (*ploidy_p & 0x80);
                        int const sample_status =
                            (setter.set_sample(i) * 0x1) + (missing * 0x2);

                        if (sample_status & 0x1)
                        {
                            setter.set_number_of_entries(
                                2, 4, ePerPhasedHaplotypePerAllele,
                                eProbability);
                        }


                        if (!valueConsumer.check(2)) { throw BGenError(); }
                        // Consume values and interpret them.
                        for (uint32_t hap = 0; hap < 2; ++hap)
                        {
                            double const value = valueConsumer.next();
                            switch (sample_status)
                            {
                            case eIgnore: break;
                            case eSetAsMissing:
                                setter.set_value(2 * hap + 0,
                                                 genfile::MissingValue());
                                setter.set_value(2 * hap + 1,
                                                 genfile::MissingValue());
                                break;
                            case eSetThisSample:
                                setter.set_value(2 * hap + 0, value);
                                setter.set_value(2 * hap + 1, 1.0 - value);
                                break;
                            }
                        }
                        setter.sample_completed();
                    }
                }
                else
                {
                    for (uint32_t i = 0; i < pack.numberOfSamples;
                         ++i, ++ploidy_p)
                    {
                        bool const missing = (*ploidy_p & 0x80);
                        int const sample_status =
                            (setter.set_sample(i) * 0x1) + (missing * 0x2);

                        if (sample_status & 0x1)
                        {
                            setter.set_number_of_entries(
                                2, 3, ePerUnorderedGenotype, eProbability);
                        }
                        if (!valueConsumer.check(2)) { throw BGenError(); }
                        double const value1 = valueConsumer.next();
                        double const value2 = valueConsumer.next();

                        switch (sample_status)
                        {
                        case eIgnore: break;
                        case eSetAsMissing:
                            setter.set_value(0, genfile::MissingValue());
                            setter.set_value(1, genfile::MissingValue());
                            setter.set_value(2, genfile::MissingValue());
                            break;
                        case eSetThisSample:
                            setter.set_value(0, value1);
                            setter.set_value(1, value2);
                            // Clamp the value to 0 to avoid small -ve values
                            setter.set_value(
                                2, std::max(1.0 - value1 - value2, 0.0));
                            break;
                        }
                        setter.sample_completed();
                    }
                }
            }
            call_finalise(setter);
        }

        template <typename Setter, typename BitParser>
        void parse_probability_data_general(GenotypeDataBlock const& pack,
                                            BitParser valueConsumer,
                                            Context const& context,
                                            Setter& setter)
        {
            // These values are specific bit combinations and should not be
            // changed.
            enum SampleStatus
            {
                eIgnore = 0,
                eSetThisSample = 1,
                eSetAsMissing = 3
            };

            byte_t const* ploidy_p = pack.ploidy;

#if DEBUG_BGEN_FORMAT
            std::cerr << "parse_probability_data_v12(): numberOfSamples = "
                      << numberOfSamples << ", phased = " << phased << ".\n";
            std::cerr << "parse_probability_data_v12(): *buffer: "
                      << bgen::impl::to_hex(buffer, end) << ".\n";
#endif

            setter.initialise(pack.numberOfSamples,
                              uint32_t(pack.numberOfAlleles));
            call_set_min_max_ploidy(setter, uint32_t(pack.ploidyExtent[0]),
                                    uint32_t(pack.ploidyExtent[1]),
                                    pack.numberOfAlleles, pack.phased);

            {
                if (pack.phased)
                {
                    for (uint32_t i = 0; i < pack.numberOfSamples;
                         ++i, ++ploidy_p)
                    {
                        uint32_t const ploidy = uint32_t(*ploidy_p & 0x3F);
                        bool const missing = (*ploidy_p & 0x80);
                        int const sample_status =
                            (setter.set_sample(i) * 0x1) + (missing * 0x2);

                        uint32_t const valueCount =
                            (ploidy * pack.numberOfAlleles);

                        if (sample_status & 0x1)
                        {
                            setter.set_number_of_entries(
                                ploidy, valueCount,
                                ePerPhasedHaplotypePerAllele, eProbability);
                        }

                        // Consume values and interpret them.
                        double sum = 0.0;
                        uint32_t reportedValueCount = 0;
                        if (!valueConsumer.check(ploidy
                                                 * (pack.numberOfAlleles - 1)))
                        { throw BGenError(); }
                        for (uint32_t hap = 0; hap < ploidy; ++hap)
                        {
                            for (uint32_t allele = 0;
                                 allele < (pack.numberOfAlleles - 1); ++allele)
                            {
                                double const value = valueConsumer.next();
                                switch (sample_status)
                                {
                                case eIgnore: break;
                                case eSetAsMissing:
                                    setter.set_value(reportedValueCount++,
                                                     genfile::MissingValue());
                                    break;
                                case eSetThisSample:
                                    setter.set_value(reportedValueCount++,
                                                     value);
                                    sum += value;
                                    break;
                                }
                            }

                            // set value for kth allele
                            switch (sample_status)
                            {
                            case eIgnore: break;
                            case eSetAsMissing:
                                setter.set_value(reportedValueCount++,
                                                 genfile::MissingValue());
                                break;
                            case eSetThisSample:
                                setter.set_value(reportedValueCount++, 1 - sum);
                                sum = 0.0;
                                break;
                            }
                        }
                        setter.sample_completed();
                    }
                }
                else
                {
                    for (uint32_t i = 0; i < pack.numberOfSamples;
                         ++i, ++ploidy_p)
                    {
                        uint32_t const ploidy = uint32_t(*ploidy_p & 0x3F);
                        bool const missing = (*ploidy_p & 0x80);
                        int const sample_status =
                            (setter.set_sample(i) * 0x1) + (missing * 0x2);

                        uint32_t const valueCount =
                            genfile::bgen::impl::n_choose_k(
                                uint32_t(ploidy + pack.numberOfAlleles - 1),
                                uint32_t(pack.numberOfAlleles - 1));
                        uint32_t const storedValueCount = valueCount - 1;

                        if (sample_status & 0x1)
                        {
                            setter.set_number_of_entries(ploidy, valueCount,
                                                         ePerUnorderedGenotype,
                                                         eProbability);
                        }

                        if (!valueConsumer.check(storedValueCount))
                        { throw BGenError(); }

                        double sum = 0.0;
                        uint32_t reportedValueCount = 0;
                        for (uint32_t h = 0; h < storedValueCount; ++h)
                        {
                            double const value = valueConsumer.next();
                            switch (sample_status)
                            {
                            case eIgnore: break;
                            case eSetAsMissing:
                                setter.set_value(reportedValueCount++,
                                                 genfile::MissingValue());
                                break;
                            case eSetThisSample:
                                setter.set_value(reportedValueCount++, value);
                                sum += value;
                                break;
                            }
                        }

                        // set final value
                        switch (sample_status)
                        {
                        case eIgnore: break;
                        case eSetAsMissing:
                            setter.set_value(reportedValueCount++,
                                             genfile::MissingValue());
                            break;
                        case eSetThisSample:
                            setter.set_value(reportedValueCount++, 1.0 - sum);
                            break;
                        }
                        setter.sample_completed();
                    }
                }
            }
            call_finalise(setter);
        }

        struct ProbabilityDataWriter
            : public genfile::bgen::impl::ProbabilityDataWriterBase
        {
            enum
            {
                eMinPloidyByte = 6,
                eMaxPloidyByte = 7,
                ePloidyBytes = 8
            };
            enum Missing
            {
                eNotSet = 0,
                eMissing = 1,
                eNotMissing = 2
            };
            enum State
            {
                eUninitialised = 0,
                eInitialised = 1,
                eSampleSet = 2,
                eNumberOfEntriesSet = 3,
                eValueSet = 4,
                eBaked = 5,
                eFinalised = 6
            };

            ~ProbabilityDataWriter() {}

            ProbabilityDataWriter(
                uint8_t const number_of_bits,
                double const max_rounding_error_per_prob = 0.0005)
                : m_number_of_bits(number_of_bits)
                ,
                // Actually likely rounding error is rounding error from limit
                // precision number, plus error in floating point representation
                // (which is at most epsilon).
                m_max_error_per_prob(max_rounding_error_per_prob
                                     + std::numeric_limits<double>::epsilon())
                , m_state(eUninitialised)
                , m_order_type(eUnknownOrderType)
                , m_number_of_samples(0)
                , m_number_of_alleles(0)
                , m_ploidy(0)
                , m_sample_i(0)
                , m_number_of_entries(0)
                , m_entries_per_bake(0)
                , m_entry_i(0)
                , m_missing(eNotSet)
                , m_sum(0.0)
                , m_data(0)
            {
                m_ploidyExtent[0] = 63;
                m_ploidyExtent[1] = 0;
            }

            void initialise(uint32_t nSamples, uint16_t nAlleles,
                            byte_t* buffer, byte_t* const end)
            {
                assert(m_state == eUninitialised);
                m_p = m_buffer = buffer;
                m_end = end;

                // Write whatever fieds we can now.
                m_p = write_little_endian_integer(m_p, m_end, nSamples);
                m_p = write_little_endian_integer(m_p, m_end, nAlleles);
                // skip the ploidy and other non-data bytes, which we'll write
                // later
                m_number_of_samples = nSamples;
                m_number_of_alleles = nAlleles;
                m_p += nSamples + 4;
                m_data = 0;
                m_offset = 0;
                if (m_number_of_samples == 0)
                {
                    m_ploidyExtent[0] = 0;
                    m_order_type = ePerPhasedHaplotypePerAllele;
                }
                // We set data as phased until we learn otherwise
                // In the case of samples with ploidy < 2, this means we
                // nominally interpret data the data as phased until we see an
                // informative sample.
                m_buffer[8 + m_number_of_samples] = 1;
                m_buffer[9 + m_number_of_samples] = m_number_of_bits;
                m_state = eInitialised;
            }

            bool set_sample(std::size_t i)
            {
                assert(m_state == eInitialised || m_state == eBaked
                       || m_state == eSampleSet);
                // ensure samples are visited in order.
                assert((m_sample_i == 0 && i == 0) || (i == m_sample_i + 1));
                if (m_state == eSampleSet)
                {
                    // Last sample was completely missing, mark as 0 ploid &
                    // missing
                    m_buffer[ePloidyBytes + m_sample_i] = 0x80;
                }
                m_sample_i = i;
                m_state = eSampleSet;
                return true;
            }

            void set_number_of_entries(uint32_t ploidy,
                                       uint32_t number_of_entries,
                                       OrderType const order_type,
                                       ValueType const value_type)
            {
                assert(m_state == eSampleSet);

                assert(ploidy < 64);
                m_ploidy = ploidy;
                uint8_t ploidyByte(ploidy & 0xFF);
                m_buffer[ePloidyBytes + m_sample_i] = ploidyByte;
                m_ploidyExtent[0] = std::min(m_ploidyExtent[0], ploidyByte);
                m_ploidyExtent[1] = std::max(m_ploidyExtent[1], ploidyByte);

                assert(order_type == ePerUnorderedGenotype
                       || order_type == ePerPhasedHaplotypePerAllele);
                if (ploidy > 1)
                {
                    if (m_order_type == eUnknownOrderType)
                    { m_order_type = order_type; }
                    else
                    {
                        if (m_order_type != order_type) { throw BGenError(); }
                    }
                }
                if (value_type != eProbability) { throw BGenError(); }
                m_number_of_entries = number_of_entries;
                m_entries_per_bake =
                    (m_order_type == ePerUnorderedGenotype || m_ploidy == 0)
                        ? m_number_of_entries
                        : m_number_of_alleles;
                m_entry_i = 0;
                m_missing = eNotSet;
                m_state = eNumberOfEntriesSet;
                m_sum = 0.0;
            }

            void set_value(uint32_t entry_i, genfile::MissingValue const value)
            {
                assert(m_state == eNumberOfEntriesSet || m_state == eValueSet
                       || m_state == eBaked);
                assert(m_entry_i < m_number_of_entries);
                assert(m_entry_i == 0 || m_missing == eMissing);
                assert(m_entry_i == (entry_i % m_entries_per_bake));
                m_values[m_entry_i++] = 0.0;
                m_missing = eMissing;
                if (m_entry_i == m_entries_per_bake)
                {
                    bake(&m_values[0], m_entry_i, m_sum);
                    m_entry_i = 0;
                    m_sum = 0.0;
                    m_state = eBaked;
                }
                else
                {
                    m_state = eValueSet;
                }
            }

            void set_value(uint32_t entry_i, double const value)
            {
                assert(m_state == eNumberOfEntriesSet || m_state == eValueSet
                       || m_state == eBaked);
                assert(m_missing == eNotSet || m_missing == eNotMissing);
                assert(m_entry_i < m_number_of_entries);
                assert(m_entry_i == (entry_i % m_entries_per_bake));
                m_values[m_entry_i++] = value;
                m_sum += value;

#if DEBUG_BGEN_FORMAT
                std::cerr << "set_value( " << entry_i << ", " << value
                          << "); m_entry_i = " << m_entry_i << "\n";
#endif
                if (value != value || value < 0.0
                    || value > (1.0 + m_max_error_per_prob))
                {
                    std::cerr << "Sample " << m_sample_i << ", value "
                              << entry_i << " is " << std::setprecision(17)
                              << value << ", expected within bounds 0 - "
                              << (1.0 + m_max_error_per_prob) << ".\n";
                    throw BGenError();
                }
                if (value != 0.0) { m_missing = eNotMissing; }
                if (m_entry_i == m_entries_per_bake)
                {
                    bake(&m_values[0], m_entry_i, m_sum);
                    m_entry_i = 0;
                    m_sum = 0.0;
                    m_state = eBaked;
                }
                else
                {
                    m_state = eValueSet;
                }
            }

            void finalise()
            {
                if (m_order_type == eUnknownOrderType)
                { m_order_type = ePerPhasedHaplotypePerAllele; }
                m_buffer[8 + m_number_of_samples] =
                    (m_order_type == ePerPhasedHaplotypePerAllele) ? 1 : 0;
                m_buffer[9 + m_number_of_samples] = m_number_of_bits;

                if (!((m_number_of_samples == 0 && m_state == eInitialised)
                      || (m_number_of_entries == 0
                          && m_state == eNumberOfEntriesSet)
                      || m_state == eBaked))
                {
                    std::cerr
                        << "genfile::bgen::v12::ProbabilityDataWriter::"
                           "finalise(): m_number_of_samples = "
                        << m_number_of_samples << " m_state = " << m_state
                        << " m_number_of_entries = " << m_number_of_entries
                        << " m_entries_per_bake = " << m_entries_per_bake
                        << " m_entry_i = " << m_entry_i
                        << " m_order_type = " << m_order_type
                        << " m_ploidy = " << m_ploidy << ".\n";
                    throw BGenError();
                }
                // Write any remaining data
                if (m_offset > 0)
                {
                    int const nBytes = (m_offset + 7) / 8;
#if DEBUG_BGEN_FORMAT
                    std::cerr
                        << "genfile::bgen::impl::v12::ProbabilityDataWriter::"
                           "finalise()(): final offset = "
                        << m_offset
                        << ", number of bits = " << int(m_number_of_bits)
                        << ", writing " << nBytes
                        << " final bytes (space = " << (m_end - m_p)
                        << ", ploidy extent = " << uint32_t(m_ploidyExtent[0])
                        << " " << uint32_t(m_ploidyExtent[0]) << "\n";
#endif
                    assert((m_p + nBytes) <= m_end);
                    m_p = std::copy(
                        reinterpret_cast<byte_t const*>(&m_data),
                        reinterpret_cast<byte_t const*>(&m_data) + nBytes, m_p);
                }
                m_buffer[eMinPloidyByte] = m_ploidyExtent[0];
                m_buffer[eMaxPloidyByte] = m_ploidyExtent[1];
                m_state = eFinalised;
            }

            std::pair<byte_t const*, byte_t const*> repr() const
            {
                return std::make_pair(m_buffer, m_p);
            }

        private:
            byte_t* m_buffer;
            byte_t* m_p;
            byte_t* m_end;
            uint8_t const m_number_of_bits;
            double const m_max_error_per_prob;
            State m_state;
            uint8_t m_ploidyExtent[2];
            OrderType m_order_type;
            uint32_t m_number_of_samples;
            uint32_t m_number_of_alleles;
            uint32_t m_ploidy;
            std::size_t m_sample_i;
            std::size_t m_number_of_entries;
            std::size_t m_entries_per_bake;
            std::size_t m_entry_i;
            Missing m_missing;
            double m_values[100];
            double m_sum;
            std::size_t index[100];
            uint64_t m_data;
            std::size_t m_offset;

        private:
            void bake(double* values, std::size_t count, double const sum)
            {
                if (m_missing == eMissing || m_missing == eNotSet)
                {
                    // Have never seen a non-missing value.
                    m_p = impl::write_scaled_probs(&m_data, &m_offset, values,
                                                   count, m_number_of_bits, m_p,
                                                   m_end);
#if DEBUG_BGEN_FORMAT
                    std::cerr
                        << "genfile::bgen::impl::v12::ProbabilityDataWriter::"
                           "bake()(): setting sample "
                        << m_sample_i << " to missing (byte "
                        << (ePloidyBytes + m_sample_i) << ").\n";
#endif
                    // flag this sample as missing.
                    m_buffer[ePloidyBytes + m_sample_i] |= 0x80;
                }
                else
                {
                    double const max_error_in_sum =
                        (count * m_max_error_per_prob);
                    if ((sum != sum) || (sum > (1.0 + max_error_in_sum))
                        || (sum < (1.0 - max_error_in_sum)))
                    {
                        std::cerr << "These " << count << " values sum to "
                                  << std::fixed << std::setprecision(17) << sum
                                  << ", "
                                  << "I expected the sum to be in the range "
                                  << (1.0 - max_error_in_sum) << " - "
                                  << (1.0 + max_error_in_sum) << ".\n";
                        std::cerr << "Values are:\n";
                        for (std::size_t i = 0; i < count; ++i)
                        { std::cerr << values[i] << "\n"; }
                        throw BGenError();
                    }
                    // We project onto the unit simplex before computing the
                    // approximation.
                    for (std::size_t i = 0; i < count; ++i)
                    { values[i] /= sum; }
                    impl::compute_approximate_probabilities(
                        values, &index[0], count, m_number_of_bits);
                    m_p = impl::write_scaled_probs(&m_data, &m_offset, values,
                                                   count, m_number_of_bits, m_p,
                                                   m_end);
                }
            }
        };
    }

    template <typename Setter>
    void parse_probability_data(byte_t const* buffer, byte_t const* const end,
                                Context const& context, Setter& setter)
    {
        if ((context.flags & e_Layout) == e_Layout0
            || (context.flags & e_Layout) == e_Layout1)
        { v11::parse_probability_data(buffer, end, context, setter); }
        else
        {
            v12::parse_probability_data(buffer, end, context, setter);
        }
    }

    template <typename Setter>
    void read_and_parse_genotype_data_block(std::istream& aStream,
                                            Context const& context,
                                            Setter& setter,
                                            std::vector<byte_t>* buffer1,
                                            std::vector<byte_t>* buffer2)
    {
        read_genotype_data_block(aStream, context, buffer1);
        uncompress_probability_data(context, *buffer1, buffer2);
        parse_probability_data(&(*buffer2)[0], &(*buffer2)[0] + buffer2->size(),
                               context, setter);
    }
    template <typename Setter>
    void read_and_parse_genotype_data_block(
        FileRead& aStream, const std::string& file_name, Context const& context,
        Setter& setter, std::vector<byte_t>* buffer1,
        std::vector<byte_t>* buffer2, const std::streampos idx)
    {
        read_genotype_data_block(aStream, file_name, context, buffer1, idx);
        uncompress_probability_data(context, *buffer1, buffer2);
        parse_probability_data(&(*buffer2)[0], &(*buffer2)[0] + buffer2->size(),
                               context, setter);
    }
    // Write identifying data fields for the given variant.
    template <typename AlleleGetter>
    byte_t* write_snp_identifying_data(
        std::vector<byte_t>* buffer, Context const& context, std::string SNPID,
        std::string RSID, std::string chromosome, uint32_t position,
        uint16_t const number_of_alleles, AlleleGetter get_allele)
    {
        uint32_t const layout = context.flags & e_Layout;
        assert(layout == e_Layout1 || layout == e_Layout2);

        // Make sure we have space
        std::size_t size = 10 + RSID.size() + SNPID.size() + chromosome.size()
                           + ((layout == e_Layout1) ? 4 : 2);
        for (uint16_t allele_i = 0; allele_i < number_of_alleles; ++allele_i)
        { size += 4 + get_allele(allele_i).size(); }
        buffer->resize(size);

        // Write the data to the buffer
        byte_t* p = &(buffer->operator[](0));
        byte_t* const end = p + size;

        if (layout == e_Layout1)
        {
            p = write_little_endian_integer(p, end, context.number_of_samples);
            // otherwise number of samples appears in the probability data block
            // below.
        }

        std::size_t const max_id_length = std::numeric_limits<uint16_t>::max();
        assert(SNPID.size() <= static_cast<std::size_t>(max_id_length));
        assert(RSID.size() <= static_cast<std::size_t>(max_id_length));
        p = write_length_followed_by_data(p, end, uint16_t(SNPID.size()),
                                          SNPID.data());
        p = write_length_followed_by_data(p, end, uint16_t(RSID.size()),
                                          RSID.data());
        p = write_length_followed_by_data(p, end, uint16_t(chromosome.size()),
                                          chromosome);
        p = write_little_endian_integer(p, end, position);

        if (layout == e_Layout2)
        {
            // v12 has an explicit allele count
            p = write_little_endian_integer(p, end, number_of_alleles);
        }
        else if (number_of_alleles != 2u)
        {
            throw BGenError();
        }

        std::size_t const max_allele_length =
            std::numeric_limits<uint32_t>::max();
        for (uint16_t allele_i = 0; allele_i < number_of_alleles; ++allele_i)
        {
            std::string const& allele = get_allele(allele_i);
            if (allele.size() > static_cast<std::size_t>(max_allele_length))
            { throw BGenError(); }
            p = write_length_followed_by_data(p, end, uint32_t(allele.size()),
                                              allele.data());
        }
        assert(p == end);
        return p;
    }

    struct GenotypeDataBlockWriter
    {
        GenotypeDataBlockWriter(std::vector<byte_t>* buffer1,
                                std::vector<byte_t>* buffer2,
                                Context const& context,
                                int const number_of_bits,
                                double permitted_rounding_error = 0.0005)
            : m_buffer1(buffer1)
            , m_buffer2(buffer2)
            , m_context(context)
            , m_layout(m_context.flags & e_Layout)
            , m_number_of_bits(number_of_bits)
            , m_layout1_writer()
            , m_layout2_writer(number_of_bits, permitted_rounding_error)
            , m_writer(0)
        {
            assert(m_buffer1 != 0 && m_buffer2 != 0);
            assert(m_layout == e_Layout1 || m_layout == e_Layout2);
            if (m_layout == e_Layout1) { m_writer = &m_layout1_writer; }
            else if (m_layout == e_Layout2)
            {
                m_writer = &m_layout2_writer;
            }
            else
            {
                assert(0);
            }
        }

        void initialise(uint32_t nSamples, uint16_t nAlleles)
        {
            assert(nSamples == m_context.number_of_samples);
            std::size_t const max_ploidy = 15;
            std::size_t const buffer_size =
                (m_layout == e_Layout1)
                    ? (6 * nSamples)
                    : (10 + nSamples
                       + (((nSamples
                            * (impl::n_choose_k(max_ploidy + nAlleles - 1,
                                                std::size_t(nAlleles) - 1)
                               - 1)
                            * m_number_of_bits)
                           + 7)
                          / 8));
            ;
            m_buffer1->resize(buffer_size);
            m_writer->initialise(nSamples, nAlleles, &(*m_buffer1)[0],
                                 &(*m_buffer1)[0] + buffer_size);
        }

        bool set_sample(std::size_t i) { return m_writer->set_sample(i); }

        void set_number_of_entries(uint32_t ploidy, uint32_t number_of_entries,
                                   OrderType const order_type,
                                   ValueType const value_type)
        {
            m_writer->set_number_of_entries(ploidy, number_of_entries,
                                            order_type, value_type);
        }

        void set_value(uint32_t entry_i, genfile::MissingValue const value)
        {
            m_writer->set_value(entry_i, value);
        }

        void set_value(uint32_t entry_i, double const value)
        {
            m_writer->set_value(entry_i, value);
        }

        void finalise()
        {
            m_writer->finalise();
            // Sanity check: did we get the size right?
            assert(m_writer->repr().first == &(*m_buffer1)[0]);
            assert(
                (m_writer->repr().second >= m_writer->repr().first)
                && std::size_t(m_writer->repr().second - m_writer->repr().first)
                       <= m_buffer1->size());
            uLongf const uncompressed_data_size =
                (m_writer->repr().second - m_writer->repr().first);

#if DEBUG_BGEN_FORMAT
            std::cerr << (m_writer->repr().first) << "  :"
                      << m_writer->repr().second << ", diff = "
                      << (m_writer->repr().second - m_writer->repr().first)
                      << "\n";
            std::cerr << "expected " << uncompressed_data_size << "\n";
#endif
            uint32_t const compressionType =
                (m_context.flags & e_CompressedSNPBlocks);
            if (compressionType != e_NoCompression)
            {
#if HAVE_ZLIB
                std::size_t offset = (m_layout == e_Layout2) ? 8 : 4;
                if (compressionType == e_ZlibCompression)
                {
                    zlib_compress(&(*m_buffer1)[0],
                                  &(*m_buffer1)[0] + uncompressed_data_size,
                                  m_buffer2, offset,
                                  9 // highest compression setting.
                    );
                }
                else if (compressionType == e_ZstdCompression)
                {
                    assert(0);
                    /*
                    zstd_compress(
                        &(*m_buffer1)[0],
                        &(*m_buffer1)[0] + uncompressed_data_size, m_buffer2,
                        offset,
                        17 // reasonable balance between speed and compression.
                    );
                    */
                }
                else
                {
                    assert(0);
                }
                // compression_buffer_size is now the compressed length of the
                // data. Now write total compressed data size to the start of
                // the buffer, including the uncompressed data size if we are in
                // layout 1.2.
                if (m_layout == e_Layout2)
                {
                    write_little_endian_integer(
                        &(*m_buffer2)[0], &(*m_buffer2)[0] + 4,
                        uint32_t(m_buffer2->size()) - 4);
                    write_little_endian_integer(
                        &(*m_buffer2)[0] + 4, &(*m_buffer2)[0] + 8,
                        uint32_t(uncompressed_data_size));
                }
                else
                {
                    write_little_endian_integer(
                        &(*m_buffer2)[0], &(*m_buffer2)[0] + 4,
                        uint32_t(m_buffer2->size()) - 4);
                }

                m_result = std::make_pair(&(*m_buffer2)[0],
                                          &(*m_buffer2)[0] + m_buffer2->size());
#else
                assert(0);
#endif
            }
            else
            {
                // Copy uncompressed data to compression buffer
                // This is inefficient but is not expected to be used much, so
                // not important.
                std::size_t offset = (m_layout == e_Layout2) ? 4 : 0;
                m_buffer2->resize(m_buffer1->size() + offset);
                if (m_layout == e_Layout2)
                {
                    write_little_endian_integer(
                        &(*m_buffer2)[0], &(*m_buffer2)[0] + 4,
                        uint32_t(uncompressed_data_size));
                }
                std::copy(&(*m_buffer1)[0],
                          &(*m_buffer1)[0] + uncompressed_data_size,
                          &(*m_buffer2)[0] + offset);
                m_result = std::make_pair(
                    &(*m_buffer2)[0],
                    &(*m_buffer2)[0] + uncompressed_data_size + offset);
            }
        }

        std::pair<byte_t const*, byte_t const*> repr() const
        {
            return m_result;
        }

    private:
        std::vector<byte_t>* m_buffer1;
        std::vector<byte_t>* m_buffer2;
        Context const& m_context;
        uint32_t const m_layout;
        std::size_t m_number_of_bits;
        v11::ProbabilityDataWriter m_layout1_writer;
        v12::ProbabilityDataWriter m_layout2_writer;
        impl::ProbabilityDataWriterBase* m_writer;
        std::pair<byte_t const*, byte_t const*> m_result;
    };
}
}

#endif
