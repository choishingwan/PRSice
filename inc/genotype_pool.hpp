#ifndef GenotypePool_HPP
#define GenotypePool_HPP
#include <cstddef>
#include <cstdint>
#include <plink_common.hpp>
#include <vector>
// modified based on https://thinkingeek.com/2017/11/19/simple-memory-pool/

class IndividualGenotype
{
private:
    IndividualGenotype* m_next;
    uintptr_t* m_geno_start;

public:
    IndividualGenotype() {}
    IndividualGenotype(uintptr_t* current)
    {
        m_geno_start = current;
        m_next = nullptr;
    }
    IndividualGenotype* get_next_item() const { return m_next; }
    void set_next_item(IndividualGenotype* n) { m_next = n; }
    void set_start_location(uintptr_t* i) { m_geno_start = i; }
    uintptr_t* get_geno() { return m_geno_start; }
    // Methods for the storage of the item.
};


class GenotypePool
{
private:
    class MemoryPool
    {
        // actual memory storage
        std::unique_ptr<uintptr_t[]> m_genotypes alignas(alignof(uintptr_t));
        // std::vector<uintptr_t> m_genotypes alignas(alignof(uintptr_t));
        // storing the pointer location as a linked list, this allow us to find
        // the next empty space
        std::unique_ptr<IndividualGenotype[]> m_genotype_list;
        // pointer to the next memory storage
        std::unique_ptr<MemoryPool> m_next;

    public:
        ~MemoryPool() {}
        MemoryPool(const size_t num_snps, const size_t memory_per_snp)
            : m_genotypes(new uintptr_t[num_snps * memory_per_snp])
            , m_genotype_list(new IndividualGenotype[num_snps])

        {
            // m_genotypes.resize(num_snps * memory_per_snp, 0);
            // initialize the linked list
            m_genotype_list[0].set_start_location(m_genotypes.get());
            // m_genotype_list[0].set_start_location(m_genotypes.data());
            for (size_t i = 1; i < num_snps; ++i)
            {
                m_genotype_list[i - 1].set_next_item(&m_genotype_list[i]);
                m_genotype_list[i].set_start_location(m_genotypes.get()
                                                      + i * memory_per_snp);
            }
            m_genotype_list[num_snps - 1].set_next_item(nullptr);
        }
        // return pointer to the start of the linked list, can be used to
        // appoint to the free list pointer
        IndividualGenotype* get_genotype_location()
        {
            return m_genotype_list.get();
        }
        uintptr_t* get_genotype_storage() { return m_genotypes.get(); }
        // connect the new memory pool to the current memory pool
        void set_next_collection(std::unique_ptr<MemoryPool>&& n)
        {
            assert(!m_next);
            m_next.reset(n.release());
        }
    };
    size_t m_num_snps;
    size_t m_memory_per_snp;
    std::unique_ptr<MemoryPool> m_memory_pool;
    IndividualGenotype* m_free_list;

public:
    GenotypePool() {};
    GenotypePool(size_t num_snps, size_t memory_per_snp)
        : m_num_snps(num_snps)
        , m_memory_per_snp(round_up_pow2(memory_per_snp, CACHELINE))
        , m_memory_pool(new MemoryPool(m_num_snps, m_memory_per_snp))
        , m_free_list(m_memory_pool->get_genotype_location())
    {
    }
    IndividualGenotype* alloc()
    {
        if (m_free_list == nullptr)
        {
            std::unique_ptr<MemoryPool> new_pool(
                new MemoryPool(m_num_snps, m_memory_per_snp));
            new_pool->set_next_collection(std::move(m_memory_pool));
            m_memory_pool.reset(new_pool.release());
            m_free_list = m_memory_pool->get_genotype_location();
        }
        IndividualGenotype* current_item = m_free_list;
        // Update the free list to the next free item.
        m_free_list = current_item->get_next_item();
        // Get the storage for T.
        assert(current_item != nullptr);
        return current_item;
    }
    void free(IndividualGenotype* t)
    {
        // Convert this pointer to T to its enclosing pointer of an item of the
        // arena.
        if (t == nullptr)
        { throw std::runtime_error("Error: Can't free null pointer!"); }
        IndividualGenotype* current_item = std::move(t);
        t = nullptr;
        // Add the item at the beginning of the free list.
        current_item->set_next_item(m_free_list);
        m_free_list = current_item;
    }
};
#endif // GenotypePool_HPP
