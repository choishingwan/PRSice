#ifndef GenotypePool_HPP
#define GenotypePool_HPP
#include <cstddef>
#include <cstdint>
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
    // Methods for the storage of the item.
    uintptr_t* get_storage() { return m_geno_start; }
};


class GenotypePool
{
private:
    class MemoryPool
    {
        // actual memory storage
        std::vector<uintptr_t> m_genotype_collections;
        // storing the pointer location as a linked list, this allow us to find
        // the next empty space
        std::unique_ptr<IndividualGenotype[]> m_genotype_list;
        // pointer to the next memory storage
        std::unique_ptr<MemoryPool> m_next;

    public:
        MemoryPool(const size_t num_snps, const size_t memory_per_snp)
            : m_genotype_list(new IndividualGenotype[num_snps])

        {
            m_genotype_collections.resize(num_snps * memory_per_snp, 0);
            // initialize the linked list
            for (size_t i = 1; i < num_snps; ++i)
            { m_genotype_list[i - 1].set_next_item(&m_genotype_list[i]); }
            m_genotype_list[num_snps - 1].set_next_item(nullptr);
        }
        // return pointer to the start of the linked list, can be used to
        // appoint to the free list pointer
        IndividualGenotype* get_genotypes() { return m_genotype_list.get(); }
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
    GenotypePool(size_t num_snps, size_t memory_per_snp)
        : m_num_snps(num_snps)
        , m_memory_per_snp(memory_per_snp)
        , m_memory_pool(new MemoryPool(m_num_snps, m_memory_per_snp))
        , m_free_list(m_memory_pool->get_genotypes())
    {
    }
    uintptr_t* alloc()
    {
        if (m_free_list == nullptr)
        {
            std::unique_ptr<MemoryPool> new_pool(
                new MemoryPool(m_num_snps, m_memory_per_snp));
            new_pool->set_next_collection(std::move(m_memory_pool));
            m_memory_pool.reset(new_pool.release());
            m_free_list = m_memory_pool->get_genotypes();
        }
        IndividualGenotype* current_item = m_free_list;
        // Update the free list to the next free item.
        m_free_list = current_item->get_next_item();
        // Get the storage for T.
        uintptr_t* result = current_item->get_storage();
        return result;
    }
    void free(IndividualGenotype* t)
    {
        // Convert this pointer to T to its enclosing pointer of an item of the
        // arena.
        IndividualGenotype* current_item = t;
        // Add the item at the beginning of the free list.
        current_item->set_next_item(m_free_list);
        m_free_list = current_item;
    }
};
#endif // GenotypePool_HPP
