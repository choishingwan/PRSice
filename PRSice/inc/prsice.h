#ifndef PRSICE_H
#define PRSICE_H

#include <string>
#include <fstream>
#include <armadillo>
#include <vector>

class PRSice
{
    public:
        PRSice();
        virtual ~PRSice();
        void run();
    protected:
    private:
        void prs_clump();
        void prs_prune();
};

#endif // PRSICE_H
