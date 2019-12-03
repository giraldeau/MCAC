#include "aggregats/aggregat_list.hpp"
#include "tools.hpp"
#include <iostream>

using namespace std;
namespace MCAC {
void ListAggregat::setpointers() {
    auto newdeb((*Storage)[0].begin());
    auto newfin((*Storage)[0].end());
    if ((newdeb == ptr_deb) && (newfin == ptr_fin)) {
        return;
    }
    for (Aggregate *aggregate : list) {
        aggregate->setpointers();
    }
    ptr_deb = newdeb;
    ptr_fin = newfin;
}
ListAggregat::ListAggregat(PhysicalModel *the_physical_model) noexcept:
    storage_list<AggregatesFields::AGGREGAT_NFIELDS, Aggregate>(),
    physicalmodel(the_physical_model),
    maxradius(0.),
    avg_npp(1),
    max_time_step(0.),
    index_sorted_time_steps(),
    cumulative_time_steps(),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    writer(new ThreadedIO(*physicalmodel, physicalmodel->N)),
    last_saved(0),
    spheres(*the_physical_model, physicalmodel->N),
    verlet() {
    verlet.Init(the_physical_model->GridDiv, the_physical_model->L);
    storage_list<AggregatesFields::AGGREGAT_NFIELDS, Aggregate>::Init(physicalmodel->N, *this);
    setpointers();

    //Initialize the aggregates

    size_t testmem = 0;
    setpointers();
    for (size_t i = 0; i < size(); i++) {

        //random size
        double x = Random();
        double dp = 0;
        if (physicalmodel->Mode == 1) {
            dp = physicalmodel->Dpm + sqrt(2.0) * physicalmodel->sigmaDpm * inverf(2.0 * x - 1.0); //Loi normale
        } else {
            dp = exp(log(physicalmodel->Dpm)
                     + sqrt(2.0) * log(physicalmodel->sigmaDpm) * inverf(2.0 * x - 1.0)); //Loi log-normale
        }
        dp = dp / 1E9;
        if (dp <= 0) {
            dp = physicalmodel->Dpm * 1E-9;
        }
        bool placed = false;
        while (!placed) {
            //random position
            array<double, 3> newpos{{Random() * physicalmodel->L,
                                     Random() * physicalmodel->L,
                                     Random() * physicalmodel->L}};

            //++++++++++++ Test de superposition des sphérules lors de leur génération aléatoire ++++++++++++
            if (test_free_space(newpos, dp)) {
                list[i]->init(*physicalmodel, verlet, newpos, i, spheres, dp);
                placed = true;
            } else {
                i--;
                testmem++;
            }
            if (testmem > size()) {
                cout << "Impossible de générer tous les monomères sans superposition." << endl;
                cout << "La fraction volumique doit être diminuée." << endl;
                exit(0);
            }
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        }
    }
}
ListAggregat::~ListAggregat() noexcept {
    delete writer;

    //#pragma omp simd
    for (Aggregate *aggregate : list) {
        aggregate->unset_verlet();
    }
}
}// namespace MCAC

