#include "aggregat.h"

Aggregat::Aggregat(Sphere _mysphere)
{
    mysphere = &_mysphere;
    parents[0] = NULL;
    parents[1] = NULL;
    son = NULL;
    physicalmodel = mysphere->physicalmodel;
    creation_date = physicalmodel->temps;
}

Aggregat::Aggregat(Aggregat Agg1, Aggregat Agg2)
{
    mysphere = new Sphere(*(Agg1.physicalmodel));

    parents[0] = &Agg1;
    parents[1] = &Agg2;
    son = NULL;
    parents[0]->son = this;
    parents[1]->son = this;
    physicalmodel = mysphere->physicalmodel;
    creation_date = physicalmodel->temps;
}
