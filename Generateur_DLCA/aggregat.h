#ifndef AGGREGAT_H
#define AGGREGAT_H

#include "Sphere.h"

class Aggregat
{
    public:
        Aggregat(Sphere);
        Aggregat(Aggregat Agg1, Aggregat Agg2);
    private:
        Sphere* mysphere;
        Aggregat* parents[2];
        Aggregat* son;
        PhysicalModel* physicalmodel;
        double creation_date;
};




class ListAggregat
{
    friend class Aggregat;

    /* Generic */

    private:
        int N;
        std::vector < Aggregat* > Aggregats;
        ListSphere spheres;
        PhysicalModel* physicalmodel;

        ListAggregat(void);
        ListAggregat(PhysicalModel& _physicalmodel, const int N);

        ~ListAggregat(void);
        Aggregat& operator[](const int);

        int size() const;

        /* Storage specific */

        std::array< std::vector<double>, 7> Storage_dble;
        std::array< std::vector<int>, 7> Storage_int;
        std::vector < int > index;

        void setpointers();

};


#endif // AGGREGAT_H
