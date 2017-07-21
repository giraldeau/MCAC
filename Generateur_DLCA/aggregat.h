#ifndef AGGREGAT_H
#define AGGREGAT_H

#include "Sphere.h"

class Aggregat;
class ListAggregat;

class Aggregat
{
    private:
        Sphere* InclusiveSphere;
        Aggregat* parents[2];
        Aggregat* son;
        PhysicalModel* physicalmodel;

        double creation_date;
        int Label;
        int Nc;     //Number of contacts
        int Np;     //Number of spheres

        double *rg;  //Gyration Radius
        double *dm;  //Mobility Diameter
        double *lpm; //Mean Free Path
        double *time_step;
        double *rmax;                   //Radius of the sphere containing Agg
        double *volAgregat;             //Etimation of the Aggregate's volume
        double *surfAgregat;            //Estimation of the sufrace of the aggregate
        double *Tv;                     //Taux de recouvrement volumique
        double *volAgregat_without_cov; //Volume of the aggregate without considering the spheres covering each other
        double *cov;                    //Covering Parameter
        double *ratio_surf_vol;         //Ratio surface / volume
        double *free_surface;           //Free surface of the aggregate (without covering)

        double *x,*y,*z; // position of the gravity center

        ListSphere myspheres;
    public:
        Aggregat(Sphere);
        Aggregat(Aggregat Agg1, Aggregat Agg2);

        void Init(void);
        const std::array<double, 4> Position(void);


    /* Storage specific */

    public:
        Aggregat(void);
        Aggregat(ListSphere& Storage, const int id);
        Aggregat(PhysicalModel&);
        ~Aggregat(void);

        Aggregat(PhysicalModel&, const double x, const double y, const double z, const double r);
        Aggregat(PhysicalModel&, const double* position, const double r);
        Aggregat(Sphere&);

    private:
        std::array< std::vector<double>, 16>* Storage;
        ListAggregat* external_storage;

        void setpointers(void);
        void add(void);
        double operator[](const int);


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

        std::array< std::vector<double>, 16> Storage;
        std::vector < int > index;

        void setpointers();

};


#endif // AGGREGAT_H
