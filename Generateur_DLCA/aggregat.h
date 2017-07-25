#ifndef AGGREGAT_H
#define AGGREGAT_H

#include "Sphere.h"
#include <list>

using namespace std;

class Aggregate;
class ListAggregat;
class Verlet;

class Aggregate
{
    private:
        Sphere* InclusiveSphere;
        Aggregate* parents[2];
        Aggregate* son;
        PhysicalModel* physicalmodel;
        Verlet* verlet;
        array<int, 4> IndexVerlet;

        double creation_date;
        int Label;
        int Nc;     //Number of contacts
        int Np;     //Number of spheres
        bool InVerlet;

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
        Aggregate(Sphere);
        Aggregate(Aggregate Agg1, Aggregate Agg2);

        void Init(void);
        const array<double, 4> Position(void);
        void Position(const array<double, 4> position);
        void Position(const double x,const double y,const double z);
        void Translate(const array<double, 4> vector);
        void Translate(const double* vector);
        array<int, 4> VerletIndex();
        void Init(PhysicalModel& ,Verlet&,const array<double, 4> position ,const int _label, ListSphere&,double r);
        void UpdatesSpheres(ListSphere&, int* index);


    /* Storage specific */

    public:
        Aggregate(void);
        Aggregate(ListSphere& Storage, const int id);
        Aggregate(PhysicalModel&);
        ~Aggregate(void);

        Aggregate(PhysicalModel&, const double x, const double y, const double z, const double r);
        Aggregate(PhysicalModel&, const double* position, const double r);
        Aggregate(Sphere&);

        void AfficheVerlet();


    private:
        array< vector<double>, 16>* Storage;
        ListAggregat* external_storage;

        void setpointers(void);
        void add(void);
        double operator[](const int);


};


/*

class ListAggregat
{
    friend class Aggregate;

    private:
        int N;
        vector < Aggregate* > Aggregats;
        ListSphere spheres;
        PhysicalModel* physicalmodel;

        ListAggregat(void);
        ListAggregat(PhysicalModel& _physicalmodel, const int N);

        ~ListAggregat(void);
        Aggregate& operator[](const int);

        int size() const;

        array< vector<double>, 16> Storage;
        vector < int > index;

        void setpointers();

};

*/

class Verlet
{
public:
    void Remove(const int id,const array<int, 4> Index);
    list<int>* GetCell(const int i,const int j,const int k)const;
    void Init(const int GridDiv);

    private:
    list<int>**** verletlist;

};

#endif // AGGREGAT_H
