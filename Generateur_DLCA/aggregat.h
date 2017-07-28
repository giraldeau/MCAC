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
        PhysicalModel* physicalmodel;

        Sphere* InclusiveSphere;
        ListSphere myspheres;

        array<Aggregate*,2> parents;
        Aggregate* son;

        Verlet* verlet;
        array<int, 4> IndexVerlet;

        array< vector<double>, 16>* Storage;
        ListAggregat* external_storage;

        double creation_date;
        double nctmp,nptmp;

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

        int Label;
        int Nc;     //Number of contacts
        int Np;     //Number of spheres

        bool InVerlet;

    public:
        Aggregate(Sphere);
        Aggregate(Aggregate Agg1, Aggregate Agg2);

        void Init(void);
        const array<double, 4> GetPosition(void) const;
        void SetPosition(const array<double, 4> position);
        void SetPosition(const double x,const double y,const double z);
        void Translate(const array<double, 4> vector);
        void Translate(const double vector[]);
        array<int, 4> VerletIndex();
        void Init(PhysicalModel& ,Verlet&,const array<double, 4> position ,const int _label, ListSphere&,double r);
        void UpdatesSpheres(ListSphere&, int index[]);
        void ReplacePosi();
        void RayonGiration(void);
        double Distance_Aggregate(Aggregate&, array<double,4> vectorOther, array<double,4> Vectdir);
        void Update();


    /* Storage specific */

    public:
        Aggregate(void);
        Aggregate(ListSphere& Storage, const int id);
        Aggregate(PhysicalModel&);
        ~Aggregate(void);

        Aggregate(PhysicalModel&, const double x, const double y, const double z, const double r);
        Aggregate(PhysicalModel&, const double position[], const double r);
        Aggregate(Sphere&);

        void AfficheVerlet();
        double& operator[](const int);


    private:

        void setpointers(void);
        void add(void);
        Aggregate& operator=(const Aggregate& other);
        Aggregate& operator=(const Aggregate&& other) noexcept;
        Aggregate(const Aggregate& other);


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
    ~Verlet(void);
    void destroy(void);

    private:
    list<int>**** verletlist;
    int GridDiv;
};

#endif // AGGREGAT_H
