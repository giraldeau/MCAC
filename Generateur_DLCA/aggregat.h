#ifndef AGGREGAT_H
#define AGGREGAT_H

#include "Sphere.h"
#include <storage.h>
#include <list>

using namespace std;

class Aggregate;
class ListAggregat;
class Verlet;

class Aggregate : public storage_elem<16,ListAggregat>
{

    friend class ListAggregat;

    private:
        PhysicalModel* physicalmodel;

        Sphere* InclusiveSphere;
        ListSphere myspheres;

        array<Aggregate*,2> parents;
        Aggregate* son;

        Verlet* verlet;
        array<int, 4> IndexVerlet;

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

        int Nc;     //Number of contacts
        int Np;     //Number of spheres

        bool InVerlet;

    public:

        void Init(void);
        void Init(PhysicalModel& ,Verlet&,const array<double, 4> position ,const int _label, ListSphere&,double r);

        void SetPosition(const array<double, 4> position);
        void SetPosition(const double x,const double y,const double z);
        void Translate(const array<double, 4> vector);
        void Translate(const double vector[]);

        const array<double, 4> GetPosition(void) const;
        array<int, 4> VerletIndex();

        void Update();
        void UpdatesSpheres(ListSphere&, int indexInStorage[]);
        void RayonGiration(void);
        double Distance_Aggregate(Aggregate&, array<double,4> vectorOther, array<double,4> Vectdir);

        void AfficheVerlet();

        double& operator[](const int var);

        /* Storage specific */
    private:
        void setpointers(void);

    public:
        /** Default constructor in local storage */
        Aggregate(void);
        Aggregate(PhysicalModel&);

        /** Constructor in local storage with initialization */
        Aggregate(PhysicalModel&, const double x, const double y, const double z, const double r);
        Aggregate(PhysicalModel&, const array<double, 4> position, const double r);
        Aggregate(PhysicalModel&, const double position[], const double r);
        Aggregate(const Sphere&);


        /** Constructor with external storage */
        Aggregate(ListSphere& Storage, const int id);
        Aggregate(Sphere);
        Aggregate(Aggregate Agg1, Aggregate Agg2);
        Aggregate(ListAggregat&, const int label);

        /** Copy constructor */
        Aggregate(const Aggregate&);

        /** Move constructor */
        Aggregate (Aggregate&&) noexcept; /* noexcept needed to enable optimizations in containers */

        /** Destructor */
        ~Aggregate(void) noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

        /** Copy assignment operator */
        Aggregate& operator= (const Aggregate& other);

        /** Move assignment operator */
        Aggregate& operator= (Aggregate&& other) noexcept;

};


class Verlet
{
private:
    list<int>**** verletlist;
    int GridDiv;

public:
    void Remove(const int id,const array<int, 4> Index);
    list<int>* GetCell(const int i,const int j,const int k)const;
    void Init(const int GridDiv);
    void destroy(void);

public:
    /** Default constructor */
    Verlet(void);

    /** Copy constructor */
    Verlet(const Verlet&);

    /** Move constructor */
    Verlet (Verlet&&) noexcept; /* noexcept needed to enable optimizations in containers */

    /** Destructor */
    ~Verlet(void) noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

    /** Copy assignment operator */
    Verlet& operator= (const Verlet& other);

    /** Move assignment operator */
    Verlet& operator= (Verlet&& other) noexcept;
};



class ListAggregat : public storage_list<16,Aggregate>
{
    friend class Aggregate;

    private:
        PhysicalModel* physicalmodel;

    public:
        ListSphere spheres;
        Verlet verlet;
        void Init(PhysicalModel&, const int size);


        /* Storage specific */
    private:
        void setpointers();

    public:
        /** Default constructor in local storage */
        ListAggregat(void);
        ListAggregat(PhysicalModel& _physicalmodel, const int size);

        /** Copy constructor */
        ListAggregat(const ListAggregat& other);

        /** Move constructor */
        ListAggregat (ListAggregat&&) noexcept; /* noexcept needed to enable optimizations in containers */

        /** Destructor */
        ~ListAggregat(void) noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

        /** Copy assignment operator */
        ListAggregat& operator= (const ListAggregat& other);

        /** Move assignment operator */
        ListAggregat& operator= (ListAggregat&& other) noexcept;
};



#endif // AGGREGAT_H
