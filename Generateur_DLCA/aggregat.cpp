#include "aggregat.h"
#include <math.h>
#include <cmath>

Aggregate::Aggregate(void)
{
    Storage = new array< vector<double>, 16>;
    for (int j=0;j<=15;j++)
        (*Storage)[j].assign(1, 0.);
    external_storage=NULL;

    InclusiveSphere = new Sphere;
    parents[0] = NULL;
    parents[1] = NULL;
    son = NULL;
    physicalmodel = NULL;
    Label = 0;
    creation_date = 0.;
    Nc = 0.;

    Init();
}

Aggregate::Aggregate(PhysicalModel& _physicalmodel)
{
    Storage = new array< vector<double>, 16>;
    for (int j=0;j<=15;j++)
        (*Storage)[j].assign(1, 0.);
    external_storage=NULL;

    InclusiveSphere = new Sphere;
    parents[0] = NULL;
    parents[1] = NULL;
    son = NULL;
    physicalmodel = &_physicalmodel;
    Label = 0;
    creation_date = 0.;
    Nc = 0.;

    Init();
}

void Aggregate::Init(void)
{
    setpointers();

    *rg=0.;  //Gyration Radius
    *dm=0.;  //Mobility Diameter
    *lpm=0.; //Mean Free Path
    *time_step=0.;
    *rmax=0.;                   //Radius of the sphere containing Agg
    *volAgregat=0.;             //Etimation of the Aggregate's volume
    *surfAgregat=0.;            //Estimation of the sufrace of the aggregate
    *Tv=0.;                     //Taux de recouvrement volumique
    *volAgregat_without_cov=0.; //Volume of the aggregate without considering the spheres covering each other
    *cov=0.;                    //Covering Parameter
    *ratio_surf_vol=0.;         //Ratio surface / volume
    *free_surface=0.;           //Free surface of the aggregate (without covering)

    *x=0.;
    *y=0.;
    *z=0.;

    IndexVerlet = {0,0,0,0};
}


void Aggregate::setpointers(void)
{
    rg=&(*Storage)[0][Label];  //Gyration Radius
    dm=&(*Storage)[1][Label];  //Mobility Diameter
    lpm=&(*Storage)[2][Label]; //Mean Free Path
    time_step=&(*Storage)[3][Label];
    rmax=&(*Storage)[4][Label];                   //Radius of the sphere containing Agg
    volAgregat=&(*Storage)[5][Label];             //Etimation of the Aggregate's volume
    surfAgregat=&(*Storage)[6][Label];            //Estimation of the sufrace of the aggregate
    Tv=&(*Storage)[7][Label];                     //Taux de recouvrement volumique
    volAgregat_without_cov=&(*Storage)[8][Label]; //Volume of the aggregate without considering the spheres covering each other
    cov=&(*Storage)[9][Label];                    //Covering Parameter
    ratio_surf_vol=&(*Storage)[10][Label];         //Ratio surface / volume
    free_surface=&(*Storage)[11][Label];           //Free surface of the aggregate (without covering)
    x=&(*Storage)[12][Label];
    y=&(*Storage)[13][Label];
    z=&(*Storage)[14][Label];
}

void Aggregate::Init(PhysicalModel& _physicalmodel,list<int>****& _Verlet,const array<double, 4> position ,const int _label)
{
    physicalmodel = &_physicalmodel;
    Verlet = &_Verlet;
    Label = _label;

    if (physicalmodel->use_verlet)
    {
        //(*Verlet[IndexVerlet[1]][IndexVerlet[2]][IndexVerlet[3]])->push_front(Label);
    }
    Position(position);
}

const array<double, 4> Aggregate::Position(void)
{
    array<double, 4> mypos;
    mypos[1]=*x;
    mypos[2]=*y;
    mypos[3]=*z;
    return mypos;
}

void Aggregate::Position(const double newx,const double newy,const double newz)
{
    *x = newx;
    *y = newy;
    *z = newz;
    /*
    if (physicalmodel->use_verlet)
    {
        //$ Update Verlet
        array<int, 4> newindexVerlet = VerletIndex();

        if (VerletIndex() != IndexVerlet)
        {
            (*Verlet[IndexVerlet[1]][IndexVerlet[2]][IndexVerlet[3]])->remove(Label);
            IndexVerlet = newindexVerlet;
            (*Verlet[IndexVerlet[1]][IndexVerlet[2]][IndexVerlet[3]])->push_front(Label);
        }
    }
    */
}

void Aggregate::Position(const array<double, 4> position)
{
    Position(position[1], position[2], position[3]);
}

void Aggregate::Translate(const array<double, 4> vector)
{
    Position(*x +vector[1], *y + vector[2], *z +vector[3]);
}

void Aggregate::Translate(const double* vector)
{
    Position(*x +vector[1], *y + vector[2], *z +vector[3]);
}

Aggregate::~Aggregate(void)
{
    if (physicalmodel->use_verlet)
    {
        //(*Verlet[IndexVerlet[1]][IndexVerlet[2]][IndexVerlet[3]])->remove(Label);
    }

    if (external_storage==NULL)
    {
        delete Storage;
    }

    external_storage=NULL;

    delete InclusiveSphere;
    parents[0] = NULL;
    parents[1] = NULL;
    son = NULL;
    physicalmodel = NULL;
    Label = 0;
    creation_date = 0.;
    Nc = 0.;

    rg=NULL;
    dm=NULL;
    lpm=NULL;
    time_step=NULL;
    rmax=NULL;
    volAgregat=NULL;
    surfAgregat=NULL;
    Tv=NULL;
    volAgregat_without_cov=NULL;
    cov=NULL;
    ratio_surf_vol=NULL;
    free_surface=NULL;
    x=NULL;
    y=NULL;
    z=NULL;



}


array<int, 4> Aggregate::VerletIndex()
{
    array<int, 4>  index= {0.,0.,0.,0.};
    double step = physicalmodel->GridDiv/physicalmodel->L;
    int origin = physicalmodel->GridDiv+1;

    index[1]=floor((*x)*step)+origin;
    index[2]=floor((*y)*step)+origin;
    index[3]=floor((*z)*step)+origin;
    return index;
}

/*
Aggregat::Aggregat(Sphere _mysphere)
{
    InclusiveSphere = &_mysphere;
    parents[0] = NULL;
    parents[1] = NULL;
    son = NULL;
    physicalmodel = InclusiveSphere->physicalmodel;
    creation_date = physicalmodel->temps;
}

Aggregat::Aggregat(Aggregat Agg1, Aggregat Agg2)
{
    InclusiveSphere = new Sphere(*(Agg1.physicalmodel));

    parents[0] = &Agg1;
    parents[1] = &Agg2;
    son = NULL;
    parents[0]->son = this;
    parents[1]->son = this;
    physicalmodel = InclusiveSphere->physicalmodel;
    creation_date = physicalmodel->temps;
}


Sphere::Sphere(PhysicalModel& _physicalmodel)
{
    Storage = new array< vector<double>, 15>;

    for (int j=0;j<=6;j++)
        (*Storage)[j].assign(1, 0.);

    external_storage=NULL;
    physicalmodel = &_physicalmodel;
    AggLabel = 0;
    SphereLabel = 0;

    Init();
}

Sphere::Sphere(ListSphere& aggregat,const int id)
{
    external_storage =&aggregat;

    Storage = aggregat.Storage;
    physicalmodel = aggregat.physicalmodel;
    AggLabel = 0;
    SphereLabel = id;

    Init();

    external_storage->setpointers();
}


Sphere::Sphere(PhysicalModel& _physicalmodel, const double newx,const double newy,const double newz,const double newr)
{
    Storage = new array< vector<double>, 7>;

    for (int j=0;j<=6;j++)
        (*Storage)[j].assign(1, 0.);

    external_storage=NULL;
    physicalmodel = &_physicalmodel;
    AggLabel = 0;
    SphereLabel = 0;

    Init();

    *x = newx;
    *y = newy;
    *z = newz;
    *r = newr;
    UpdateVolAndSurf();
}

Sphere::Sphere(PhysicalModel& _physicalmodel, const double* newp,const double newr) : Sphere(_physicalmodel,newp[1],newp[2],newp[3],newr){}

Sphere::Sphere(Sphere& c) : Sphere(*(c.physicalmodel), *c.x,*c.y,*c.z,*c.r){}

void Sphere::add(void)
{
    for (int j=0;j<=6;j++)
        (*Storage)[j].push_back(0.);
    SphereLabel = (*Storage)[0].size()-1;

    Init();
}

void Sphere::Init(void)
{
    setpointers();

    *x = 0.;
    *y = 0;
    *z = 0;
    *r = 0;
    *volume = 0.;
    *surface = 0.;
}


Sphere::~Sphere(void)
{
    for (int j=0;j<=6;j++)
        (*Storage)[j].erase((*Storage)[j].begin() + SphereLabel);

    //setpointers();
    if(external_storage!=NULL)
        external_storage->setpointers();
    else if (Storage!=NULL)
        delete Storage;

    external_storage=NULL;
    SphereLabel = 0;
    AggLabel = 0;

    x = NULL;
    y = NULL;
    z = NULL;
    r = NULL;
    volume = NULL;
    surface = NULL;
    physicalmodel = NULL;

}
*/
