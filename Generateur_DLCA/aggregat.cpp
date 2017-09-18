#include "aggregat.h"
#include <math.h>
#include <cmath>
#include <iomanip>


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW2(a) ((a)*(a))
#define POW3(a) ((a)*(a)*(a))

using namespace std;

const double PI = atan(1.0)*4;

Aggregate::Aggregate(void):
    storage_elem<13,ListAggregat>(),
    physicalmodel(new PhysicalModel),
    myspheres(),
    verlet(nullptr),
    IndexVerlet({{0,0,0}}),
    distances(),
    volumes(),
    surfaces(),
    rg(nullptr),
    dm(nullptr),
    lpm(nullptr),
    time_step(nullptr),
    rmax(nullptr),
    volAgregat(nullptr),
    surfAgregat(nullptr),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    Np(0),
    InVerlet(false)
{
    Init();
}

Aggregate::Aggregate(PhysicalModel& _physicalmodel) :
    storage_elem<13,ListAggregat>(),
    physicalmodel(&_physicalmodel),
    myspheres(),
    verlet(nullptr),
    IndexVerlet({{0,0,0}}),
    distances(),
    volumes(),
    surfaces(),
    rg(nullptr),
    dm(nullptr),
    lpm(nullptr),
    time_step(nullptr),
    rmax(nullptr),
    volAgregat(nullptr),
    surfAgregat(nullptr),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    Np(0),
    InVerlet(false)
{
    Init();
}


Aggregate::Aggregate(ListAggregat& _storage, const int _N):
    storage_elem<13,ListAggregat>(_storage, _N),
    physicalmodel(_storage.physicalmodel),
    myspheres(),
    verlet(nullptr),
    IndexVerlet({{0,0,0}}),
    distances(),
    volumes(),
    surfaces(),
    rg(nullptr),
    dm(nullptr),
    lpm(nullptr),
    time_step(nullptr),
    rmax(nullptr),
    volAgregat(nullptr),
    surfAgregat(nullptr),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    Np(0),
    InVerlet(false)
{
    Init();
}

void Aggregate::Init(void)
{
    if(!(external_storage==nullptr))
        external_storage->setpointers();
    setpointers();

    *rg=0.;  //Gyration Radius
    *dm=0.;  //Mobility Diameter
    *lpm=0.; //Mean Free Path
    *time_step=0.;
    *rmax=0.;                   //Radius of the sphere containing Agg
    *volAgregat=0.;             //Etimation of the Aggregate's volume
    *surfAgregat=0.;            //Estimation of the sufrace of the aggregate

    *rx = 0;
    *ry = 0;
    *rz = 0;

    IndexVerlet = {{0,0,0}};

    SetPosition(0.,0.,0.);
    UpdateDistances();
}


void Aggregate::setpointers(void)
{
    rg=&(*Storage)[0][indexInStorage];  //Gyration Radius
    dm=&(*Storage)[1][indexInStorage];  //Mobility Diameter
    lpm=&(*Storage)[2][indexInStorage]; //Mean Free Path
    time_step=&(*Storage)[3][indexInStorage];
    rmax=&(*Storage)[4][indexInStorage];                   //Radius of the sphere containing Agg
    volAgregat=&(*Storage)[5][indexInStorage];             //Etimation of the Aggregate's volume
    surfAgregat=&(*Storage)[6][indexInStorage];            //Estimation of the sufrace of the aggregate
    x=&(*Storage)[7][indexInStorage];
    y=&(*Storage)[8][indexInStorage];
    z=&(*Storage)[9][indexInStorage];
    rx=&(*Storage)[10][indexInStorage];
    ry=&(*Storage)[11][indexInStorage];
    rz=&(*Storage)[12][indexInStorage];
}

void Aggregate::Init(PhysicalModel& _physicalmodel,Verlet& _verlet,const array<double, 3> position ,const int _label, ListSphere& spheres,const double Dp)
{
    if(physicalmodel->toBeDestroyed)
        delete physicalmodel;
    physicalmodel = &_physicalmodel;

    verlet = &_verlet;
    indexInStorage = _label;

    if(!(external_storage==nullptr))
        external_storage->setpointers();
    setpointers();

    if (physicalmodel->use_verlet)
    {
        (*verlet)[IndexVerlet[0]][IndexVerlet[1]][IndexVerlet[2]].push_front(_label);
        InVerlet=true;
    }
    SetPosition(position);
    spheres[_label].SetLabel(_label);
    spheres[_label].InitVal(position, Dp/2);
    myspheres = ListSphere(spheres,{_label});
    Np = myspheres.size();

    UpdateDistances();
    *dm = Dp;                   //Diamètre de mobilité
    Update();

}


//############################################# Calcul de la distance inter-agrégats ############################################
bool Aggregate::Contact(Aggregate& other) const noexcept
{
/*
    Sphere SphereMe(*this);
    Sphere SphereOther(other);

    if (! SphereMe.Contact(SphereOther))
        return false;
*/
    //$ Loop on all the spheres of the other aggregate
    for (const Sphere* othersphere : other.myspheres)
    {
        //$ For every sphere in the aggregate :
        for (const Sphere* mysphere: myspheres)
        {
            if (mysphere->Contact(*othersphere))
                return true;
        }
    }
    return false;
}

double Aggregate::Distance(Aggregate& other,array<double,3> Vectdir) const
{
    double mindist(*lpm);
    bool contact(false);

    //$ Loop on all the spheres of the other aggregate
    for (const Sphere* othersphere : other.myspheres)
    {
        //$ For every sphere in the aggregate :
        for (const Sphere* mysphere: myspheres)
        {
            double dist=mysphere->Collision(*othersphere, Vectdir);
            if (0. <= dist && dist <= mindist)
            {
                mindist = dist;
                contact = true;
            }
        }
    }
    if (contact)
        return mindist;
    else
        return -1.;
}
//###############################################################################################################################

const array<double, 3> Aggregate::GetPosition(void) const noexcept
{
    array<double, 3> mypos;
    mypos[0]=*x;
    mypos[1]=*y;
    mypos[2]=*z;
    return mypos;
}

void Aggregate::SetPosition(const double newx,const double newy,const double newz) noexcept
{
    IndexVerlet = GetVerletIndex();

    *x = periodicPosition(newx,physicalmodel->L);
    *y = periodicPosition(newy,physicalmodel->L);
    *z = periodicPosition(newz,physicalmodel->L);

    if (InVerlet && physicalmodel->use_verlet)
    {
        //$ Update Verlet
        array<int, 3> newindexVerlet = GetVerletIndex();

        if (newindexVerlet != IndexVerlet)
        {
            (*verlet)[IndexVerlet[0]][IndexVerlet[1]][IndexVerlet[2]].remove(indexInStorage);
            IndexVerlet = newindexVerlet;
            (*verlet)[IndexVerlet[0]][IndexVerlet[1]][IndexVerlet[2]].push_front(indexInStorage);
        }
    }

}

void Aggregate::SetPosition(const array<double, 3> position) noexcept
{
    SetPosition(position[0], position[1], position[2]);
}

void Aggregate::Translate(const array<double, 3> vector) noexcept
{
    // keep track of the previous position
    array<double,3> oldpos = GetPosition();

    // move the aggregate
    SetPosition(*x +vector[0], *y + vector[1], *z +vector[2]);

    // Real movement of the aggregate (periodic)
    array<double,3> newpos = GetPosition();
    newpos[0] -= oldpos[0];
    newpos[1] -= oldpos[1];
    newpos[2] -= oldpos[2];

    // move the first sphere taking care of the periodicity
    *myspheres[0].x = periodicPosition(*myspheres[0].x + newpos[0],physicalmodel->L);
    *myspheres[0].y = periodicPosition(*myspheres[0].y + newpos[1],physicalmodel->L);
    *myspheres[0].z = periodicPosition(*myspheres[0].z + newpos[2],physicalmodel->L);

    // move all the other sphere relatively to the first
    for (Sphere* mysphere : myspheres)
    {
        mysphere->SetPosition(*myspheres[0].x+*mysphere->rx,
                              *myspheres[0].y+*mysphere->ry,
                              *myspheres[0].z+*mysphere->rz);
    }
}

Aggregate::~Aggregate(void) noexcept
{
    if(physicalmodel->toBeDestroyed)
        delete physicalmodel;

    if (InVerlet && physicalmodel->use_verlet)
    {
        (*verlet)[IndexVerlet[0]][IndexVerlet[1]][IndexVerlet[2]].remove(indexInStorage);
    }
}

 __attribute__((pure)) double Aggregate::GetLpm() const noexcept
{
    return *lpm;
}

 __attribute__((pure)) double Aggregate::GetVolAgregat() const noexcept
{
    return *volAgregat;
}

 __attribute__((pure)) int Aggregate::GetLabel() const noexcept
{
    return indexInStorage;
}

//######### Mise à jour des paramètres physiques d'un agrégat (rayon de giration, masse, nombre de sphérules primaires) #########
void Aggregate::Update()
{
    // This function will update the parameter of Agg
    Volume();

    MassCenter();

    CalcRadius();

    //$ Determination of the Radius of gyration of Agg using RayonGiration()
    RayonGiration();

    //$ Determination of Dm using ConvertRg2Dm
    *dm = physicalmodel->ConvertRg2DmFromStart(Np,*rg,*dm);

    if (physicalmodel->ActiveModulephysique == 1)
    {
        double masse = physicalmodel->Rho*(*volAgregat);//Determination of the real mass of Agg
        double vit =  physicalmodel->velocity(masse);
        double diff = physicalmodel->diffusivity(*dm);
        *lpm = 8*diff/PI/vit;
        *time_step = *lpm/vit; //Displacement duration
    }
    else
    {
        *lpm = physicalmodel->Dpm*1E-9;
        *time_step = 1E-6;
    }

    if(!(external_storage==nullptr))
    {
        if (*rmax > external_storage->maxradius)
            external_storage->maxradius = *rmax;
    }
}



//############# Calculation of the volume, surface, center of mass and Giration radius of gyration of an aggregate ##############
void Aggregate::Volume(void)
{
    *volAgregat = *surfAgregat = 0.0; // Volume and surface of Agg Id

    //$ Initialisation of the arrays of volume, surface of each sphere, and the center of mass
    volumes.resize(Np);
    surfaces.resize(Np);

    //$ For the Spheres i in Agg Id
    for (int i = 0; i < Np; i++)
    {
        //$ Calculation of the volume and surface of monomere i of Agg id
        volumes[i] = *myspheres[i].volume; //Calculation of the volume of i
        surfaces[i] = *myspheres[i].surface;    //Calculation of the surface of i
    }

    for (int i = 0; i < Np; i++)
    {
        for (int j = i+1; j < Np; j++) //for the j spheres composing Aggregate n°id
        {

            double voli, volj, surfi, surfj;
            voli = volj = surfi = surfj = 0.;

            //$ Calculation of the intersection between the spheres i and j
            myspheres[i].Intersection(myspheres[j],distances[i][j], voli,volj,surfi,surfj);

            //$ The volume and surface covered by j is substracted from those of i
            volumes[i] = volumes[i] - voli;    //Calcul du volume de la sphérule i moins le volume de
                                             //la calotte due à la sphérule j
            surfaces[i] = surfaces[i] - surfi; //Calcul de la surface de la sphérule i moins la surface de
                                             //la calotte due à la sphérule j

            //$ The volume and surface covered by i is substracted from those of j
            volumes[j] = volumes[j] - volj;    //Calcul du volume de la sphérule j moins le volume de
                                             //la calotte due à la sphérule i
            surfaces[j] = surfaces[j] - surfj; //Calcul de la surface de la sphérule j moins la surface de
                                             //la calotte due à la sphérule i
        }
        //$ Calculation of the total volume and surface of the aggregate
        *volAgregat = *volAgregat + volumes[i];    //Total Volume of Agg id
        *surfAgregat = *surfAgregat + surfaces[i]; //Total Surface of Agg id
    }
}

void Aggregate::MassCenter(void)
{
    array<double, 3> newpos({{0.,0.,0.}});

    //$ For the Spheres i in Agg Id
    for (int i = 0; i < Np; i++)
    {
        //$ Calculation of the total volume and surface of the aggregate
        double dx = *myspheres[i].rx-*rx;
        double dy = *myspheres[i].ry-*ry;
        double dz = *myspheres[i].rz-*rz;

        //$ Calculation of the position of the center of mass
        newpos[0] += dx * volumes[i];
        newpos[1] += dy * volumes[i];
        newpos[2] += dz * volumes[i];
    }
    //$ Filling of PosiGravite

    for (int k = 0; k < 3; k++)
        newpos[k] = newpos[k]/(*volAgregat);

    *rx = newpos[0]+*rx;
    *ry = newpos[1]+*ry;
    *rz = newpos[2]+*rz;

    SetPosition(*myspheres[0].x+*rx,
                *myspheres[0].y+*ry,
                *myspheres[0].z+*rz);

    const int loopsize(Np);
    const double _rx(*rx);
    const double _ry(*ry);
    const double _rz(*rz);
    for (int i = 0; i < loopsize; i++)
    {
        double dx = *myspheres[i].rx - _rx;
        double dy = *myspheres[i].ry - _ry;
        double dz = *myspheres[i].rz - _rz;
        distances[i][Np]=sqrt(POW2(dx)+POW2(dy)+POW2(dz));
    }
}

void Aggregate::CalcRadius(void)
{

    *rmax = 0.0; // Maximum radius of the aggregate, this corresponds to the distance between the center of mass of the aggregate and the edge of the furthest ball from said center.
                // It is used to assimilate the aggregate to a sphere when checking for intersections

    for (int i = 0; i < Np; i++)
        *rmax=MAX(*rmax,*myspheres[i].r+distances[i][Np]);
}

void Aggregate::RayonGiration(void)
{
    // This function determines the Gyration Radius of the Aggregate Id.

    // These correspond to the sum of the volumes of each spheres multiplied by their respective coefficient, they are used  used in the final formula of the Radius of Gyration
    double Arg(0.);
    double Brg(0.);

    for (int i = 0; i < Np; i++)
    {
        //$ Calculation of Rg
        Arg = Arg + volumes[i]*POW2(distances[i][Np]); // distance to the gravity center
        Brg = Brg + volumes[i]*POW2(*myspheres[i].r);
    }

    *rg = sqrt(fabs((Arg+3.0/5.0*Brg)/(*volAgregat)));
    *volAgregat=fabs(*volAgregat);
}
//###############################################################################################################################


array<int, 3> Aggregate::GetVerletIndex() noexcept
{
    array<int, 3>  index({{0,0,0}});
    double step = physicalmodel->GridDiv/physicalmodel->L;

    index[0]=int(floor((*x)*step));
    index[1]=int(floor((*y)*step));
    index[2]=int(floor((*z)*step));

    return index;
}


Sphere::Sphere(const Aggregate& Agg) : Sphere()
{
    if(physicalmodel->toBeDestroyed)
        delete physicalmodel;
    physicalmodel = Agg.physicalmodel;
    InitVal(*Agg.x,*Agg.y,*Agg.z,*Agg.rmax);
}



void Aggregate::Merge(Aggregate& other)
{

    //$ Update of the labels of the spheres that were in the deleted aggregate
    //$ And their new relative position
    double dx = periodicDistance(*other.myspheres[0].x-*myspheres[0].x,physicalmodel->L);
    double dy = periodicDistance(*other.myspheres[0].y-*myspheres[0].y,physicalmodel->L);
    double dz = periodicDistance(*other.myspheres[0].z-*myspheres[0].z,physicalmodel->L);

    // For all the spheres that were in the deleted aggregate
    for (Sphere* othersphere : other.myspheres)
    {
        // change the Label to the new owner
        othersphere->SetLabel(indexInStorage);

        // change the relative position to the new aggregate
        othersphere->RelativeTranslate(dx,dy,dz);

        // Move them accordingly (periodicity)
        othersphere->SetPosition(*myspheres[0].x+*othersphere->rx,
                                 *myspheres[0].y+*othersphere->ry,
                                 *myspheres[0].z+*othersphere->rz);
    }

    // Merge the spheresLists
    myspheres.merge(other.myspheres);
    Np = myspheres.size();

    UpdateDistances();
    Update();
}

void Aggregate::DecreaseLabel(void) noexcept
{
    if (InVerlet)
        verlet->Remove(indexInStorage,GetVerletIndex());

    indexInStorage--;
    myspheres.DecreaseLabel();

    if (InVerlet)
        verlet->Add(indexInStorage,GetVerletIndex());
}

void Aggregate::UpdateDistances(void) noexcept
{
    distances.resize(Np);
    for (int i = 0; i < Np; i++)
    {
        // The last index is the distance to the mass center
        distances[i].resize(Np+1);
    }

    for (int i = 0; i < Np; i++)
    {
        // a sphere is that close to itself
        distances[i][i]=0;

        for (int j = i+1; j < Np; j++)
        {
            // Compute the distance between sphere i and j without taking periodicity into account
            distances[i][j] = myspheres[i].RelativeDistance(myspheres[j]);

            // distances are symetric !
            distances[j][i] = distances[i][j];
        }
    }
}

void Aggregate::check(void)
{
    int k =0;
    for (const Sphere* mySphere : myspheres)
    {

        const array<double, 3> pos = mySphere->Position();
        cout << k << "\t";
        for (int j = 0; j < 3; j++)
            cout << pos[j] << "\t";
        cout << mySphere->Radius() << endl;
        k++;
    }
    for (int i=0; i<Np;i++)
    {
        for (int j=i+1; j<Np;j++)
            cout << i << "\t"
                 << j << "\t"
                 << myspheres[i].Distance(myspheres[j]) << "\t"
                 << myspheres[i].Contact(myspheres[j]) << endl;
    }
}
