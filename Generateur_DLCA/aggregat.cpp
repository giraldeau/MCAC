#include "aggregat.h"
#include <math.h>
#include <cmath>
#include <iomanip>


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW2(a) ((a)*(a))
#define POW3(a) ((a)*(a)*(a))



const double PI = atan(1.0)*4;

Aggregate::Aggregate(void):
    storage_elem<15,ListAggregat>(),
    physicalmodel(new PhysicalModel),
    myspheres(),
    verlet(nullptr),
    IndexVerlet({{0,0,0}}),
    distances(),
    volumes(),
    surfaces(),
    nctmp(0.),
    nptmp(1.),
    rg(nullptr),
    dm(nullptr),
    lpm(nullptr),
    time_step(nullptr),
    rmax(nullptr),
    volAgregat(nullptr),
    surfAgregat(nullptr),
    Tv(nullptr),
    volAgregat_without_cov(nullptr),
    cov(nullptr),
    ratio_surf_vol(nullptr),
    free_surface(nullptr),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    Nc(0),
    Np(0),
    InVerlet(false)
{
    Init();
}

Aggregate::Aggregate(PhysicalModel& _physicalmodel) :
    storage_elem<15,ListAggregat>(),
    physicalmodel(&_physicalmodel),
    myspheres(),
    verlet(nullptr),
    IndexVerlet({{0,0,0}}),
    distances(),
    volumes(),
    surfaces(),
    nctmp(0.),
    nptmp(1.),
    rg(nullptr),
    dm(nullptr),
    lpm(nullptr),
    time_step(nullptr),
    rmax(nullptr),
    volAgregat(nullptr),
    surfAgregat(nullptr),
    Tv(nullptr),
    volAgregat_without_cov(nullptr),
    cov(nullptr),
    ratio_surf_vol(nullptr),
    free_surface(nullptr),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    Nc(0),
    Np(0),
    InVerlet(false)
{
    Init();
}


Aggregate::Aggregate(ListAggregat& _storage, const int _N):
    storage_elem<15,ListAggregat>(_storage, _N),
    physicalmodel(_storage.physicalmodel),
    myspheres(),
    verlet(nullptr),
    IndexVerlet({{0,0,0}}),
    distances(),
    volumes(),
    surfaces(),
    nctmp(0.),
    nptmp(1.),
    rg(nullptr),
    dm(nullptr),
    lpm(nullptr),
    time_step(nullptr),
    rmax(nullptr),
    volAgregat(nullptr),
    surfAgregat(nullptr),
    Tv(nullptr),
    volAgregat_without_cov(nullptr),
    cov(nullptr),
    ratio_surf_vol(nullptr),
    free_surface(nullptr),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    Nc(0),
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
    *Tv=0.;                     //Taux de recouvrement volumique
    *volAgregat_without_cov=0.; //Volume of the aggregate without considering the spheres covering each other
    *cov=0.;                    //Covering Parameter
    *ratio_surf_vol=0.;         //Ratio surface / volume
    *free_surface=0.;           //Free surface of the aggregate (without covering)

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
    Tv=&(*Storage)[7][indexInStorage];                     //Taux de recouvrement volumique
    volAgregat_without_cov=&(*Storage)[8][indexInStorage]; //Volume of the aggregate without considering the spheres covering each other
    cov=&(*Storage)[9][indexInStorage];                    //Covering Parameter
    ratio_surf_vol=&(*Storage)[10][indexInStorage];         //Ratio surface / volume
    free_surface=&(*Storage)[11][indexInStorage];           //Free surface of the aggregate (without covering)
    x=&(*Storage)[12][indexInStorage];
    y=&(*Storage)[13][indexInStorage];
    z=&(*Storage)[14][indexInStorage];
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
            double dist=mysphere->Collision(*othersphere, Vectdir);;
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
    for (Sphere* mysphere: myspheres)
    {
        mysphere->Translate(vector);
    }

    SetPosition(*x +vector[0], *y + vector[1], *z +vector[2]);
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


double& Aggregate::operator[](const int var)
{
    switch (var)
    {
    case 0:
        return *rg;
    case 1:
        nptmp=double(Np);
        return nptmp;
    case 3:
        return *dm;
    case 4:
        return *lpm;
    case 5:
        return *time_step;
    case 6:
        return *rmax;
    case 7:
        return *volAgregat;
    case 8:
        return *surfAgregat;
    case 9:
        return *volAgregat_without_cov;
    case 11:
        return *free_surface;
    }

    *cov = 0.0;// Coeficient of mean covering
    Nc = 0; // Number of contacts
    double terme = 0.0; // Volume and surface of Agg Id
    *Tv = 0.0;

    //$ For the Spheres i in Agg Id
    for (int i = 0; i < Np; i++)
    {
        for (int j = i+1; j < Np; j++) //for the j spheres composing Aggregate n°id
        {
            double dist = distances[i][j];
            double rpmoy = (*myspheres[i].r + *myspheres[j].r)/2.0; //Mean Radius between i and j monomeres
            double dbordabord = ((dist-2.0*rpmoy)/(2.0*rpmoy))*1E6; //distance between the two particles
            //$ Check if i is covering j
            //$ [dbordabord <= 1]
            if (dbordabord <= 1.)
            {
                //$ Calculation of the Number of contacts
                *cov = *cov - dbordabord; //Coefficient of total Covering of Agg id
                Nc += 1; //Nombre de coordination (nombre de points de contact entre deux sphérules)
            }
        }
        terme = terme + volumes[i]/(*myspheres[i].volume);

    }
    *Tv = 1 - terme /Np;

    Nc = Nc/2;//and determine the coefficient of mean covering of Agg Id
    //$ Check if there is more than one monomere in the aggregate
    //$ [nmonoi == 1]
    if (Np == 1 || Nc == 0)
    {
        //$ Cov = 0
        *cov = 0;
    }
    else
    {
        //$ Determination of the coefficient of mean covering, using the one determined in the precedent loop
        *cov = *cov/(double(Nc))/2.0;
    }

    switch (var)
    {
    case 2:
        nctmp=double(Nc);
        return nctmp;
    case 10:
        return *cov;
    default:
        cout << "Unknown element " << var <<endl;
    }
    exit(5);
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
    volumes.assign(Np, 0.);
    surfaces.assign(Np, 0.);

    //$ For the Spheres i in Agg Id
    for (int i = 0; i < Np; i++)
    {
        //$ Calculation of the volume and surface of monomere i of Agg id
        volumes[i] += *myspheres[i].volume; //Calculation of the volume of i
        surfaces[i] += *myspheres[i].surface;    //Calculation of the surface of i

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

        double dx = periodicDistance(*myspheres[i].x-*x,physicalmodel->L);
        double dy = periodicDistance(*myspheres[i].y-*y,physicalmodel->L);
        double dz = periodicDistance(*myspheres[i].z-*z,physicalmodel->L);

        //$ Calculation of the position of the center of mass
        newpos[0] += dx * volumes[i];
        newpos[1] += dy * volumes[i];
        newpos[2] += dz * volumes[i];
    }
    //$ Filling of PosiGravite

    for (int k = 0; k < 3; k++)
        newpos[k] = newpos[k]/(*volAgregat);

    //cout << Position()<< " "<< newpos <<endl;
    SetPosition(newpos[0]+*x,newpos[1]+*y,newpos[2]+*z);

    const int loopsize(Np);
    for (int i = 0; i < loopsize; i++)
        distances[i][0]=myspheres[i].Distance(*x,*y,*z);

}

void Aggregate::CalcRadius(void)
{

    *rmax = 0.0; // Maximum radius of the aggregate, this corresponds to the distance between the center of mass of the aggregate and the edge of the furthest ball from said center.
                // It is used to assimilate the aggregate to a sphere when checking for intersections

    for (int i = 0; i < Np; i++)
        *rmax=MAX(*rmax,*myspheres[i].r+distances[i][0]);
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
        Arg = Arg + volumes[i]*POW2(distances[i][0]); // distance to the gravity center
        Brg = Brg + volumes[i]*POW2(*myspheres[i].r);
    }

    *rg = sqrt(fabs((Arg+3.0/5.0*Brg)/(*volAgregat)));
    *volAgregat=fabs(*volAgregat);
}
//###############################################################################################################################

Sphere Aggregate::GetInclusiveSphere(void) const
{
    Sphere InclusiveSphere1(*physicalmodel, *x, *y, *z, *rmax);
    return InclusiveSphere1;
}


array<int, 3> Aggregate::GetVerletIndex() noexcept
{
    array<int, 3>  index({{0,0,0}});
    double step = physicalmodel->GridDiv/physicalmodel->L;

    index[0]=int(floor((*x)*step));
    index[1]=int(floor((*y)*step));
    index[2]=int(floor((*z)*step));

    return index;
}


void Aggregate::AfficheVerlet() const
{
    _List_iterator<int> i;

    if (physicalmodel->use_verlet && InVerlet)
    {
        list<int> cell = (*verlet)[IndexVerlet[0]][IndexVerlet[1]][IndexVerlet[2]];
        cout<<"Coordinates: "<< IndexVerlet[0] << " " << IndexVerlet[1] << " " << IndexVerlet[2] << endl;
        cout<<"Taille: "<<cell.size()<< endl;
        cout << " friends :"<< endl;
        for(const int member : cell)
        {
            cout << " " << member <<endl;
        }
    }
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
    for (Sphere* othersphere : other.myspheres)
    {
        othersphere->SetLabel(indexInStorage);
    }
    myspheres.merge(other.myspheres);
    Np = myspheres.size();
    UpdateDistances();
}

void Aggregate::DecreaseLabel(void) noexcept
{
    indexInStorage--;
    myspheres.DecreaseLabel();
}

void Aggregate::UpdateDistances(void) noexcept
{
    distances.resize(Np);

    //$ For the Spheres i in Agg Id
    for (int i = 0; i < Np; i++)
        distances[i].resize(Np)
                ;
    //$ For the Spheres i in Agg Id
    for (int i = 0; i < Np; i++)
    {
        distances[i][i]=0;
        for (int j = i+1; j < Np; j++) //for the j spheres composing Aggregate n°id
        {
            distances[i][j]=myspheres[i].Distance(myspheres[j]);
            distances[j][i] = distances[i][j];
        }
    }
}
