#include "aggregat.h"
#include <math.h>
#include <cmath>
#include <iomanip>


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

const double PI = atan(1.0)*4;

Aggregate::Aggregate(void):
    physicalmodel(nullptr),
    InclusiveSphere(new Sphere()),
    myspheres(),
    parents(),
    son(nullptr),
    verlet(nullptr),
    IndexVerlet({0,0,0}),
    Storage(new array< vector<double>, 16>()),
    external_storage(nullptr),
    creation_date(0.),
    nctmp(0.),
    nptmp(0.),
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
    Label(0),
    Nc(0),
    Np(0),
    InVerlet(false)
{
    for (int j=0;j<=15;j++)
        (*Storage)[j].assign(1, 0.);
    Init();
}
Aggregate::Aggregate(PhysicalModel& _physicalmodel) : Aggregate(){ physicalmodel = &_physicalmodel; }

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

    IndexVerlet = {0,0,0,0};

    Position(0.,0.,0.);
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

void Aggregate::Init(PhysicalModel& _physicalmodel,Verlet& _verlet,const array<double, 4> position ,const int _label, ListSphere& spheres,const double Dp)
{
    physicalmodel = &_physicalmodel;
    verlet = &_verlet;
    Label = _label;

    if (physicalmodel->use_verlet)
    {
        verlet->GetCell(IndexVerlet[1],IndexVerlet[2],IndexVerlet[3])->push_front(Label);
        InVerlet=true;
    }
    Position(position);
    int listlabel[2] = {1,Label};
    spheres[Label].SetLabel(Label);
    spheres[Label].Init(position, Dp/2);
    myspheres = ListSphere(spheres,listlabel);


    *rmax = Dp/2;                 //Rayon de la sphère d'enveloppe de l'agrégat réunifié
    *dm = Dp;                   //Diamètre de mobilité
    *volAgregat = PI*pow(Dp,3)/6;
    *surfAgregat = PI*pow(Dp,2);

    Nc = 0;                  //Nombre de coordination Nc
    *cov =0;
    *rg = sqrt(3.0/5.0)*Dp/2;

    *volAgregat_without_cov = PI*pow(Dp, 3)/6;      //Volume de l'agrégat réunifié sans recouvrement (Avant c'était 0; : Taux de recouvrement volumique)
    *free_surface = PI*pow(Dp, 2);       //Surface de l'agrégat réunifié sans recouvrement (Avant c'était surface/(masse/Rho); : Surface/volume de l'agrégat réunifié)

    if (physicalmodel->ActiveModulephysique==1)
    {
        double Diff = physicalmodel->diffusivity(Dp);
        double masse = physicalmodel->Rho*PI*pow(Dp,3)/6;
        double Vit =  physicalmodel->velocity(masse);
        double tmp = 8*Diff/PI/Vit;
        *lpm = tmp;
        *time_step = tmp/Vit;              //Durée du déplacement
    }
    else
    {
        *lpm = physicalmodel->Dpm*1E-9;    //Libre parcours moyen
        *time_step = 1e-6;              //Durée du déplacement
    }

}


//############################################# Calcul de la distance inter-agrégats ############################################
double Aggregate::Distance_Aggregate(Aggregate& other,array<double,4> vectorOther,array<double,4> Vectdir)
{
    double dist, ret;
    ret = 1.0;

    //Liste des sphérules constituant l'agrégat n°Agg

    int nmonoi = myspheres.size();
    int np = other.myspheres.size();

    //$ Loop on all the spheres of the other aggregate
    for (int j = 1; j <= np; j++)
    {
        //$ spheredecale is used to replace the sphere into the corresponding box
        Sphere spheredecale(other.myspheres[j]);
        spheredecale.Translate(vectorOther);

        //$ For every sphere in the aggregate :
        for (int knum = 1; knum <= nmonoi; knum++)
        {
            //$ Check if j and k are contact and if not, what is the distance between them
            dist=myspheres[knum].Collision(spheredecale, Vectdir, *lpm);
            if (dist < ret)
            {
                ret = dist;
            }

        }
    }
    return ret;
}
//###############################################################################################################################


//######### Mise à jour des paramètres physiques d'un agrégat (rayon de giration, masse, nombre de sphérules primaires) #########
void Aggregate::Update()
{// This function will update the parameter of Agg

    //$ Determination of the Radius of gyration of Agg using RayonGiration()
    RayonGiration();

    //$ Determination of the mean radius of degrees 1,2,3 of Agg
    int np = myspheres.size();
    double rpmoy, rpmoy2, rpmoy3;
    rpmoy = rpmoy2 = rpmoy3 = 0.0;
    for (int i = 1; i <= np; i++)
    {
        rpmoy = rpmoy + myspheres[i].Radius(); //Sum of the radius of each sphere in Agg
        rpmoy2 = rpmoy2 + pow(myspheres[i].Radius(), 2);
        rpmoy3 = rpmoy3 + pow(myspheres[i].Radius(), 3);
    }
    rpmoy = rpmoy/(double(np));
    rpmoy2 = rpmoy2/(double(np));
    rpmoy3 = rpmoy3/(double(np));


    //$ Determination of Dm using ConvertRg2Dm
    *dm = physicalmodel->ConvertRg2DmFromStart(np,*rg,*dm/2);

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

    *volAgregat_without_cov = np*4.0*PI*rpmoy3/3.0; //Volume of the aggregate without considering the spheres covering each other (Avant c'était Tv : Taux de recouvrement volumique)
    *surfAgregat = np*4.0*PI*rpmoy2;  //Free surface of the aggregate (without covering)(Avant c'était surfAgregat/volAgregat : Estimation du rapport surface/volume de l'agrégat)
}



void Aggregate::UpdatesSpheres(ListSphere& spheres,int* index)
{
    myspheres = ListSphere(spheres,index);
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

    if (InVerlet && physicalmodel->use_verlet)
    {
        //$ Update Verlet
        array<int, 4> newindexVerlet = VerletIndex();

        if (newindexVerlet != IndexVerlet)
        {

            verlet->GetCell(IndexVerlet[1],IndexVerlet[2],IndexVerlet[3])->remove(Label);
            IndexVerlet = newindexVerlet;
            verlet->GetCell(IndexVerlet[1],IndexVerlet[2],IndexVerlet[3])->push_front(Label);
        }
    }

    //ReplacePosi();
}

void Aggregate::Position(const array<double, 4> position)
{
    Position(position[1], position[2], position[3]);
}

void Aggregate::Translate(const array<double, 4> vector)
{

    int nmonoi = myspheres.size();
    for (int i = 1; i <= nmonoi; i++)
    {
        myspheres[i].Translate(vector);
    }

    Position(*x +vector[1], *y + vector[2], *z +vector[3]);
}

void Aggregate::Translate(const double* vector)
{

    int nmonoi = myspheres.size();
    for (int i = 1; i <= nmonoi; i++)
    {
        myspheres[i].Translate(vector);
    }

    Position(*x +vector[1], *y + vector[2], *z +vector[3]);
}

//############################################# Conditions aux limites périodiques ##############################################

void Aggregate::ReplacePosi()
{
    // This function will relocate an aggregate when it gets out of the box limiting the space

    if (physicalmodel == nullptr)
        return;

    array<double, 4> trans;
    bool move=false;

    const array<double, 4> pos = Position();

    //$ for every dimension
    for (int i = 1; i <= 3; i++)
    {
        //$ Check if it is getting out
        if (pos[i] > physicalmodel->L)
        {
            trans[i] = - physicalmodel->L;
            move = true;
        }
        else if (pos[i] < 0)
        {
            trans[i] = physicalmodel->L;
            move = true;
        }
        else
        {
            trans[i] = 0;
        }
    }

    //$ If it is getting out
    if (move)
    {
        //$ Update the position of aggregate
        Translate(trans);
    }
}

Aggregate::~Aggregate(void)
{
    if (InVerlet && physicalmodel->use_verlet)
    {
        verlet->GetCell(IndexVerlet[1],IndexVerlet[2],IndexVerlet[3])->remove(Label);
        InVerlet=false;
    }

    if (external_storage==nullptr)
    {
        delete Storage;
    }

    external_storage=nullptr;

    delete InclusiveSphere;
    parents[0] = nullptr;
    parents[1] = nullptr;
    son = nullptr;
    physicalmodel = nullptr;
    Label = 0;
    creation_date = 0.;
    Nc = 0.;

    rg=nullptr;
    dm=nullptr;
    lpm=nullptr;
    time_step=nullptr;
    rmax=nullptr;
    volAgregat=nullptr;
    surfAgregat=nullptr;
    Tv=nullptr;
    volAgregat_without_cov=nullptr;
    cov=nullptr;
    ratio_surf_vol=nullptr;
    free_surface=nullptr;
    x=nullptr;
    y=nullptr;
    z=nullptr;



}


double& Aggregate::operator[](const int i)
{
    switch (i)
    {
    case 0:
        return *rg;
    case 1:
        nptmp=double(myspheres.size());
        return nptmp;
    case 2:
        nctmp=double(Nc);
        return nctmp;
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
    case 10:
        return *cov;
    case 11:
        return *free_surface;
    default:
        cout << "Unknown element " << i <<endl;
    }
    exit(5);
}


//############# Calculation of the volume, surface, center of mass and Giration radius of gyration of an aggregate ##############

void Aggregate::RayonGiration(void)
{// This function determines the Gyration Radius of the Aggregate Id.
    double dist, rpmoy, dbordabord, li, r, Arg, Brg, terme;
    int nmonoi;
    vector<double> tabVol;
    vector<double> tabSurf;

    *cov = 0.0;// Coeficient of mean covering
    Nc = 0; // Number of contacts
    *volAgregat = *surfAgregat = terme = 0.0; // Volume and surface of Agg Id
    *rmax = 0.0; // Maximum radius of the aggregate, this corresponds to the distance between the center of mass of the aggregate and the edge of the furthest ball from said center.
                // It is used to assimilate the aggregate to a sphere when checking for intersections
    *Tv = 0.0;

    nmonoi = myspheres.size();

    tabVol.assign(nmonoi+1, 0.);
    tabSurf.assign(nmonoi+1, 0.);

    //$ Initialisation of the arrays of volume, surface of each sphere, and the center of mass
    array<double, 4> newpos={0.,0.,0.,0.};

    //$ For the Spheres i in Agg Id
    for (int i = 1; i <= nmonoi; i++)
    {
        //$ Calculation of the volume and surface of monomere i of Agg id
        tabVol[i] += myspheres[i].Volume(); //Calculation of the volume of i
        tabSurf[i] += myspheres[i].Surface();    //Calculation of the surface of i

        for (int j = i+1; j <= nmonoi; j++) //for the j spheres composing Aggregate n°id
        {

            double voli, volj, surfi, surfj;
            voli = volj = surfi = surfj = 0.;

            //$ Calculation of the intersection between the spheres i and j
            dist = myspheres[i].Intersection(myspheres[j],voli,volj,surfi,surfj);

            rpmoy = (myspheres[i].Radius()+myspheres[j].Radius())/2.0; //Mean Radius between i and j monomeres
            dbordabord = ((dist-2.0*rpmoy)/(2.0*rpmoy))*1E6; //distance between the two particles
            //$ Check if i is covering j
            //$ [dbordabord <= 1]
            if (dbordabord <= 1.)
            {
                //$ Calculation of the Number of contacts
                *cov = *cov - dbordabord; //Coefficient of total Covering of Agg id
                Nc += 1; //Nombre de coordination (nombre de points de contact entre deux sphérules)
            }

            //$ The volume and surface covered by j is substracted from those of i
            tabVol[i] = tabVol[i] - voli;    //Calcul du volume de la sphérule i moins le volume de
                                             //la calotte due à la sphérule j
            tabSurf[i] = tabSurf[i] - surfi; //Calcul de la surface de la sphérule i moins la surface de
                                             //la calotte due à la sphérule j

            //$ The volume and surface covered by i is substracted from those of j
            tabVol[j] = tabVol[j] - volj;    //Calcul du volume de la sphérule j moins le volume de
                                             //la calotte due à la sphérule i
            tabSurf[j] = tabSurf[j] - surfj; //Calcul de la surface de la sphérule j moins la surface de
                                             //la calotte due à la sphérule i
        }
        //$ Calculation of the total volume and surface of the aggregate
        *volAgregat = *volAgregat + tabVol[i];    //Total Volume of Agg id
        *surfAgregat = *surfAgregat + tabSurf[i]; //Total Surface of Agg id

        terme = terme + tabVol[i]/myspheres[i].Volume();

        //$ Calculation of the position of the center of mass
        const array<double, 4> pos = myspheres[i].Position();
        for (int k = 1; k <= 3; k++)
        {
            newpos[k] += pos[k]*tabVol[i];
        }
    }
    *Tv = 1 - terme /nmonoi;

    Nc = Nc/2;//and determine the coefficient of mean covering of Agg Id
    //$ Check if there is more than one monomere in the aggregate
    //$ [nmonoi == 1]
    if (nmonoi == 1 || Nc == 0)
    {
        //$ Cov = 0
        *cov = 0;
    }
    else
    {
        //$ Determination of the coefficient of mean covering, using the one determined in the precedent loop
        *cov = *cov/(double(Nc))/2.0;
    }
    //$ Filling of PosiGravite

    for (int k = 1; k <= 3; k++)
    {
        newpos[k] = newpos[k]/(*volAgregat);
    } //Centre of mass of Agg Id

    //cout << Position()<< " "<< newpos <<endl;
    Position(newpos);
    ReplacePosi();

    //$ Determination of the maximal radius of Agg Id and the Radius of gyration
    Arg = Brg = 0.0; // These correspond to the sum of the volumes of each spheres multiplied by their respective coefficient, they are used  used in the final formula of the Radius of Gyration

    for (int i = 1; i <= nmonoi; i++)
    {
        //$ Determination of the distance between each monomere and the center of mass of Agg Id
        li = myspheres[i].Distance(Position()); //Distance entre les centres de masse de la sphérule i et de l'agrégat n°id

        r = li + myspheres[i].Radius();

        //$ Calculation of rmax
        *rmax=MAX(*rmax,r);

        //$ Calculation of Rg
        Arg = Arg + tabVol[i]*pow(li, 2);
        Brg = Brg + tabVol[i]*pow(myspheres[i].Radius(), 2);
    }

    *rg = sqrt(fabs((Arg+3.0/5.0*Brg)/(*volAgregat)));
    *volAgregat=fabs(*volAgregat);
}
//###############################################################################################################################




array<int, 4> Aggregate::VerletIndex()
{
    array<int, 4>  index= {0,0,0,0};
    double step = physicalmodel->GridDiv/physicalmodel->L;
    int origin = physicalmodel->GridDiv+1;

    index[1]=int(floor((*x)*step)+origin);
    index[2]=int(floor((*y)*step)+origin);
    index[3]=int(floor((*z)*step)+origin);
    return index;
}


void Aggregate::AfficheVerlet()
{   _List_iterator<int> i;

    if (physicalmodel->use_verlet && InVerlet)
    {
        list<int>* cell = (*verlet).GetCell(IndexVerlet[1],IndexVerlet[2],IndexVerlet[3]);
        cout<<"Coordinates: "<< IndexVerlet[1] << " " << IndexVerlet[2] << " " << IndexVerlet[3] << endl;
        cout<<"Taille: "<<cell->size()<< endl;
        cout << " friends :"<< endl;
        for(i = cell->begin();
            i!= cell->end();
            i++)
        {
            cout << " " << *i<<endl;
        }
    }
}



void Verlet::Remove(const int id,const array<int, 4> Index)
{
    verletlist[Index[1]][Index[2]][Index[3]]->remove(id);
}

void Verlet::Init(const int GridDiv)
{
    verletlist=new list<int>***[3*GridDiv+1];
    for(int i=0;i<=3*GridDiv;i++)
    {
        verletlist[i]=new list<int>**[3*GridDiv+1];

        for(int j=0;j<=3*GridDiv;j++)
        {
            verletlist[i][j]=new list<int>*[3*GridDiv+1];

            for(int k=0;k<=3*GridDiv;k++)
            {
                verletlist[i][j][k]= new list<int>;
            }
        }

    }
}

__attribute((pure)) list<int>* Verlet::GetCell(const int i,const int j,const int k)const
{
    return verletlist[i][j][k];
}

