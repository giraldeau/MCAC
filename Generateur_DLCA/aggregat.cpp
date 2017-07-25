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
    InVerlet = false;
    verlet = NULL;

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
    InVerlet = false;
    verlet = NULL;

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

void Aggregate::Init(PhysicalModel& _physicalmodel,Verlet& _verlet,const array<double, 4> position ,const int _label, ListSphere& spheres,const double _r)
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
    spheres[Label].Init(position, _r);
    myspheres = ListSphere(spheres,listlabel);

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

Aggregate::~Aggregate(void)
{
    if (InVerlet && physicalmodel->use_verlet)
    {
        verlet->GetCell(IndexVerlet[1],IndexVerlet[2],IndexVerlet[3])->remove(Label);
        InVerlet=false;
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

list<int>* Verlet::GetCell(const int i,const int j,const int k)const
{
    return verletlist[i][j][k];
}

