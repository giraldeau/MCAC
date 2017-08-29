#include "aggregat.h"
#include <math.h>
#include <cmath>
#include <iomanip>


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW2(a) ((a)*(a))
#define POW3(a) ((a)*(a)*(a))

ListAggregat::ListAggregat(void):
    physicalmodel(nullptr),
    maxradius(0.),
    maxtime_step(0.),
    spheres(),
    verlet()
{}

double ListAggregat::GetMaxTimeStep()
{
    double m = *list[1]->time_step;
    for (int i = 0;i < size; i++)
        m =max(*list[i]->time_step, m);

    return m;
}

void ListAggregat::Init(PhysicalModel& _physicalmodel,const int _N)
{
    physicalmodel=&_physicalmodel;
    spheres.Init(_physicalmodel, _N);
    verlet.Init(_physicalmodel.GridDiv,_physicalmodel.L);

    storage_list<16,Aggregate>::Init(_N,*this);
}


//########################################## Determination of the contacts between agrgates ##########################################
vector<int> ListAggregat::PotentialContacts(int AggMe,array<double,4> Vectdir, vector<int> SearchSpace)
{
    vector<int> listOfPotentialContacts;

    Sphere SphereMe(*list[AggMe-1]);
    Sphere SphereOther(*list[AggMe-1]);

    //$ [For all other agregate]
    for (unsigned i = 0;i < SearchSpace.size(); i++)
    {
        int AggOther = SearchSpace[i];

        if (AggOther != AggMe)
        {
            SphereOther.InitVal(*list[AggOther-1]->x,
                                *list[AggOther-1]->y,
                                *list[AggOther-1]->z,
                                *list[AggOther-1]->rmax);

            double distForContact;

            if (SphereMe.Contact(SphereOther))
                distForContact = 0;
            else
                distForContact = SphereMe.Collision(SphereOther,Vectdir);

            if(0. <= distForContact && distForContact <= *list[AggMe-1]->lpm )
                listOfPotentialContacts.push_back(AggOther);
        }
    }
    return listOfPotentialContacts;
}
//###############################################################################################################################


vector<int> ListAggregat::GetSearchSpace(const int source, const array<double,4> Vectdir)
{
    if (!physicalmodel->use_verlet)
        return indexInStorage;
    else
    {
        double lpm ( *list[source-1]->lpm );
        double mindist ( *list[source-1]->rmax + maxradius );
        array<double, 4> sourceposition = list[source-1]->GetPosition();

        array<double, 4> Vector;
        Vector[0] = 0;
        Vector[1] = lpm * Vectdir[1];
        Vector[2] = lpm * Vectdir[2];
        Vector[3] = lpm * Vectdir[3];

        return verlet.GetSearchSpace(sourceposition , mindist, Vector);
    }
}



//########################################## Determination of the contacts between agrgates ##########################################
int ListAggregat::DistanceToNextContact(const int source, const array<double,4> Vectdir, double &distmin)
{
    // Use Verlet to reduce search
    vector<int> SearchSpace(GetSearchSpace(source,Vectdir));

    // Assimilate Aggregate as sphere to drasticly speed-up search
    vector<int> PotentialContact(PotentialContacts(source,Vectdir,SearchSpace));

    const int npossible(PotentialContact.size());
    int aggcontact(0);
    distmin = *list[source-1]->lpm;

    //$ loop on the agregates potentially in contact
    for (int s = 0; s < npossible; s++) //For every aggregate that could be in contact
    {
        int agg = PotentialContact[s];

        if (list[source-1]->Contact(*list[agg-1]))
        {
            cout << "Already contact !! " << endl;
            cout << " ignoring" << endl;
        }
        else
        {
            double dist = list[source-1]->Distance(*list[agg-1],Vectdir);
            if (dist >= 0 && dist <= distmin)
            {
                distmin = dist;
                aggcontact = agg; //Prise en compte de l'image par translation d'un agrégat cible
            }
        }
    }

    return aggcontact;
}
//###############################################################################################################################



void ListAggregat::setpointers()
{
    //#pragma omp for simd
    for (int i = 0; i < size; i++)
    {
        list[i]->setpointers();
    }
}


ListAggregat::~ListAggregat(void) noexcept
{
    //#pragma omp for simd
    for (int i = 0; i < size; i++)
    {
        list[i]->InVerlet=false;
    }
}
