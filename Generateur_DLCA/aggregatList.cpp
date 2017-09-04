#include "aggregat.h"
#include <math.h>
#include <cmath>
#include <iomanip>
#include <algorithm>


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW2(a) ((a)*(a))
#define POW3(a) ((a)*(a)*(a))

ListAggregat::ListAggregat(void):
    storage_list<16,Aggregate>(),
    physicalmodel(new PhysicalModel),
    maxradius(0.),
    indexSortedTimeSteps(),
    CumulativeTimeSteps(),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    spheres(),
    verlet()
{
}

  __attribute__((pure)) double ListAggregat::GetMaxTimeStep() const
{

    double m = *list[0]->time_step;
    for (const Aggregate* Agg : list)
        m =max(*Agg->time_step, m);

    return m;
}

void ListAggregat::Init(PhysicalModel& _physicalmodel,const int _N)
{
    if(physicalmodel->toBeDestroyed)
        delete physicalmodel;

    physicalmodel=&_physicalmodel;
    spheres.Init(_physicalmodel, _N);
    verlet.Init(_physicalmodel.GridDiv,_physicalmodel.L);

    storage_list<16,Aggregate>::Init(_N,*this);
    setpointers();
}


//########################################## Determination of the contacts between agrgates ##########################################
vector<int> ListAggregat::PotentialContacts(int AggMe,array<double,4> Vectdir, vector<int> SearchSpace) const
{
    vector<int> listOfPotentialContacts;

    Sphere SphereMe(*list[AggMe-1]);
    Sphere SphereOther(*list[AggMe-1]);

    //$ [For all other agregate]
    for (const int& AggOther : SearchSpace)
    {
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


vector<int> ListAggregat::GetSearchSpace(const int source, const array<double,4> Vectdir) const
{
    if (!physicalmodel->use_verlet)
    {
        vector < int > SearchSpace(indexInStorage);
        for (int& idx : SearchSpace)
            idx++;
        return SearchSpace;
    }
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
int ListAggregat::DistanceToNextContact(const int source, const array<double,4> Vectdir, double &distmin) const
{
    // Use Verlet to reduce search
    vector<int> SearchSpace(GetSearchSpace(source,Vectdir));

    // Assimilate Aggregate as sphere to drasticly speed-up search
    vector<int> PotentialContact(PotentialContacts(source,Vectdir,SearchSpace));

    int aggcontact(0);
    distmin = *list[source-1]->lpm;

    //$ loop on the agregates potentially in contact
    for (const int & agg : PotentialContact) //For every aggregate that could be in contact
    {

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
    vector<double>::iterator newdeb((*Storage)[1].begin());
    vector<double>::iterator newfin((*Storage)[1].end());
    if ((newdeb == ptr_deb) && (newfin == ptr_fin))
        return;
    for (Aggregate* Agg : list)
    {
        Agg->setpointers();
    }
    ptr_deb = newdeb;
    ptr_fin = newfin;
}


ListAggregat::~ListAggregat(void) noexcept
{
    if(physicalmodel->toBeDestroyed)
        delete physicalmodel;

    //#pragma omp for simd
    for (Aggregate* Agg : list)
    {
        Agg->InVerlet=false;
    }
}



int ListAggregat::Merge(const int first, const int second)
{
    const int keeped(MIN(first,second)-1);
    const int removed(MAX(first,second)-1);

    list[keeped]->Merge(*list[removed]);

    if (list[removed]->InVerlet)
        verlet.Remove(removed+1,list[removed]->GetVerletIndex());
    for (int i=removed+1;i<_size;i++)
    {
        if (list[i]->InVerlet)
        {
            verlet.Remove(i+1,list[i]->IndexVerlet);
            verlet.Add(i,list[i]->IndexVerlet);
        }
        list[i]->DecreaseLabel();
    }

    remove(*list[removed]);

    setpointers();

    return keeped+1;
}


template <typename T>
vector<int> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<int> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](int i1, int i2) {return v[i1] < v[i2];});

  return idx;
}

void ListAggregat::SortTimeSteps(double factor)
{
    vector<double> TpT(_size);

    #pragma omp for simd
    for (int i=0; i < _size; i++)
        TpT[i] = factor/(*list[i]->time_step);

    indexSortedTimeSteps = sort_indexes(TpT);

    CumulativeTimeSteps.resize(_size);

    //$ Accumulate the timesteps
    CumulativeTimeSteps[0] = TpT[indexSortedTimeSteps[0]];
    for (int i=1; i < _size; i++)
    {
        CumulativeTimeSteps[i] = CumulativeTimeSteps[i-1]+TpT[indexSortedTimeSteps[i]];
    }
}

int ListAggregat::RandomPick(double &deltatemps, const double random)
{
    //$ Pick a random sphere
    double valAlea=random*CumulativeTimeSteps[_size-1];
    long int n = lower_bound(CumulativeTimeSteps.begin(), CumulativeTimeSteps.end(), valAlea) - CumulativeTimeSteps.begin();

    deltatemps = CumulativeTimeSteps[_size-1];

    return indexSortedTimeSteps[n]+1;
}

