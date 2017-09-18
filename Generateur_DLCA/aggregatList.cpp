#include "aggregat.h"
#include <math.h>
#include <cmath>
#include <iomanip>
#include <algorithm>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW2(a) ((a)*(a))
#define POW3(a) ((a)*(a)*(a))

using namespace std;

ListAggregat::ListAggregat(void):
    storage_list<13,Aggregate>(),
    physicalmodel(new PhysicalModel),
    maxradius(0.),
    indexSortedTimeSteps(),
    CumulativeTimeSteps(),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    Writer(new ThreadedIO(*physicalmodel, 0)),
    spheres(),
    verlet()
{
}

  __attribute__((pure)) double ListAggregat::GetMaxTimeStep() const
{

    double m = *list[0]->time_step;
    for (const Aggregate* Agg : list)
        m =MAX(*Agg->time_step, m);

    return m;
}

void ListAggregat::Init(PhysicalModel& _physicalmodel,const int _N)
{
    if(physicalmodel->toBeDestroyed)
    {
        delete physicalmodel;
        delete Writer;
    }

    physicalmodel=&_physicalmodel;
    Writer = new ThreadedIO(_physicalmodel, _N),
    spheres.Init(_physicalmodel, _N);
    verlet.Init(_physicalmodel.GridDiv,_physicalmodel.L);

    storage_list<13,Aggregate>::Init(_N,*this);
    setpointers();
}


//########################################## Determination of the contacts between agrgates ##########################################
vector<int> ListAggregat::PotentialCollision(int AggMe,array<double,3> Vectdir, vector<int> SearchSpace) const
{
    vector<int> listOfPotentialCollision;

    Sphere SphereMe(*list[AggMe]);
    Sphere SphereOther(*list[AggMe]);

    //$ [For all other agregate]
    for (const int& AggOther : SearchSpace)
    {
        if (AggOther != AggMe)
        {
            SphereOther.InitVal(*list[AggOther]->x,
                                *list[AggOther]->y,
                                *list[AggOther]->z,
                                *list[AggOther]->rmax);

            double distForCollision;

            if (SphereMe.Contact(SphereOther))
                distForCollision = 0;
            else
                distForCollision = SphereMe.Collision(SphereOther,Vectdir);

            if(0. <= distForCollision && distForCollision <= *list[AggMe]->lpm )
                listOfPotentialCollision.push_back(AggOther);
        }
    }
    return listOfPotentialCollision;
}

vector<int> ListAggregat::GetSearchSpace(const int source, const array<double,3> Vectdir) const
{

    vector < int > SearchSpace;
    if (!physicalmodel->use_verlet)
    {
        // The full aggregat list index
        SearchSpace.resize(size());
        iota(SearchSpace.begin(), SearchSpace.end(), 0);

        // Except me
        SearchSpace.erase(SearchSpace.begin()+source);

        return SearchSpace;
    }
    else
    {
        // Extract from verlet

        double lpm ( *list[source]->lpm );
        double mindist ( *list[source]->rmax + maxradius );
        array<double, 3> sourceposition = list[source]->GetPosition();

        array<double, 3> Vector;
        Vector[0] = lpm * Vectdir[0];
        Vector[1] = lpm * Vectdir[1];
        Vector[2] = lpm * Vectdir[2];

        SearchSpace = verlet.GetSearchSpace(sourceposition , mindist, Vector);

        // Remove me
        for(size_t i =0;i<SearchSpace.size();i++)
            if (SearchSpace[i]==source)
            {
                SearchSpace.erase(SearchSpace.begin()+i);
                return SearchSpace;
            }
        cout << "I'm not on the verlet list ???"<<endl;
        cout << "This is an error"<<endl;
        exit(68);
        return SearchSpace;
    }
}



//########################################## Determination of the contacts between agrgates ##########################################
int ListAggregat::DistanceToNextContact(const int source, const array<double,3> Vectdir, double &distmin) const
{
    // Use Verlet to reduce search
    vector<int> SearchSpace(GetSearchSpace(source,Vectdir));

    // Assimilate Aggregate as sphere to drasticly speed-up search
    vector<int> Potential(PotentialCollision(source,Vectdir,SearchSpace));
    //vector<int> PotentialContact(SearchSpace);

    int aggcontact(-1);
    distmin = *list[source]->lpm;

    //$ loop on the agregates potentially in contact
    for (const int & agg : Potential) //For every aggregate that could be in contact
    {
        // If two aggragates are already in contact due to surface growing
        if (list[source]->Contact(*list[agg]))
        {
            distmin = 0;
            return agg;
        }
        else
        {
            double dist = list[source]->Distance(*list[agg],Vectdir);
            if (dist >= 0 && dist <= distmin)
            {
                distmin = dist;
                aggcontact = agg; //Prise en compte de l'image par translation d'un agrégat cible
            }
        }
    }

    return aggcontact;
}

void ListAggregat::setpointers()
{
    vector<double>::iterator newdeb((*Storage)[0].begin());
    vector<double>::iterator newfin((*Storage)[0].end());
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
    {
        delete physicalmodel;
        delete Writer;
    }

    //#pragma omp for simd
    for (Aggregate* Agg : list)
    {
        Agg->InVerlet=false;
    }
}



int ListAggregat::Merge(const int first, const int second)
{
    const int keeped(MIN(first,second));
    const int removed(MAX(first,second));

    list[keeped]->Merge(*list[removed]);

    if (list[removed]->InVerlet)
    {
        verlet.Remove(removed,list[removed]->GetVerletIndex());
        list[removed]->InVerlet = false;
    }

    for (int i=removed+1;i<size();i++)
    {
        list[i]->DecreaseLabel();
    }

    remove(*list[removed]);

    setpointers();

    return keeped;
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
    vector<double> TpT(size());

    #pragma omp for simd
    for (int i=0; i < size(); i++)
        TpT[i] = factor/(*list[i]->time_step);

    indexSortedTimeSteps = sort_indexes(TpT);

    CumulativeTimeSteps.resize(size());

    //$ Accumulate the timesteps
    CumulativeTimeSteps[0] = TpT[indexSortedTimeSteps[0]];
    for (int i=1; i < size(); i++)
    {
        CumulativeTimeSteps[i] = CumulativeTimeSteps[i-1]+TpT[indexSortedTimeSteps[i]];
    }
}

int ListAggregat::RandomPick(double &deltatemps, const double random)
{
    //$ Pick a random sphere
    double valAlea=random*CumulativeTimeSteps[size()-1];
    long int n = lower_bound(CumulativeTimeSteps.begin(), CumulativeTimeSteps.end(), valAlea) - CumulativeTimeSteps.begin();

    deltatemps = CumulativeTimeSteps[size()-1];

    return indexSortedTimeSteps[n];
}

