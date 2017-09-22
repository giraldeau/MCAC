#include "aggregatList.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW2(a) ((a)*(a))
#define POW3(a) ((a)*(a)*(a))

using namespace std;

namespace DLCA{


ListAggregat::ListAggregat():
    storage_list<15,Aggregate>(),
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
    {
        m =MAX(*Agg->time_step, m);
    }

    return m;
}

void ListAggregat::Init(PhysicalModel& _physicalmodel,const size_t _size)
{
    if(physicalmodel->toBeDestroyed)
    {
        delete physicalmodel;
        delete Writer;
    }

    physicalmodel=&_physicalmodel;
    Writer = new ThreadedIO(_physicalmodel, _size),
    spheres.Init(_physicalmodel, _size);
    verlet.Init(_physicalmodel.GridDiv,_physicalmodel.L);

    storage_list<15,Aggregate>::Init(_size,*this);
    setpointers();
}


//########################################## Determination of the contacts between agrgates ##########################################
vector<size_t> ListAggregat::PotentialCollision(size_t MovingAgg,array<double,3> Vectdir, vector<size_t> SearchSpace) const
{
    vector<size_t> listOfPotentialCollision;

    Sphere SphereMe(*list[MovingAgg]);
    Sphere SphereOther(*list[MovingAgg]);

    //$ [For all other agregate]
    for (const size_t& AggOther : SearchSpace)
    {
        if (AggOther != MovingAgg)
        {
            SphereOther.InitVal(*list[AggOther]->x,
                                *list[AggOther]->y,
                                *list[AggOther]->z,
                                *list[AggOther]->rmax);

            double distForCollision;

            if (SphereMe.Contact(SphereOther))
            {
                distForCollision = 0;
            }
            else
            {
                distForCollision = SphereMe.Collision(SphereOther,Vectdir);
            }

            if(0. <= distForCollision && distForCollision <= *list[MovingAgg]->lpm )
            {
                listOfPotentialCollision.push_back(AggOther);
            }
        }
    }
    return listOfPotentialCollision;
}

vector<size_t> ListAggregat::GetSearchSpace(const size_t source, const array<double,3> Vectdir) const
{

    vector < size_t > SearchSpace;
    if (physicalmodel->use_verlet)
    {   
        // Extract from verlet

        double lpm ( *list[source]->lpm );
        double mindist ( *list[source]->rmax + maxradius );
        array<double, 3> sourceposition = list[source]->GetPosition();

        array<double, 3> Vector({lpm * Vectdir[0],
                                 lpm * Vectdir[1],
                                 lpm * Vectdir[2]});

        SearchSpace = verlet.GetSearchSpace(sourceposition , mindist, Vector);

        // Remove me
        for(size_t i =0;i<SearchSpace.size();i++)
        {
            if (SearchSpace[i]==source)
            {
                SearchSpace.erase(SearchSpace.begin()+i);
                return SearchSpace;
            }
        }
        cout << "I'm not on the verlet list ???"<<endl;
        cout << "This is an error"<<endl;
        exit(68);
        return SearchSpace;
    }

    // The full aggregat list index
    SearchSpace.resize(size());
    iota(SearchSpace.begin(), SearchSpace.end(), 0);

    // Except me
    SearchSpace.erase(SearchSpace.begin()+source);

    return SearchSpace;
}



//########################################## Determination of the contacts between agrgates ##########################################
int ListAggregat::DistanceToNextContact(const size_t source, const array<double,3> Vectdir, double &distmin) const
{
    // Use Verlet to reduce search
    vector<size_t> SearchSpace(GetSearchSpace(source,Vectdir));

    // Assimilate Aggregate as sphere to drasticly speed-up search
    vector<size_t> Potential(PotentialCollision(source,Vectdir,SearchSpace));
    //vector<size_t> PotentialContact(SearchSpace);

    int aggcontact(-1);
    distmin = *list[source]->lpm;

    //$ loop on the agregates potentially in contact
    for (const size_t & agg : Potential) //For every aggregate that could be in contact
    {
        // If two aggragates are already in contact due to surface growing
        if (list[source]->Contact(*list[agg]))
        {
            distmin = 0;
            cout << "Surprise there is already a contact !" << endl;
            cout << "This is probably due to surface growing" << endl;
            return int(agg);
        }
        double dist = list[source]->Distance(*list[agg],Vectdir);
        if (dist >= 0 && dist <= distmin)
        {
            distmin = dist;
            aggcontact = int(agg); //Prise en compte de l'image par translation d'un agrégat cible
        }
    }

    return aggcontact;
}

void ListAggregat::setpointers()
{
    auto newdeb((*Storage)[0].begin());
    auto newfin((*Storage)[0].end());
    if ((newdeb == ptr_deb) && (newfin == ptr_fin))
    {
        return;
    }
    for (Aggregate* Agg : list)
    {
        Agg->setpointers();
    }
    ptr_deb = newdeb;
    ptr_fin = newfin;
}


ListAggregat::~ListAggregat() noexcept
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



size_t ListAggregat::Merge(const size_t first, const size_t second)
{
    const size_t keeped(MIN(first,second));
    const size_t removed(MAX(first,second));

    list[keeped]->Merge(*list[removed]);

    if (list[removed]->InVerlet)
    {
        verlet.Remove(removed,list[removed]->GetVerletIndex());
        list[removed]->InVerlet = false;
    }

    for (size_t i=removed+1;i<size();i++)
    {
        list[i]->DecreaseLabel();
    }

    remove(*list[removed]);

    setpointers();

    return keeped;
}


template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

void ListAggregat::SortTimeSteps(double factor)
{
    vector<double> TpT(size());

    #pragma omp for simd
    for (size_t i=0; i < size(); i++)
    {
        TpT[i] = factor/(*list[i]->time_step);
    }

    indexSortedTimeSteps = sort_indexes(TpT);

    CumulativeTimeSteps.resize(size());

    //$ Accumulate the timesteps
    CumulativeTimeSteps[0] = TpT[indexSortedTimeSteps[0]];
    for (size_t i=1; i < size(); i++)
    {
        CumulativeTimeSteps[i] = CumulativeTimeSteps[i-1]+TpT[indexSortedTimeSteps[i]];
    }
}

size_t ListAggregat::RandomPick(double &deltatemps, const double random)
{
    //$ Pick a random sphere
    double valAlea=random*CumulativeTimeSteps[size()-1];
    size_t n = lower_bound(CumulativeTimeSteps.begin(), CumulativeTimeSteps.end(), valAlea) - CumulativeTimeSteps.begin();

    deltatemps = CumulativeTimeSteps[size()-1];

    return indexSortedTimeSteps[n];
}

}// namespace DLCA

