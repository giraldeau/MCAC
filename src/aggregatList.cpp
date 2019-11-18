#include "aggregatList.hpp"
#include <algorithm>
#include <numeric>
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
    lastSaved(0),
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
    Writer = new ThreadedIO(_physicalmodel, _size);
    spheres.Init(_physicalmodel, _size);
    verlet.Init(_physicalmodel.GridDiv,_physicalmodel.L);

    storage_list<15,Aggregate>::Init(_size,*this);
    setpointers();
}


//########################################## Determination of the contacts between agrgates ##########################################
vector< pair<size_t,double> > ListAggregat::SortSearchSpace(size_t MovingAgg,array<double,3> Vectdir, vector<size_t> SearchSpace) const
{
    vector< pair<size_t,double> > SortedSearchSpace;

    Sphere SphereMe(*list[MovingAgg]);
    Sphere SphereOther(*list[MovingAgg]);

    //$ [For all other agregate]
    for (const size_t& AggOther : SearchSpace)
    {
        SphereOther.InitVal(*list[AggOther]->x,
                            *list[AggOther]->y,
                            *list[AggOther]->z,
                            *list[AggOther]->rmax);

        auto pos = SortedSearchSpace.begin();
        double dist = 0.;
        bool collision = SphereMe.Contact(SphereOther);

        if (!collision)
        {
            pair<bool,double> Collision = SphereMe.Collision(SphereOther,Vectdir);
            if (Collision.first && Collision.second <= *list[MovingAgg]->lpm )
            {
                collision=Collision.first;
                dist = Collision.second;
                pos = lower_bound (SortedSearchSpace.begin(), SortedSearchSpace.end(), dist,
                                   [](pair<size_t,double> i1, double d){return i1.second < d;});
            }
        }
        if (collision)
        {
            pair<size_t,double> suspect = {AggOther,dist};
            SortedSearchSpace.insert(pos,suspect);
        }
    }
    return SortedSearchSpace;
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

        array<double, 3> Vector{{lpm * Vectdir[0],
                                 lpm * Vectdir[1],
                                 lpm * Vectdir[2]}};

        SearchSpace = verlet.GetSearchSpace(sourceposition , mindist, Vector);

        // Remove me
        for(size_t i =0;i<SearchSpace.size();i++)
        {
            if (SearchSpace[i]==source)
            {
                SearchSpace.erase(SearchSpace.begin()+long(i));
                return SearchSpace;
            }
        }
        cout << "I'm not on the verlet list ???"<<endl;
        cout << "This is an error"<<endl;
        exit(68);
        //return SearchSpace;
    }

    // The full aggregat list index
    SearchSpace.resize(size());
    iota(SearchSpace.begin(), SearchSpace.end(), 0);

    // Except me
    SearchSpace.erase(SearchSpace.begin()+long(source));

    return SearchSpace;
}


//########################################## Determination of the contacts between agrgates ##########################################
bool ListAggregat::TestFreeSpace(const std::array<double,3> pos, const double diameter) const
{
    // Use Verlet to reduce search
    double mindist(diameter*0.5 + maxradius);

    vector<size_t> SearchSpace = verlet.GetSearchSpace(pos , mindist);

    double mindist2 = POW2(mindist);

    //$ loop on the agregates potentially in contact
    for (const size_t & suspect : SearchSpace) //For every aggregate that could be in contact
    {
        //$ Loop on all the spheres of the aggregate
        for (const Sphere* sphere : list[suspect]->myspheres)
        {
            if (sphere->Distance2(pos[0], pos[1], pos[2]) <= mindist2)
            {
                return true;
            }
        }
    }
    return false;
}


//########################################## Determination of the contacts between agrgates ##########################################
int ListAggregat::DistanceToNextContact(const size_t source, const array<double,3> Vectdir, double &distmin) const
{
    // Use Verlet to reduce search
    vector<size_t> SearchSpace(GetSearchSpace(source,Vectdir));

    // Assimilate Aggregate as sphere to drasticly speed-up search
    vector< pair<size_t,double> > Potential(SortSearchSpace(source,Vectdir,SearchSpace));

    int aggcontact(-1);
    distmin = *list[source]->lpm;
    double mindistagg(0.);

    //$ loop on the agregates potentially in contact
    for (const pair<size_t,double> & suspect : Potential) //For every aggregate that could be in contact
    {
        size_t agg=suspect.first;
        double distagg=suspect.second;

        // We already found the closest one
        if (aggcontact>=0)
        {
            double secu = 2*(*list[size_t(aggcontact)]->rmax + *list[source]->rmax);
            if (distagg - mindistagg > secu)
            {
                return aggcontact;
            }
        }
        /*
        // If two aggragates are already in contact due to surface growing
        if (list[source]->Contact(*list[agg]))
        {
            distmin = 0;
            cout << "Surprise there is already a contact !" << endl;
            cout << "This is probably due to surface growing" << endl;
            return int(agg);
        }
        */
        double dist = list[source]->Distance(*list[agg],Vectdir);
        if (dist >= 0 && dist <= distmin)
        {
            distmin = dist;
            aggcontact = int(agg); //Prise en compte de l'image par translation d'un agrégat cible
            mindistagg = distagg;
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
    }

    delete Writer;

    //#pragma omp simd
    for (Aggregate* Agg : list)
    {
        Agg->InVerlet=false;
    }
}

size_t ListAggregat::Merge(const size_t first, const size_t second)
{
    const size_t keeped(MIN(first,second));
    const size_t removed(MAX(first,second));

    // compute proper time of the final aggregate
    // keeping global time constant
    double newtime = double(size() - 1) * (*list[keeped]->time + *list[removed]->time) / double(size()) - physicalmodel->Time;

    // Merge the two aggregate but do not remove the deleted one
    list[keeped]->Merge(*list[removed]);

    remove(removed);

    setpointers();

    *list[keeped]->time = newtime;

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

    #pragma omp simd
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

__attribute__((pure)) size_t ListAggregat::RandomPick(const double random)
{
    //$ Pick a random sphere
    double valAlea = random*CumulativeTimeSteps[size()-1];
    long n = lower_bound(CumulativeTimeSteps.begin(), CumulativeTimeSteps.end(), valAlea) - CumulativeTimeSteps.begin();

    size_t NumAgg = indexSortedTimeSteps[size_t(n)];

    return NumAgg;
}

__attribute__((pure)) double ListAggregat::GetTimeStep(const double max) const
{
    double deltatemps = max / CumulativeTimeSteps[size()-1];

    return deltatemps ;
}

__attribute__((pure)) size_t ListAggregat::PickLast()
{
      double time = *list[0]->time;
      size_t latest = 0;
      for (const Aggregate* Agg : list)
      {
          if (*Agg->time < time)
          {
              time = *Agg->time;
              latest = Agg->Label;
          }
      }

      return latest;
}

ListAggregat::ListAggregat(PhysicalModel& _physicalmodel, const size_t _size):
    storage_list<15,Aggregate>(),
    physicalmodel(&_physicalmodel),
    maxradius(0.),
    indexSortedTimeSteps(),
    CumulativeTimeSteps(),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    Writer(new ThreadedIO(*physicalmodel, _size)),
    lastSaved(0),
    spheres(),
    verlet()
{
    spheres.Init(_physicalmodel, _size);
    verlet.Init(_physicalmodel.GridDiv,_physicalmodel.L);

    storage_list<15,Aggregate>::Init(_size,*this);
    setpointers();
}


Aggregate* ListAggregat::add(const Aggregate& oldAgg)
{
    Aggregate* newAgg = storage_list<15,Aggregate>::add(oldAgg,*this);

    setpointers();
    for (Aggregate* Agg : list)
    {
        Agg->setpointers();
        for (Sphere* Sph : Agg->myspheres)
        {
            Sph->setpointers();
        }
    }

    return newAgg;
}

void ListAggregat::Duplication()
{
    size_t oldNAgg = size();
    double oldL = physicalmodel->L;
    physicalmodel->L *= 2;
    physicalmodel->N *= 8;

    // TODO Rework this in order not to to it aggregate by aggregate

    for (size_t iagg=0;iagg<oldNAgg; iagg++)
    {
        for (int i=0;i<=1; i++){
            for (int j=0;j<=1; j++){
                for (int k=0;k<=1; k++){
                    if (i!=0 or j!=0 or k!=0){

                        Aggregate* newAgg = storage_list<15,Aggregate>::add(*list[iagg],*this);
                        newAgg->Label = size() - 1;
                        newAgg->verlet = &verlet;
                        newAgg->InVerlet = true;

                        setpointers();
                        for (Aggregate* Agg : list)
                        {
                            Agg->setpointers();
                            for (Sphere* Sph : Agg->myspheres)
                            {
                                Sph->setpointers();
                            }
                        }
                        newAgg->verlet = &verlet;
                        newAgg->InVerlet = true;

                        array<double, 3> vec_move={i*oldL,
                                                   j*oldL,
                                                   k*oldL};
                        newAgg->Translate(vec_move);
                    }
                }
            }
        }
    }
    //$ Update Verlet
    verlet.Init(physicalmodel->GridDiv, physicalmodel->L);
    for (Aggregate* Agg : list){
        Agg->IndexVerlet = Agg->GetVerletIndex();
        verlet.Add(Agg->GetLabel(), Agg->IndexVerlet);
    }
}

}// namespace DLCA

