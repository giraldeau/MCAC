#include "statistics.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW2(a) ((a)*(a))
#define POW3(a) ((a)*(a)*(a))

using namespace std;

namespace DLCA{

const double PI = atan(1.0)*4;


void Aggregate::Statistics()
{
    *Dp=0.;
    //$ For the Spheres i in Agg Id
    for (size_t i = 0; i < Np; i++)
    {
        *Dp += myspheres[i].Radius();
    }
    *Dp = 2*(*Dp) / double(Np);
    *DgOverDp=2*(*rg)/(*Dp);
}



AnalizedAggregate::AnalizedAggregate(const Aggregate& source):
    Aggregate::Aggregate(source),
    Nc(0),
    Dp3(0.),
    DmOverDp(0.),
    DgeoOverDp(0.),
    SurfaceOverVolume(0.),
    cov(0.),
    Npe(0.),
    Tv(0.),
    Ts(0.)
{
    Analyze();
}

void AnalizedAggregate::Analyze()
{
    //$ For the Spheres i in Agg Id
    for (size_t i = 0; i < Np; i++)
    {
        for (size_t j = i+1; j < Np; j++) //for the j spheres composing Aggregate n°id
        {
            double dist = distances[i][j];
            double rpmoy = (myspheres[i].Radius() + myspheres[j].Radius())*0.5; //Mean Radius between i and j monomeres
            double dbordabord = dist-2.0*rpmoy; //distance between the two particles
            //$ Check if i is covering j
            //$ [dbordabord <= 1]
            if (dbordabord < 0.)
            {
                //$ Calculation of the Number of contacts
                cov = cov - dbordabord/(2.0*rpmoy); //Coefficient of total Covering of Agg id
                Nc += 1; //Nombre de coordination (nombre de points de contact entre deux sphérules)
            }
        }
        Tv = Tv + volumes[i]/(myspheres[i].Volume());
        Ts = Ts + surfaces[i]/(myspheres[i].Surface());
        Dp3 += POW3(2*myspheres[i].Radius());
    }

    if (Nc > 0)
    {
        cov = cov/double(Nc);
    }

    Tv = 1 - Tv / double(Np);
    Ts = 1 - Ts / double(Np);
    Dp3 = pow(Dp3 / double(Np) , 1./3.);

    DmOverDp=*dm/(*Dp);
    DgeoOverDp=2*(*rmax)/(*Dp);

    SurfaceOverVolume = *surfAgregat / *volAgregat;
    Npe = *volAgregat/(PI * POW3(Dp3)/6.0);
}


__attribute__((pure)) bool Aggregate::operator <(const Aggregate& a) const
{
   double d = *DgOverDp - *a.DgOverDp;
   if (abs(d)<0.1)
   {
       d=0;
   }
   return d<0;
}


bool StatisticStorage::InsertIfNew(const Aggregate& Agg)
{

    if (Agg.Np <10)
    {
        return false;
    }

    if (! targetNp[Agg.Np] ||
        ! binary_search (SavedAggregates.begin(), SavedAggregates.end(), Agg)
       )
    {
        AnalizedAggregate Analyzed(Agg);

        size_t n(0);

        // Search for the first bigger DgOverDp
        n = lower_bound(SavedAggregates.begin(), SavedAggregates.end(),Agg)
                    - SavedAggregates.begin();

        // If this is the last I don't know which is bigger
        if (n == SavedAggregates.size()-1 && !SavedAggregates.empty() && SavedAggregates[n]<Agg)
        {
                SavedAggregates.push_back(Analyzed);
        }
        else
        {
            SavedAggregates.insert(SavedAggregates.begin()+n,Analyzed);
        }
        targetNp[Agg.Np]=true;
        return true;
    }
    return false;
}


StatisticStorage::StatisticStorage(PhysicalModel& _physicalmodel):
    targetNp(),
    SavedAggregates(),
    physicalmodel(&_physicalmodel)
{
}

void StatisticStorage::Init()
{
    targetNp.assign(physicalmodel->N,false);
}

void StatisticStorage::Analyze(const ListAggregat& current)
{
//    return;
    for (const Aggregate* Agg : current)
    {
        InsertIfNew(*Agg);
    }
}



void StatisticStorage::print() const
{
    ofstream myfile;
    myfile.open ("testStats.dat", ios::out | ios::trunc);
    myfile << "Current number of saved aggregates : " << SavedAggregates.size() << endl;
    for (AnalizedAggregate Agg : SavedAggregates)
    {
        myfile << *Agg.DgOverDp << "\t" << Agg.Np << endl;
    }
    myfile.close();

}

/** Destructor */
StatisticStorage::~StatisticStorage() noexcept /* explicitly specified destructors should be annotated noexcept as best-practice */
{
    if(physicalmodel->toBeDestroyed)
    {
        delete physicalmodel;
    }
}

}  // namespace DLCA

