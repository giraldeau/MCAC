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

    auto ret = FractalLaw[Agg.Np].emplace(*Agg.DgOverDp);

    return ret.second;

    /*
    if (SavedAggregates[Agg.Np].empty())
    {
        AnalizedAggregate Analyzed(Agg);
        SavedAggregates[Agg.Np].push_back(Analyzed);
        return true;
    }

    if ( ! binary_search (SavedAggregates[Agg.Np].begin(), SavedAggregates[Agg.Np].end(), Agg)  )
    {
        AnalizedAggregate Analyzed(Agg);

        size_t n(0);

        // Search for the first bigger DgOverDp
        n = lower_bound(SavedAggregates[Agg.Np].begin(), SavedAggregates[Agg.Np].end(),Agg)
                    - SavedAggregates[Agg.Np].begin();

        // If this is the last I don't know which is bigger
        if (n == SavedAggregates[Agg.Np].size()-1 && !SavedAggregates[Agg.Np].empty() && SavedAggregates[Agg.Np][n]<Agg)
        {
                SavedAggregates[Agg.Np].push_back(Analyzed);
        }
        else
        {
            SavedAggregates[Agg.Np].insert(SavedAggregates[Agg.Np].begin()+n,Analyzed);
        }
        return true;
    }
    return false;
    */
}

tuple<bool,double,double,double> linreg(const vector<double>& x, const vector<double>& y)
{
    double   sumx = 0.0;                        /* sum of x                      */
    double   sumx2 = 0.0;                       /* sum of x**2                   */
    double   sumxy = 0.0;                       /* sum of x * y                  */
    double   sumy = 0.0;                        /* sum of y                      */
    double   sumy2 = 0.0;                       /* sum of y**2                   */

    size_t n = x.size();


   for (size_t i=0;i<n;i++)
   {
      double X(log(x[i]));
      double Y(log(y[i]));
      sumx  += X;
      sumx2 += POW2(X);
      sumxy += X * Y;
      sumy  += Y;
      sumy2 += POW2(Y);
   }



   double denom = (n * sumx2 - POW2(sumx));
   if (n==0 || abs(denom) < 1e-9) {
       // singular matrix. can't solve the problem.
       tuple<bool,double,double,double> res{false,0.,0.,0.};
       return res;
   }

   double a = (n * sumxy  -  sumx * sumy) / denom;
   double b = (sumy * sumx2  -  sumx * sumxy) / denom;

   /* compute correlation coeff     */
   double r = (sumxy - sumx * sumy / n) /
            POW2((sumx2 - POW2(sumx)/n) *
            (sumy2 - POW2(sumy)/n));

   tuple<bool,double,double,double> res{true,a,b,r};
   return res;
}


StatisticStorage::StatisticStorage(PhysicalModel& _physicalmodel):
    SavedAggregates(),
    FractalLaw(),
    physicalmodel(&_physicalmodel)
{
}

void StatisticStorage::Init()
{
    FractalLaw.resize(physicalmodel->N);
}

void StatisticStorage::Analyze(const ListAggregat& current)
{
//    return;

    vector<double> Nps;
    vector<double> DgOverDps;
    Nps.reserve(current.size());
    DgOverDps.reserve(current.size());

    for (const Aggregate* Agg : current)
    {
        InsertIfNew(*Agg);
        Nps.push_back(double(Agg->Np));
        DgOverDps.push_back(*Agg->DgOverDp);
    }


    tuple<bool,double,double,double> InstantaneousFractalLaw = linreg(DgOverDps,Nps);
    if(get<0>(InstantaneousFractalLaw))
    {
        cout << exp(get<2>(InstantaneousFractalLaw))
             << " * x^ "
             << get<1>(InstantaneousFractalLaw)
             << "  --- r= "
             << get<3>(InstantaneousFractalLaw) << endl;
    }


}



void StatisticStorage::print() const
{
    ofstream myfile;
    myfile.open ("testStats.dat", ios::out | ios::trunc);

    size_t total(0);
    /*
    for (const set <AnalizedAggregate>& ListAgg : SavedAggregates)
    {
        total += ListAgg.size();
        for (const AnalizedAggregate& Agg : ListAgg)
        {
            myfile << *Agg.DgOverDp << "\t" << Agg.Np << endl;
        }
    }
    */
    for (size_t Np=0;Np<FractalLaw.size();Np++)
    {
        total += FractalLaw[Np].size();
        for (const double& DgOverDp : FractalLaw[Np])
        {
            myfile << DgOverDp << "\t" << Np << endl;
        }
    }
    myfile.close();
    cout << "Total number of saved aggregates : " << total << endl;


    vector<double> Nps;
    vector<double> DgOverDps;
    Nps.reserve(total);
    DgOverDps.reserve(total);
    for (size_t Np=0;Np<FractalLaw.size();Np++)
    {
        for (const double& DgOverDp : FractalLaw[Np])
        {
            Nps.push_back(double(Np));
            DgOverDps.push_back(DgOverDp);
        }
    }

    tuple<bool,double,double,double> CompleteFractalLaw = linreg(DgOverDps,Nps);
    if(get<0>(CompleteFractalLaw))
    {
        cout << exp(get<2>(CompleteFractalLaw))
             << " * x^ "
             << get<1>(CompleteFractalLaw)
             << "  --- r= "
             << get<3>(CompleteFractalLaw) << endl;
    }

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

