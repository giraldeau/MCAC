#include "statistics.hpp"
#include "aggregat.hpp"

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


void Aggregate::partialStatistics()
{
    // This routine will be called for each update of the aggregate (and before full statistics)
    fullStatistics();
}


void Aggregate::fullStatistics()
{
    // This routine will be called after each aggregation

    Dp=0.;
    //$ For the Spheres i in Agg Id
    for (size_t i = 0; i < Np; i++)
    {
    }

    //$ For the Spheres i in Agg Id
    for (size_t i = 0; i < Np; i++)
    {
        /*
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
        }*/
        //Tv = Tv + volumes[i]/(myspheres[i].Volume());
        //Ts = Ts + surfaces[i]/(myspheres[i].Surface());
        Dp += myspheres[i].Radius();
        //Dp3 += POW3(2*myspheres[i].Radius());
    }

    //if (Nc > 0)
    //{
    //    cov = cov/double(Nc);
    //}

    //Tv = 1 - Tv / double(Np);
    //Ts = 1 - Ts / double(Np);
    Dp = 2 * Dp / double(Np);
    //Dp3 = pow(Dp3 / double(Np) , 1./3.);

    DgOverDp=2*(*rg)/Dp;
    //DmOverDp=*dm/Dp;
    //DgeoOverDp=2*(*rmax)/Dp;

    //SurfaceOverVolume = *surfAgregat / *volAgregat;
    //Npe = *volAgregat/(PI * POW3(Dp3)/6.0);
}


 __attribute__((pure)) bool StatcmpAgg::operator() (const Aggregate& lhs, const Aggregate& rhs) const
{
    Statcmpdouble cmp;
    return cmp(lhs.DgOverDp, rhs.DgOverDp);
}

 __attribute__((pure)) bool Statcmpdouble::operator() (const double& lhs, const double& rhs) const
{
    double d = lhs/rhs-1;
    if (abs(d)<0.01)
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

    auto ret = FractalLaw[Agg.Np - 1].emplace(Agg.DgOverDp);
    if(ret.second)
    {
        Aggregate* newagg = SavedAggregates->add(Agg);
        times.push_back(newagg->physicalmodel->Time);
//        cout << "before" << endl;
//        Agg.print();
//        cout << endl << "after" << endl;
//        newagg->print();
    }

    return ret.second;
}

tuple<bool,double,double,double> linreg(const vector<double>& x, const vector<double>& y)
{
    double   sumx = 0.0;                        /* sum of x                      */
    double   sumx2 = 0.0;                       /* sum of x**2                   */
    double   sumxy = 0.0;                       /* sum of x * y                  */
    double   sumy = 0.0;                        /* sum of y                      */
    double   sumy2 = 0.0;                       /* sum of y**2                   */

    size_t n = x.size();
    auto N = double(n);

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



   double denom = (N * sumx2 - POW2(sumx));
   if (n==0 || abs(denom) < 1e-9) {
       // singular matrix. can't solve the problem.
       tuple<bool,double,double,double> res{false,0.,0.,0.};
       return res;
   }

   double a = (N * sumxy  -  sumx * sumy) / denom;
   double b = (sumy * sumx2  -  sumx * sumxy) / denom;

   /* compute correlation coeff     */
   double r = (sumxy - sumx * sumy / N) /
            POW2((sumx2 - POW2(sumx)/N) *
            (sumy2 - POW2(sumy)/N));

   tuple<bool,double,double,double> res{true,a,b,r};
   return res;
}


StatisticStorage::StatisticStorage(PhysicalModel& _physicalmodel):
    physicalmodel(&_physicalmodel),
    SavedAggregates(new ListAggregat(*physicalmodel,0)),
    FractalLaw(),
    times(),
    WriterAgg(),
    WriterSph()
{
}

void StatisticStorage::Init()
{
    FractalLaw.resize(physicalmodel->N);
    WriterAgg = new ThreadedIO(*physicalmodel, physicalmodel->N);
    WriterSph = new ThreadedIO(*physicalmodel, physicalmodel->N);

}

void StatisticStorage::Analyze(const ListAggregat& current)
{
    // Store aggregates for statistical analysis
    for (const Aggregate* Agg : current)
    {
        InsertIfNew(*Agg);
    }
}

tuple<bool,double,double,double> StatisticStorage::getInstantaneousFractalLaw(const ListAggregat& current) const
{
    vector<double> Nps;
    vector<double> DgOverDps;
    Nps.reserve(current.size());
    DgOverDps.reserve(current.size());

    for (const Aggregate* Agg : current)
    {
        Nps.push_back(double(Agg->Np));
        DgOverDps.push_back(Agg->DgOverDp);
    }


    return linreg(DgOverDps,Nps);
}

tuple<bool,double,double,double> StatisticStorage::getCompleteFractalLaw() const
{
    size_t total(0);
    for (size_t Np=1;Np<=FractalLaw.size();Np++)
    {
        total += FractalLaw[Np - 1].size();
    }
    vector<double> Nps;
    vector<double> DgOverDps;
    Nps.reserve(total);
    DgOverDps.reserve(total);

    for (size_t Np=1;Np<=FractalLaw.size();Np++)
    {
        for (const double& DgOverDp : FractalLaw[Np - 1])
        {
            Nps.push_back(double(Np));
            DgOverDps.push_back(DgOverDp);
        }
    }

    return linreg(DgOverDps,Nps);
}

void StatisticStorage::print() const
{

    /*
    ofstream myfile;
    myfile.open ("testStats.dat", ios::out | ios::trunc);

    for (size_t Np=1;Np<=FractalLaw.size();Np++)
    {
        for (const double& DgOverDp : FractalLaw[Np])
        {
            myfile << DgOverDp << "\t" << Np << endl;
        }
    }

    myfile.close();
    cout << "Total number of saved aggregates : " << total << endl;
    */

    size_t total(0);
    for (size_t Np=1;Np<=FractalLaw.size();Np++)
    {
        total += FractalLaw[Np - 1].size();
    }
    cout << "Total number of saved aggregates : " << total << endl;

    auto CompleteFractalLaw = getCompleteFractalLaw();
    if(get<0>(CompleteFractalLaw))
    {
        cout << "Estimation of the fractal law : " << endl << "  "
             << exp(get<2>(CompleteFractalLaw))
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
    delete WriterAgg;
    delete WriterSph;

}


void StatisicsData::partialStatistics(){}
void StatisicsData::fullStatistics(){}


}  // namespace DLCA

