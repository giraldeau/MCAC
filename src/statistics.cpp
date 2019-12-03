#include "statistics.hpp"
#include "aggregats/aggregat.hpp"
#include "aggregats/aggregat_list.hpp"
#include "spheres/sphere.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW_2(a) ((a)*(a))
#define POW_3(a) ((a)*(a)*(a))

using namespace std;

namespace MCAC{

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
    for (size_t i = 0; i < n_spheres; i++)
    {
    }

    //$ For the Spheres i in Agg Id
    for (size_t i = 0; i < n_spheres; i++)
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
        //Tv = Tv + volumes[i]/(myspheres[i].get_volume());
        //Ts = Ts + surfaces[i]/(myspheres[i].get_surface());
        Dp += myspheres[i].get_radius();
        //Dp3 += POW3(2*myspheres[i].get_radius());
    }

    //if (Nc > 0)
    //{
    //    cov = cov/double(Nc);
    //}

    //Tv = 1 - Tv / double(Np);
    //Ts = 1 - Ts / double(Np);
    Dp = 2 * Dp / double(n_spheres);
    //Dp3 = pow(Dp3 / double(Np) , 1./3.);

    DgOverDp=2*(*rg)/Dp;
    //DgeoOverDp=2*(*rmax)/Dp;

    //SurfaceOverVolume = *surfAgregat / *volAgregat;
    //Npe = *volAgregat/(PI * POW3(Dp3)/6.0);
}

tuple<bool,double,double,double> ListAggregat::get_instantaneous_fractal_law() const
{
    vector<double> Nps;
    vector<double> DgOverDps;
    Nps.reserve(size());
    DgOverDps.reserve(size());

    for (const Aggregate* Agg : list)
    {
        Nps.push_back(double(Agg->size()));
        DgOverDps.push_back(Agg->DgOverDp);
    }


    return linreg(DgOverDps,Nps);
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

void StatisicsData::partialStatistics(){}
void StatisicsData::fullStatistics(){}


}  // namespace MCAC

