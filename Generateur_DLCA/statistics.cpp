#include "statistics.h"
#include <math.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW2(a) ((a)*(a))
#define POW3(a) ((a)*(a)*(a))

using namespace std;

namespace DLCA{

const double PI = atan(1.0)*4;

AnalizeAggregate::AnalizeAggregate(const Aggregate& source):
    Aggregate::Aggregate(source)
{
    cov = 0.0;// Coeficient of mean covering
    Nc = 0; // Number of contacts
    Tv = 0.0;
    Ts = 0.0;
    Dp = 0.0;
    Dp3 = 0.0;

    //$ For the Spheres i in Agg Id
    for (int i = 0; i < Np; i++)
    {
        for (int j = i+1; j < Np; j++) //for the j spheres composing Aggregate n°id
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
        Dp += 2*myspheres[i].Radius();
        Dp3 += POW3(2*myspheres[i].Radius());
    }

    if (Nc > 0)
    {
        cov = cov/Nc;
    }

    Tv = 1 - Tv / Np;
    Ts = 1 - Ts / Np;
    Dp = Dp / Np;
    Dp3 = pow(Dp3 / Np , 1./3.);

    DgOverDp=2*(*rg)/Dp;
    DmOverDp=2*(*dm)/Dp;
    DgeoOverDp=*rmax/Dp;
    SurfaceOverVolume = *surfAgregat / *volAgregat;
    Npe = *volAgregat/(PI * POW3(Dp3)/6.0);
}


bool Statistics::keep(const Aggregate& Agg) const
{
    return !target[Agg.Np-1];
}


Statistics::Statistics(PhysicalModel& _physicalmodel):
    physicalmodel(&_physicalmodel)
{
}

void Statistics::Init(void)
{
    target.assign(physicalmodel->N,false);
}



void Statistics::Analyze(const ListAggregat& current)
{
    for (const Aggregate* Agg : current)
    {
        if(keep(*Agg))
        {
            AnalizeAggregate Analized(*Agg);
            append(Analized);
        }
    }
}


void Statistics::append(const AnalizeAggregate& Agg)
{
    SavedAggregates.push_back(Agg);
    target[Agg.Np-1] = true;
}


void Statistics::print(void) const
{
    ofstream myfile;
    myfile.open ("testStats.dat", ios::out | ios::trunc);
    myfile << "Current number of saved aggregates : " << SavedAggregates.size() << endl;
    for (size_t i=0 ; i<SavedAggregates.size();i++)
    {
        myfile << SavedAggregates[i].DgOverDp << "\t" << SavedAggregates[i].Np << endl;
    }
    myfile.close();

}
}
