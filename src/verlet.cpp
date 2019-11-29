#include "verlet.hpp"
#include "physical_model.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

namespace MCAC{


void Verlet::Remove(const size_t id,const array<size_t, 3> Index)
{
    Cell[Index[0]][Index[1]][Index[2]].remove(id);
}
void Verlet::Add(const size_t id,const array<size_t, 3> Index)
{
    Cell[Index[0]][Index[1]][Index[2]].push_front(id);
}

void Verlet::Init(const size_t _GridDiv, const double _L)
{

    Cell.clear();

    GridDiv = _GridDiv;
    L=_L;

    Cell.resize(GridDiv);
    for(vector< vector< list< size_t > > >& Vx : Cell)
    {
        Vx.resize(GridDiv);

        for(vector< list< size_t > >& Vy : Vx)
        {
            Vy.resize(GridDiv);
        }
    }

}


vector<size_t> Verlet::GetSearchSpace(const array<double, 3> sourceposition , const double width) const
{

    double xp{sourceposition[0]+width};
    double xm{sourceposition[0]-width};
    double yp{sourceposition[1]+width};
    double ym{sourceposition[1]-width};
    double zp{sourceposition[2]+width};
    double zm{sourceposition[2]-width};

    auto bornei1 {static_cast<int>(floor(xm*static_cast<double>(GridDiv)/L))};
    auto bornei2 {static_cast<int>(floor(xp*static_cast<double>(GridDiv)/L)+1)};
    auto bornej1 {static_cast<int>(floor(ym*static_cast<double>(GridDiv)/L))};
    auto bornej2 {static_cast<int>(floor(yp*static_cast<double>(GridDiv)/L)+1)};
    auto bornek1 {static_cast<int>(floor(zm*static_cast<double>(GridDiv)/L))};
    auto bornek2 {static_cast<int>(floor(zp*static_cast<double>(GridDiv)/L)+1)};

    if (bornei2-bornei1>=static_cast<int>(GridDiv))
    {
        bornei1=0;
        bornei2=static_cast<int>(GridDiv)-1;
    }
    if (bornej2-bornej1>=static_cast<int>(GridDiv))
    {
        bornej1=0;
        bornej2=static_cast<int>(GridDiv)-1;
    }
    if (bornek2-bornek1>=static_cast<int>(GridDiv))
    {
        bornek1=0;
        bornek2=static_cast<int>(GridDiv)-1;
    }


    list<size_t> tmpSearchSpace;

    // ///////
    for (int i=bornei1;i<=bornei2;i++)
    {
        for (int j=bornej1;j<=bornej2;j++)
        {
            for (int k=bornek1;k<=bornek2;k++)
            {
                // periodic
                auto ii = static_cast<size_t>(periodicPosition(i,static_cast<int>(GridDiv)));
                auto jj = static_cast<size_t>(periodicPosition(j,static_cast<int>(GridDiv)));
                auto kk = static_cast<size_t>(periodicPosition(k,static_cast<int>(GridDiv)));
                tmpSearchSpace.insert(tmpSearchSpace.end(),
                                      Cell[ii][jj][kk].begin(),
                                      Cell[ii][jj][kk].end());

            }
        }
    }
 
    vector<size_t> SearchSpace{ make_move_iterator(tmpSearchSpace.begin()),
                             make_move_iterator(tmpSearchSpace.end()) };
    return SearchSpace;
}


vector<size_t> Verlet::GetSearchSpace(const array<double, 3> sourceposition,
    const double width,
    const array<double, 3> Vector) const
{

    double xp{sourceposition[0]+width + MAX(Vector[0],0)};
    double xm{sourceposition[0]-width + MIN(Vector[0],0)};
    double yp{sourceposition[1]+width + MAX(Vector[1],0)};
    double ym{sourceposition[1]-width + MIN(Vector[1],0)};
    double zp{sourceposition[2]+width + MAX(Vector[2],0)};
    double zm{sourceposition[2]-width + MIN(Vector[2],0)};

    auto bornei1 {static_cast<int>(floor(xm*static_cast<double>(GridDiv)/L))};
    auto bornei2 {static_cast<int>(floor(xp*static_cast<double>(GridDiv)/L)+1)};
    auto bornej1 {static_cast<int>(floor(ym*static_cast<double>(GridDiv)/L))};
    auto bornej2 {static_cast<int>(floor(yp*static_cast<double>(GridDiv)/L)+1)};
    auto bornek1 {static_cast<int>(floor(zm*static_cast<double>(GridDiv)/L))};
    auto bornek2 {static_cast<int>(floor(zp*static_cast<double>(GridDiv)/L)+1)};

    if (bornei2-bornei1>=static_cast<int>(GridDiv))
    {
        bornei1=0;
        bornei2=static_cast<int>(GridDiv)-1;
    }
    if (bornej2-bornej1>=static_cast<int>(GridDiv))
    {
        bornej1=0;
        bornej2=static_cast<int>(GridDiv)-1;
    }
    if (bornek2-bornek1>=static_cast<int>(GridDiv))
    {
        bornek1=0;
        bornek2=static_cast<int>(GridDiv)-1;
    }


    list<size_t> tmpSearchSpace;

    // ///////
    for (int i=bornei1;i<=bornei2;i++)
    {
        for (int j=bornej1;j<=bornej2;j++)
        {
            for (int k=bornek1;k<=bornek2;k++)
            {
                // periodic
                auto ii = static_cast<size_t>(periodicPosition(i,static_cast<int>(GridDiv)));
                auto jj = static_cast<size_t>(periodicPosition(j,static_cast<int>(GridDiv)));
                auto kk = static_cast<size_t>(periodicPosition(k,static_cast<int>(GridDiv)));
                tmpSearchSpace.insert(tmpSearchSpace.end(),
                                      Cell[ii][jj][kk].begin(),
                                      Cell[ii][jj][kk].end());

            }
        }
    }

    vector<size_t> SearchSpace{ make_move_iterator(tmpSearchSpace.begin()),
                             make_move_iterator(tmpSearchSpace.end()) };
    return SearchSpace;
}


void Verlet::print(const std::array<size_t, 3> Index) const
{
    cout << "list of agg registered in cell " << Index[0] << " "
                                              << Index[1] << " "
                                              << Index[2] << " " << endl;

    for(const size_t agg : Cell[Index[0]][Index[1]][Index[2]])
    {
        cout << agg<<endl;
    }
}

void Verlet::search(const size_t id) const
{
    cout << "Searching " << id << " in verlet"<<endl;
    for (size_t i=0;i<GridDiv;i++)
    {
        for (size_t j=0;j<GridDiv;j++)
        {
            for (size_t k=0;k<GridDiv;k++)
            {
                for(const size_t agg : Cell[i][j][k])
                {
                    if (agg==id)
                    {
                        cout << i << " "
                             << j << " "
                             << k << endl;
                    }
                }
            }
        }
    }
}

/* Default constructor */
Verlet::Verlet():
    Cell(),
    GridDiv(0),
    L(0)
{}

}  // namespace MCAC
