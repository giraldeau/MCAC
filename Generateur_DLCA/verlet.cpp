#include "verlet.h"
#include <math.h>
#include <cmath>
#include <iomanip>


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW2(a) ((a)*(a))
#define POW3(a) ((a)*(a)*(a))


void Verlet::Remove(const int id,const array<int, 4> Index)
{
    (*this)[Index[1]][Index[2]][Index[3]].remove(id);
}
void Verlet::Add(const int id,const array<int, 4> Index)
{
    (*this)[Index[1]][Index[2]][Index[3]].push_front(id);
}

void Verlet::Init(const int _GridDiv, const double _L)
{
    (*this).clear();

    GridDiv = _GridDiv;
    L=_L;

    (*this).resize(GridDiv);
    for(vector< vector< list< int > > >& Vx : (*this))
    {
        Vx.resize(GridDiv);

        for(vector< list< int > >& Vy : Vx)
        {
            Vy.resize(GridDiv);
        }

    }
}


vector<int> Verlet::GetSearchSpace(const array<double, 4> sourceposition , const double witdh, const array<double, 4> Vector) const
{

    double xp(sourceposition[1]+witdh+MAX(Vector[1],0));
    double xm(sourceposition[1]-witdh+MIN(Vector[1],0));
    double yp(sourceposition[2]+witdh+MAX(Vector[2],0));
    double ym(sourceposition[2]-witdh+MIN(Vector[2],0));
    double zp(sourceposition[3]+witdh+MAX(Vector[3],0));
    double zm(sourceposition[3]-witdh+MIN(Vector[3],0));

    int bornei1 (int(floor(xm*GridDiv/L)));
    int bornei2 (int(floor(xp*GridDiv/L)+1));
    int bornej1 (int(floor(ym*GridDiv/L)));
    int bornej2 (int(floor(yp*GridDiv/L)+1));
    int bornek1 (int(floor(zm*GridDiv/L)));
    int bornek2 (int(floor(zp*GridDiv/L)+1));

    if (bornei2-bornei1>=GridDiv)
        {bornei1=0 ; bornei2=GridDiv-1;}
    if (bornej2-bornej1>=GridDiv)
        {bornej1=0 ; bornej2=GridDiv-1;}
    if (bornek2-bornek1>=GridDiv)
        {bornek1=0 ; bornek2=GridDiv-1;}


    list<int> tmpSearchSpace;

    // ///////
    for (int i=bornei1;i<=bornei2;i++)
    {
        for (int j=bornej1;j<=bornej2;j++)
        {
            for (int k=bornek1;k<=bornek2;k++)
            {
                int ii(i),jj(j),kk(k);

                // periodic
                while (ii<0)
                    ii += GridDiv;
                while (jj<0)
                    jj += GridDiv;
                while (kk<0)
                    kk += GridDiv;
                while (ii>=GridDiv)
                    ii -= GridDiv;
                while (jj>=GridDiv)
                    jj -= GridDiv;
                while (kk>=GridDiv)
                    kk -= GridDiv;

                tmpSearchSpace.insert(tmpSearchSpace.end(),
                                      (*this)[ii][jj][kk].begin(),
                                      (*this)[ii][jj][kk].end());
            }
        }
    }
    vector<int> SearchSpace{ make_move_iterator(tmpSearchSpace.begin()),
                             make_move_iterator(tmpSearchSpace.end()) };

    return SearchSpace;
}


