#include "math.h"
#include "agg.h"
#include <QTime>
#include <stdio.h>

const double PI = atan(1.0)*4;

float VolCalotte(float h, float R)
{
    float VolCal;
    VolCal = PI*(pow(R,2.0)*h-pow(R,3.0)/3.0+1.0/3.0*pow(R-h,3.0));
    return VolCal;
}

float Volintercept(float xi,float yi,float zi,float Rp,float R)
{
    float Volint, d, hp, h, Rmin, Rmax;

    if (Rp < R)     {Rmin = Rp; Rmax = R;}
    else        {Rmin = R; Rmax = Rp;}

    //Volume intersection entre une sphère de rayon R centrée en (0,0,0) et une autre de rayon Rp centrée en (xi,yi,zi)
    d = sqrt(pow(xi,2.0)+pow(yi,2.0)+pow(zi,2.0));

    if (d >= (Rmin+Rmax))
        Volint = 0;
    else if ((d > Rmin && d < Rmax) || d == 0.0)
        Volint = 4.0*PI*pow(Rmin, 3.0)/3.0;
    else
    {
        hp = (pow(Rmax,2.0)-pow(Rmin-d,2.0))/2.0/d;
        h = (pow(Rmin,2.0)-pow(Rmax-d,2.0))/2.0/d;
        Volint = VolCalotte(hp, Rmin)+VolCalotte(h, Rmax);
    }

    return Volint;
}

float VolinterceptComplet(int Np1, float* X1, float* Y1, float* Z1, float* R1, int Np2, float* X2, float* Y2, float* Z2, float* R2) //Calcul le volume intersection entre deux agregats
{
    int i, j;
    float Volint, xi, yi, zi, dist;

    Volint = 0.0;

    for (i = 0; i < Np1; i++)
        for (j = 0; j < Np2; j++)
        {
            xi = X2[j]-X1[i];
            yi = Y2[j]-Y1[i];
            zi = Z2[j]-Z1[i];
            dist = sqrt(pow(xi,2.0) + pow(yi,2.0) + pow(zi,2.0));
            if (dist <= (R1[i]+R2[j])) Volint = Volint + Volintercept(xi, yi, zi, R2[j], R1[i]);
        }

    return Volint;
}

int Agg::CalculAutocorrelation3D(int Nbr, int Nborientations, float* RayonTab, float* Volint, float* tmpfloat)
{
    float Rmin, Rmax, alea, phi, tetha, Volmoyenne, rayon, VolumeTotal;
    int  test, i, j, k;
    float* Volume;

    //Tableaux représentant le clone de l'agrégat étudié
    float* CloneR;
    float* CloneX;
    float* CloneY;
    float* CloneZ;
    CloneR = new float[Np3D];
    CloneX = new float[Np3D];
    CloneY = new float[Np3D];
    CloneZ = new float[Np3D];

    Volume = new float[Nborientations]; //Tableau pour le calcul du volume moyen
    Rmin = Dpmoy/10;
    Rmax = Rgeo3D*3;

    VolumeTotal = VolinterceptComplet(Np3D, X, Y, Z, R, Np3D, X, Y, Z, R); //Volume total de l'agrégat

    tmpfloat[1] = (float)VolumeTotal;

    //On déplace l'agregat par rapport à lui-même pour le calcul de son l'autocorrelation
    test = 1;
    k = 0;
    while (k < Nbr && test != 0)
    {
        //printf('avancement en rayon %d\n',(k-1)/(Nbr-1)*100)
        rayon = Rmin*pow(Rmax/Rmin,(double)(k)/(double)(Nbr-1)); //le décalage suit une loi log en vu de l'affichage final en log-log
        //rayon=(float)(k)/(float)(Nbr-1)*(Rmax-Rmin)+Rmin;  //décalage linéaire

        for (j = 0; j < Nborientations; j++)
        {
            alea = (float)qrand() / (float)RAND_MAX;
            tetha = alea*2*PI;
            alea = (float)qrand() / (float)RAND_MAX;
            phi = acos(1-2*alea);
            for (i = 0; i < Np3D; i++)
            {
                CloneX[i]=X[i]+rayon*sin(phi)*cos(tetha);
                CloneY[i]=Y[i]+rayon*sin(phi)*sin(tetha);
                CloneZ[i]=Z[i]+rayon*cos(phi);
                CloneR[i]=R[i];
            }
            Volume[j] = VolinterceptComplet(Np3D, X, Y, Z, R, Np3D, CloneX, CloneY, CloneZ, CloneR);
        }

        for (j = 0, Volmoyenne = 0; j < Nborientations; j++) Volmoyenne = Volmoyenne + Volume[j];
        Volmoyenne = Volmoyenne/(float)Nborientations;

        if (Volmoyenne != 0)
        {
            RayonTab[k] = rayon;
            Volint[k] = Volmoyenne/VolumeTotal; //On normalise à 1 l'autocorrelation
            k = k+1;
            test = 1;
        }
        else test = 0;
    }
    return k;

    delete[] CloneR;
    delete[] CloneX;
    delete[] CloneY;
    delete[] CloneZ;
}

