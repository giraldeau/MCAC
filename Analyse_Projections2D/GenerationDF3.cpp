#include "math.h"
#include "agg.h"
#include <stdio.h>

#include <iostream>
#include <vector>
#include <string.h>
#include <fstream>

//Utilisation de librairies d'algèbre linéaire  :
#include <Dense>
#include <Eigenvalues>
using Eigen::Matrix2f;
using Eigen::Vector3cf;

float RechercheMaxTab(int Nbelem, float* Tab) //Calcule la valeur max d'un tableau
{
   int i;
   float max;

   max = 0.0;

   for (i = 0; i < Nbelem; i++)
   {
       if (Tab[i] > max)
           max = Tab[i];
   }

   return max;
}

float RechercheMinTab(int Nbelem, float* Tab) //Calcule la valeur max d'un tableau
{
   int i;
   float min;

   min = 1E20;

   for (i = 0; i < Nbelem; i++)
   {
       if (Tab[i] < min)
           min = Tab[i];
   }

   return min;
}

float CalculTabMoy(int Nbelem, float* Tab) //Calcul la valeur moyenne d'un tableau
{
   int i;
   float moy;

   moy = 0.0;

   for (i = 0; i < Nbelem; i++)
   {
       moy = moy + Tab[i];
   }

   return moy/((float)Nbelem);
}

float f3Dval(int Np, float x, float y, float z, float* TabX, float* TabY , float* TabZ, float* TabR, float alphagangue) //fonction LevelSet permettant l'évaluation 3D de la gangue
{
   int i;
   float d, val;

   val = 0.0;

   for (i = 0; i < Np; i++)
   {
       d = sqrt(pow(x-TabX[i],2.0)+pow(y-TabY[i],2.0)+pow(z-TabZ[i],2.0));
       val = val+0.5+0.5*erf(-(d/TabR[i]-1)/alphagangue);
   }

   return val;
}

void Agg::ProgGenerationDF3(float DensiteDipolaire, float Cov, float alphagangue, char* commentaires, char* NomFichier, int* NbDipoles, float* Rg, float* Anisotropie, float* aeff, float* volumeagg, int* sizeDF3, float* ResoDF3)
{
    float gama, Rmax, decal, Xmin, Ymin, Zmin, Xmax, Ymax, Zmax, Rmoyen, val, valmax, volumevoxel;
    int i, j, k, NbsliceX, NbsliceY, NbsliceZ, NbCells;
    size_t tmp;
    char NomFinal[500];
    float sumx, sumy, sumz, A, B, C, D, E, F, x, y, z;
    FILE* g;
    FILE* g2;

    int NbDipolestmp;
    float Rgtmp, Anisotropietmp, aefftmp, volumeaggtmp;

    //Tableaux représentant le clone de l'agrégat étudié
    float* CloneR;
    float* CloneX;
    float* CloneY;
    float* CloneZ;
    float* tabx;
    float* taby;
    float* tabz;
    int* tabDipoleX;
    int* tabDipoleY;
    int* tabDipoleZ;
    int* tabTypeDipole;

    CloneR = new float[Np3D];
    CloneX = new float[Np3D];
    CloneY = new float[Np3D];
    CloneZ = new float[Np3D];

    using namespace Eigen;
    Matrix3f mat;
    std::complex<float> valpropre;
    float norm1, norm2, norm3, normmax, normmin;

    gama = 1/(1-Cov); //Taux d'acroissement des rayons utilisé pour générer le recouvrement


    for (i = 0; i < Np3D; i++)
    {
        CloneX[i] = X[i]/gama;
        CloneY[i] = Y[i]/gama;
        CloneZ[i] = Z[i]/gama;
        CloneR[i] = R[i];
    }


    Rmax = RechercheMaxTab(Np3D, CloneR); //Rayon maximum (utilisé pour calculer le décalage entre l'agrégat et le bord du domaine)
    Xmin = RechercheMinTab(Np3D, CloneX); //coordonnées minimale en X
    Ymin = RechercheMinTab(Np3D, CloneY); //coordonnées minimale en Y
    Zmin = RechercheMinTab(Np3D, CloneZ); //coordonnées minimale en Z
    Rmoyen = CalculTabMoy(Np3D, CloneR);
    decal = 2*Rmax;

    //On décale l'agrégat pour le recentrer dans un domaine de coordonnées positives
    //Ecriture du fichier shape.txt pour autres exports
    sprintf(NomFinal, "%s_shape.txt", NomFichier);
    g2=fopen(NomFinal, "w"); // ouverture fichier pour generer un fichier de donnees pour représenter les agrégats sans necking
    fprintf(g2,"X\tY\tZ\tR   (agregat sans necking)\n");
    for (i = 0; i < Np3D; i++)
    {
        CloneX[i] = CloneX[i] - Xmin + decal;
        CloneY[i] = CloneY[i] - Ymin + decal;
        CloneZ[i] = CloneZ[i] - Zmin + decal;        
        fprintf(g2,"%f\t%f\t%f\t%f\n",CloneX[i],CloneY[i],CloneZ[i],CloneR[i]);
    }
    //Recherche des coordonnées max
    Xmax = RechercheMaxTab(Np3D, CloneX) + decal;
    Ymax = RechercheMaxTab(Np3D, CloneY) + decal;
    Zmax = RechercheMaxTab(Np3D, CloneZ) + decal;

    //On définie la dimension en X,, Y et Z de l'image 3D (DF3)
    NbsliceX = int(DensiteDipolaire*Xmax/2/Rmoyen)+1; sizeDF3[0] = NbsliceX;
    NbsliceY = int(DensiteDipolaire*Ymax/2/Rmoyen)+1; sizeDF3[1] = NbsliceY;
    NbsliceZ = int(DensiteDipolaire*Zmax/2/Rmoyen)+1; sizeDF3[2] = NbsliceZ;

    //On définit des tableaux de correspondance spatiale de chaque tranche de l'image 3D
    tabx = new float[NbsliceX];
    taby = new float[NbsliceY];
    tabz = new float[NbsliceZ];
    for (i = 0; i < NbsliceX; i++) tabx[i] = (float)i/(float)(NbsliceX)*Xmax;
    for (i = 0; i < NbsliceY; i++) taby[i] = (float)i/(float)(NbsliceY)*Ymax;
    for (i = 0; i < NbsliceZ; i++) tabz[i] = (float)i/(float)(NbsliceZ)*Zmax;
    volumevoxel = Xmax*Ymax*Zmax/(float)(NbsliceX)/(float)(NbsliceY)/(float)(NbsliceZ);
    ResoDF3[0] = pow(volumevoxel, -(float)1/(float)3); //Calcul de la résolution des images DF3 en voxel/nm

     //Génération de l'entête du fichier DF3
    NbCells = NbsliceX*NbsliceY*NbsliceZ;
    std::vector<char> data(6, 0);
    data[0] = NbsliceX/256; data[1] = NbsliceX%256;
    data[2] = NbsliceY/256; data[3] = NbsliceY%256;
    data[4] = NbsliceZ/256; data[5] = NbsliceZ%256;

    sprintf(NomFinal,"%s_Pov.df3",NomFichier);
    std::ofstream os(NomFinal, std::ios::binary);
    os.write(&data[0], data.size());

    //Préparation des tableaux nécéssaires à l'export DDSCAT
    tabDipoleX = new int[NbCells];
    tabDipoleY = new int[NbCells];
    tabDipoleZ = new int[NbCells];
    tabTypeDipole = new int[NbCells];
    NbDipolestmp = 0; valmax = 0;

    //On calcule pour chaque point de l'image 3D la valeur de la fonction LevelSet
    for (i = 0; i < NbsliceZ; i++)
        for (j = 0; j < NbsliceY; j++)
            for (k = 0; k < NbsliceX; k++)
            {
                val = f3Dval(Np3D, tabx[k], taby[j], tabz[i], CloneX, CloneY, CloneZ, CloneR, alphagangue);
                if (val>valmax) valmax = val;

                //Export fichier DF3 : Cas de l'export en 16bit
                if (val>1) tmp = (size_t)65535;
                else tmp = (size_t)(65535.0*val);
                data[0] = tmp/256; data[1] = tmp%256;
                os.write(&data[0],2*sizeof(char));

                //Export vers DDSCAT
                if (val >= 0.5) //intérieur de l'agrégat
                {
                    tabDipoleX[NbDipolestmp] = k;
                    tabDipoleY[NbDipolestmp] = j;
                    tabDipoleZ[NbDipolestmp] = i;
                    tabTypeDipole[NbDipolestmp] = 1;
                    NbDipolestmp++;
                }
//              else if (val>=0.3 && val<0.5) //Couche extérieure de l'agrégat
//              {
//                  tabDipoleX[NbDipolestmp] = k;
//                  tabDipoleY[NbDipolestmp] = j;
//                  tabDipoleZ[NbDipolestmp] = i;
//                  tabTypeDipole[NbDipolestmp] = 2;
//                  NbDipolestmp++;
//              }
            }
    os.close();

    volumeaggtmp = volumevoxel*(float)NbDipolestmp;
    aefftmp = pow(3.0/4.0/3.1415927*(double)volumeaggtmp, (double)1.0/(double)3.0);

    //On reprend le résultat afin de calculer la matrice d'inertie de l'agrégat
    sumx = 0; sumy = 0; sumz = 0; A = 0; B = 0; C = 0; D = 0; E = 0; F = 0;
    for (i = 0; i < NbDipolestmp; i++)
    {
        x = tabx[tabDipoleX[i]];
        y = taby[tabDipoleY[i]];
        z = tabz[tabDipoleZ[i]];
        sumx = sumx + x;
        sumy = sumy + y;
        sumz = sumz + z;
    }
    for (i = 0; i < NbDipolestmp; i++)
    {
        x = tabx[tabDipoleX[i]] - sumx/((float)NbDipolestmp);
        y = taby[tabDipoleY[i]] - sumy/((float)NbDipolestmp);
        z = tabz[tabDipoleZ[i]] - sumz/((float)NbDipolestmp);
        A = A + pow((double)y,2.0) + pow((double)z,2.0);
        B = B + pow((double)x,2.0) + pow((double)z,2.0);
        C = C + pow((double)x,2.0) + pow((double)y,2.0);
        D = D + y*z;
        E = E + x*z;
        F = F + x*y;
    }

    //On défini la matrice d'inertie
    mat << A/((float)NbDipolestmp), -F/((float)NbDipolestmp), -E/((float)NbDipolestmp), -F/((float)NbDipolestmp), B/((float)NbDipolestmp), -D/((float)NbDipolestmp), -E/((float)NbDipolestmp), -D/((float)NbDipolestmp), C/((float)NbDipolestmp);

    fprintf(g2,"Elements Matrice Inertie\n");
    fprintf(g2,"%f\t%f\t%f\n",A/((float)NbDipolestmp),-F/((float)NbDipolestmp),-E/((float)NbDipolestmp));
    fprintf(g2,"%f\t%f\t%f\n",-F/((float)NbDipolestmp),B/((float)NbDipolestmp),-D/((float)NbDipolestmp));
    fprintf(g2,"%f\t%f\t%f\n",-E/((float)NbDipolestmp),-D/((float)NbDipolestmp),C/((float)NbDipolestmp));
    fclose(g2);

    Vector3cf eivals = mat.eigenvalues(); //On calcul les valeurs propres
    valpropre = eivals[0];
    norm1 = sqrt(valpropre.real()*valpropre.real()+valpropre.imag()*valpropre.imag());
    valpropre = eivals[1];
    norm2 = sqrt(valpropre.real()*valpropre.real()+valpropre.imag()*valpropre.imag());
    valpropre = eivals[2];
    norm3 = sqrt(valpropre.real()*valpropre.real()+valpropre.imag()*valpropre.imag());
    Rgtmp = sqrt((norm1+norm2+norm3)/2);



    //Recherche des valeurs propres min et max : calcul de l'anisotropie
    normmax = 0; normmin = 1E20;
    if (norm1 > normmax) normmax = norm1;
    if (norm2 > normmax) normmax = norm2;
    if (norm3 > normmax) normmax = norm3;
    if (norm1 < normmin) normmin = norm1;
    if (norm2 < normmin) normmin = norm2;
    if (norm3 < normmin) normmin = norm3;
    Anisotropietmp = normmax/normmin;

    sprintf(commentaires,"Param alpha : %e, Cov : %e, Rg (nm): %e, Anisotropie : %e,  aeff (nm) : %e", alphagangue, Cov, Rgtmp, Anisotropietmp, aefftmp); //Commentaire dans le fichier shape.dat

    //Ecriture du fichier shape.dat pour DDSCAT
    sprintf(NomFinal, "%s_shape.dat", NomFichier);
    g=fopen(NomFinal, "w"); // ouverture fichier pour generer un fichier de donnees DDSCAT "FROM_FILE"
    fprintf(g,"%s\n", commentaires);
    fprintf(g,"%d\t= NAT\n", NbDipolestmp);
    fprintf(g,"1.000000 0.000000 0.000000 = targer vector a1 (in TF)\n");
    fprintf(g,"0.000000 1.000000 0.000000 = targer vector a2 (in TF)\n");
    fprintf(g,"%f\t%f\t%f\t= lattice spacings (d_x,d_y,d_z)/d\n", 1.0, 1.0, 1.0);
    fprintf(g,"%f\t%f\t%f\t= lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d for dipole 0 0 0\n", -(float)NbsliceX/2.0, -(float)NbsliceY/2.0, -(float)NbsliceZ/2.0);
    fprintf(g,"JA IX IY IZ ICOMP(x,y,z)\n");
    for (i = 0; i < NbDipolestmp; i++)
    {
          fprintf(g,"%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i, tabDipoleX[i], tabDipoleY[i], tabDipoleZ[i], tabTypeDipole[i], tabTypeDipole[i], tabTypeDipole[i]);
    }
    fclose(g);



    NbDipoles[0] = NbDipolestmp;
    Rg[0] = Rgtmp;
    Anisotropie[0] = Anisotropietmp;
    aeff[0] = aefftmp;
    volumeagg[0] = volumeaggtmp;

    sprintf(commentaires, "Np : %d   Cov : %f   Gangue : %f   DD : %f   \nNb Dipoles : %d  Rg: %f nm   Anisotropie:%f   aeff:%f nm  VolumeAgg %f nm3\n", Np3D, Cov, alphagangue, DensiteDipolaire, NbDipolestmp, Rgtmp, Anisotropietmp, aefftmp, volumeaggtmp);

    delete[] CloneR;
    delete[] CloneX;
    delete[] CloneY;
    delete[] CloneZ;
    delete[] tabx;
    delete[] taby;
    delete[] tabz;
    delete[] tabDipoleX;
    delete[] tabDipoleY;
    delete[] tabDipoleZ;
    delete[] tabTypeDipole;
}
