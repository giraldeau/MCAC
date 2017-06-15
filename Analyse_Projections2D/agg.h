#ifndef AGG_H
#define AGG_H

#include <cv.h>
#include <cxcore.h>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

class Agg
{
private:
public:
    float* R;
    float* X;
    float* Y;
    float* Z;
    int nsphere;
    int numfichier;
    int idagg;
    float coef;

    IplImage* img;
    IplImage* imgBC ;
    IplImage* imgAffichage;
    IplImage* img_Inv;
    IplImage* img_EDMshow;
    IplImage* img_EDM;
    IplImage* img_SurfInv;
    IplImage* img_Surf;
    IplImage* img_EDMdilate;
    IplImage* img_SurfInv8U;
    IplImage* img_Surfdilate;

    CvMemStorage* storage;
    CvSeq* contours;
    CvSeq* poly;
    CvBox2D boxe;

    bool ok_calcul;
    bool ok_fractal;

    void Calculmoments(IplImage* img);
    void CVDrawBox2D(IplImage* img, CvBox2D b, CvScalar co);

    CvPoint g2D;
    float W, H, L, Pa, Patotal, Aa, Rg2D, Rg3D, Rgeo3D, Dm3D, g2Dx, g2Dy, g3Dx, g3Dy, g3Dz, Dpmoy,VolDLCA,VolWO,SurfWO;
    int Np3D, TailleImageMax,Nc;

    char* CheminFichierSphere;

    Agg(int nspheres, float coef);
    void Infos();
    void Libere();
    int AnalyseEDM(int maxelem, float* TabDiam, float* Surf);
    int Fractal(int maxelem, float* BOITE, float* N);
    void LectureAgg(int numfichier, int idagg, float tetha, float phi, double SizeStructElem, double Cov2D, int typeAnalyse);
    void Traitement(int typeAnalyse);
    void FinTraitement(float coefmem, bool activeEDM);
    void Affiche(bool activeEDM, char* nomimage);
    void Reoriente(float tetha, float phi);
    //Fonctions associées à l'export 3D des agrégats (DF3)
    void ProgGenerationDF3(float DensiteDipolaire, float Cov, float alphagangue, char* commentaires, char* NomFichier, int* NbDipoles, float* Rg, float* Anisotropie, float* aeff, float* volumeagg, int* sizeDF3, float* ResoDF3);
    //Fonctions associées à l'autocorrelation3D
    int CalculAutocorrelation3D(int Nbr, int Nborientations, float* RayonTab, float* Volint, float* tmpfloat);
};

#endif // AGG_H

