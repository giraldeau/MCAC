#include "agg.h"
#include <stdio.h>
#include <iostream>

using namespace cv;
using namespace std;

Agg::Agg(int nspheres, float coef)
{
    nsphere = nspheres;
    this->coef = coef;
    R = new float[nsphere];
    X = new float[nsphere];
    Y = new float[nsphere];
    Z = new float[nsphere];
    CheminFichierSphere = new char[500];
}

void Agg::Libere()
{
    delete[] R;
    delete[] X;
    delete[] Y;
    delete[] Z;
}

void Agg::Infos()
{
    printf("Infos: numfichier:%d idagg=%d\n", numfichier, idagg);
    if (ok_calcul)
        printf("Rg2D:%6.1f Pa:%6.1f Aa:%6.1f G:%d x %d W:%6.1f H:%6.1f L:%6.1f\n", Rg2D, Pa, Aa, g2D.x, g2D.y, W, H, L);
}

void Agg::CVDrawBox2D(IplImage* img, CvBox2D b, CvScalar co)
{
    CvPoint p1, p2, p3, p4;
    CvPoint2D32f pt[4];

    cvBoxPoints(b, pt);

    p1.x = (int)pt[0].x;
    p2.x = (int)pt[1].x;
    p3.x = (int)pt[2].x;
    p4.x = (int)pt[3].x;

    p1.y = (int)pt[0].y;
    p2.y = (int)pt[1].y;
    p3.y = (int)pt[2].y;
    p4.y = (int)pt[3].y;

    cvLine(img, p1, p2, co, 1, 8, 0);
    cvLine(img, p2, p3, co, 1, 8, 0);
    cvLine(img, p3, p4, co, 1, 8, 0);
    cvLine(img, p4, p1, co, 1, 8, 0);
}

void Agg::Affiche(bool activeEDM, char *nomimage)
{
    CvPoint pa, pb;
    CvScalar ca;
    CvScalar cb;

    IplImage* img_Inv = NULL;

    cvNamedWindow("imgAffichage", 0);
    ca = cvScalar(0, 0, 255, 0);
    cb = cvScalar(255, 0, 0, 0);
    //cvDrawContours(imgAffichage, contours, cb, cb, 0, 1, 8);
    //cvDrawContours(imgAffichage, poly, ca, ca, 0, 1, 8);
    //CVDrawBox2D(imgAffichage, boxe, ca);
    pa.x = g2D.x-10;
    pb.x = g2D.x+10;
    pa.y = g2D.y;
    pb.y = pa.y;
    //cvLine(imgAffichage, pa, pb, ca, 1, 8, 0);
    pa.x = g2D.x;
    pb.x = pa.x;
    pa.y = g2D.y-10;
    pb.y = g2D.y+10;
    //cvLine(imgAffichage, pa, pb, ca, 1, 8, 0);
    cvShowImage("imgAffichage", imgAffichage);

    if (activeEDM == 1)
    {
        cvNamedWindow("DistanceEuclidienne", 0);
        cvShowImage("DistanceEuclidienne", img_EDMshow);
    }

    //Inversion de l'image imgAffichage pour avoir noir sur fond blanc
    img_Inv = cvCreateImage(cvGetSize(imgAffichage), IPL_DEPTH_8U, 1);
    cvThreshold(imgAffichage, img_Inv, 155, 255, CV_THRESH_BINARY_INV);
    cvSaveImage(nomimage,img_Inv);

}

void Agg::LectureAgg(int numfichier, int idagg, float tetha, float phi, double nmSizeStructElem, double Cov2D, int typeAnalyse)
{
    char* NomFichier;
    FILE* f;
    int i, label, np, dimx, dimy, tmpint, SizeStructElem;
    char t[200];
    float r, x, y, z, d, xmin, ymin, xmax, ymax, tmp, tmp1, tmp2, tmp3;

    np = 0;
    this->numfichier = numfichier;
    this->idagg = idagg;

    //Lecture des caractéristiques 3D de l'agrégat considéré dans le fichier Agg
    NomFichier = new char[500];
    sprintf(NomFichier, "%s/Agg%05d.txt", CheminFichierSphere, numfichier);
    f = fopen(NomFichier, "r");
    delete[] NomFichier;
    for (i = 0; i < 29; i++) fscanf(f, "%20s", t); //On saute l'entête
    for (i = 0; i < idagg-1; i++) fscanf(f, "%f%d%d%f%f%f%f%f%f%f%f%f%f%f%f", &tmp, &tmpint, &tmpint, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp); //On saute les lignes inutiles
    fscanf(f, "%f%d%d%f%f%f%f%f%f%f%f%f%f%f%f", &Rg3D, &Np3D, &Nc, &Dm3D, &tmp, &tmp, &Rgeo3D, &tmp1, &tmp2, &tmp3, &VolDLCA, &tmp, &VolWO, &tmp, &SurfWO); //Lecture des paramètres de l'agrégat
    fclose(f);
    g3Dx = tmp1;
    g3Dy = tmp2;
    g3Dz = tmp3;
    VolDLCA=VolDLCA*1E-25;
    VolWO=VolWO*1E-25;
    SurfWO=SurfWO*1E-16;

    //Lecture des positions 3D des sphérules de l'agrégat considéré dans le fichier sphère
    NomFichier = new char[500];
    sprintf(NomFichier, "%s/Sphere%05d.txt", CheminFichierSphere, numfichier);
    f = fopen(NomFichier, "r");
    delete[] NomFichier;
    Dpmoy = 0.0;

    for (i = 0; i < 19; i++)    fscanf(f, "%20s", t);
    for (i = 0; i < nsphere; i++)
    {
        fscanf(f, "%d%f%f%f%f", &label, &r, &x, &y, &z);
        if (label == idagg)
        {
            R[np] = r; X[np] = x*(1-Cov2D); Y[np] = y*(1-Cov2D); Z[np] = z*(1-Cov2D); np++;
            Dpmoy = Dpmoy + 2.0*r;
        }
    }
    fclose(f);

    Dpmoy = Dpmoy/(float)np;

    Reoriente(tetha, phi); //Réoriente de façon aléatoire l'agrégat

    xmin = ymin = 1E23;
    xmax = ymax = -xmin;

    for (i = 0; i < np; i++)
    {
        d = 2*R[i]; x = X[i]; y = Y[i];
        if (x-d < xmin) xmin = x-d;
        if (y-d < ymin) ymin = y-d;
        if (x+d > xmax) xmax = x+d;
        if (y+d > ymax) ymax = y+d;
    }

    xmax = xmax-xmin;
    ymax = ymax-ymin;

    //Adaptation de la résolution pour ne pas dépasser une taille maximale autorisée d'image
    if ((xmax > ymax) && ((int)(xmax*coef) > TailleImageMax))
    {
        coef = ((float)TailleImageMax)/xmax;
    }

    if ((ymax > xmax) && ((int)(ymax*coef) > TailleImageMax))
    {
        coef = ((float)TailleImageMax)/ymax;
    }

    dimx = (int)(xmax*coef);
    dimy = (int)(ymax*coef);

    //On recherche le côté le plus grand
    if (dimx < dimy)    dimx = dimy;
    else dimy = dimx;

    for (i = 0; i < np; i++)
    {
        X[i] = (X[i]-xmin)*coef;
        Y[i] = (Y[i]-ymin)*coef;
        R[i] = R[i]*coef;
    }

    img = cvCreateImage(cvSize(dimx,dimy), IPL_DEPTH_8U, 1);
    imgAffichage = cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
    imgBC = cvCreateImage(cvSize(dimx,dimy), IPL_DEPTH_8U, 1);

    cvZero(img);
    cvZero(imgBC);
    cvZero(imgAffichage);

    SizeStructElem = (int)(nmSizeStructElem*coef);

    if (SizeStructElem < 1)    {
        for (i = 0; i < np; i++)
        {
            cvCircle(img, cvPoint(X[i],Y[i]), R[i], cvScalarAll(255), -1, 8);
            cvCircle(imgAffichage, cvPoint(X[i],Y[i]), R[i], cvScalarAll(255), -1, 8);
        }
    }
    else   {
        for (i = 0; i < np; i++)
        {
            cvCircle(img,cvPoint(X[i],Y[i]),R[i]+(float)(SizeStructElem),cvScalarAll(255),-1,8);  //On crée un recouvrement des sphérules en grossissant les diamètres
            cvCircle(imgAffichage,cvPoint(X[i],Y[i]),R[i]+(float)(SizeStructElem),cvScalarAll(255),-1,8);
        }
        //On procède à une érosion des images d'agrégats de façon à supprimer la surface ajoutée par le grossissement des sphérules tout en conservant la gangue de connection des sphérules
        //CHANGEMENT
        IplConvKernel *kernel;
        kernel = cvCreateStructuringElementEx(2*SizeStructElem+1, 2*SizeStructElem+1, SizeStructElem+1, SizeStructElem+1, CV_SHAPE_ELLIPSE);
        cvErode(img, img, kernel, 1);
        cvErode(imgAffichage, imgAffichage, kernel, 1);
    }



    if (typeAnalyse == 0)   cvCopy(img, imgBC);

    //On ré-exprime les tableaux X, Y, Z et R en nm (et non plus en pixel)
    for (i = 0; i < np; i++)
    {
        X[i] = (X[i])/coef + xmin;
        Y[i] = (Y[i])/coef + ymin;
        R[i] = R[i]/coef;
    }

    ok_calcul = true;
    ok_fractal = false;
}

void Agg::Calculmoments(IplImage* img)
{
    int i, j, w, h;
    float x1, y1, x2, y2, n;
    unsigned char* b;

    w = img->width;
    h = img->height;
    n = x1 = x2 = y1 = y2 = 0.0;

    for (j = 0; j < h; j++)
    {
        b = (unsigned char*)(img->imageData+j*img->widthStep);
        for (i = 0; i < w; i++)
            if (b[i] > 0)
            {
                n = n + 1.0;
                x1 = x1 + (float)i;
                y1 = y1 + (float)j;
                x2 = x2 + pow((double)i, 2.0);
                y2 = y2 + pow((double)j, 2.0);
            }
    }

    g2Dx = x1/n;
    g2Dy = y1/n;
    g2D.x = g2Dx;
    g2D.y = g2Dy;
    Rg2D = sqrt((x2+y2)/n-pow(g2Dx, 2.0)-pow(g2Dy, 2.0))/coef;
    Aa = n/coef/coef;
}

void Agg::Traitement(int typeAnalyse)
{
    int n;
    float x, y;

    CvRect r;
    CvSeq* cc;
    CvPoint* p1;
    CvPoint* p2;
    CvScalar ca;
    CvScalar cb;

    storage = cvCreateMemStorage(0);
    contours = 0;

    ca = cvScalar(0, 0, 255, 0);
    cb = cvScalar(255, 0, 0, 0);

    //Calcul du contour complet (intérieur plus extérieur)
    Patotal = 0.0;
    n = cvFindContours(img, storage, &contours, sizeof(CvContour), CV_RETR_LIST);
    printf("contour complet : n=%d\n", n);

    if (typeAnalyse == 2)   cvDrawContours(imgBC, contours, cb, cb, 2, 1, 8);

    cc = NULL;

    for (CvSeq* c = contours; c != NULL; c = c->h_next)
    {
        Patotal = Patotal + cvContourPerimeter(c);
    }

    Patotal = Patotal/coef;

    //Calcul du contour extérieur
    n = cvFindContours(img, storage, &contours, sizeof(CvContour), CV_RETR_EXTERNAL);
    printf("contour exterieur : n=%d\n", n);

    if (typeAnalyse == 1)   cvDrawContours(imgBC, contours, cb, cb, 1, 1, 8);

    if (contours != NULL)
    {
         Calculmoments(img);
         Pa = cvContourPerimeter(contours)/coef;
         boxe = cvMinAreaRect2(contours);
         r = cvBoundingRect(contours);
         poly = cvApproxPoly(contours, sizeof(CvContour), storage, CV_POLY_APPROX_DP,
                            r.width<r.height?r.width:r.height);
         p1 = (CvPoint*)cvGetSeqElem(poly, 0);
         p2 = (CvPoint*)cvGetSeqElem(poly, 1);


         W = boxe.size.width/coef;
         H = boxe.size.height/coef;
         x = p1->x-p2->x;
         y = p1->y-p2->y;
         L = sqrt(pow(x, 2.0)+pow(y, 2.0))/coef;
    }

    ok_calcul = true;
}

int Agg::Fractal(int maxelem, float* BOITE, float* N)
{
    int id, n, w, h;
    IplImage* imgAffichage = NULL;
    n = cvCountNonZero(imgBC);
    w = imgBC->width;
    h = imgBC->height;
    float w0 = w;
    float largeurBoite = w0/w/coef;
    id = 0;

    while ((w > 1) && (h > 1) && (id < maxelem))
    {
        N[id] = n/coef/coef; //Cas de l'analyse surfacique représente la surface en nm² recouverte par les boites
        BOITE[id] = largeurBoite; //Représente la largeur de la boite en nm

        w = w/1.2;
        h = h/1.2;
        largeurBoite = w0/w/coef;

        id++;
        if (imgAffichage != NULL) cvReleaseImage(&imgAffichage);
        if ((w > 2) && (h > 2))
        {
            imgAffichage = cvCreateImage(cvSize(w, h), 8, 1);
            cvResize(imgBC, imgAffichage, CV_INTER_AREA);
            n = cvCountNonZero(imgAffichage);
            cvReleaseImage(&imgAffichage);
        }
    }

    ok_fractal = true;

    return id;
}

int Agg::AnalyseEDM(int maxelem, float* TabDiam, float* Surf)
{
    int k = 0;
    double min_val, max_val, coeffdl, max_valdilate, min_valdilate;
    //double coeff;
    char buffer[300];

    IplImage* img_EDM = cvCreateImage(cvGetSize(imgBC), IPL_DEPTH_32F, 1);
    IplImage* img_SurfInv = NULL;
    IplImage* img_Surf = NULL;
    IplImage* img_EDMdilate = NULL;
    IplImage* img_SurfInv8U = NULL;
    IplImage* img_Surfdilate = NULL;
    IplImage* img_EDMshow = cvCreateImage(cvGetSize(imgBC), IPL_DEPTH_8U, 1);

    cvDistTransform(imgBC, img_EDM, CV_DIST_WELSCH, NULL ,NULL, NULL); //Calcul De l'image EDM
    cvMinMaxLoc(img_EDM, &min_val, &max_val); //Distances max et min de la fonction EDM en pixels
    //coeff=255/max_val;
    //cvConvertScale(img_EDM, img_EDMshow, coeff, 0);

    for(k = 0; k <= (int)max_val && (k < maxelem); k++)
    {
        img_Surf = cvCreateImage(cvGetSize(img_EDMshow), IPL_DEPTH_32F, 1);
        cvThreshold(img_EDM, img_Surf, k, 255, CV_THRESH_BINARY);

        //Inversion de l'image img_Surf pour avoir noir sur fond blanc
        img_SurfInv = cvCreateImage(cvGetSize(img_EDM), IPL_DEPTH_32F, 1);
        cvThreshold(img_Surf, img_SurfInv, k, 255, CV_THRESH_BINARY_INV);

        // Application de cvDistTransform sur l'image inversé ;
        img_EDMdilate = cvCreateImage(cvGetSize(img_EDM), IPL_DEPTH_32F, 1);
        img_SurfInv8U = cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
        //je dois convertir en 8 bit imgSurfInv
        cvMinMaxLoc(img_SurfInv, &min_valdilate, &max_valdilate);
        coeffdl = 255/max_valdilate;
        cvConvertScale(img_SurfInv, img_SurfInv8U, coeffdl, 0);
        cvDistTransform(img_SurfInv8U, img_EDMdilate, CV_DIST_WELSCH, NULL, NULL, NULL);

        //binarisation de l'image EDM suivant k pour la dilatation
        img_Surfdilate = cvCreateImage(cvGetSize(img_EDMdilate), IPL_DEPTH_32F, 1);
        cvThreshold(img_EDMdilate, img_Surfdilate, k, 255, CV_THRESH_BINARY_INV);

        sprintf(buffer, "DistEuclSilhouette_%d.tif", k);
        IplImage* img_SurfVal5 = cvCreateImage(cvGetSize(img_Surfdilate), IPL_DEPTH_8U, 1);
        cvConvertScale(img_Surfdilate, img_SurfVal5, coeffdl, 0);
        cvSaveImage(buffer, img_SurfVal5);
        cvReleaseImage(&img_SurfVal5);

        Surf[k] = (float)cvCountNonZero(img_Surfdilate); //Nb de pixel dans l'image EDM seuillée
        TabDiam[k] = (float)k;
        cvReleaseImage(&img_Surf);
        cvReleaseImage(&img_SurfInv);
        cvReleaseImage(&img_EDMdilate);
        cvReleaseImage(&img_SurfInv8U);
        cvReleaseImage(&img_Surfdilate);
    }

    return max_val;

    cvReleaseImage(&img_Surf);
    cvReleaseImage(&img_SurfInv);
    cvReleaseImage(&img_EDMdilate);
    cvReleaseImage(&img_SurfInv8U);
    cvReleaseImage(&img_Surfdilate);
    cvReleaseImage(&img_EDM);
    cvReleaseImage(&img_EDMshow);
}

void Agg::FinTraitement(float coefmem, bool activeEDM)
{
    cvReleaseImage(&img);
    cvReleaseImage(&imgAffichage);
    cvReleaseImage(&imgBC);

    if (activeEDM == 1)
    {
        cvReleaseImage(&img_EDMshow);
    }

    cvReleaseMemStorage(&storage);
    coef = coefmem; //On propose à nouveau la résolution proposée par l'utilisateur
}

void Agg::Reoriente(float tetha, float phi) //Oriente de façon aléatoire l'agrégat considéré
{
    int i;
    float Xtemp[Np3D+1];
    float Ytemp[Np3D+1];
    float Ztemp[Np3D+1];

    //On calcule la matrice rotation aléatoire de l'agrégat
    for (i = 0; i < Np3D; i++)
    {
        Xtemp[i] = g3Dx+sin(phi)*cos(tetha)*(X[i]-g3Dx)-sin(tetha)*(Y[i]-g3Dy)-cos(phi)*cos(tetha)*(Z[i]-g3Dz);
        Ytemp[i] = g3Dy+sin(phi)*sin(tetha)*(X[i]-g3Dx)+cos(tetha)*(Y[i]-g3Dy)-cos(phi)*sin(tetha)*(Z[i]-g3Dz);
        Ztemp[i] = g3Dz+cos(phi)*(X[i]-g3Dx)+sin(tetha)*(Z[i]-g3Dz);
    }
    //On remplace dans le tableau les coordonnées des sphérules par les nouvelles
    for (i = 0; i < Np3D; i++)
    {
        X[i] = Xtemp[i];
        Y[i] = Ytemp[i];
        Z[i] = Ztemp[i];
    }
}
