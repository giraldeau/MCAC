#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "agg.h"
#include <QFileDialog>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string.h>
#include <dirent.h>
#include <list>
#include <QTime>
#include <stdio.h>

#include <cv.h>
#include <cxcore.h>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

QString FichierParam = "___";
QString pathParam = "___";
char commentaires[500];
char CheminDATA[500];
char CheminSauve[500];
int PasAnalyse;

const double PI = atan(1.0)*4;

float coef = 20;               //en pixel/nm
int typeAnalyse = 0;           //Spécifie si l'analyse BOX counting s'effectue sur la surface pleine (=0), sur le périmètre extérieur (=1) ou périmètre intérieur et extérieur (=2)
int maxelem = 500;             //Dimension tableau BoxCounting
int Nborientations = 1;        //Nombre d'orientations considérées pour chaque agrégat
double SizeStructElem = 0.0;   //Recouvrement en pixel des sphérules (corrigé par la suite par un opérateur morphologique d'érosion) >2 recommandé
double Cov2D=0.0;              //Cov à ajouter en post traitement pour l'affichage 2D
int Nbr = 1;                   //Nombre de point dans le graph d'autocorrelation 3D
int DefTailleImageMax = 5000;  //Largeur maximale des images générées
double DensiteDipolaire = 1.0; //Densité dipolaire utilisée pour l'export 3D (DF3)
double Cov = 0.0;              //Paramètre de recouvrement utilisé pour l'export 3D (DF3)
double alphagangue = 0.0;      //Paramètre de gangue utilisé pour l'export 3D (DF3)

int NDatas;
int* NumFichier;
int* NumAgg;
int* Np;
float* Bc;
float* NBc;
float* DEDM;
float* SEDM;
Agg* agg;
float* RayonTab;
float* Volint;

float coefmem;
char path[500], sauve[500], NomComplet[500], NomComplet2[500], NomComplet3[500];
int N, Mode, DeltaSauve, Nagg;
double T, Dpm, sigmaDpm, FV, P, Rho, temps;
double* TabLabel;
double* TabNp;
double* TabDp;
double* Tabposx;
double* Tabposy;
double* Tabposz;
double* TabRg;
double* TabDm;
double* Tablpm;
double* Tabdeltat;
double* TabRgeo;
double* TabXg;
double* TabYg;
double* TabZg;
double* nTab;
FILE* f;
FILE* f1;
FILE* f2;
bool with_dots;


void test_locale()
{
    double testfloat = 1.5;
    char* teststr1 = "1.5";
    char* teststr2 = "1,5";
    double test1=atof(teststr1);
    double test2=atof(teststr2);

    if (fabs(test1-testfloat)<1e-3)
        with_dots = true;
    else if (fabs(test2-testfloat)<1e-3)
            with_dots = false;
    else
    {
        printf("What locale are you using ?\n");
        exit(1);
    }
}

double latof(const char* string)
{
    std::string mystring = string;
    if (!with_dots)
    {
        int f = mystring.find(".");
        if (f>0)
            mystring.replace(f, 1, ",");
    }
    return atof(mystring.c_str());
}


void LectureParam()
{   
    char t1[500], com[500];
    f = fopen(qPrintable(FichierParam), "rt"); //Lecture du fichier des paramètres

    test_locale();

    fgets(com, 500, f);
    sscanf(com, "%s  %s", t1, com);
    N = atoi(t1);
    fgets(com, 500, f);
    sscanf(com, "%s  %s", t1, com);
    T = latof(t1);
    fgets(com, 500, f);
    sscanf(com, "%s  %s", t1, com);
    Dpm = latof(t1);
    fgets(com,500,f);
    sscanf(com, "%s  %s", t1, com);
    sigmaDpm = latof(t1);
    fgets(com, 500 ,f);
    sscanf(com, "%s  %s", t1, com);
    FV = latof(t1)*1E-6;
    fgets(com, 500, f);
    sscanf(com, "%s  %s", t1, com);
    P = latof(t1);
    fgets(com, 500, f);
    sscanf(com, "%s  %s", t1, com);
    Rho = latof(t1);
    fgets(com, 500, f);
    sscanf(com, "%s  %s", t1, com);
    Mode = atoi(t1);
    fgets(com, 500, f);
    sscanf(com, "%s  %s", t1, com);
    DeltaSauve = atoi(t1);
    fgets(com, 500, f);
    sscanf(com, "%s  %s", sauve, com);
    fclose(f);

    strcat(CheminDATA, qPrintable(pathParam));
    strcat(CheminDATA, "//");
    strcat(CheminDATA, sauve);
    strcat(CheminSauve, qPrintable(pathParam));
    strcat(CheminSauve, "//");
    sprintf(commentaires, "N=%d\n T=%1.3f\n Dpm=%1.3f\n sigmaDpm=%1.3f\n FV=%1.3e\n P=%1.3f\n Mode=%d\n Rho=%1.3f\n DeltaSauve=%d\n CheminDATA=%s\n", N, T, Dpm, sigmaDpm, FV, P, Mode, Rho, DeltaSauve, CheminDATA);
}

void LectureFichierDat()
{
    FILE* ff;
    int ok, i;
    float f;
    char com[500];
    ff = fopen(NomComplet, "rt");
    fgets(com, 500, ff);

    ok = 1;
    NDatas = 0;
    while (ok > 0)
    {
        ok = fscanf(ff, "%d\t%d\t%d\t%e\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e", &i, &i, &i, &f, &i, &f, &f, &f, &f, &f, &f, &f, &f, &f, &f);

        if (ok > 0) NDatas++;
    }

    fclose(ff);

    NumFichier = new int[NDatas];
    NumAgg = new int[NDatas];
    Np = new int[NDatas];
    ff = fopen(NomComplet, "r");
    fgets(com, 500, ff);
    ok = 1;

    NDatas = 0;
    while (ok > 0)
    {
        ok = fscanf(ff, "%d\t%d\t%d\t%e\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e",
                    &NumFichier[NDatas],
                    &NumAgg[NDatas],
                    &Np[NDatas], &f, &i, &f, &f, &f, &f, &f, &f, &f, &f, &f, &f);
        if (ok > 0) NDatas++;        
    }

    fclose(ff);
}

void MainWindow::ProgAnalyse()
{
    int i, j, k, indexmin, indexmax;
    int nf, nf2, nf3;
    FILE* ff;
    FILE* fBC;
    FILE* fEDM;
    FILE* fAutoco;
    FILE* fDf3;
    FILE* fEstimVol;
    float tetha, phi, alea;
    float* tmpfloat;
    int* NbDipoles;
    float* RgDf3;
    float* Anisotropie;
    float* aeff;
    float* volumeaggDf3;
    int* sizeDF3;
    float* ResoDF3;
    float TabMoy[8];

    NbDipoles = new int[1];
    tmpfloat = new float[1];
    RgDf3 = new float[1];
    Anisotropie = new float[1];
    aeff = new float[1];
    ResoDF3 = new float[1];
    volumeaggDf3 = new float[1];
    sizeDF3 = new int[3];
    Bc = new float[maxelem];
    NBc = new float[maxelem];
    DEDM = new float[maxelem];
    SEDM = new float[maxelem];
    RayonTab = new float[Nbr]; //Tableau des distances associées
    Volint = new float[Nbr];   //Tableau donnant la valeur de l'autocorrelation

    agg = new Agg(N, coef);
    agg->CheminFichierSphere = CheminDATA;

    sprintf(NomComplet2, "%s/Analyse2D.txt", CheminSauve);
    ff = fopen(NomComplet2, "w");
    sprintf(commentaires,"Génération du fichier : Analyse2D.txt\n");
    ui->AfficheurRep->append(commentaires);
    fprintf(ff, "NumFichier\tNumAgg\tDpmoy\tAa\tRg2D\tPa\tPatotal\tW\tH\tL\tNp3D\tRg3D\tRgeo3D\tDm3D\tResolution--Nb Orientations : %d\n", Nborientations);

    if (ui->checkBoxBC->isChecked())
    {
        if (ui->radioButtonPerExt->isChecked()) sprintf(NomComplet2,"%s/BoxCounting_Ext_Perim.txt", CheminSauve);
        else if (ui->radioButtonPerimExtInt->isChecked()) sprintf(NomComplet2,"%s/BoxCounting_Ext-Int_Perim.txt", CheminSauve);
        else sprintf(NomComplet2,"%s/BoxCounting_Plain_Surf.txt", CheminSauve);

        fBC = fopen(NomComplet2, "w");
        sprintf(commentaires, "Génération du fichier : BoxCounting.txt\n");
        ui->AfficheurRep->append(commentaires);
        fprintf(fBC, "NumFichier\tNumAgg\tNbOrientations\tLmax\tDpmoy\tNbLignes\tReso [pixel/nm]\nLBox\tNbBox\n");
    }

    if (ui->checkBoxEDM->isChecked())
    {
        if (ui->radioButtonPerExt->isChecked()) sprintf(NomComplet3, "%s/EDM_Ext_Perim.txt", CheminSauve);
        else if (ui->radioButtonPerimExtInt->isChecked()) sprintf(NomComplet3, "%s/EDM_Ext-Int_Perim.txt", CheminSauve);
        else sprintf(NomComplet3, "%s/EDM_Plain_Surf.txt", CheminSauve);

        fEDM = fopen(NomComplet3, "w");
        sprintf(commentaires, "Génération du fichier : EDM.txt\n");
        ui->AfficheurRep->append(commentaires);
        fprintf(fEDM, "NumFichier\tNumAgg\tNbOrientations\tLmax\tDpmoy\tNbLignes\tReso [pixel/nm]\nDpix\tSurfpix\n");
    }

    if (ui->checkBoxAutoco->isChecked())
    {
        sprintf(NomComplet2, "%s/Autocorrelation.txt", CheminSauve);
        fAutoco = fopen(NomComplet2, "w");
        sprintf(commentaires, "Génération du fichier : Autocorrelation.txt\n");
        ui->AfficheurRep->append(commentaires);
        fprintf(fAutoco, "NumFichier\tNumAgg\tRg3D\tRgeo3D\tDpmoy\tNbpointsGraph\nr\tAutoco\n");
    }

    if (ui->checkBoxDF3->isChecked())
    {
        sprintf(NomComplet2, "%s/Grandeurs_Realistes.txt", CheminSauve);
        fDf3 = fopen(NomComplet2, "a");
        fprintf(fDf3, "Generation d'agregats realistes avec Cov : %f   Param-Gangue : %f    DensitéDipolaire : %f\n", Cov, alphagangue, DensiteDipolaire); //Ligne d'entête du fichier DF3
        fprintf(fDf3, "Rang\tNp\tDpmoy\tNbDipoles\tRg(nm)\tAnisotropie\taeff(nm)\tVolumeAgg(nm3)\tNbX\tNbY\tNbZ\tResolution_DF3\n");//Ligne d'entête du fichier DF3
        sprintf(commentaires, "Un fichier DF3 va être généré pour chaque agrégat étudié ainsi qu'un fichier pour DDSCAT\nLes grandeurs associées sont reportées dans le fichier Grandeurs_Realistes.txt");
        ui->AfficheurRep->append(commentaires);

        sprintf(NomComplet2, "%s/Estimation_Volumes.txt", CheminSauve);
        fEstimVol = fopen(NomComplet2, "w");
        fprintf(fEstimVol, "Nc_Np\tVolDLCA_VolWO\tVolReel_VolWO\n");//Ligne d'entête du fichier d'analyse de l'estimation du volume par le programme DLCA
    }

    indexmin = ui->spinBoxChoixAggregatmin->value();
    indexmax = ui->spinBoxChoixAggregatmax->value();
    qsrand(QTime::currentTime().msec()); //time(NULL)

    for (i = indexmin-1; i <= indexmax-1; i++) //i<NDatas
    {
        ui->lcdNumberAggCourant->display(i-indexmin+2);

        for (k = 0; k < 8; k++)     TabMoy[k] = 0; //Initialisation du tableau pour le calcul des valeurs moyennes

        for (j = 0; j < Nborientations; j++)       //On va effectuer l'analyse pour différentes orientations des agrégats
        {
            ui->progressBar->setValue(((double)j)/((double)Nborientations)*100);

            alea = (float)qrand()/(float)RAND_MAX;
            tetha = alea*2*PI;
            alea = (float)qrand()/(float)RAND_MAX;
            phi = acos(1-2*alea);


            agg->coef = coefmem;
            agg->TailleImageMax = DefTailleImageMax;

            agg->LectureAgg(NumFichier[i], NumAgg[i], tetha, phi, SizeStructElem, Cov2D, typeAnalyse);
            //agg->LectureAgg(NumFichier[i], NumAgg[i], 0.0, 0.0, SizeStructElem, Cov2D, typeAnalyse);
            agg->Traitement(typeAnalyse);

            sprintf(commentaires, "coef:%f  W:%d  L:%f\n", agg->coef, agg->imgAffichage->width, agg->L);
            ui->AfficheurRep->append(commentaires);

            if (ui->checkBoxBC->isChecked())
            {
                nf = agg->Fractal(maxelem, Bc, NBc);
                fprintf(fBC, "%d\t%d\t%d\t%e\t%e\t%d\t%e\n", NumFichier[i], NumAgg[i], j+1, agg->L, agg->Dpmoy, nf, agg->coef);//Début d'une ligne pour l'analyse BOX Counting

                if (typeAnalyse == 1 || typeAnalyse == 2) //cas des analyses périmètriques : Le nombre de boites décrit un périmètres : une longueur
                    for (k = 0; k < nf; k++) fprintf(fBC, "%e\t%e\n", Bc[k], NBc[k]*agg->coef);  //description résultat BC
                else
                    for (k = 0; k < nf; k++) fprintf(fBC, "%e\t%e\n", Bc[k], NBc[k]);  //cas des analyses surfaciques : Le nombre de boites décrit une surface
            }

            if (ui->checkBoxEDM->isChecked())
            {
                nf3 = agg->AnalyseEDM(maxelem,DEDM,SEDM);
                fprintf(fEDM, "%d\t%d\t%d\t%e\t%e\t%d\t%e\n", NumFichier[i], NumAgg[i], j+1, agg->L, agg->Dpmoy, nf3, agg->coef);//Début d'une ligne pour l'analyse EDM
                if (typeAnalyse == 1 || typeAnalyse == 2) //cas des analyses périmètriques : Le nombre de boites décrit un périmètres : une longueur
                    for (k = 0; k < nf3; k++) fprintf(fBC, "%e\t%e\n", DEDM[k], SEDM[k]);  //description résultat BC
                else
                    for (k = 0; k < nf3; k++) fprintf(fEDM, "%e\t%e\n", DEDM[k], SEDM[k]);  //cas des analyses surfaciques : Le nombre de boites décrit une surface
            }

            agg->Infos();

            if (ui->radioButtonAfficheimage->isChecked())
            {
                sprintf(NomComplet2, "%s/Agg_%6.5d.tif", CheminSauve, i+1); //Nom du Fichier image généré
                agg->Affiche(ui->checkBoxEDM->isChecked(),NomComplet2);
            }

            QCoreApplication::processEvents();

            TabMoy[0] = TabMoy[0] + agg->Aa;
            TabMoy[1] = TabMoy[1] + agg->Rg2D;
            TabMoy[2] = TabMoy[2] + agg->Pa;
            TabMoy[3] = TabMoy[3] + agg->Patotal;
            TabMoy[4] = TabMoy[4] + agg->W;
            TabMoy[5] = TabMoy[5] + agg->H;
            TabMoy[6] = TabMoy[6] + agg->L;
            TabMoy[7] = TabMoy[7] + agg->coef;

            agg->FinTraitement(coefmem, ui->checkBoxEDM->isChecked());
        }

        agg->Aa = TabMoy[0]/((float)Nborientations);
        agg->Rg2D = TabMoy[1]/((float)Nborientations);
        agg->Pa = TabMoy[2]/((float)Nborientations);
        agg->Patotal = TabMoy[3]/((float)Nborientations);
        agg->W = TabMoy[4]/((float)Nborientations);
        agg->H = TabMoy[5]/((float)Nborientations);
        agg->L = TabMoy[6]/((float)Nborientations);
        fprintf(ff, "%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%d\t%e\t%e\t%e\t%e\n", NumFichier[i], NumAgg[i], agg->Dpmoy, agg->Aa, agg->Rg2D, agg->Pa, agg->Patotal, agg->W, agg->H, agg->L, agg->Np3D, agg->Rg3D, agg->Rgeo3D, agg->Dm3D, TabMoy[7]/((float)Nborientations));

        //Calcul d'autocorrelation3D
        if (ui->checkBoxAutoco->isChecked())
        {
            agg->LectureAgg(NumFichier[i], NumAgg[i], 0.0, 0.0, SizeStructElem, Cov2D, typeAnalyse); //Lecture des paramètres de l'agrégat sans changer son orientation
            nf2 = agg->CalculAutocorrelation3D(Nbr, Nborientations, RayonTab, Volint, tmpfloat);
            fprintf(fAutoco, "%d\t%d\t%e\t%e\t%e\t%d\n", NumFichier[i], NumAgg[i], agg->Rg3D, agg->Rgeo3D, agg->Dpmoy, nf2);//Début d'une ligne pour l'analyse d'Autocorrelation
            for (k = 0; k < nf2; k++) fprintf(fAutoco, "%e\t%e\n", RayonTab[k], Volint[k]);  //description résultat BC
        }

        //Génération du descriptif 3D de l'agrégat (format DF3)
        if (ui->checkBoxDF3->isChecked())
        {
            agg->LectureAgg(NumFichier[i], NumAgg[i], 0.0, 0.0, SizeStructElem, Cov2D, typeAnalyse); //Lecture des paramètres de l'agrégat sans changer son orientation
            sprintf(NomComplet2, "%s/Agg_%6.5d", CheminSauve, i+1); //Nom du Fichier Df3 généré
            agg->ProgGenerationDF3(DensiteDipolaire, Cov, alphagangue, commentaires, NomComplet2, NbDipoles, RgDf3, Anisotropie, aeff, volumeaggDf3, sizeDF3, ResoDF3);
            fprintf(fDf3, "%d\t%d\t%e\t%d\t%e\t%e\t%e\t%e\t%d\t%d\t%d\t%e\n", i+1, agg->Np3D, agg->Dpmoy, NbDipoles[0], RgDf3[0], Anisotropie[0], aeff[0], volumeaggDf3[0], sizeDF3[0], sizeDF3[1], sizeDF3[2], ResoDF3[0]); //Ajout d'une ligne pour les résultats associés au DF3
            ui->AfficheurRep->append(commentaires);
            fprintf(fEstimVol, "%e\t%e\t%e\n",(double)(agg->Nc)/(double)(agg->Np3D),agg->VolDLCA/agg->VolWO,volumeaggDf3[0]*1E-27/agg->VolWO);
        }
    }

    fclose(ff);
    if (ui->checkBoxBC->isChecked()) fclose(fBC);
    if (ui->checkBoxEDM->isChecked()) fclose(fEDM);
    if (ui->checkBoxAutoco->isChecked()) fclose(fAutoco);
    if (ui->checkBoxDF3->isChecked()) {fclose(fDf3);fclose(fEstimVol);}
    delete[] Volint;
    delete[] RayonTab;
    delete[] Bc;
    delete[] NBc;
    agg->Libere();
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->progressBar->setVisible(false);
    ui->lcdNumberAggCourant->setVisible(false);
    ui->lcdNumber_AggMax->setVisible(false);
    ui->line->setVisible(false);
    ui->label_3->setVisible(false);
    ui->label_4->setVisible(false);

    connect(ui->pushButton,SIGNAL(clicked()),this,SLOT(BoutonQuitter()));
    connect(ui->BoutonRep,SIGNAL(clicked()),this,SLOT(BoutonRechercheParam()));
    connect(ui->BoutonExecAnalyse,SIGNAL(clicked()),this,SLOT(ExecuterAnalyse()));
    connect(ui->checkBoxBC,SIGNAL(clicked()),this,SLOT(GestionBC()));
    connect(ui->checkBoxEDM,SIGNAL(clicked()),this,SLOT(GestionEDM()));
    connect(ui->checkBoxDF3,SIGNAL(clicked()),this,SLOT(GestionDF3()));
    connect(ui->checkBoxAutoco,SIGNAL(clicked()),this,SLOT(GestionAutoco()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::GestionBC() //Gère les options liées aux calculs box counting
{
    if (ui->checkBoxBC->isChecked() || ui->checkBoxEDM->isChecked())
    {
        ui->radioButtonBCPlain->setEnabled(true);
        ui->radioButtonPerExt->setEnabled(true);
        ui->radioButtonPerimExtInt->setEnabled(true);
    }
    else
    {
        ui->radioButtonBCPlain->setEnabled(false);
        ui->radioButtonPerExt->setEnabled(false);
        ui->radioButtonPerimExtInt->setEnabled(false);
    }
}

void MainWindow::GestionEDM() //Gère les options liées à l'analyse EDM
{
    if (ui->checkBoxEDM->isChecked() || ui->checkBoxBC->isChecked())
    {
        ui->radioButtonBCPlain->setEnabled(true);
        ui->radioButtonPerExt->setEnabled(true);
        ui->radioButtonPerimExtInt->setEnabled(true);
    }
    else
    {
        ui->radioButtonBCPlain->setEnabled(false);
        ui->radioButtonPerExt->setEnabled(false);
        ui->radioButtonPerimExtInt->setEnabled(false);
    }
}

void MainWindow::GestionDF3() //Gère les options liées aux calculs box counting
{
    if (ui->checkBoxDF3->isChecked())
    {
        ui->doubleSpinBoxCov->setEnabled(true);
        ui->label_15->setEnabled(true);
        ui->doubleSpinBoxGangue->setEnabled(true);
        ui->label_16->setEnabled(true);
        ui->doubleSpinBoxDensiteDipolaire->setEnabled(true);
        ui->label_17->setEnabled(true);
    }
    else
    {
        ui->doubleSpinBoxCov->setEnabled(false);
        ui->label_15->setEnabled(false);
        ui->doubleSpinBoxGangue->setEnabled(false);
        ui->label_16->setEnabled(false);
        ui->label_17->setEnabled(false);
        ui->doubleSpinBoxDensiteDipolaire->setEnabled(false);
    }
}

void MainWindow::GestionAutoco() //Gère les options liées aux calculs box counting
{
    if (ui->checkBoxAutoco->isChecked())
    {
        ui->spinBoxNbPAutoco->setEnabled(true);
        ui->label_10->setEnabled(true);
    }
    else
    {
        ui->spinBoxNbPAutoco->setEnabled(false);
        ui->label_10->setEnabled(false);
    }
}

void MainWindow::BoutonQuitter()
{
    this->close();
    exit(0);
}

void MainWindow::BoutonRechercheParam()
{
    ui->progressBar->setVisible(false);
    FichierParam = QFileDialog::getOpenFileName(this,"Sélectionner le fichier de données DLCA",
                                                     "C:/Users/dlca/Desktop/DLCA_sous_Qt",
                                                     "Fichier de paramètres (*.par *.txt *.dat)");
    QFileInfo tmp2 = FichierParam;
    pathParam = tmp2.absolutePath(); //Cette variable ne retient que le chemin du fichier param

    //Affiche le chemin dans la zone de texte
    ui->AfficheurRep->append(FichierParam);
    LectureParam();
    ui->AfficheurRep->append(commentaires);

    sprintf(NomComplet, "%s/NumFichier_NumAgg_Np_Npe_Nc_Dg-Dp_Dm-Dp_Dgeo-Dp_Dpmoy_Dpmoy3_Vol_Surf_Tv_vraicov_SurfsurVol.dat", CheminSauve); //Fichier contenant la liste d'agrégats à traiter
    f = fopen(NomComplet, "rt"); //Fichier contenant la liste d'agrégats à traiter
    fclose(f);
    if (f == NULL)  //Le fichier n'a pas été généré
    {
        sprintf(commentaires,"----Attention le fichier NumFichier_NumAgg_Np_Npe_Nc_Dg-Dp_Dm-Dp_Dgeo-Dp_Dpmoy_Dpmoy3_Vol_Surf_Tv_vraicov_SurfsurVol.dat n'a pas été généré----\n");
        ui->AfficheurRep->append(commentaires);
        sprintf(commentaires,"----Le programme Analyse_DLCA doit être exécuté au préalable----\n");
        ui->AfficheurRep->append(commentaires);
    }
    else
    {
        sprintf(commentaires,"----Le fichier NumFichier_NumAgg_Np_Npe_Nc_Dg-Dp_Dm-Dp_Dgeo-Dp_Dpmoy_Dpmoy3_Vol_Surf_Tv_vraicov_SurfsurVol.dat a été trouvé !----\n");
        ui->AfficheurRep->append(commentaires);
        LectureFichierDat();
        sprintf(commentaires, "Il y a %d agrégats sélectionnés à analyser\n", NDatas);
        ui->AfficheurRep->append(commentaires);
        ui->label_11->setEnabled(true);
        ui->label_12->setEnabled(true);
        ui->label_13->setEnabled(true);
        ui->BoutonExecAnalyse->setEnabled(true);
        ui->spinBoxChoixAggregatmin->setEnabled(true);
        ui->spinBoxChoixAggregatmin->setMaximum(NDatas);
        ui->spinBoxChoixAggregatmax->setMaximum(NDatas);
        ui->spinBoxChoixAggregatmax->setValue(NDatas);
        ui->spinBoxChoixAggregatmax->setEnabled(true);
    }
}

void MainWindow::ExecuterAnalyse()
{
    Nbr = ui->spinBoxNbPAutoco->value(); //On lit le nombre de points souhaités pour le graphe d'autocorrelation

    Nborientations = ui->spinBoxNbOrientations->value();

    coefmem = ui->doubleSpinBoxCoef->value();
    Cov2D = ui->doubleSpinBoxCov2D->value(); //On lit le Cov 2D pour post Traitement
    SizeStructElem = ui->doubleSpinBoxSizeSE->value(); //On lit la taille de l'élément structurant utilisé pour gérer l'effet sintering
    ui->lcdNumber_AggMax->display(ui->spinBoxChoixAggregatmax->value()-ui->spinBoxChoixAggregatmin->value()+1); //On affiche le nombre d'agrégats étudiés
    sprintf(commentaires, "Nombre d'orientations considérées : %d\n", Nborientations);
    ui->AfficheurRep->append(commentaires);

    DensiteDipolaire = ui->doubleSpinBoxDensiteDipolaire->value(); //Lecture Densité dipolaire utilisée pour l'export 3D (DF3)
    Cov = ui->doubleSpinBoxCov->value();                           //Lecture du paramètre de recouvrement utilisé pour l'export 3D (DF3)
    alphagangue = ui->doubleSpinBoxGangue->value();                //Lecture Paramètre de gangue utilisé pour l'export 3D (DF3)

    ui->progressBar->setMaximum(100);
    ui->progressBar->setVisible(true);
    ui->progressBar->setValue(0);
    ui->lcdNumberAggCourant->setVisible(true);
    ui->lcdNumber_AggMax->setVisible(true);
    ui->line->setVisible((true));
    ui->label_3->setVisible(true);
    ui->label_4->setVisible(true);
    if (ui->radioButtonPerExt->isChecked()) typeAnalyse = 1;
    else if (ui->radioButtonPerimExtInt->isChecked()) typeAnalyse = 2;
    else typeAnalyse = 0;
    ProgAnalyse();
    ui->progressBar->setVisible(false);
    ui->lcdNumberAggCourant->setVisible(false);
    ui->lcdNumber_AggMax->setVisible(false);
    ui->line->setVisible(false);
    ui->label_3->setVisible(false);
    ui->label_4->setVisible(false);
    sprintf(commentaires,"Fin de l'analyse ... \n");
    ui->AfficheurRep->append(commentaires);
}
