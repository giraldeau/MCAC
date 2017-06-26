#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string.h>
#include <dirent.h>
#include <list>

using namespace std;

QString FichierParam = "___";
QString pathParam = "___";
char commentaires[500];
char CheminDATA[500];
char CheminSauve[500];
int PasAnalyse = 0;

char path[500], sauve[500], NomComplet[500];
int N, Mode, DeltaSauve, Nagg;
const double PI = atan(1.0)*4;
double T, Dpm, sigmaDpm, FV, P, Rho, X, temps, L;
double* compteur;
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
double* TabNc;
double* TabVolume;
double* TabSurface;
double* TabTv;
double* Tabcov;
//double* TabTs;
double* TabSurfsurVol;
double* nTab;

FILE* f;
FILE* f1;
FILE* f2;
FILE* f3;
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



struct TabTri
{
    int NumFichier;
    int NumAgg;
    int Np;
    double Npe;
    int Nc;
    double DgsurDp;
    double DmsurDp;
    double DgeosurDp;
    double Dpmoy;
    double Dpmoy3;
    double Volume;
    double Surface;
    double Tv;
    double cov;
    double SurfsurVol;
};

void LectureParam()
{
    char t1[500], com[500];

    test_locale();


    f = fopen(qPrintable(FichierParam), "rt"); //Lecture du fichier de paramètres
    fgets(com, 500 ,f);
    sscanf(com, "%s  %s", t1, com);
    N = atoi(t1);
    fgets(com, 500, f);
    sscanf(com, "%s  %s", t1, com);
    T = atof(t1);
    fgets(com, 500, f);
    sscanf(com, "%s  %s", t1, com);
    Dpm = atof(t1);
    fgets(com, 500, f);
    sscanf(com, "%s  %s", t1, com);
    sigmaDpm = atof(t1);
    fgets(com, 500, f);
    sscanf(com, "%s  %s", t1, com);
    FV = atof(t1)*1E-6;
    fgets(com, 500, f);
    sscanf(com, "%s  %s", t1, com);
    P = atof(t1);
    fgets(com, 500, f);
    sscanf(com,"%s  %s", t1, com);
    Rho = atof(t1);
    fgets(com, 500, f);
    sscanf(com, "%s  %s", t1, com);
    Mode = atoi(t1);
    fgets(com, 500, f);
    sscanf(com, "%s  %s", t1, com);
    DeltaSauve = atoi(t1);
    fgets(com,500,f);
    sscanf(com, "%s  %s", sauve, com);
    fclose(f);

    strcat(CheminDATA, qPrintable(pathParam));
    strcat(CheminDATA, "//");
    strcat(CheminDATA, sauve);

    strcat(CheminSauve, qPrintable(pathParam));
    strcat(CheminSauve, "//");

    sprintf(commentaires,"N=%d\nT=%1.3f\nDpm=%1.3f\nsigmaDpm=%1.3f\nFV=%1.3e\nP=%1.3f\nMode=%d\nRho=%1.3f\nDeltaSauve=%d\nCheminDATA=%s\n",N,T,Dpm,sigmaDpm,FV,P,Mode,Rho,DeltaSauve,CheminDATA);
}

void LectureSphere(int numfichier)
{
    int i;
    char t1[500], t2[500], t3[500], t4[500], t5[500], com[500];

    //Lecture du fichier Sphere
    sprintf(NomComplet,"%s/Sphere%05d.txt", CheminDATA, numfichier);
    f = fopen(NomComplet, "r");
    fgets(com, 500, f);
    fgets(com, 500, f);
    fgets(com, 500, f);
    sscanf(com, "%s  %s", t1, com);
    X = latof(t1); //Lecture du paramètre X
    L = X*Dpm*1E-9; //Largeur de la boite en m
    fgets(com, 500, f);
    fgets(com, 500, f);
    fgets(com, 500, f);
    sscanf(com, "%s  %s", t1, com);
    Nagg = atoi(t1); //Lecture du nombre d'agrégats
    fgets(com, 500, f);
    sscanf(com, "%s  %s", t1, com);
    temps = latof(t1)*1E-6; //Instant considéré
    fgets(com, 500, f);

    for (i = 1; i <= N; i++)
    {
        fgets(com,500,f);
        sscanf(com, "%s %s %s %s %s", t1, t2, t3, t4, t5);
        TabLabel[i] = latof(t1);
        TabDp[i] = latof(t2)*2E-9;
        Tabposx[i] = latof(t3)*1E-9;
        Tabposy[i] = latof(t4)*1E-9;
        Tabposz[i] = latof(t5)*1E-9;
    }

    fclose(f);
}

void LectureAggregat(int numfichier)
{
    int i;
    char t1[500], t2[500], t3[500], t4[500], t5[500], t6[500], t7[500], t8[500], t9[500], t10[500], t11[500], t12[500], t13[500], t14[500], t15[500], com[500];

    //Lecture du fichier Aggregate
    sprintf(NomComplet, "%s/Agg%05d.txt", CheminDATA, numfichier);

    f = fopen(NomComplet, "r");
    fgets(com, 500, f);
    fgets(com, 500, f);
    fgets(com, 500, f);
    fgets(com, 500, f);
    fgets(com, 500, f);
    fgets(com, 500, f);
    fgets(com, 500, f);
    fgets(com, 500, f);

    for (i = 1; i <= Nagg; i++)
    {
        fgets(com, 500, f);
        sscanf(com, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s", t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15);
        TabRg[i] = latof(t1)*1E-9;
        TabNp[i] = latof(t2);
        TabNc[i] = latof(t3);
        TabDm[i] = latof(t4)*1E-9;
        Tablpm[i] = latof(t5)*1E-9;
        Tabdeltat[i] = latof(t6)*1E-6;
        TabRgeo[i] = latof(t7)*1E-9;
        TabXg[i] = latof(t8)*1E-9;
        TabYg[i] = latof(t9)*1E-9;
        TabZg[i] = latof(t10)*1E-9;
        TabVolume[i] = latof(t11)*1E-25;
        TabSurface[i] = latof(t12)*1E-16;
        TabTv[i] = latof(t13);
        Tabcov[i] = latof(t14);
        TabSurfsurVol[i] = latof(t15)*1E9;
    }

    fclose(f);
}

void EffacerFichier(int numfichier)
{
    bool valid;

    //Lecture du fichier Sphere
    sprintf(NomComplet,"%s/Sphere%05d.txt", CheminDATA, numfichier);
    valid = QFile::remove(NomComplet);

    //Lecture du fichier Aggregat
    sprintf(NomComplet, "%s/Agg%05d.txt", CheminDATA, numfichier);
    valid = QFile::remove(NomComplet);
}

int compterFichier(DIR* dir) //Permet de calculer le nombre de fichiers dans un répertoire
{
    int nbr = 0;
    struct dirent* ent = NULL;

    while ((ent = readdir(dir)) != NULL)
    {
        if (strcmp(ent->d_name, ".") != 0 &&  /* Si le fichier lu n'est pas . */
            strcmp(ent->d_name, "..") != 0 && /*  Et n'est pas .. non plus */
            strchr(ent->d_name, '.') != NULL) /* et n'est pas un sous repertoire*/
            nbr++; /* Alors on incrémente le compteur */
    }

    return nbr;
}

/*
Calcule la valeur moyenne du tableau Tab pour les nbelem premiers éléments si critere=0 ou
si les éléments satisfont le critère reporté pour chaque élément dans le tableau TabCritere
Pond permet l'évaluation d'un moment d'ordre pond
*/
double Calculmoycritere(double* Tab, int nbelem, double* TabCritere, double critere, double* compteur, double pond)
{
    double moy;
    int i;

    moy = 0.0;
    *compteur = 0.0;

    for (i = 1; i <= nbelem; i++)
    {
        if (TabCritere[i] == critere || critere == 0.0)
        {
            moy = moy + pow(Tab[i], pond);
            *compteur = *compteur + 1.0;
        }
    }

    moy  = moy/(*compteur);

    return moy;
}

//Calcule le tableau de nombre d'agrégats ayant Np monomères
void Calculdensimonomeres(double* nTab, int Nagg, int N)
{
    int i;
    for (i = 1; i <= N; i++)      nTab[i] = 0.0; //Initialisation du tableau de densité
    for (i = 1; i <= Nagg; i++)     nTab[(int)(TabNp[i])] = nTab[(int)(TabNp[i])]+1.0;
}

void  MainWindow::ProgAnalyse()
{
    DIR* rep = NULL;
    int NbFichier, i, j, q;
    double lpmmoy, Rgeomoy, Dmmoy, lambda, invNagg, lambdasurdist, Knudsen, DpAggmoy, Dpmoy,Dpmoy2,Dpmoy3, Npmoy, Rgmoy, Ncmoy, Volmoy, Surfmoy, Tvmoy, covmoy, SurfsurVolmoy, crit;

    list<TabTri> mylist;
    list<TabTri>::iterator it;
    TabTri TT1;

    compteur = new double[1];
    TabLabel = new double[N+1];
    TabDp = new double[N+1];
    Tabposx = new double[N+1];
    Tabposy = new double[N+1];
    Tabposz = new double[N+1];
    TabRg = new double[N+1];
    TabNp = new double[N+1];
    TabNc = new double[N+1];
    TabDm = new double[N+1];
    Tablpm = new double[N+1];
    Tabdeltat = new double[N+1];
    TabRgeo = new double[N+1];
    TabXg = new double[N+1];
    TabYg = new double[N+1];
    TabZg = new double[N+1];
    TabVolume = new double[N+1];
    TabSurface = new double[N+1];
    TabTv = new double[N+1];
    Tabcov = new double[N+1];
    //TabTs = new double[N+1];
    TabSurfsurVol = new double[N+1];
    nTab = new double[N+1];

    std::locale::global(std::locale("C"));

    sprintf(commentaires,"________________________\n");
    ui->AfficheurRep->append(commentaires);

    rep = opendir(CheminDATA); // Ouverture du dossier contenant les données
    NbFichier = compterFichier(rep);	NbFichier = NbFichier/2; //Nombre de données sauvées

    //On commence par établir un fichier contenant la liste des diamètres des sphérules
    LectureSphere(0);
    Dpmoy = Calculmoycritere(TabDp, N, TabDp, 0.0, compteur, 1.0); //Calcule le Dp moyen

    sprintf(NomComplet, "%s/Distrib_Dp(nm).txt", CheminSauve);
    f = fopen(NomComplet, "w");
    sprintf(commentaires, "Génération du fichier : Distrib_Dp(nm).txt\n");
    ui->AfficheurRep->append(commentaires);

    for (i = 1; i <= N; i++)
    {
        fprintf(f, "%e\n", TabDp[i]*1E9);
    }
    fclose(f);

    //Analyse nécessitant de lire tous les fichiers temporels
    lambda = 66.5E-9*(101300/P)*(T/293.15)*(1+110/293.15)/(1+110/T);

    sprintf(NomComplet,"%s/Analyse_Temporelle.dat", CheminSauve);
    f1 = fopen(NomComplet, "w");
    fprintf(f1, "t\tNagg\tinvNagg\tlpm-dist\tKnudsen\tNpmoy\tRgmoy-Rpmoy\tRgeomoy-Rpmoy\tRmmoy-Rpmoy\tNcmoy\tVolmoy\tSurfmoy\tTvmoy\tvraicovmoy\tSurfsurVolmoy\n");
    sprintf(commentaires, "Génération du fichier : Analyse_Temporelle.dat\n");
    ui->AfficheurRep->append(commentaires);

    sprintf(NomComplet,"%s/AnalyseDistrib_t_Nagg_Npmoy_Np_X_Y.dat", CheminSauve);
    f2 = fopen(NomComplet, "w");
    fprintf(f2, "t\tNagg\tNpmoy\tNp\tX\tY\n");
    sprintf(commentaires, "Génération du fichier : AnalyseDistrib_t_Nagg_Npmoy_Np_X_Y.dat\n");
    ui->AfficheurRep->append(commentaires);

    for (i = 1, q = 0; i <= NbFichier; i = i + PasAnalyse)
    {
        //printf("Avancement %1.2f %\n",((double)i)/((double)NbFichier)*100);
        ui->progressBar->setValue(((double)i)/((double)NbFichier)*100);
        QCoreApplication::processEvents();
        LectureSphere(i-1);
        LectureAggregat(i-1);

        Dpmoy = Calculmoycritere(TabDp, N, TabDp, 0.0, compteur, 1.0);  //Calcule le Dp moyen
        Dpmoy2 = Calculmoycritere(TabDp, N, TabDp, 0.0, compteur, 2.0);  //Calcule le Dp moyen
        Dpmoy3 = Calculmoycritere(TabDp, N, TabDp, 0.0, compteur, 3.0); //Calcule le Dp moyen d'ordre 3

        invNagg = pow(L,3.0)*(1.0/((double)Nagg)-1.0/((double)N));
        lpmmoy = Calculmoycritere(Tablpm, Nagg, Tablpm, 0.0,compteur, 1.0);     //Calcule le libre parcours moyen moyenné sur l'ensemble des agrégats
        Rgeomoy = Calculmoycritere(TabRgeo, Nagg, TabRgeo, 0.0, compteur, 1.0); //Calcule le rayon d'enveloppe moyen
        lambdasurdist = lpmmoy/(L*pow((double)Nagg,(-1.0/3.0))-2*Rgeomoy);      //Renvoie le rapport lpmmoy/dist
        Dmmoy = Calculmoycritere(TabDm,Nagg,TabDm,0.0,compteur,1.0);            //Calcule le diamètre de mobilité moyen
        Knudsen = lambda/Dmmoy*2;	//Knudsen
        Npmoy = Calculmoycritere(TabNp, Nagg, TabNp, 0.0, compteur, 1.0); //Calcule le nombre de sphérules par agrégat moyen
        Rgmoy = Calculmoycritere(TabRg, Nagg, TabRg, 0.0, compteur, 1.0); //Calcule le rayon de giration moyen
        //+++++++++++++++
        Ncmoy = Calculmoycritere(TabNc, Nagg, TabNc, 0.0, compteur, 1.0);                         //Calcule le nombre de coordination moyen moyenné sur l'ensemble des agrégats
        Volmoy = Calculmoycritere(TabVolume, Nagg, TabVolume, 0.0, compteur, 1.0);                //Calcule le volume moyen moyenné sur l'ensemble des agrégats
        Surfmoy = Calculmoycritere(TabSurface, Nagg, TabSurface, 0.0, compteur, 1.0);             //Calcule la surface moyenne moyenné sur l'ensemble des agrégats
        Tvmoy = Calculmoycritere(TabTv, Nagg, TabTv, 0.0, compteur, 1.0);                         //Calcule le taux de recouvrement volumique moyen moyenné sur l'ensemble des agrégats
        covmoy = Calculmoycritere(Tabcov, Nagg, Tabcov, 0.0, compteur, 1.0);                      //Calcule le taux de recouvrement surfacique moyen moyenné sur l'ensemble des agrégats
        SurfsurVolmoy = Surfmoy/PI/Dpmoy2; //Calculmoycritere(TabSurfsurVol, Nagg, TabSurfsurVol, 0.0, compteur, 1.0); //Calcule le rapport surface/volume moyen moyenné sur l'ensemble des agrégats
        //+++++++++++++++
        Calculdensimonomeres(nTab,Nagg,N); //Calcule la fonction de densité de nombre de sphérules

        fprintf(f1, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", temps, (double)Nagg, invNagg, lambdasurdist, Knudsen, Npmoy, Rgmoy*2/Dpmoy, Rgeomoy*2/Dpmoy ,Dmmoy/Dpmoy, Ncmoy, Volmoy, Surfmoy, Tvmoy, covmoy, SurfsurVolmoy);

        //On génère un fichier contenant les densités de nombre de sphérules primaires (loi self-preserving)
        for (j = 1; j <= N; j++)
            if (nTab[j] > 0.0)
                fprintf(f2,"%e\t%d\t%e\t%d\t%e\t%e\n",temps,Nagg,Npmoy,j,(double)j/Npmoy,nTab[j]*pow(Npmoy,2.0)/((double)N));

        //Analyse sur l'ensemble des agrégats à chaque instant
        int testeffacement=1;
        for(j = 1; j <= Nagg; j++, q++)
        {
            DpAggmoy = Calculmoycritere(TabDp, N, TabLabel, (double)j, compteur, 1.0); //Calcule le Dp moyen de l'agrégat considéré
            TT1.DgsurDp = TabRg[j]*2/DpAggmoy;
            TT1.NumFichier = i-1;
            TT1.NumAgg = j;
            TT1.Np = TabNp[j];
            TT1.DmsurDp = TabDm[j]/DpAggmoy;
            TT1.DgeosurDp = TabRgeo[j]*2/DpAggmoy;
            TT1.Dpmoy = DpAggmoy;
            TT1.Dpmoy3 = pow(Calculmoycritere(TabDp, N, TabLabel, (double)j, compteur, 3.0), 1./3.); //Calcule le Dp moyen d'ordre 3 de l'agrégat considéré
            TT1.Nc = TabNc[j];
            TT1.Volume = TabVolume[j]; //Calcule le volume de l'agrégat considéré
            TT1.Surface = TabSurface[j]; //Calcule la surface de l'agrégat considéré
            TT1.Tv = TabTv[j]; //Calcule le recouvrement volumique de l'agrégat considéré
            TT1.cov = Tabcov[j];
            //TT1.Ts=TabTs[j]; //Calcule le recouvrement surfacique de l'agrégat considéré
            TT1.SurfsurVol = TabSurfsurVol[j]; //Calcule le rapport surface/volume de l'agrégat considéré

            TT1.Npe = TabVolume[j]/(PI*pow(TT1.Dpmoy3, 3.0)/6.0);

            //On recherche dans la liste si l'agrégat a déjà été vu
            it=mylist.begin();

            while (it!=mylist.end() && (*it).DgsurDp<TT1.DgsurDp)       ++it;
            if (mylist.size()>0)        crit=(int)(((*it).DgsurDp-TT1.DgsurDp)*100);
            if (it==mylist.end())       {mylist.push_back (TT1);testeffacement=0;}
            else
            {
                if (crit>0)       {mylist.insert (it,TT1);testeffacement=0;}
            }
        }
        if (testeffacement==1 && ui->radioButtoneffacmement->isChecked())  EffacerFichier(i-1);
    }
        fclose(f1);
        fclose(f2);

        //Analyse Fractale
        sprintf(NomComplet,"%s/NumFichier_NumAgg_Np_Npe_Nc_Dg-Dp_Dm-Dp_Dgeo-Dp_Dpmoy_Dpmoy3_Vol_Surf_Tv_vraicov_SurfsurVol.dat",CheminSauve);
        f3=fopen(NomComplet,"w");
        fprintf(f3,"NumFichier\tNumAgg\tNp\tNpe\tNc\tDg-Dp\tDm-Dp\tDgeo-Dp\tDpmoy\tDpmoy3\tVol\tSurf\tTv\tvraicov\tSurfsurVol\n");
        sprintf(commentaires,"Génération du fichier : NumFichier_NumAgg_Np_Npe_Nc_Dg-Dp_Dm-Dp_Dgeo-Dp_Dpmoy_Dpmoy3_Vol_Surf_Tv_vraicov_SurfsurVol.dat\n");
        ui->AfficheurRep->append(commentaires);

        for (it=mylist.begin(); it!=mylist.end(); it++)
        {
                fprintf(f3,"%d\t%d\t%d\t%e\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",(*it).NumFichier,(*it).NumAgg,(*it).Np,(*it).Npe,(*it).Nc,(*it).DgsurDp,(*it).DmsurDp,(*it).DgeosurDp,(*it).Dpmoy,(*it).Dpmoy3,(*it).Volume,(*it).Surface,(*it).Tv,(*it).cov,(*it).SurfsurVol);
        }
        fclose(f3);

        printf("Nb de cas enregistres : %d sur %d\n",mylist.size(),q);
        sprintf(commentaires,"Nb de cas enregistres : %d sur %d\n",mylist.size(),q);
        ui->AfficheurRep->append(commentaires);
}

void Fermeture()
{
    delete[] compteur;
    delete[] TabLabel;
    delete[] TabDp;
    delete[] Tabposx;
    delete[] Tabposy;
    delete[] Tabposz;
    delete[] TabRg;
    delete[] TabNp;
    delete[] TabNc;
    delete[] TabDm;
    delete[] Tablpm;
    delete[] Tabdeltat;
    delete[] TabRgeo;
    delete[] TabXg;
    delete[] TabYg;
    delete[] TabZg;
    delete[] TabVolume;
    delete[] TabSurface;
    delete[] TabTv;
    //delete[] TabTs;
    delete[] Tabcov;
    delete[] TabSurfsurVol;
    delete[] nTab;

}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->pushButton,SIGNAL(clicked()),this,SLOT(BoutonQuitter()));
    connect(ui->BoutonRep,SIGNAL(clicked()),this,SLOT(BoutonRechercheParam()));
    connect(ui->BoutonExecAnalyse,SIGNAL(clicked()),this,SLOT(ExecuterAnalyse()));
    connect(ui->pushButton_2,SIGNAL(clicked()),this,SLOT(LecturePasAnalyse()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::BoutonQuitter()
{
    this->close();
}

void MainWindow::BoutonRechercheParam()
{
    ui->progressBar->setVisible(false);
    FichierParam = QFileDialog::getOpenFileName(this, "Sélectionner le fichier de données DLCA",
                                                      "C:/Users/dlca/Desktop/DLCA_sous_Qt/",
                                                      "Fichier de paramètres (*.par *.txt *.dat)");
    QFileInfo tmp2 = FichierParam;
    pathParam = tmp2.absolutePath(); //Cette variable ne retient que le chemin du fichier param

    //Affiche le chemin dans la zone de texte
    ui->AfficheurRep->append(FichierParam);
    LectureParam();
    ui->AfficheurRep->append(commentaires);
    ui->label_2->setEnabled(true);
    ui->spinBoxPasLecture->setEnabled(true);
    ui->pushButton_2->setEnabled(true);
    ui->AfficheurRep->append(CheminDATA);
    ui->radioButtoneffacmement->setEnabled(true);
}

void MainWindow::LecturePasAnalyse()
{
    PasAnalyse=ui->spinBoxPasLecture->value();
    sprintf(commentaires,"Choix du pas d'analyse : %d  \n",PasAnalyse);
    ui->AfficheurRep->append(commentaires);
    ui->BoutonExecAnalyse->setEnabled(true);
}

void MainWindow::ExecuterAnalyse()
{
    ui->progressBar->setMaximum(100);
    ui->progressBar->setVisible(true);
    ui->progressBar->setValue(0);
    ProgAnalyse();
    Fermeture();
    ui->progressBar->setVisible(false);
    sprintf(commentaires,"Fin de l'analyse ... Quatre fichiers ont été générés\n");
    ui->AfficheurRep->append(commentaires);
}
