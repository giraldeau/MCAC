#include "physical_model.hpp"
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/stat.h>


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW_2(a) ((a)*(a))
#define POW_3(a) ((a)*(a)*(a))

namespace fs = std::experimental::filesystem;
using namespace std;

namespace MCAC{

const double PI = atan(1.0)*4;

PhysicalModel::PhysicalModel() :
    Asurfgrowth(0),
    dfe(1.4),
    kfe(1.8),
    xsurfgrowth(2),
    coeffB(0),
    lambda(0),
    Dpeqmass(0),
    rpeqmass(0),
    gamma_(0),
    P(101300),
    T(1500),
    Mu(0),
    K(1.38066E-23),
    Rho(1.8e-3),
    Dpm(30),
    sigmaDpm(1.25),
    Time(0),
    X(0),
    FV(3e-3),
    L(DBL_MAX),
    precision(1e-5),
    FactorModelBeta(0),
    CPUStart(0),
    CPULimit(-1),
    PHY_Time_limit(-1),
    NPP_avg_limit(0),
    GridDiv(10),
    N(2500),
    AggMin(1),
    DeltaSauve(10),
    root_method(2),
    Mode(1),
    Wait(0),
    WaitLimit(-1),
    ActiveModulephysique(true),
    ActiveVariationTempo(false),
    use_verlet(true),
    toBeDestroyed(true),
    CheminSauve()
{}


PhysicalModel::PhysicalModel(const string& FichierParam) : PhysicalModel()
{
    FILE* f;
    char sauve[500];

    f = fopen(FichierParam.c_str(), "rt");

    if (!f)
    {
        cout << FichierParam.c_str() << " not found" << endl;
        exit(1);
    }
    char commentaires[500];
    char com[500]; // Char array used in the ASCII Save

    if( !fgets(com, 500, f))
    {
        cout << "I need the N parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        N=size_t(atoi(commentaires));
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the FV parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        FV=latof(commentaires)*1E-6;
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the Dpm parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        Dpm=latof(commentaires);
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the Mode parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        Mode=atoi(commentaires);
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the sigmaDpm parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        sigmaDpm=latof(commentaires);
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the ActiveModulephysique parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        ActiveModulephysique=atoi(commentaires);
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the T parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        T=latof(commentaires);
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the P parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        P=latof(commentaires);
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the Rho parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        Rho=latof(commentaires);
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the kfe parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        kfe=atof(commentaires);
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the dfe parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        dfe=atof(commentaires);
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the coeffB parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        coeffB=atof(commentaires);
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the xsurfgrowth parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        xsurfgrowth=atof(commentaires);
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the ActiveVariationTempo parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        ActiveVariationTempo=atoi(commentaires);
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the AggMin parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        AggMin=size_t(atoi(commentaires));
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the WaitLimit parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        WaitLimit=atoi(commentaires);
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the CPULimit parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        CPULimit=atoi(commentaires);
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the PHY_Time_limit parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        PHY_Time_limit=atof(commentaires);
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the NPP_avg_limit parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        NPP_avg_limit=size_t(atoi(commentaires));
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the DeltaSauve parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",commentaires,com);
        DeltaSauve=size_t(atoi(commentaires));
    }
    if( !fgets(com, 500, f))
    {
        cout << "I need the output_dir parameter" << endl;
        exit(1);
    }
    else
    {
        sscanf(com,"%s  %s",sauve,com);
    }
    fclose(f);

    if (Mode == 1)
    {
        X = pow(double(N)*PI/6.0/FV*(1.0+3.0*sigmaDpm*sigmaDpm/Dpm/Dpm),1.0/3.0); //Loi normale
    }
    else
    {
        X = pow(double(N)*PI/6.0/FV*exp(9.0/2.0*log(sigmaDpm)*log(sigmaDpm)),1.0/3.0); //Loi log-normale
    }

    fs::path pathParam = extractPath(FichierParam);
    CheminSauve = pathParam / sauve;

    if( ! fs::exists(CheminSauve))
    {
        if( ! fs::create_directory(CheminSauve))
        {
            cout << "Error creating directory " << CheminSauve << endl;
            exit(1);
        }
    }
    else
    {
        if( ! fs::is_directory(CheminSauve))
        {
            cout << "Error not a directory " << CheminSauve << endl;
            exit(1);
        }
    }
}


void PhysicalModel::Init()
{

    L = X*Dpm*1E-9;
    Time=0;

    K = 1.38066E-23;
    lambda = 66.5E-9*(101300/P)*(T/293.15)*(1+110/293.15)/(1+110/T);
    Dpeqmass = Dpm*exp(1.5*log(sigmaDpm)*log(sigmaDpm));        //Diamètre équivalent massique moyen des monomères
                                                                //donné par l'équation de Hatch-Choate
    rpeqmass = (Dpeqmass*1E-9)/2.0;                             //Rayon équivalent massique moyen des monomères
    gamma_ = 1.378*(0.5+0.5*erf(((lambda/rpeqmass)+4.454)/10.628));
    Mu = 18.203E-6*(293.15+110)/(T+110)*pow(T/293.15,1.5);

    Asurfgrowth = coeffB*1E-3;

    // Bool used to chose if the script will run a Verlet list, significantly reducing the cost of Calcul_Distance
    use_verlet = true;
    GridDiv = 10;               // Number of Divisions of the box

    AggMin = MAX(AggMin, 1);
    CPUStart = clock();
    SetPrecision(1e-5);
    UseSecante();

    print();

    toBeDestroyed = false;
}

bool PhysicalModel::Finished(const size_t Nagg, const size_t NPP_avg) const
{
    if (Nagg <= AggMin)
    {
        cout << "We reach the AggMin condition" << endl << endl;
        return true;
    }
    if (WaitLimit > 0 && Wait >= WaitLimit)
    {
        cout << "We reach the WaitLimit condition" << endl << endl;
        return true;
    }

    if (CPULimit > 0)
    {
        clock_t currentCPU = clock();
        double elapse = double(currentCPU - CPUStart) / CLOCKS_PER_SEC;
        if (elapse >= CPULimit)
        {
            cout << "We reach the CPULimit condition" << endl << endl;
            return true;
        }
    }

    if (PHY_Time_limit>0 && Time >= PHY_Time_limit)
    {
        cout << "We reach the Maximum physical time condition " << Time << "/" << PHY_Time_limit  << endl;
        return true;
    }

    if (NPP_avg >= NPP_avg_limit)
    {
        cout << "We reach the NPP_avg_limit condition " << NPP_avg << "/" << NPP_avg_limit << endl;
        return true;
    }

    return false;
}

void PhysicalModel::SetPrecision(const double _precision)
{
    precision = _precision;
}

void PhysicalModel::UseDichotomia()
{
    root_method = 0;
}

void PhysicalModel::UseBrent()
{
    root_method = 1;
}

void PhysicalModel::UseSecante()
{
    root_method = 2;
}

void PhysicalModel::print() const
{

    string RootMethod;
    switch(root_method)
    {
    case 0 :
        RootMethod = "Dichotomia";
        break;
    case 1:
        RootMethod = "Brent";
        break;
    default:
        cout << "Root method unknown, using secante" << endl;
       /* FALLTHRU */
    [[clang::fallthrough]]; case 2 :
        RootMethod = "Secante";
    }


    cout << "Physical parameters:"              << endl
         << " Initial Nagg : " << N             << endl
         << " Box size     : " << L             << endl
         << " X            : " << X             << endl
         << " Pressure     : " << P             << endl
         << " Temperature  : " << T             << endl
         << " diffusivity  : " << Mu            << endl
         << " K            : " << K             << endl
         << " FV           : " << FV            << endl
         << " density      : " << Rho           << endl
         << " Dpm          : " << Dpm           << endl
         << " sigmaDpm     : " << sigmaDpm      << endl
         << " Asurfgrowth  : " << Asurfgrowth   << endl
         << " xsurfgrowth  : " << xsurfgrowth   << endl
         << " coeffB       : " << coeffB        << endl
         << " dfe          : " << dfe           << endl
         << " kfe          : " << kfe           << endl
         << " lambda       : " << lambda        << endl
         << " Dpeqmass     : " << Dpeqmass      << endl
         << " rpeqmass     : " << rpeqmass      << endl
         << " gamma_       : " << gamma_        << endl
         << endl
         << "Options for Pysical model: "       << endl
         << " precision   : " << precision      << endl
         << " root method : " << RootMethod     << endl
         << " Mode : " << Mode                  << endl
         << endl
         << "Ending calcul when:"               << endl
         << " - There is " << AggMin << " aggregats left or less" << endl;
    if (WaitLimit > 0)
    {
        cout << " - It has been " << WaitLimit << " iterations without collision" <<endl;
    }
    if (CPULimit > 0)
    {
        cout << " - The simulations is running for more than " << CPULimit << " seconds" << endl;
    }
    if (PHY_Time_limit > 0)
    {
        cout << " - The residence time is larger than " << PHY_Time_limit << " seconds" << endl;
    }
    if (NPP_avg_limit > 0)
    {
        cout << " - The average Npp per aggregate reach " << NPP_avg_limit << " monomers" << endl;
    }
    cout << endl;
}


//#####################################################################################################################

 __attribute__((pure)) double PhysicalModel::Cunningham(const double R) const //Facteur correctif de Cunningham
{
    double A = 1.142;
    double B = 0.558;
    double C = 0.999;
    return 1.0+A*lambda/R+B*lambda/R*exp(-C*R/lambda);
}


//############################# Fonctions pour le calcul du diamètre de mobilité ################################

/*
 Fonction permettant de retrouver le rayon de mobilité en régime transitoire
 On obtient le bon rayon de mobilité lorsque la fonction retourne 0
*/

double PhysicalModel::ModeleBeta(const double rm) const
{
    if (rm<0.) {throw 42;}
    return Cunningham(rm) - FactorModelBeta*rm;
}

double PhysicalModel::Dichotomie(const double x0) const
{
    double 	rmin, rmax, rmed, frmed, frmin, frmax;
    rmin = x0/100;   //pow(np/1.5,1/1.8)*rp/40; //borne inférieure de rm
    rmax = 2*x0; //pow(np/1.5,1/1.8)*rp*40; //bornes de recherche de rm

    frmin = ModeleBeta(rmin);
    frmax = ModeleBeta(rmax);

    rmed = (rmin+rmax)/2 ;
    frmed = ModeleBeta(rmed);

    while (fabs(frmed)>precision)
    {
        if (frmin*frmax>=0)
        {
            cout << "Intervalle incorrect : " << frmin << " " << frmax << endl;
            return -1;
        }
        if (frmin*frmed < 0)
        {
            rmax = rmed;
            frmax = frmed;
        }
        else
        {
            rmin = rmed;
            frmin = frmed;
        }
        rmed = (rmin+rmax)/2 ;
        frmed = ModeleBeta(rmed);
    }
    return rmed;
}

double PhysicalModel::brentq(const double x0) const
{
    double xa,xb,xtol,rtol;
    int iter=500;

    xa = x0/100;    //pow(np/1.5,1/1.8)*rp/40; //borne inférieure de rm
    xb = 2*x0;      //pow(np/1.5,1/1.8)*rp*40; //bornes de recherche de rm
    xtol = rtol = precision;

    double xpre = xa, xcur = xb;
    double xblk = 0., fpre, fcur, fblk = 0., spre = 0., scur = 0.;
    double stry, dpre, dblk;
    int i;

    fpre = ModeleBeta(xpre);
    if (fabs(fpre) < precision) {
        return xpre;
    }
    fcur = ModeleBeta(xcur);
    if (fabs(fcur) < precision) {
        return xcur;
    }

    if (fpre*fcur > 0)
    {
        cout << "Intervalle incorrect : " << fpre << " " << fcur << endl;
        return -1;
    }

    for (i = 0; i < iter; i++) {
        if (fpre*fcur < 0) {
            xblk = xpre;
            fblk = fpre;
            spre = scur = xcur - xpre;
        }
        if (fabs(fblk) < fabs(fcur)) {
            xpre = xcur;
            xcur = xblk;
            xblk = xpre;

            fpre = fcur;
            fcur = fblk;
            fblk = fpre;
        }

        double delta = (xtol + rtol*fabs(xcur))/2;
        double sbis = (xblk - xcur)/2;
        //if (fcur == 0 || fabs(sbis) < delta) {
        if (fabs(fcur) < precision) {
            return xcur;
        }

        if (fabs(spre) > delta && fabs(fcur) < fabs(fpre)) {
            if (abs(xpre - xblk) <= 1e-16) {
                /* interpolate */
                stry = -fcur*(xcur - xpre)/(fcur - fpre);
            }
            else {
                /* extrapolate */
                dpre = (fpre - fcur)/(xpre - xcur);
                dblk = (fblk - fcur)/(xblk - xcur);
                stry = -fcur*(fblk*dblk - fpre*dpre)
                    /(dblk*dpre*(fblk - fpre));
            }
            if (2*fabs(stry) < MIN(fabs(spre), 3*fabs(sbis) - delta)) {
                /* good short step */
                spre = scur;
                scur = stry;
            } else {
                /* bisect */
                spre = sbis;
                scur = sbis;
            }
        }
        else {
            /* bisect */
            spre = sbis;
            scur = sbis;
        }

        xpre = xcur; fpre = fcur;
        if (fabs(scur) > delta) {
            xcur += scur;
        }
        else {
            xcur += (sbis > 0 ? delta : -delta);
        }

        fcur = ModeleBeta(xcur);
    }
    return xcur;
}

double PhysicalModel::Secante(const double x0) const
{
    double x[2];
    double fx[2];
    bool i=true;

    x[0]=x0+4.e-13;
    fx[0] = ModeleBeta(x[0]);

    x[1]=x0;

    try
    {
        while ((fabs(fx[!i])>precision))
        {
            fx[i]=ModeleBeta(x[i]);
            x[!i]=x[i]-fx[i]*(x[i]-x[!i])/(fx[i]-fx[!i]);
            i = !i;
        }
    }
    catch (int)
    {
        return brentq(x0);
    }
    return x[i];
}

__attribute((pure)) double PhysicalModel::ConvertRg2Dm(const size_t np, const double rg,const double rmoy)
{
    double start = 2*rmoy*pow(np,gamma_/dfe);
    return ConvertRg2DmFromStart(np,rg,start);
}

double PhysicalModel::ConvertRg2DmFromStart(const size_t np, const double rg, const double start)
{
    FactorModelBeta = pow(kfe,-1/dfe)*pow(np,(1-gamma_)/dfe)*Cunningham(rpeqmass)/rg;
/*
    double npeqmass = kfe*pow(rg/rpeqmass,dfe);
    double nptmp = pow(npeqmass*pow(np,gamma_-1),1./gamma_);
    FactorModelBeta = 1./(rpeqmass*pow(nptmp,gamma_/dfe)/Cunningham(rpeqmass));
*/
    switch (root_method)
    {
    case 2:
        return Secante(start*0.5)*2;
    case 1:
        return brentq(start*0.5)*2;
    default:
        return Dichotomie(start*0.5)*2;
    }
}

 __attribute__((pure)) double PhysicalModel::Grow(const double R,const double dt) const
{
    return R + Asurfgrowth*pow(R, xsurfgrowth-2)*dt;
}

 __attribute((pure)) double PhysicalModel::friction_coeff(const size_t npp) const
 {
     double cc = Cunningham(rpeqmass);
     return (6.0*PI*Mu*rpeqmass/cc)*pow(npp,gamma_/dfe); //to be expressed in terms of V_agg/V_pp
 }

 __attribute((pure)) double PhysicalModel::friction_coeff2(const double rgg) const
 {
     double rmm = rgg*1.3;
     double cc = Cunningham(rmm);
     return (6.0*PI*Mu*rmm/cc);
 }

 __attribute((pure)) double PhysicalModel::diffusivity(const double f_agg) const
 {
     return K*T/f_agg;
 }

 __attribute__((pure)) double PhysicalModel::velocity(const double masse) const
{
    return sqrt(8.0*K*T/PI/masse);
}

 __attribute__((pure)) double PhysicalModel::relax_time(const double masse, const double f_agg) const
{
    return masse/f_agg;
}

__attribute((const)) double inverfc(const double p)
{
        double x,t,pp;
        if (p >= 2.0){ return -100.;}
        if (p <= 0.0){ return 100.;}
        pp = (p < 1.0)? p : 2. - p;
        t = sqrt(-2.*log(pp/2.));
        x = -0.70711*((2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t);
        for (int j=0;j<2;j++) {
                double err = erfc(x) - pp;
                x += err/(1.12837916709551257*exp(-(x*x))-x*err);
        }
        return (p < 1.0? x : -x);
}

__attribute((const)) double inverf(const double p) {
    return inverfc(1.-p);
}

bool locale_with_dots()
{
    static bool tested = false;
    static bool with_dots;

    if (tested)
    {
        return with_dots;
    }
    double testfloat = 1.5;
    string teststr1 = "1.5";
    string teststr2 = "1,5";
    double test1=atof(teststr1.c_str());
    double test2=atof(teststr2.c_str());

    if (fabs(test1-testfloat)<1e-3)
    {
        with_dots = true;
    }
    else if (fabs(test2-testfloat)<1e-3)
    {
        with_dots = false;
    }
    else
    {
        cout << "What locale are you using ?" << endl;
        exit(1);
    }
    return with_dots;
}

double latof(const char* _char)
{
    string mystring = _char;
    if (!locale_with_dots())
    {
        size_t f = mystring.find('.');
        if (f>0)
        {
            mystring.replace(f, 1, ",");
        }
    }
    return atof(mystring.c_str());
}

fs::path extractPath(const string& filename)
{
    fs::path path = filename;
    fs::path fullpath = fs::absolute(path);

    if(! fs::exists(fullpath))
    {
        cout << "File does not exist\n" << endl;
        exit(1);
    }

    fs::path parentpath = fullpath.parent_path();

    return parentpath; //Cette variable ne retient que le chemin du fichier param
}

bool dirExists(const char *path)
{
    struct stat info{};

    if(stat( path, &info ) != 0)
    {
        return false;
    }
    if(info.st_mode & S_IFDIR)
    {
        return true;
    }
    return false;
}


tuple<bool,double,double,double> linreg(const vector<double>& x, const vector<double>& y)
{
    double   sumx = 0.0;                        /* sum of x                      */
    double   sumx2 = 0.0;                       /* sum of x**2                   */
    double   sumxy = 0.0;                       /* sum of x * y                  */
    double   sumy = 0.0;                        /* sum of y                      */
    double   sumy2 = 0.0;                       /* sum of y**2                   */

    size_t n = x.size();
    auto N = double(n);

   for (size_t i=0;i<n;i++)
   {
      double X(log(x[i]));
      double Y(log(y[i]));
      sumx  += X;
      sumx2 += POW_2(X);
      sumxy += X * Y;
      sumy  += Y;
      sumy2 += POW_2(Y);
   }



   double denom = (N * sumx2 - POW_2(sumx));
   if (n==0 || abs(denom) < 1e-9) {
       // singular matrix. can't solve the problem.
       tuple<bool,double,double,double> res{false,0.,0.,0.};
       return res;
   }

   double a = (N * sumxy  -  sumx * sumy) / denom;
   double b = (sumy * sumx2  -  sumx * sumxy) / denom;

   /* compute correlation coeff     */
   double r = (sumxy - sumx * sumy / N) /
              POW_2((sumx2 - POW_2(sumx) / N) *
                    (sumy2 - POW_2(sumy) / N));

   tuple<bool,double,double,double> res{true,a,b,r};
   return res;
}


}  // namespace MCAC

