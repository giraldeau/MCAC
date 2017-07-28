#ifndef MAINWINDOW_H
#define MAINWINDOW_H


#include <vector>
#include <array>

using namespace std;

#ifdef WITH_GUI

#include <QMainWindow>


namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    void print(char* str);
    void progress(int value);


private:
    Ui::MainWindow *ui;

private slots:
    void BoutonRechercheParam();
    void ModulePhysique();
    void BoutonQuitter();
    void ExecuterDLCA();
    void BoutonRechercheSuiviTempo();
};
#else
#ifdef WITH_QT
#include <QMainWindow>

class STUB : public QMainWindow
{
    Q_OBJECT;
};
#endif
#endif


int No_GUI(int argc, char *argv[]);
void SauveASCII(int value, int id);
void InitRandom();
double Random();
double Maxi2D(int colonne, int nmax);
double MinEtIndex(double* tableau, int size, int& position);
void SupprimeLigne(int ligne);
void InsertionSort(int n, double arr[], int index[]);
void quickSort(double arr[], vector<int>& index, int left, int right);
void quickSort(int n, double arr[], vector<int>& index);
void MonTri(int n, double arr[], vector<int>& index);
int Probabilite(bool trier,double &deltatemps);
void CalculDistance(int id, double &distmin, int &aggcontact);
int Reunit(int AggI, int AggJ, int &err);
int LectureSuiviTempo();
int rechercheValTab();
void Init();
void Fermeture();
bool locale_with_dots();
double latof(const char* _char);
void LectureParams();
void Calcul();
void print(char* str);
int dirExists(const char *path);


#endif // MAINWINDOW_H
