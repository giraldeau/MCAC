#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "physical_model.h"
#include "aggregatList.h"

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
void CalculDistance(int id, double &distmin, int &aggcontact);
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
