#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "physical_model.h"
#include "statistics.h"

#include <vector>
#include <array>

#ifdef WITH_GUI

#include <QMainWindow>
namespace DLCA{

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
}
#else
#ifdef WITH_QT
#include <QMainWindow>

namespace DLCA{

class STUB : public QMainWindow
{
    Q_OBJECT;
};
}
#endif
#endif

namespace DLCA{

int No_GUI(int argc, char *argv[]);
void SauveASCII(int value, int id);
void InitRandom();
double Random();
void CalculDistance(int id, double &distmin, int &aggcontact);
int LectureSuiviTempo();
int rechercheValTab();
int Init(PhysicalModel&,Statistics&, ListAggregat&);
void Fermeture();
bool locale_with_dots();
double latof(const char* _char);
PhysicalModel LectureParams(const std::string& FichierParam);
void Calcul(PhysicalModel&);
void print(std::string str);
int dirExists(const char *path);
std::string extractPath(const std::string& file);

}

#endif // MAINWINDOW_H
