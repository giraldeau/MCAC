#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "physical_model.hpp"
#include "aggregatList.hpp"
#include "statistics.hpp"

#include <array>
#include <vector>
#include <experimental/filesystem>


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
void Init(PhysicalModel&,StatisticStorage&, ListAggregat&);
void Fermeture();
bool locale_with_dots();
double latof(const char* _char);
PhysicalModel LectureParams(const std::string& FichierParam);
void Calcul(PhysicalModel&);
void print(std::string str);
bool dirExists(const char *path);
std::experimental::filesystem::path extractPath(const std::string& file);

}  // namespace DLCA


#endif // MAINWINDOW_H
