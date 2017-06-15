#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;

private slots:
    void BoutonRechercheParam();
    void ModulePhysique();
    void BoutonQuitter();
    void ExecuterDLCA();
    void Calcul();
    void BoutonRechercheSuiviTempo();
};

#endif // MAINWINDOW_H
