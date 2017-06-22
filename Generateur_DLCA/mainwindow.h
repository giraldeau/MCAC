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

int No_GUI(int argc, char *argv[]);
#endif // MAINWINDOW_H
