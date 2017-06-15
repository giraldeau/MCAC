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
    void BoutonQuitter();
    void ExecuterAnalyse();
    void ProgAnalyse();
    void GestionBC();
    void GestionEDM();
    void GestionDF3();
    void GestionAutoco();
};

#endif // MAINWINDOW_H
