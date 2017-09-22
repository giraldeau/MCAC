#ifdef WITH_GUI
#include <QApplication>
#endif
#include "mainwindow.hpp"
#include <cstdio>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
    if(argc>1)
    {
        return DLCA::No_GUI(argc, argv);
    }
#ifdef WITH_GUI
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
#else
    cout << "Missing argument : param file" << endl;
    return 1;
#endif

}
