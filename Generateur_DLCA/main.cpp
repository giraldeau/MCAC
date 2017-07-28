#ifdef WITH_GUI
#include <QApplication>
#endif
#include "mainwindow.h"
#include <stdio.h>
#include <iostream>

int main(int argc, char *argv[])
{
    if(argc>1)
    {
        return No_GUI(argc, argv);
    }
    else
    {
#ifdef WITH_GUI
        QApplication a(argc, argv);
        MainWindow w;
        w.show();

        return a.exec();
#else
        cout << "Missing argument : param file" << endl;
#endif

    }

    return 1;
}
