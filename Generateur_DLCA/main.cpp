#include <QApplication>
#include "mainwindow.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
    if(argc>1)
    {
        return No_GUI(argc, argv);
    }
    else
    {
        QApplication a(argc, argv);
        MainWindow w;
        w.show();

        return a.exec();
    }

    return 1;
}
