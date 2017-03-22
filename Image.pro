#-------------------------------------------------
#
# Project created by QtCreator 2015-12-28T16:55:42
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Test
TEMPLATE = app


SOURCES += main.cpp\
        MainWindow.cpp \
    MyView.cpp \
    Image.cpp

HEADERS  += MainWindow.h \
    MyView.h \
    Image.h

FORMS    += MainWindow.ui

RESOURCES += \
    resources.qrc
