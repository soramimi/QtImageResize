#include "MainWindow.h"
#include "ui_MainWindow.h"

MainWindow *the_mainwindow;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    the_mainwindow = this;
	ui->setupUi(this);
	setWindowFlags(Qt::Window | Qt::CustomizeWindowHint);

	label = new QLabel(this);
	label->setGeometry(10, 10, 50, 20);
	label->setText("Hello, world");
	label->setStyleSheet("QLabel { background: #ffffff; }");
}

MainWindow::~MainWindow()
{
	delete ui;
}

void MainWindow::setLabelText(QString const &text)
{
	label->setText(text);
}
