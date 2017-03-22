#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QLabel>
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

	void setLabelText(QString const &text);
private:
	Ui::MainWindow *ui;
	QLabel *label;
};

extern MainWindow *the_mainwindow;

#endif // MAINWINDOW_H
