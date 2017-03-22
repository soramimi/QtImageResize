#include "MainWindow.h"
#include "MyView.h"
#include <QPainter>
#include <QTime>
#include "Image.h"

MyView::MyView(QWidget *parent) : QWidget(parent)
{
	image = QImage(":/image/checker2.png");
	image.toPixelFormat(QImage::Format_ARGB32);
	QPixmap pm(":/image/transparent.png");
	background_brush = QBrush(pm);
}

void MyView::paintEvent(QPaintEvent *)
{
	QTime time;
	time.start();

	QPainter pr(this);
	int w = width();
	int h = height();

	pr.fillRect(0, 0, w, h, background_brush);

	QImage newimg = resizeImage(image, w, h, EnlargeMethod::Bicubic, true);
//	QImage newimg = image.scaled(w, h, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);

	int ms = time.elapsed();
	char tmp[100];
	sprintf(tmp, "%dms", ms);
	the_mainwindow->setLabelText(tmp);

	pr.drawImage(0, 0, newimg, 0, 0, w, h);
}

