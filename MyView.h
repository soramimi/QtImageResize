#ifndef MYVIEW_H
#define MYVIEW_H

#include <QWidget>
#include <QImage>

class MyView : public QWidget
{
	Q_OBJECT
private:
	QImage image;
	QBrush background_brush;
public:
	explicit MyView(QWidget *parent = 0);

signals:

public slots:

	// QWidget interface
protected:
	void paintEvent(QPaintEvent *);
};

#endif // MYVIEW_H
