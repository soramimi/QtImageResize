#ifndef IMAGE_H
#define IMAGE_H

class QImage;

enum class EnlargeMethod {
	Nearest,
	Bilinear,
	Bicubic,
};

QImage resizeImage(QImage image, int dst_w, int dst_h, EnlargeMethod method = EnlargeMethod::Bilinear, bool alphachannel = true);

#endif // IMAGE_H
