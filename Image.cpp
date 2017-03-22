#include "Image.h"
#include <QImage>
#include <math.h>

static double bicubic(double t)
{
	if (t < 0) t = -t;
	double tt = t * t;
	double ttt = t * t * t;
	const double a = -0.5;
	if (t < 1) return (a + 2) * ttt - (a + 3) * tt + 1;
	if (t < 2) return a * ttt - 5 * a * tt + 8 * a * t - 4 * a;
	return 0;
}

class FloatRGB {
public:
	float r;
	float g;
	float b;
	FloatRGB()
		: r(0)
		, g(0)
		, b(0)
	{
	}
	FloatRGB(float r, float g, float b)
		: r(r)
		, g(g)
		, b(b)
	{
	}
	FloatRGB(QRgb const &src)
		: r(qRed(src))
		, g(qGreen(src))
		, b(qBlue(src))
	{
	}
	FloatRGB operator + (FloatRGB const &right) const
	{
		return FloatRGB(r + right.r, g + right.g, b + right.b);
	}
	FloatRGB operator * (float t) const
	{
		return FloatRGB(r * t, g * t, b * t);
	}
	void operator += (FloatRGB const &o)
	{
		r += o.r;
		g += o.g;
		b += o.b;
	}
	void add(FloatRGB const &p, float v)
	{
		r += p.r * v;
		g += p.g * v;
		b += p.b * v;
	}

	void operator *= (float t)
	{
		r *= t;
		g *= t;
		b *= t;
	}
	int r8() const
	{
		if (r <= 0) return 0;
		if (r >= 255) return 255;
		return (int)r;
	}
	int g8() const
	{
		if (g < 0) return 0;
		if (g > 255) return 255;
		return (int)g;
	}
	int b8() const
	{
		if (b < 0) return 0;
		if (b > 255) return 255;
		return (int)b;
	}
	QRgb toQPixel(float amount) const
	{
		if (amount == 1) {
			return qRgb(r8(), g8(), b8());
		} else if (amount == 0) {
			return qRgb(0, 0, 0);
		}
		float m = 1 / amount;
		FloatRGB p = *this * m;
		return qRgb(p.r8(), p.g8(), p.b8());
	}
	QRgb toQRgb() const
	{
		return qRgb(r8(), g8(), b8());
	}
};

class FloatRGBA {
public:
	float r;
	float g;
	float b;
	float a;
	FloatRGBA()
		: r(0)
		, g(0)
		, b(0)
		, a(0)
	{
	}
	FloatRGBA(float r, float g, float b,  float a = 255)
		: r(r)
		, g(g)
		, b(b)
		, a(a)
	{
	}
	FloatRGBA(QRgb const &src)
		: r(qRed(src))
		, g(qGreen(src))
		, b(qBlue(src))
		, a(qAlpha(src))
	{
	}
	FloatRGBA operator + (FloatRGBA const &right) const
	{
		return FloatRGBA(r + right.r, g + right.g, b + right.b);
	}
	FloatRGBA operator * (float t) const
	{
		return FloatRGBA(r * t, g * t, b * t);
	}
	void operator += (FloatRGBA const &o)
	{
		r += o.r;
		g += o.g;
		b += o.b;
	}
	void operator *= (float t)
	{
		r *= t;
		g *= t;
		b *= t;
	}
	void add(FloatRGBA const &p, float v)
	{
		v *= p.a;
		a += v;
		v /= 255;
		r += p.r * v;
		g += p.g * v;
		b += p.b * v;
	}
	int r8() const
	{
		if (r <= 0) return 0;
		if (r >= 255) return 255;
		return (int)r;
	}
	int g8() const
	{
		if (g < 0) return 0;
		if (g > 255) return 255;
		return (int)g;
	}
	int b8() const
	{
		if (b < 0) return 0;
		if (b > 255) return 255;
		return (int)b;
	}
	int a8() const
	{
		if (a < 0) return 0;
		if (a > 255) return 255;
		return (int)a;
	}
	QRgb toQPixel(float amount) const
	{
		if (amount == 0) {
			return qRgba(0, 0, 0, 0);
		}
		FloatRGBA pixel(*this);
		pixel *= (255.0f / pixel.a);
		pixel.a = pixel.a / amount;
		if (pixel.a < 0) pixel.a = 0;
		else if (pixel.a > 255) pixel.a = 255;
		return qRgba(pixel.r8(), pixel.g8(), pixel.b8(), pixel.a8());
	}
	QRgb toQRgba(float amount) const
	{
		return toQPixel(amount);
	}
};

QImage resizeNearestNeighbor(QImage const &image, int dst_w, int dst_h)
{
	const int src_w = image.width();
	const int src_h = image.height();
	QImage newimg(dst_w, dst_h, QImage::Format_ARGB32);
	for (int y = 0; y < dst_h; y++) {
		float fy = (float)y * src_h / dst_h;
		QRgb *dst = (QRgb *)newimg.scanLine(y);
		QRgb const *src = (QRgb const *)image.scanLine((int)fy);
		float mul = (float)src_w / dst_w;
		for (int x = 0; x < dst_w; x++) {
			float fx = (float)x * mul;
			dst[x] = src[(int)fx];
		}
	}
	return std::move(newimg);
}

template <typename PIXEL>
QImage resizeAveragingT(QImage const &image, int dst_w, int dst_h)
{
	const int src_w = image.width();
	const int src_h = image.height();
	QImage newimg(dst_w, dst_h, QImage::Format_ARGB32);
	for (int y = 0; y < dst_h; y++) {
		float lo_y = (float)y * src_h / dst_h;
		float hi_y = (float)(y + 1) * src_h / dst_h;
		QRgb *dst = (QRgb *)newimg.scanLine(y);
		float mul = (float)src_w / dst_w;
		for (int x = 0; x < dst_w; x++) {
			float lo_x = (float)x * mul;
			float hi_x = (float)(x + 1) * mul;
			int lo_iy = (int)lo_y;
			int hi_iy = (int)hi_y;
			int lo_ix = (int)lo_x;
			int hi_ix = (int)hi_x;
			PIXEL pixel;
			float volume = 0;
			for (int sy = lo_iy; sy <= hi_iy; sy++) {
				float vy = 1;
				if (sy < src_h) {
					if (lo_iy == hi_iy) {
						vy = hi_y - lo_y;
					} else if (sy == lo_iy) {
						vy = 1 - (lo_y - sy);
					} else if (sy == hi_iy) {
						vy = hi_y - sy;
					}
				}
				QRgb const *src = (QRgb const *)image.scanLine(sy < src_h ? sy : (src_h - 1));
				for (int sx = lo_ix; sx <= hi_ix; sx++) {
					PIXEL p = src[sx < src_w ? sx : (src_w - 1)];
					float vx = 1;
					if (sx < src_w) {
						if (lo_ix == hi_ix) {
							vx = hi_x - lo_x;
						} else if (sx == lo_ix) {
							vx = 1 - (lo_x - sx);
						} else if (sx == hi_ix) {
							vx = hi_x - sx;
						}
					}
					float v = vy * vx;
					pixel.add(p, v);
					volume += v;
				}
			}
			dst[x] = pixel.toQPixel(volume);
		}
	}
	return std::move(newimg);
}

template <typename PIXEL>
QImage resizeAveragingHT(QImage const &image, int dst_w)
{
	const int src_w = image.width();
	const int src_h = image.height();
	QImage newimg(dst_w, src_h, QImage::Format_ARGB32);
	for (int y = 0; y < src_h; y++) {
		QRgb *dst = (QRgb *)newimg.scanLine(y);
		for (int x = 0; x < dst_w; x++) {
			float lo_x = (float)x * src_w / dst_w;
			float hi_x = (float)(x + 1) * src_w / dst_w;
			int lo_ix = (int)lo_x;
			int hi_ix = (int)hi_x;
			PIXEL pixel;
			float volume = 0;
			QRgb const *src = (QRgb const *)image.scanLine(y < src_h ? y : (src_h - 1));
			for (int sx = lo_ix; sx <= hi_ix; sx++) {
				PIXEL p = src[sx < src_w ? sx : (src_w - 1)];
				float v = 1;
				if (sx < src_w) {
					if (lo_ix == hi_ix) {
						v = hi_x - lo_x;
					} else if (sx == lo_ix) {
						v = 1 - (lo_x - sx);
					} else if (sx == hi_ix) {
						v = hi_x - sx;
					}
				}
				pixel.add(p, v);
				volume += v;
			}
			dst[x] = pixel.toQPixel(volume);
		}
	}
	return std::move(newimg);
}

template <typename PIXEL>
QImage resizeAveragingVT(QImage const &image, int dst_h)
{
	const int src_w = image.width();
	const int src_h = image.height();
	QImage newimg(src_w, dst_h, QImage::Format_ARGB32);
	for (int y = 0; y < dst_h; y++) {
		float lo_y = (float)y * src_h / dst_h;
		float hi_y = (float)(y + 1) * src_h / dst_h;
		QRgb *dst = (QRgb *)newimg.scanLine(y);
		for (int x = 0; x < src_w; x++) {
			int lo_iy = (int)lo_y;
			int hi_iy = (int)hi_y;
			PIXEL pixel;
			float volume = 0;
			for (int sy = lo_iy; sy <= hi_iy; sy++) {
				float v = 1;
				if (sy < src_h) {
					if (lo_iy == hi_iy) {
						v = hi_y - lo_y;
					} else if (sy == lo_iy) {
						v = 1 - (lo_y - sy);
					} else if (sy == hi_iy) {
						v = hi_y - sy;
					}
				}
				QRgb const *src = (QRgb const *)image.scanLine(sy < src_h ? sy : (src_h - 1));
				PIXEL p = src[x];
				pixel.add(p, v);
				volume += v;
			}
			dst[x] = pixel.toQPixel(volume);
		}
	}
	return std::move(newimg);
}

struct bilinear_t {
	int i0, i1;
	float v0, v1;
};

template <typename PIXEL>
QImage resizeBilinearT(QImage const &image, int dst_w, int dst_h)
{
	const int src_w = image.width();
	const int src_h = image.height();
	QImage newimg(dst_w, dst_h, QImage::Format_ARGB32);

	std::vector<bilinear_t> lut(dst_w);
	bilinear_t *lut_p = &lut[0];
	for (int x = 0; x < dst_w; x++) {
		float tx = (float)x * src_w / dst_w - 0.5;
		int x0, x1;
		if (tx < 0) {
			x0 = x1 = 0;
			tx = 0;
		} else {
			x0 = x1 = (int)tx;
			if (x0 + 1 < src_w) {
				x1 = x0 + 1;
				tx -= x0;
			} else {
				x0 = x1 = src_w - 1;
				tx = 0;
			}
		}
		lut_p[x].i0 = x0;
		lut_p[x].i1 = x1;
		lut_p[x].v1 = tx;
		lut_p[x].v0 = 1 - tx;
	}

	for (int y = 0; y < dst_h; y++) {
		float yt = (float)y * src_h / dst_h - 0.5;
		int y0, y1;
		if (yt < 0) {
			y0 = y1 = 0;
			yt = 0;
		} else {
			y0 = y1 = (int)yt;
			if (y0 + 1 < src_h) {
				y1 = y0 + 1;
				yt -= y0;
			} else {
				y0 = y1 = src_h - 1;
				yt = 0;
			}
		}
		float ys = 1 - yt;
		QRgb *dst = (QRgb *)newimg.scanLine(y);
		QRgb const *src1 = (QRgb const *)image.scanLine(y0);
		QRgb const *src2 = (QRgb const *)image.scanLine(y1);
		for (int x = 0; x < dst_w; x++) {
			float a11 = lut_p[x].v0 * ys;
			float a12 = lut_p[x].v1 * ys;
			float a21 = lut_p[x].v0 * yt;
			float a22 = lut_p[x].v1 * yt;
			PIXEL pixel;
			pixel.add(PIXEL(src1[lut_p[x].i0]), a11);
			pixel.add(PIXEL(src1[lut_p[x].i1]), a12);
			pixel.add(PIXEL(src2[lut_p[x].i0]), a21);
			pixel.add(PIXEL(src2[lut_p[x].i1]), a22);
			dst[x] = pixel.toQPixel(a11 + a12 + a21 + a22);
		}
	}
	return std::move(newimg);
}

template <typename PIXEL>
QImage resizeBilinearHT(QImage const &image, int dst_w)
{
	const int src_w = image.width();
	const int src_h = image.height();
	QImage newimg(dst_w, src_h, QImage::Format_ARGB32);
	for (int y = 0; y < src_h; y++) {
		QRgb *dst = (QRgb *)newimg.scanLine(y);
		QRgb const *src = (QRgb const *)image.scanLine(y);
		float mul = (float)src_w / dst_w;
		for (int x = 0; x < dst_w; x++) {
			float xt = (float)x * mul - 0.5;
			int x0, x1;
			if (xt < 0) {
				x0 = x1 = 0;
				xt = 0;
			} else {
				x0 = x1 = (int)xt;
				if (x0 + 1 < src_w) {
					x1 = x0 + 1;
					xt -= x0;
				} else {
					x0 = x1 = src_w - 1;
					xt = 0;
				}
			}
			float xs = 1 - xt;
			PIXEL p1(src[x0]);
			PIXEL p2(src[x1]);
			PIXEL p;
			p.add(p1, xs);
			p.add(p2, xt);
			dst[x] = p.toQPixel(1);
		}
	}
	return std::move(newimg);
}

template <typename PIXEL>
QImage resizeBilinearVT(QImage const &image, int dst_h)
{
	const int src_w = image.width();
	const int src_h = image.height();
	QImage newimg(src_w, dst_h, QImage::Format_ARGB32);
	for (int y = 0; y < dst_h; y++) {
		float yt = (float)y * src_h / dst_h - 0.5;
		int y0, y1;
		if (yt < 0) {
			y0 = y1 = 0;
			yt = 0;
		} else {
			y0 = y1 = (int)yt;
			if (y0 + 1 < src_h) {
				y1 = y0 + 1;
				yt -= y0;
			} else {
				y0 = y1 = src_h - 1;
				yt = 0;
			}
		}
		float ys = 1 - yt;
		QRgb *dst = (QRgb *)newimg.scanLine(y);
		QRgb const *src1 = (QRgb const *)image.scanLine(y0);
		QRgb const *src2 = (QRgb const *)image.scanLine(y1);
		for (int x = 0; x < src_w; x++) {
			PIXEL p1(src1[x]);
			PIXEL p2(src2[x]);
			PIXEL p;
			p.add(p1, ys);
			p.add(p2, yt);
			dst[x] = p.toQPixel(1);
		}
	}
	return std::move(newimg);
}

typedef float (*bicubic_lut_t)[4];

static bicubic_lut_t makeBicubicLookupTable(int src, int dst, std::vector<float> *out)
{
	out->resize(dst * 4);
	float (*lut)[4] = (float (*)[4])&(*out)[0];
	for (int x = 0; x < dst; x++) {
		float sx = (float)x * src / dst - 0.5;
		int ix = (int)floor(sx);
		float tx = sx - ix;
		for (int x2 = -1; x2 <= 2; x2++) {
			int x3 = ix + x2;
			if (x3 >= 0 && x3 < src) {
				lut[x][x2 + 1] = bicubic(x2 - tx);
			}
		}
	}
	return lut;
}

template <typename PIXEL>
QImage resizeBicubicT(QImage const &image, int dst_w, int dst_h)
{
	const int src_w = image.width();
	const int src_h = image.height();
	QImage newimg(dst_w, dst_h, QImage::Format_ARGB32);

	std::vector<float> bicubic_lut_x;
	std::vector<float> bicubic_lut_y;
	bicubic_lut_t bicubic_lut_x_p = makeBicubicLookupTable(src_w, dst_w, &bicubic_lut_x);
	bicubic_lut_t bicubic_lut_y_p = makeBicubicLookupTable(src_h, dst_h, &bicubic_lut_y);

	for (int y = 0; y < dst_h; y++) {
		QRgb *dst = (QRgb *)newimg.scanLine(y);
		float sy = (float)y * src_h / dst_h - 0.5;
		int iy = (int)floor(sy);
		for (int x = 0; x < dst_w; x++) {
			float sx = (float)x * src_w / dst_w - 0.5;
			int ix = (int)floor(sx);
			PIXEL pixel;
			float amount = 0;
			for (int y2 = -1; y2 <= 2; y2++) {
				int y3 = iy + y2;
				if (y3 >= 0 && y3 < src_h) {
					float vy = bicubic_lut_y_p[y][y2 + 1];
					QRgb const *src = (QRgb const *)image.scanLine(y3);
					for (int x2 = -1; x2 <= 2; x2++) {
						int x3 = ix + x2;
						if (x3 >= 0 && x3 < src_w) {
							float vx = bicubic_lut_x_p[x][x2 + 1];
							PIXEL p = src[x3];
							float v = vx * vy;
							pixel.add(p, v);
							amount += v;
						}
					}
				}
			}
			dst[x] = pixel.toQPixel(amount);
		}
	}
	return std::move(newimg);
}

template <typename PIXEL>
QImage resizeBicubicHT(QImage const &image, int dst_w)
{
	const int src_w = image.width();
	const int src_h = image.height();
	QImage newimg(dst_w, src_h, QImage::Format_ARGB32);

	std::vector<float> bicubic_lut_x;
	bicubic_lut_t bicubic_lut_x_p = makeBicubicLookupTable(src_w, dst_w, &bicubic_lut_x);

	for (int y = 0; y < src_h; y++) {
		QRgb *dst = (QRgb *)newimg.scanLine(y);
		QRgb const *src = (QRgb const *)image.scanLine(y);
		for (int x = 0; x < dst_w; x++) {
			PIXEL pixel;
			float volume = 0;
			float sx = (float)x * src_w / dst_w - 0.5;
			int ix = (int)floor(sx);
			for (int x2 = -1; x2 <= 2; x2++) {
				int x3 = ix + x2;
				if (x3 >= 0 && x3 < src_w) {
					float v = bicubic_lut_x_p[x][x2 + 1];
					PIXEL p = src[x3];
					pixel.add(p, v);
					volume += v;
				}
			}
			dst[x] = pixel.toQPixel(volume);
		}
	}
	return std::move(newimg);
}

template <typename PIXEL>
QImage resizeBicubicVT(QImage const &image, int dst_h)
{
	const int src_w = image.width();
	const int src_h = image.height();
	QImage newimg(src_w, dst_h, QImage::Format_ARGB32);

	std::vector<float> bicubic_lut_y;
	bicubic_lut_t bicubic_lut_y_p = makeBicubicLookupTable(src_h, dst_h, &bicubic_lut_y);

	for (int x = 0; x < src_w; x++) {
		for (int y = 0; y < dst_h; y++) {
			QRgb *dst = (QRgb *)newimg.scanLine(y) + x;
			PIXEL pixel;
			float volume = 0;
			float sy = (float)y * src_h / dst_h - 0.5;
			int iy = (int)floor(sy);
			for (int y2 = -1; y2 <= 2; y2++) {
				int y3 = iy + y2;
				if (y3 >= 0 && y3 < src_h) {
					float v = bicubic_lut_y_p[y][y2 + 1];
					PIXEL p = ((QRgb const *)image.scanLine(y3))[x];
					pixel.add(p, v);
					volume += v;
				}
			}
			*dst = pixel.toQPixel(volume);
		}
	}
	return std::move(newimg);
}

QImage resizeImage(QImage image, int dst_w, int dst_h, EnlargeMethod method, bool alphachannel)
{
	if (dst_w > 0 && dst_h > 0) {
		int w, h;
		w = image.width();
		h = image.height();
		if (w != dst_w || h != dst_h) {
			if (dst_w < w || dst_h < h) {
				if (dst_w < w && dst_h < h) {
					if (alphachannel) {
						image = resizeAveragingT<FloatRGBA>(image, dst_w, dst_h);
					} else {
						image = resizeAveragingT<FloatRGB>(image, dst_w, dst_h);
					}
				} else if (dst_w < w) {
					if (alphachannel) {
						image = resizeAveragingHT<FloatRGBA>(image, dst_w);
					} else {
						image = resizeAveragingHT<FloatRGB>(image, dst_w);
					}
				} else if (dst_h < h) {
					if (alphachannel) {
						image = resizeAveragingVT<FloatRGBA>(image, dst_h);
					} else {
						image = resizeAveragingVT<FloatRGB>(image, dst_h);
					}
				}
			}
			w = image.width();
			h = image.height();
			if (dst_w > w || dst_h > h) {
				if (method == EnlargeMethod::Bilinear) {
					if (dst_w > w && dst_h > h) {
						if (alphachannel) {
							image = resizeBilinearT<FloatRGBA>(image, dst_w, dst_h);
						} else {
							image = resizeBilinearT<FloatRGB>(image, dst_w, dst_h);
						}
					} else if (dst_w > w) {
						if (alphachannel) {
							image = resizeBilinearHT<FloatRGBA>(image, dst_w);
						} else {
							image = resizeBilinearHT<FloatRGB>(image, dst_w);
						}
					} else if (dst_h > h) {
						if (alphachannel) {
							image = resizeBilinearVT<FloatRGBA>(image, dst_h);
						} else {
							image = resizeBilinearVT<FloatRGB>(image, dst_h);
						}
					}
				} else if (method == EnlargeMethod::Bicubic) {
					if (dst_w > w && dst_h > h) {
						if (alphachannel) {
							image = resizeBicubicT<FloatRGBA>(image, dst_w, dst_h);
						} else {
							image = resizeBicubicT<FloatRGB>(image, dst_w, dst_h);
						}
					} else if (dst_w > w) {
						if (alphachannel) {
							image = resizeBicubicHT<FloatRGBA>(image, dst_w);
						} else {
							image = resizeBicubicHT<FloatRGB>(image, dst_w);
						}
					} else if (dst_h > h) {
						if (alphachannel) {
							image = resizeBicubicVT<FloatRGBA>(image, dst_h);
						} else {
							image = resizeBicubicVT<FloatRGB>(image, dst_h);
						}
					}
				} else {
					image = resizeNearestNeighbor(image, dst_w, dst_h);
				}
			}
		}
		return std::move(image);
	}
	return QImage();
}
