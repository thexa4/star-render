#include <iostream>
#include <fstream>
#include <cfloat>
#include <OpenEXR/OpenEXRConfig.h>
#include <OpenEXR/ImfHeader.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfEnvmapAttribute.h>
#include <math.h>
#include <cmath>
#include <string>
#include "ciedata.h"
#include "star_record.hpp"
#include "temperature_lookup.hpp"

struct COORD {
	float x;
	float y;
	float z;

	COORD() : x(0), y(0), z(0) {}
	COORD(float x, float y, float z): x(x), y(y), z(z) {}

	float dot(COORD &other) {
		return x * other.x + y * other.y + z * other.z;
	}

	COORD cross(COORD &other) {
		COORD result;
		result.x = y * other.z - z * other.y;
		result.y = z * other.x - x * other.z;
		result.z = x * other.y - y * other.x;
		return result;
	}

	COORD normalized() {
		if (x == 0 && y == 0 && z == 0)
			return *this;
		float length = sqrt(x * x + y * y + z * z);
		COORD result;
		result.x = x / length;
		result.y = y / length;
		result.z = z / length;
		return result;
	}
};

struct EMPTY {
};

template<typename T, typename META>
struct PIXEL {
	T x;
	T y;
	T z;

	META meta;

	template<typename B>
	void add(PIXEL<T, B> &other) {
		x += other.x;
		y += other.y;
		z += other.z;
	}
};

template<typename T, typename META, int WIDTH>
struct SCANLINE {
	PIXEL<T, META> pixels[WIDTH];
};

template<typename T, typename META, int WIDTH, int HEIGHT>
struct IMAGE {
	SCANLINE<T, META, WIDTH> lines[HEIGHT];
};

template<int WIDTH, int HEIGHT>
void slice_downscale(IMAGE<float, EMPTY, WIDTH, HEIGHT> &dest, float maximum) {

	for (int y = 0; y < HEIGHT; y++) {
		for (int x = 0; x < WIDTH; x++) {
			auto &pixel = dest.lines[y].pixels[x];
			pixel.x /= maximum;
			pixel.y /= maximum;
			pixel.z /= maximum;
			//if (pixel.y > 0)
			//	std::cout << "Color: " << pixel.x << " " << pixel.y << " " << pixel.z << std::endl;
		}
	}
}

template<typename META, int WIDTH, int HEIGHT>
double slice_colormap(IMAGE<double, META, WIDTH, HEIGHT> &source, IMAGE<float, EMPTY, WIDTH, HEIGHT> &dest) {
	double maximum = 0.0000001;
	for (int y = 0; y < HEIGHT; y++) {
		for (int x = 0; x < WIDTH; x++) {
			auto &source_pixel = source.lines[y].pixels[x];
			// Unsure what color space
			//double r = source_pixel.x * 2.36461385 + source_pixel.y * -0.89654057 + source_pixel.z * 0.46807328;
			//double g = source_pixel.x * 0.51516621 + source_pixel.y * 1.4264081 + source_pixel.z * 0.0887581;
			//double b = source_pixel.x * 0.0052037 + source_pixel.y * -0.01440816 + source_pixel.z * 1.00920446;
			// Best RGB - D50 http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
			//double r = source_pixel.x * 1.7552599 + source_pixel.y * -0.4836786 + source_pixel.z * -0.2530000;
			//double g = source_pixel.x *  -0.5441336 + source_pixel.y * 1.5068789 + source_pixel.z * 0.0215528;
			//double b = source_pixel.x *  0.0063467 + source_pixel.y * -0.0175761 + source_pixel.z * 1.2256959;
			// Wide Gamut RGB - D50 http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
			double r = source_pixel.x * 1.4628067 + source_pixel.y * -0.1840623 + source_pixel.z * -0.2743606;
			double g = source_pixel.x * -0.5217933 + source_pixel.y * 1.4472381 + source_pixel.z * 0.0677227;
			double b = source_pixel.x * 0.0349342 + source_pixel.y * -0.0968930 + source_pixel.z * 1.2884099;
			r = source_pixel.x;
			g = source_pixel.y;
			b = source_pixel.z;
			dest.lines[y].pixels[x].x = (float)std::max(0.0, r);
			dest.lines[y].pixels[x].y = (float)std::max(0.0, g);
			dest.lines[y].pixels[x].z = (float)std::max(0.0, b);

			maximum = std::max(maximum, r);
			maximum = std::max(maximum, g);
			maximum = std::max(maximum, b);
		}
	}
	return maximum;
}

template<int WIDTH, int HEIGHT>
void write_image_slice(const char *filename, IMAGE<float, EMPTY, WIDTH, HEIGHT> &image) {
	Imf_2_5::Header header(WIDTH, HEIGHT);
        header.channels().insert("R", Imf_2_5::Channel(Imf_2_5::FLOAT));
        header.channels().insert("G", Imf_2_5::Channel(Imf_2_5::FLOAT));
        header.channels().insert("B", Imf_2_5::Channel(Imf_2_5::FLOAT));
        header.compression() = Imf_2_5::ZIP_COMPRESSION;

	Imf_2_5::OutputFile file(filename, header);
	Imf_2_5::FrameBuffer framebuffer;
	framebuffer.insert("R", Imf_2_5::Slice(Imf_2_5::FLOAT, (char *)&image.lines[0].pixels[0].x, sizeof(image.lines[0].pixels[0]), sizeof(image.lines[0])));
	framebuffer.insert("G", Imf_2_5::Slice(Imf_2_5::FLOAT, (char *)&image.lines[0].pixels[0].y, sizeof(image.lines[0].pixels[0]), sizeof(image.lines[0])));
	framebuffer.insert("B", Imf_2_5::Slice(Imf_2_5::FLOAT, (char *)&image.lines[0].pixels[0].z, sizeof(image.lines[0].pixels[0]), sizeof(image.lines[0])));

	file.setFrameBuffer(framebuffer);
	file.writePixels(HEIGHT);
}

template<typename META, int WIDTH, int HEIGHT>
void render_image(IMAGE<double, META, WIDTH, HEIGHT> &compute_buffer, IMAGE<float, EMPTY, WIDTH, HEIGHT> &file_buffer) {
	//std::cout << "Normalizing" << std::endl;

	double maximum = 0.00001;
	maximum = std::max(maximum, slice_colormap(compute_buffer, file_buffer));

	std::cout << "Downscaling" << std::endl;
	slice_downscale(file_buffer, (float)maximum);

	std::cout << "Writing result";

	write_image_slice("canvas.tmp", file_buffer);
	rename("canvas.tmp", "canvas.exr");

	std::cout << "done" << std::endl;
}
const int TEMPERATURE_MAX = 50000;
const int TEMPERATURE_SPECTRAL_SAMPLES = 470;
const double TEMPERATURE_SPECTRUM_START = 380.0;
const double TEMPERATURE_SPECTRUM_END = 830.0;
PIXEL<double, EMPTY> temperatures[TEMPERATURE_MAX];
void precompute_temperatures() {
	for (int temperature = 0; temperature < TEMPERATURE_MAX; temperature++) {
		if (temperature < 470) {
			temperatures[temperature].x = 1;
			temperatures[temperature].y = 1;
			temperatures[temperature].z = 1;
			continue;
		}

		temperatures[temperature].x = 0;
		temperatures[temperature].y = 0;
		temperatures[temperature].z = 0;

		for (int i = 0; i < TEMPERATURE_SPECTRAL_SAMPLES; i++) {
			double percentage = ((double)i) / TEMPERATURE_SPECTRAL_SAMPLES;
			double wavelength = TEMPERATURE_SPECTRUM_START + (TEMPERATURE_SPECTRUM_END - TEMPERATURE_SPECTRUM_START) * percentage;
			double ePart = exp(6.62607015 * 2.99792458 / (wavelength * 1.380649 * temperature) * 1.0e6);
			double nominator = 2.0 * 6.62607015 * 2.99792458 * 2.99792458;
			double denominator = (wavelength * wavelength) * (wavelength * wavelength) * wavelength * (ePart - 1.0);
			double spectral_intensity = nominator * 1.0e27 / denominator;

			temperatures[temperature].x += CIE_X[i] * spectral_intensity;
			temperatures[temperature].y += CIE_Y[i] * spectral_intensity;
			temperatures[temperature].z += CIE_Z[i] * spectral_intensity;
		}

		double scale = temperatures[temperature].y;
		temperatures[temperature].x /= scale;
		temperatures[temperature].y /= scale;
		temperatures[temperature].z /= scale;
		//std::cout << "Temperature: " << temperature << "K: " << temperatures[temperature].x << " " << temperatures[temperature].y << " " << temperatures[temperature].z << std::endl;
	}
}
int lookup_temperature(double x, double y) {
	if (std::isnan(x) && std::isnan(y))
		return 3200;

	const int x_size = sizeof(temperature.x_lookup) / sizeof(temperature.x_lookup[0]);
	const int y_size = sizeof(temperature.y_lookup) / sizeof(temperature.y_lookup[0]);

	if (!std::isnan(x))
		x = (x - temperature.start_x) / (temperature.end_x - temperature.start_x) * x_size;
	if (!std::isnan(y))
		y = (y - temperature.start_y) / (temperature.end_y - temperature.start_y) * y_size;

	double y_result = nan("");
	if (!std::isnan(y)) {
		int left = std::max(0.0, std::min(y_size - 2.0, y));
		double frac = std::max(0.0, std::min(1.0, y - left));

		y_result = temperature.y_lookup[left] * (1 - frac) + temperature.y_lookup[left + 1] * frac;
	}
	double x_result = nan("");
	if (!std::isnan(x)) {
		int left = std::max(0.0, std::min(x_size - 2.0, x));
		double frac = std::max(0.0, std::min(1.0, x - left));

		x_result = temperature.x_lookup[left] * (1 - frac) + temperature.x_lookup[left + 1] * frac;
	}

	if (std::isnan(x_result))
		return (int)std::round(y_result);
	if (std::isnan(y_result))
		return (int)std::round(x_result);

	return (int)std::round((y_result + x_result) / 2);
}

struct Metadata {
	short max_sightings = 0;
	int num_samples = 0;
	long sum_sightings = 0;
};

long samples_done = 0;
//18260 × 6449
#define HEIGHT 6449
#define WIDTH 18260
#define BATCH_SIZE 4
IMAGE<double, Metadata, WIDTH, HEIGHT> &compute_buffer = *new IMAGE<double, Metadata, WIDTH, HEIGHT>{};
IMAGE<float, EMPTY, WIDTH, HEIGHT> &file_buffer = *new IMAGE<float, EMPTY, WIDTH, HEIGHT>{};


int main() {
	precompute_temperatures();

	Star buffer[BATCH_SIZE];
	size_t count = fread(&buffer, sizeof(buffer[0]), BATCH_SIZE, stdin);

	COORD fwd(-0.2, -0.6, -0.5);
	fwd = fwd.normalized();
	//COORD desired_up(-0.7, -0.25, 0.6);
	//COORD desired_up(0.196116, 0.588348, -0.784465);
	//COORD desired_up(-0.702864, 0.578829, -0.413449);
	COORD desired_up(-0.5, 0.6, 0);
	desired_up = desired_up.normalized();
	float zoom = 1.5;

	COORD right = fwd.cross(desired_up).normalized();
	COORD up = right.cross(fwd).normalized();

	std::cout << "up: " << up.x << "," << up.y << "," << up.z << std::endl;

	float up_scale = (float)HEIGHT / (float)WIDTH;
	right.x = right.x * zoom;
	right.y = right.y * zoom;
	right.z = right.z * zoom;
	up.x = up.x * zoom / up_scale;
	up.y = up.y * zoom / up_scale;
	up.z = up.z * zoom / up_scale;

	long hits = 0;
	long total_stars = 0;

	while(count > 0) {
		for (int i = 0; i < count; i++) {
			Star &s = buffer[i];

			if (std::isnan(s.teff_gspphot)) {
				s.teff_gspphot = (double)(lookup_temperature(s.phot_rp_mean_mag / s.phot_g_mean_mag, s.phot_bp_mean_mag / s.phot_g_mean_mag));
			}

			if (std::isnan(s.phot_g_mean_mag)) {
				double intensity = 0;
				int count = 0;
				if (!std::isnan(s.phot_bp_mean_mag)) {
					intensity += s.phot_bp_mean_mag;
					count++;
				}
				if (!std::isnan(s.phot_rp_mean_mag)) {
					intensity += s.phot_rp_mean_mag;
					count++;
				}
				if (count == 0)
					continue;
				s.phot_g_mean_mag = intensity / count;
			}

			float lat = s.dec * 3.14159265358979323846 / 180.0;
			float lng = s.ra * 3.14159265358979323846 / 180.0;

			float cos_lat = cos(lat);
			COORD coord (
				cos_lat * cos(lng),
				sin(lat),
				cos_lat * sin(lng)
			);

			float depth = coord.dot(fwd);
			float coord_u = coord.dot(up);
			float coord_v = coord.dot(right);

			total_stars++;
			if (depth < 0)
				continue;

			//std::cout << "off: " << offset_x << ";" << offset_y << std::endl;

			//std::cout << "c: " << coord_u << ";" << coord_v << std::endl;
			if (coord_u < -1.0 || coord_u > 1.0 || coord_v < -1.0 || coord_v > 1.0)
				continue;
			hits++;

			int pos_y = (int)round((coord_u + 1.0) * HEIGHT / 2.0);
			int pos_x = (int)round((coord_v + 1.0) * WIDTH / 2.0);
			if (pos_x < 0 || pos_x >= WIDTH)
				continue;
			if (pos_y < 0 || pos_y >= HEIGHT)
				continue;

			auto &dest = compute_buffer.lines[pos_y].pixels[pos_x];

			int temperature = (int)round(s.teff_gspphot);
			if (temperature < 0 || temperature > 50000) {
				std::cerr << "Bad temperature: " << temperature << std::endl;
				temperature = 50000 - 1;
				//continue;
			}
			PIXEL<double, EMPTY> star_color = temperatures[temperature];

			double intensity = pow(2.51188643151, s.phot_g_mean_mag);

			star_color.x *= intensity;
			star_color.y *= intensity;
			star_color.z *= intensity;
			dest.add(star_color);

			dest.meta.max_sightings = std::max(dest.meta.max_sightings, s.astrometric_n_obs_al);
			dest.meta.num_samples++;
			dest.meta.sum_sightings += s.astrometric_n_obs_al;

			samples_done++;
			//if ((samples_done % 500000000) == 0) {
			//	render_image(compute_buffer, file_buffer);
			//}
		}
		count = fread(&buffer, sizeof(buffer[0]), BATCH_SIZE, stdin);
	}
	std::cout << "Percentage in view: " << round(((float)hits / (float)total_stars) * 100) << "%" << std::endl;

	const int array_len = 1024;
	double sighting_intensity[array_len];
	long sighting_count[array_len];
	double intensity_sum[array_len];
	int intensity_count[array_len];
	for (int i = 0; i < array_len; i++) {
		sighting_intensity[i] = 0;
		sighting_count[i] = 0;
		intensity_sum[i] = 0;
		intensity_count[i] = 0;
	}

	render_image(compute_buffer, file_buffer);

	return 0;
}
