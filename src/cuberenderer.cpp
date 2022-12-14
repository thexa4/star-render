#include <iostream>
#include <fstream>
#include <cfloat>
#include <OpenEXR/OpenEXRConfig.h>
#include <OpenEXR/ImfHeader.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfEnvmapAttribute.h>
#include <math.h>
#include <string>
#include "ciedata.h"
#include "star_record.hpp"
#include "temperature_lookup.hpp"

struct COORD {
	float x;
	float y;
	float z;

	float dot(COORD &other) {
		return x * other.x + y * other.y + z * other.z;
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

template<typename T, typename META, int WIDTH, int HEIGHT>
struct CUBEMAP {
	IMAGE<T, META, WIDTH, HEIGHT> left;
	IMAGE<T, META, WIDTH, HEIGHT> front;
	IMAGE<T, META, WIDTH, HEIGHT> right;
	IMAGE<T, META, WIDTH, HEIGHT> back;
	IMAGE<T, META, WIDTH, HEIGHT> top;
	IMAGE<T, META, WIDTH, HEIGHT> bottom;

	IMAGE<T, META, WIDTH, HEIGHT> * get_face(int i) {
		if (i < 0 || i > 5)
			return NULL;
		IMAGE<T, META, WIDTH, HEIGHT> *faces[] = {&top, &bottom, &left, &right, &front, &back};
		return faces[i];
	}

	PIXEL<T, META>& lookup(float lat, float lng) {
		lat *= 3.14159265358979323846 / 180.0;
		lng *= 3.14159265358979323846 / 180.0;

		COORD lookup_coord;
		float cos_lat = cos(lat);
		lookup_coord.x = cos_lat * cos(lng);
		lookup_coord.y = sin(lat);
		lookup_coord.z = cos_lat * sin(lng);

		IMAGE<T, META, WIDTH, HEIGHT> *faces[] = {&top, &bottom, &left, &right, &front, &back};
		COORD face_coords[] = {
			(COORD){0, 0, 1},
			(COORD){0, 0, -1},
			(COORD){-1, 0, 0},
			(COORD){1, 0, 0},
			(COORD){0, -1, 0},
			(COORD){0, 1, 0},
		};
		COORD face_left_coord[] = {
			(COORD){0, 1, 0},
			(COORD){0, 1, 0},
			(COORD){0, 1, 0},
			(COORD){0, 1, 0},
			(COORD){1, 0, 0},
			(COORD){1, 0, 0},
		};
		COORD face_up_coord[] = {
			(COORD){1, 0, 0},
			(COORD){1, 0, 0},
			(COORD){0, 0, 1},
			(COORD){0, 0, 1},
			(COORD){0, 0, 1},
			(COORD){0, 0, 1},
		};

		int best_index = -1;
		float best_score = -2;
		for (int i = 0; i < 6; i++) {
			float score = lookup_coord.dot(face_coords[i]);
			if (score > best_score) {
				best_index = i;
				best_score = score;
			}
		}
		auto face = faces[best_index];

		float left = lookup_coord.dot(face_left_coord[best_index]) / best_score;
		float up = lookup_coord.dot(face_up_coord[best_index]) / best_score;
		int x = (int)floor(((left + 1.0) / 2.0) * WIDTH);
		int y = (int)floor(((up + 1.0) / 2.0) * HEIGHT);

		/*if (x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT || best_index == 5) {
			std::cout << "lookup: " << ((left + 1.0) / 2.0) << ", " << ((up + 1.0) / 2.0) << " @" << best_index << std::endl;
			std::cout << "input: " << lookup_coord.dot(face_left_coord[best_index]) << ", " << lookup_coord.dot(face_up_coord[best_index]) << ", " << best_score << std::endl;
			std::cout << "left: " << left << ", up: " << up << std::endl;
		}*/

		return face->lines[y].pixels[x];
	}
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
	file.writePixels(WIDTH);
}

template<typename META, int WIDTH, int HEIGHT>
void render_image(CUBEMAP<double, META, WIDTH, HEIGHT> &compute_buffer, CUBEMAP<float, EMPTY, WIDTH, HEIGHT> &file_buffer) {
	//std::cout << "Normalizing" << std::endl;

	double maximum = 0.00001;
	maximum = std::max(maximum, slice_colormap(compute_buffer.top, file_buffer.top));
	maximum = std::max(maximum, slice_colormap(compute_buffer.bottom, file_buffer.bottom));
	maximum = std::max(maximum, slice_colormap(compute_buffer.left, file_buffer.left));
	maximum = std::max(maximum, slice_colormap(compute_buffer.front, file_buffer.front));
	maximum = std::max(maximum, slice_colormap(compute_buffer.right, file_buffer.right));
	maximum = std::max(maximum, slice_colormap(compute_buffer.back, file_buffer.back));

	slice_downscale(file_buffer.top, (float)maximum);
	slice_downscale(file_buffer.bottom, (float)maximum);
	slice_downscale(file_buffer.left, (float)maximum);
	slice_downscale(file_buffer.front, (float)maximum);
	slice_downscale(file_buffer.right, (float)maximum);
	slice_downscale(file_buffer.back, (float)maximum);

	//std::cout << "Writing result";

	write_image_slice("skycube-top.tmp", file_buffer.top);
	//std::cout << '.';
	write_image_slice("skycube-bottom.tmp", file_buffer.bottom);
	//std::cout << '.';
	write_image_slice("skycube-left.tmp", file_buffer.left);
	//std::cout << '.';
	write_image_slice("skycube-front.tmp", file_buffer.front);
	//std::cout << '.';
	write_image_slice("skycube-right.tmp", file_buffer.right);
	//std::cout << '.';
	write_image_slice("skycube-back.tmp", file_buffer.back);
	rename("skycube-top.tmp", "skycube-top.exr");
	rename("skycube-bottom.tmp", "skycube-bottom.exr");
	rename("skycube-left.tmp", "skycube-left.exr");
	rename("skycube-front.tmp", "skycube-front.exr");
	rename("skycube-right.tmp", "skycube-right.exr");
	rename("skycube-back.tmp", "skycube-back.exr");

	//std::cout << "done" << std::endl;
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
#define OUTPUT_SIZE 2048
#define BATCH_SIZE 4
CUBEMAP<double, Metadata, OUTPUT_SIZE, OUTPUT_SIZE> compute_buffer;
CUBEMAP<float, EMPTY, OUTPUT_SIZE, OUTPUT_SIZE> file_buffer;

int main() {
	precompute_temperatures();

	Star buffer[BATCH_SIZE];
	size_t count = fread(&buffer, sizeof(buffer[0]), BATCH_SIZE, stdin);
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

			auto &dest = compute_buffer.lookup(s.dec, s.ra);

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
			if ((samples_done % 500000000) == 0) {
				render_image(compute_buffer, file_buffer);
			}
		}
		count = fread(&buffer, sizeof(buffer[0]), BATCH_SIZE, stdin);
	}

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

	for (int face = 0; face < 6; face++) {
		IMAGE<double, Metadata, 2048, 2048> &img = *compute_buffer.get_face(face);
		for (int y = 0; y < OUTPUT_SIZE; y++) {
			for (int x = 0; x < OUTPUT_SIZE; x++) {
				auto &p = img.lines[y].pixels[x];

				double avg = (p.x + p.y + p.z) / 3.0;
				if (p.meta.max_sightings < array_len && p.meta.max_sightings >= 0) {
					sighting_intensity[p.meta.max_sightings] += avg;
					sighting_count[p.meta.max_sightings]++;
				}

				p.x = avg / 100000000;
				p.y = std::min((short int)250, p.meta.max_sightings);
				p.z = 0;

				if (p.meta.num_samples > 0) {
					int avg_sightings = (int)std::round((double)p.meta.sum_sightings / (double)p.meta.num_samples);
					if (avg_sightings >= 0 && avg_sightings < array_len) {
						intensity_sum[avg_sightings] += avg;
						intensity_count[avg_sightings]++;
					}
					p.z = std::min(135, avg_sightings);
				}
			}
		}
	}
	render_image(compute_buffer, file_buffer);


	for (int i = 0; i < array_len; i++) {
		int max_intensity = 0;
		if (sighting_count[i] > 0)
			max_intensity = (int)(sighting_intensity[i] / sighting_count[i] / 1000000);

		int avg_intensity = 0;
		if (intensity_count[i] > 0)
			avg_intensity = (int)(intensity_sum[i] / intensity_count[i] / 1000000);

		std::cout << max_intensity << ", " << avg_intensity << std::endl;
	}

	return 0;
}
