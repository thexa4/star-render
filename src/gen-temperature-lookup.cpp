#include <iostream>
#include <fstream>
#include <cfloat>
#include <math.h>
#include <string>
#include "ciedata.h"
#include "star_record.hpp"

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

#define BATCH_SIZE 16
#define TEMPERATURE_X 96
#define TEMPERATURE_Y 96
IMAGE<float, EMPTY, TEMPERATURE_X, TEMPERATURE_Y> temperature_buffer;
int main() {
	double start_x = (480.0 / 512.0);
	double end_x = (700.0 / 512.0);
	double start_y = (400.0 / 512.0);
	double end_y = (530.0 / 512.0);

	Star buffer[BATCH_SIZE];
	size_t count = fread(&buffer, sizeof(buffer[0]), BATCH_SIZE, stdin);
	while(count > 0) {
		for (int i = 0; i < count; i++) {
			Star &s = buffer[i];

			if (!std::isnan(s.teff_gspphot)) {
				if (!std::isnan(s.phot_rp_mean_mag) && !std::isnan(s.phot_bp_mean_mag) && !std::isnan(s.phot_g_mean_mag)) {

					double x = (s.phot_bp_mean_mag / s.phot_g_mean_mag - start_x) / (end_x - start_x) * TEMPERATURE_X;
					x = std::max(0.0, std::min((double)TEMPERATURE_X, x));
					double y = (s.phot_rp_mean_mag / s.phot_g_mean_mag - start_y) / (end_y - start_y) * TEMPERATURE_Y;
					y = std::max(0.0, std::min((double)TEMPERATURE_Y, y));

					//std::cout << "x " << x << ", y " << y << std::endl;
					auto &pixel = temperature_buffer.lines[(int)y].pixels[(int)x];
					pixel.x += s.teff_gspphot;
					pixel.y += 1;
				}
			}
		}
		count = fread(&buffer, sizeof(buffer[0]), BATCH_SIZE, stdin);
	}

	for(int y = 0; y < TEMPERATURE_Y; y++) {
		for(int x = 0; x < TEMPERATURE_X; x++) {
			auto &pixel = temperature_buffer.lines[y].pixels[x];

			pixel.x /= pixel.y * 10000.0;
			pixel.y = pixel.x;
			pixel.z = pixel.x;
		}
	}
	double x_lookup[TEMPERATURE_X];
	int x_counts[TEMPERATURE_X];
	for (int i = 0; i < TEMPERATURE_X; i++) {
		x_lookup[i] = 0;
		x_counts[i] = 0;
	}
	double y_lookup[TEMPERATURE_Y];
	int y_counts[TEMPERATURE_X];
	for (int i = 0; i < TEMPERATURE_Y; i++) {
		y_lookup[i] = 0;
		y_counts[i] = 0;
	}

	for(int y = 0; y < TEMPERATURE_Y; y++) {
		for(int x = 0; x < TEMPERATURE_X; x++) {
			auto &pixel = temperature_buffer.lines[y].pixels[x];

			if (pixel.x > 0) {
				x_lookup[x] += pixel.x;
				x_counts[x]++;
				y_lookup[y] += pixel.y;
				y_counts[y]++;
				continue;
			}

			double left = nan("");
			for (int _x = x - 1; _x > 0; _x--) {
				auto &_pixel = temperature_buffer.lines[y].pixels[_x];

				if (_pixel.x > 0) {
					left = _pixel.x;
					break;
				}
			}
			double right = nan("");
			for (int _x = x + 1; _x < TEMPERATURE_X; _x++) {
				auto &_pixel = temperature_buffer.lines[y].pixels[_x];
				if (_pixel.x > 0) {
					right = _pixel.x;
					break;
				}
			}
			double top = nan("");
			for (int _y = y - 1; _y > 0; _y--) {
				auto &_pixel = temperature_buffer.lines[_y].pixels[x];
				if (_pixel.x > 0) {
					top = _pixel.x;
					break;
				}
			}
			double bottom = nan("");
			for (int _y = y + 1; _y < TEMPERATURE_Y; _y++) {
				auto &_pixel = temperature_buffer.lines[_y].pixels[x];
				if (_pixel.x > 0) {
					bottom = _pixel.x;
					break;
				}
			}

			double horiz = 0;
			double vert = 0;
			int horiz_count = 0;
			int vert_count = 0;
			if (!std::isnan(left)) {
				horiz += left;
				horiz_count++;
			}
			if (!std::isnan(right)) {
				horiz += right;
				horiz_count++;
			}
			if (!std::isnan(top)) {
				vert += top;
				vert_count++;
			}
			if (!std::isnan(bottom)) {
				vert += bottom;
				vert_count++;
			}

			double temperature = 0;
			int temperature_count = 0;
			if (horiz_count > 0) {
				temperature += horiz / horiz_count;
				temperature_count++;
			}
			if (vert_count > 0) {
				temperature += vert / vert_count;
				temperature_count++;
			}
			if (temperature_count > 0) {
				pixel.x = temperature / temperature_count;
				pixel.y = pixel.x;
				pixel.z = pixel.x;
			} else {
				pixel.x = 1;
				pixel.y = 0;
				pixel.z = 0;
			}
		}
	}

	std::cout << "struct Temperature {" << std::endl;
	std::cout << "\tunsigned short lookup[" << TEMPERATURE_Y << "][" << TEMPERATURE_X << "] = {" << std::endl;
	for(int y = 0; y < TEMPERATURE_Y; y++) {
		bool first = true;
		std::cout << "\t\t{ ";
		for(int x = 0; x < TEMPERATURE_X; x++) {
			if (!first)
				std::cout << ", ";
			else
				first = false;
			auto &pixel = temperature_buffer.lines[y].pixels[x];
			std::cout << std::round(pixel.z * 10000);
		}
		std::cout << "}," << std::endl;
	}
	std::cout << "\t};" << std::endl;
	std::cout << "\tdouble start_x = " << start_x << ";" << std::endl;
	std::cout << "\tdouble end_x = " << end_x << ";" << std::endl;
	std::cout << "\tdouble start_y = " << start_y << ";" << std::endl;
	std::cout << "\tdouble end_y = " << end_y << ";" << std::endl;
	std::cout << "\tdouble x_lookup[" << TEMPERATURE_X << "] = {";
	bool x_first = true;
	for (int i = 0; i < TEMPERATURE_X; i++)
	{
		if (!x_first)
			std::cout << ", ";
		else
			x_first = false;
		std::cout << x_lookup[i] / x_counts[i];
	}
	std::cout << "};" << std::endl;
	std::cout << "\tdouble y_lookup[" << TEMPERATURE_Y << "] = {";
	bool y_first = true;
	for (int i = 0; i < TEMPERATURE_Y; i++)
	{
		if (!y_first)
			std::cout << ", ";
		else
			y_first = false;
		std::cout << y_lookup[i] / y_counts[i];
	}
	std::cout << "};" << std::endl;
	std::cout << "} temperature;" << std::endl;

	return 0;
}
