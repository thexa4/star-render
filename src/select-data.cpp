#include "csv_parser.hpp"
#include "star_record.hpp"
#include <iostream>

int main() {
	const char *fields[] = {
		"ra",
		"dec",
		"astrometric_n_obs_al",
		"astrometric_n_obs_ac",
		"astrometric_n_good_obs_al",
		"astrometric_n_bad_obs_al",
		"visibility_periods_used",
		"phot_g_mean_mag",
		"phot_bp_mean_mag",
		"phot_rp_mean_mag",
		"distance_gspphot",
		"teff_gspphot",
		"hilbert_index"
	};
	int offsets[sizeof(fields) / sizeof(fields[0])];
	for (int i = 0; i < sizeof(offsets) / sizeof(offsets[0]); i++)
		offsets[i] = -1;
	
	aria::csv::CsvParser parser(std::cin);
	bool parsed_records = false;
	for (auto& row : parser) {
		if (!parsed_records) {
			parsed_records = true;
			for (int i = 0; i < sizeof(offsets) / sizeof(offsets[0]); i++) {
				int offset = 0;
				for (auto& field : row) {
					if (field == fields[i])
						offsets[i] = offset;
					offset++;
				}
				if (offsets[i] == -1) {
					std::cerr << "Missing field: " << fields[i] << std::endl;
					return 1;
				}
			}
		} else {
			Star result{};
			result.ra = std::stof(row[offsets[0]]);
			result.dec = std::stof(row[offsets[1]]);
			result.astrometric_n_obs_al = std::stoi(row[offsets[2]]);
			result.astrometric_n_obs_ac = std::stoi(row[offsets[3]]);
			result.astrometric_n_good_obs_al = std::stoi(row[offsets[4]]);
			result.astrometric_n_bad_obs_al = std::stoi(row[offsets[5]]);
			result.visibility_periods_used = std::stoi(row[offsets[6]]);

			if (row[offsets[7]] != "")
				result.phot_g_mean_mag = std::stof(row[offsets[7]]);
			if (row[offsets[8]] != "")
				result.phot_bp_mean_mag = std::stof(row[offsets[8]]);
			if (row[offsets[9]] != "")
				result.phot_rp_mean_mag = std::stof(row[offsets[9]]);
			if (row[offsets[10]] != "")
				result.distance_gspphot = std::stof(row[offsets[10]]);
			if (row[offsets[11]] != "")
				result.teff_gspphot = std::stof(row[offsets[11]]);

			result.hilbert_index = std::stol(row[offsets[12]]);

			fwrite(&result, sizeof(result), 1, stdout);
		}
	}

	return 0;
}
