#pragma once
#include <cmath>

class Star {
public:
	long hilbert_index;
	float ra;
	float dec;
	short astrometric_n_obs_al;
	short astrometric_n_obs_ac;
	short astrometric_n_good_obs_al;
	short astrometric_n_bad_obs_al;
	short visibility_periods_used;
	float phot_g_mean_mag;
	float phot_bp_mean_mag;
	float phot_rp_mean_mag;
	float distance_gspphot;
	float teff_gspphot;

	Star() {
		hilbert_index = 0;
		ra = nan("");
		dec = nan("");
		astrometric_n_obs_al = 0;
		astrometric_n_obs_ac = 0;
		astrometric_n_good_obs_al = 0;
		astrometric_n_bad_obs_al = 0;
		visibility_periods_used = 0;

		phot_g_mean_mag = nan("");
		phot_bp_mean_mag = nan("");
		phot_rp_mean_mag = nan("");
		distance_gspphot = nan("");
		teff_gspphot = nan("");
	}
};
