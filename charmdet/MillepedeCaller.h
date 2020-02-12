#ifndef CHARMDET_MILLEPEDECALLER_H_
#define CHARMDET_MILLEPEDECALLER_H_

#include "TObject.h"
//includes for GBL fitter from genfit
#include <vector>
#include "GblPoint.h"
#include "GblTrajectory.h"
#include "MilleBinary.h"
#include "Track.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include <map>
#include <unordered_map>
#include "TVector3.h"
#include "TDecompLU.h"
#include "TRotation.h"
#include "TMath.h"
#include <cstdint>
#include "MufluxSpectrometer.h"
#include "MufluxSpectrometerHit.h"
#include <iostream>

//includes for MC testing
#include <random>
#include <fstream>

typedef enum
{
	MODULE,
	SINGLE_TUBE,
	LAYER
} alignment_mode;

/**
 * A class for wrapping the millepede function call such that it can be called from
 * within a python script
 *
 * @author Stefan Bieschke
 * @date Apr. 9, 2019
 */
class MillepedeCaller//: public TObject
{
public:
	MillepedeCaller(const char* out_file_name);
	virtual ~MillepedeCaller();

	double perform_GBL_refit(const genfit::Track& track, double sigma_spatial) const;
	double MC_GBL_refit(unsigned int n_tracks, double smearing_sigma, unsigned int min_hits = 3);
	void write_resolution_function(const char* filename, const genfit::Track& track, std::vector<MufluxSpectrometerHit>* raw_hits = nullptr) const;

	ClassDef(MillepedeCaller,3);

private:
	gbl::MilleBinary* m_gbl_mille_binary;

	//random generator
	std::mt19937 m_mersenne_twister;
	std::unordered_map<int,std::string> m_tube_id_to_module;
	std::unordered_map<std::string, std::vector<int>> m_modules; //detector IDs making up a module
	std::unordered_map<std::string,TVector3> m_nominal_module_centerpos; //nominal geometric center of a drift tube module

	//helper methods
	std::vector<gbl::GblPoint> list_hits(const genfit::Track* track, double sigma_spatial) const;
	void add_measurement_info(gbl::GblPoint& point, const TVector3& closest_approach, const double measurement, const double sigma_spatial) const;
	/*
	 * Helpers for jacobian calculation
	 */
	TMatrixD* calc_jacobian(const genfit::Track* track, const unsigned int hit_id_1, const unsigned int hit_id_2) const;
	std::pair<double,TMatrixD*> single_jacobian_with_arclength(const genfit::Track& track, const unsigned int hit_id) const;

	/*
	 * Helpers for projection and residuals
	 */
	TVector3 calc_shortest_distance(const TVector3& wire_top, const TVector3& wire_bot, const TVector3& track_pos, const TVector3& track_mom,  TVector3* PCA_on_wire = nullptr, TVector3* PCA_on_track = nullptr) const;
	TRotation calc_rotation_of_vector(const TVector3& v) const;
	TMatrixD rot_to_matrix(const TRotation& rot) const;
	TMatrixD calc_projection_matrix(const TMatrixD& fit_system_base_vectors, const TMatrixD& rotation_global_to_measurement) const;
	TVector3 calc_module_centerpos(const std::pair<std::string,std::vector<int>>& module_name_id_list_pair) const;

	/*
	 * GBL model methods
	 */
	std::vector<TVector3> linear_model_wo_scatter(const genfit::Track& track) const;
	void print_model_parameters(const std::vector<TVector3>& model) const;

	/*
	 * Pede parameter optimization
	 */
	std::vector<int> labels(const alignment_mode mode, const int channel_id) const;
	std::vector<int> labels_case_module(const int channel_id) const;
	TMatrixD* calc_global_parameters(const TVector3& measurement_prediction, const std::vector<TVector3>& linear_model, const TVector3& wire_bot_to_top) const;


	/*
	 * Checks for consistency and debugging
	 */
	bool check_ordered_by_arclen(const genfit::Track& track) const;
	void print_seed_hits(const genfit::Track& track) const;
	void print_fitted_residuals(gbl::GblTrajectory& trajectory) const;
	void print_fitted_track(gbl::GblTrajectory& trajectory) const;


	/*
	 * Monte-Carlo Tracks for testing
	 */
	std::vector<gbl::GblPoint> MC_list_hits(const std::vector<TVector3>& mc_track_model, double smearing_sigma, unsigned int min_hits);
	std::vector<TVector3> MC_gen_track();
	std::vector<std::pair<int,double>> MC_gen_hits(const TVector3& start, const TVector3& direction, const std::vector<int>* shifted_det_ids = nullptr);
	TMatrixD* calc_jacobian(const TVector3& PCA_1, const TVector3& PCA_2) const;

};

#endif /* CHARMDET_MILLEPEDECALLER_H_ */
