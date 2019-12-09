#include "MillepedeCaller.h"

using namespace std;


/**
 * Constructor. Initializes the wrapper to the Mille class. All paramaters of this constructor
 * are passed to the constructor of the Mille class.
 *
 * @brief Constructor
 *
 * @author Stefan Bieschke
 * @date Apr 9, 2019
 * @version 1.0
 *
 * @param outFileName string containing the filename written by mille function.
 * @param asBinary Flag that states if file should be written as binary
 * @param writeZero Flag that states if zero values should be kept or not
 */
MillepedeCaller::MillepedeCaller(const char* out_file_name)
{
	m_mersenne_twister = mt19937(19937);
	m_gbl_mille_binary = new gbl::MilleBinary(out_file_name,true,2000);

	m_tube_ids.resize(0);

	m_modules["T1U"] = {};
	m_modules["T2V"] = {};
	m_modules["T3aX"] = {};
	m_modules["T3bX"] = {};
	m_modules["T3cX"] = {};
	m_modules["T3dX"] = {};
	m_modules["T4aX"] = {};
	m_modules["T4bX"] = {};
	m_modules["T4cX"] = {};
	m_modules["T4dX"] = {};
	for(auto vec : m_modules)
	{
		vec.second.reserve(48);
	}

	//generate list of tube ids
	//T1 and T2
	for (char station = 1; station < 3; ++station)
	{
		for (char view = 0; view < 2; ++view)
		{
			for (char plane = 0; plane < 2; ++plane)
			{
				for (char layer = 0; layer < 2; ++layer)
				{
					for (char tube = 1; tube < 13; ++tube)
					{
						m_tube_ids.push_back(station*10000000+view*1000000+plane*100000+layer*10000+2000+tube);
					}
				}
			}
		}

	}

	//T3 and T4
	for (char station = 3; station < 5; ++station)
	{
		for (char plane = 0; plane < 2; ++plane)
		{
			for (char layer = 0; layer < 2; ++layer)
			{
				for (char tube = 1; tube < 49; ++tube)
				{
					stringstream module_key;
					module_key << "T" << (int)station;
					if(tube <= 12)
					{
						module_key << "dX";
					}
					else if(tube <= 24)
					{
						module_key << "cX";
					}
					else if(tube <= 36)
					{
						module_key << "bX";
					}
					else
					{
						module_key << "aX";
					}
					string key = module_key.str();
					int id = station*10000000+plane*100000+layer*10000+2000+tube;
					m_modules[key].push_back(id);
					m_tube_ids.push_back(id);
				}
			}
		}
	}
	vector<int> t1x = {};
	int station = 1;
	int view = 1;
	for (char plane = 0; plane < 2; ++plane)
	{
		for (char layer = 0; layer < 2; ++layer)
		{
			for (char tube = 1; tube < 13; ++tube)
			{
				t1x.push_back(station * 10000000 + view * 1000000 + plane * 100000
								+ layer * 10000 + 2000 + tube);
			}
		}
	}

	m_modules["T1X"] = t1x;
	vector<int> t1u = {};
}



/**
 * Default destructor. Tears down the object
 *
 * @brief Destructor
 */
MillepedeCaller::~MillepedeCaller()
{
	delete m_gbl_mille_binary;
}

/**
 * List the hits used to fit the seed track from a fit predecessing the GBL fit as a @c std::vector<gbl::GblPoint>. The GblPoint objects will
 * also contain measurements after the call of this method. The GBL points are ordered by arclength, which is the distance on the track between
 * two consecutive hits.
 *
 * @brief List hits from seed track as vector<GblPoint>
 *
 * @author Stefan Bieschke
 * @date Aug. 06, 2019
 * @version 1.0
 *
 * @param track Seed track, in this case genfit::Track from Kalman fitter
 *
 * @return std::vector<gbl::GblPoint> containing the hits ordered by arclen with measurement added
 */
vector<gbl::GblPoint> MillepedeCaller::list_hits(const genfit::Track* track, double sigma_spatial) const
{
	vector<TVector3> linear_model = linear_model_wo_scatter(*track);
	vector<gbl::GblPoint> result = {};

	//define projection of measurement system to the fit system
	TMatrixD fit_system_base_vectors(2,3);
	fit_system_base_vectors.Zero();
	fit_system_base_vectors[0][0] = 1.0; 	//first row vector for x direction
	fit_system_base_vectors[1][1] = 1.0; 	//second row vector for y direction

	vector<genfit::TrackPoint* > points = track->getPointsWithMeasurement();
	size_t n_points = points.size();
	result.reserve(n_points);


	//sort hits by arc length
	sort(points.begin(),points.end(),[](const genfit::TrackPoint* p1, const genfit::TrackPoint* p2){
		//sort for z-coordinates of bottom wire end - assumes no major differences in rotation around x axis for one layer
		return p1->getRawMeasurement()->getRawHitCoords()[2] < p2->getRawMeasurement()->getRawHitCoords()[2];
	});

	TVector3 PCA_track(0,0,0);
	//#pragma omp parallel for
	for(size_t i = 0; i < n_points; ++i)
	{
		TVector3 PCA_wire;
		TVector3 PCA_track_last = PCA_track;
		genfit::TrackPoint* point = points[i];
		//get coords of upper and lower end of hit tube
		genfit::AbsMeasurement* raw_measurement = point->getRawMeasurement();
		int det_id = raw_measurement->getDetId();
		TVectorD raw = raw_measurement->getRawHitCoords();
		TVector3 vbot(raw[0],raw[1],raw[2]);
		TVector3 vtop(raw[3],raw[4],raw[5]);
		double measurement = raw[6]; //rt distance [cm]

		TVector3 closest_approach = calc_shortest_distance(vtop,vbot,linear_model[0],linear_model[1],&PCA_wire,&PCA_track);

		TMatrixD* jacobian;
		if (i != 0)
		{
			jacobian = calc_jacobian(PCA_track_last, PCA_track);
		}
		else
		{
			jacobian = new TMatrixD(5,5);
			jacobian->UnitMatrix();
		}
		result.push_back(gbl::GblPoint(*jacobian));
		TRotation rot = calc_rotation_of_vector(closest_approach);
		TMatrixD rot_mat = rot_to_matrix(rot);
		TMatrixD projection_matrix = calc_projection_matrix(fit_system_base_vectors,rot_mat);
		TVectorD rotated_residual(2);
		rotated_residual[0] = closest_approach.Mag() - measurement;
		rotated_residual[1] = 0;
		//write ascii output for testing

		TVectorD precision(rotated_residual);
		precision[0] = 1.0 / TMath::Power(sigma_spatial,2);
		result.back().addMeasurement(projection_matrix,rotated_residual,precision);

		//calculate labels and global derivatives for hit
		vector<int> label = labels(MODULE,det_id);
		TMatrixD* globals = calc_global_parameters(PCA_track,linear_model);
		result.back().addGlobals(label, *globals);
		delete globals;
		delete jacobian;
	}

	return result;
}


//TODO implement
vector<int> MillepedeCaller::labels(const alignment_mode mode, const int channel_id) const
{
	switch(mode)
	{
	case MODULE:
		return labels_case_module(channel_id);
		break;
	case SINGLE_TUBE:
		return {0,0,0,0,0,0};
		break;
	case LAYER:
		return {0,0,0,0,0,0};
		break;
	default:
		return {0,0,0,0,0,0};
		break;
	}
}

vector<int> MillepedeCaller::labels_case_module(const int channel_id) const
{
	int station = -1;
	int layer = -1;
	int view = -1;
	int tube = -1;
	int pnb = -1;

	MufluxSpectrometer::TubeDecode(channel_id, station, view, pnb, layer, tube);
	vector<int> labels(6);
	int base_label = 1000 * station;
	if(station > 2)
	{
		//TODO make modules close to beam axis first in if checks
		if(tube >= 37 && tube <= 48)
		{
			base_label += 10;
		}
		else if(tube >= 25 && tube <= 36)
		{
			base_label += 20;
		}
		else if(tube >= 13 && tube <= 24)
		{
			base_label += 30;
		}
	}
	else
	{
		base_label += 100 * view;
	}

	for(uint8_t i = 0; i < labels.size(); ++i)
	{
		labels[i] = base_label + i + 1;
	}

	return labels;
}

/**
 * Calculate the global parameters for pede alignment. This is a matrix with dimensions 3x6 for a rigid body alignment.
 * A measurement \f$m\f$ (being a three dimensional vector, e.g x,y,z), the distorted measurement \f$\tilde{m}\f$ depends on the six\
 * (global - thus affecting ALL measurements) parameters \f$ g = {\Delta x, \Delta y, \Delta z, \alpha, \beta, \gamma}\f$.
 * Then it can be expressed as:
 * \f[
 * \tilde{m} = m + \begin{pmatrix}
 * \Delta x \\ \Delta y \\ \Delta z
 * \end{pmatrix} + \alpha \begin{pmatrix}
 *	0 \\ -z_p \\ y_p
 * \end{pmatrix} + \beta \begin{pmatrix}
 *	z_p \\ 0 \\ -x_p
 * \end{pmatrix} + \gamma \begin{pmatrix}
 * -y_p \\ x_p \\ 0
 * \end{pmatrix}
 * \f]
 *
 * @brief Calculate the global parameters for alignment
 *
 * @author Stefan Bieschke
 * @date Dec. 2, 2019
 * @version 1.0
 *
 * @param measurement_prediction predicted vector from the wire to the seed track
 *
 * @result Matrix (3x6) containing the derivatives of the measurement w.r.t all parameters
 */
TMatrixD* MillepedeCaller::calc_global_parameters(const TVector3& measurement_prediction, const vector<TVector3>& linear_model) const
{
	TMatrixD dmdg(3,6);
	TMatrixD* result = new TMatrixD(3,6);
	TMatrixD drdm(3,3);
	drdm.UnitMatrix();
	TVector3 track_direction = linear_model[1];
	TVector3 nominal_measurementplane_normal(0,0,1);
	double scalar_prod = track_direction.Dot(nominal_measurementplane_normal);

	for(short i = 0; i < drdm.GetNrows(); ++i)
	{
		for(short j = 0; j < drdm.GetNcols(); ++j)
		{
			drdm[i][j] -= track_direction[i] * nominal_measurementplane_normal[j] / scalar_prod;
		}
	}

	dmdg.Zero();
	result->Zero();

	dmdg[0][0] = dmdg[1][1] = dmdg[2][2] = 1;
	dmdg[0][4] = measurement_prediction[2];
	dmdg[0][5] = - measurement_prediction[1];
	dmdg[1][3] = - measurement_prediction[2];
	dmdg[1][5] = measurement_prediction[0];
	dmdg[2][3] = measurement_prediction[1];
	dmdg[2][4] = -measurement_prediction[0];

	result->Mult(drdm, dmdg);

	return result;
}


//TODO add mathematical description of matrix parameter determination to doc comment
/**
 * Calculate the jacobian between two consecutive hits on a predetermined track. This track (in this case) is the result of a predecessing
 * fit with a KalmanFitter from genfit. It can be rewritten to be any track, even a "guessed" one from any kind of pattern recognition.
 *
 * @brief Calculate jacobian between two consecutive hits
 *
 * @author Stefan Bieschke
 * @version 1.0
 * @date July 30, 2019
 *
 * @param track genfit::Track object, of whom the jacobian for two consecutive hits is calculated
 * @param hit_id_1 ID of the first hit on the track for that the jacobian should be calculated
 * @param hit_id_2 ID of the second hit on the track for that the jacobian should be calculated
 *
 * @return Pointer to a heap object of TMatrixD type with dimensions 5x5
 *
 * @warning Requires hit_id_1 = hit_id_2 - 1
 * @warning Requires hit_id_1 < track.getNumPointsWithMeasurement() && hit_id_2 < track.getNumPointsWithMeasurement()
 * @warning Heap object without auto deletion
 */
TMatrixD* MillepedeCaller::calc_jacobian(const genfit::Track* track, const unsigned int hit_id_1, const unsigned int hit_id_2) const
{
	TMatrixD* jacobian = new TMatrixD(5,5);

	// 1.) init unity matrix
	jacobian->UnitMatrix();

	//2.) enter non-zero partial differentials
	//2.1) get the two points on track where reconstruction happened
	genfit::MeasuredStateOnPlane state_at_id_1 = track->getFittedState(hit_id_1);
	genfit::MeasuredStateOnPlane state_at_id_2 = track->getFittedState(hit_id_2);

	TVector3 pos1 = state_at_id_1.getPos();
	TVector3 pos2 = state_at_id_2.getPos();

	double dz = pos2.Z() - pos1.Z();

	//2.2) enter dx and dy to jacobian
	(*jacobian)[3][1] = dz;
	(*jacobian)[4][2] = dz;

	return jacobian;
}

/**
 * Calculate the jacobian matrix for the transport of track parameters from the previous hit to a hit labelled by @c hit_id. Alongside
 * the jacobian matrix, an arclength is calculated, which is the distance on the track between the previous hit and the hit specified by
 * @c hit_id. Note that as the transport of parameters from the previous hit is calculated, the parameter @c hit_id is required to be
 * greater than zero. See the requirements below for more details.
 * The result is returned as std::pair<double,TMatrixD*>. The TMatrixD* pointer points to heap memory holding the 5x5 jacobian matrix.
 * For more info on the jacobian see the documentation of the method
 * @c calc_jacobian(const genfit::Track* track, const unsigned int hit_id_1, const unsigned int hit_id_2) const.
 *
 * @brief Calculate a single pair of jacobian with arclength for a given @c hit_id
 *
 * @author Stefan Bieschke
 * @date Aug. 06, 2019
 * @version 1.0
 *
 * @param track Seed track from a fit predecessing the GBL refit, e.g the genfit Kalman fitter
 * @param hit_id ID (number) of the hit for that the pair shall be computed. See requirements below
 *
 * @require 0 < hit_id < track.getNumPointsWithMeasurement()
 *
 * @return std::pair<double,TMatrixD*> containing the the arclength and the jacobian matrix
 *
 * @note Only needed, when hits on seed track are not ordered correctly by arclength, ca be checked with method check_ordered_by_arclen(const genfit::Track&)
 * @warning TMatrixD* is heap memory and undeleted. Must be deleted by caller when no longer needed.
 */
pair<double,TMatrixD*> MillepedeCaller::single_jacobian_with_arclength(const genfit::Track& track, const unsigned int hit_id) const
{
	//arc length is distance from very first point on the track
	TVector3 fitted_pos_1 = track.getFittedState(0).getPos();
	TVector3 fitted_pos_2 = track.getFittedState(hit_id).getPos();
	TVector3 between_hits = fitted_pos_2 - fitted_pos_1;

	double distance = between_hits.Mag();

	TMatrixD* jacobian = calc_jacobian(&track, hit_id - 1, hit_id);

	return make_pair(distance, jacobian);
}

//TODO document
/**
 *
 */
double MillepedeCaller::perform_GBL_refit(const genfit::Track& track, double sigma_spatial) const
{
	try
	{
		vector < gbl::GblPoint > points = list_hits(&track, sigma_spatial);
		gbl::GblTrajectory traj(points,false); //param false for B=0

		traj.milleOut(*m_gbl_mille_binary);
		//check track validity
		if(!traj.isValid())
		{
			cout << "Error, GBL trajectory is invalid." << endl;
			cerr << "Error, GBL trajectory is invalid." << endl;
			return -1;
		}

		int rc, ndf;
		double chi2, lostWeight;

		cout << "------------performing refit--------------" << endl;
		cout << "Seed track chi2: " << track.getFitStatus()->getChi2() << " Ndf: " << track.getFitStatus()->getNdf() << endl;


		rc = traj.fit(chi2,ndf,lostWeight);
		cout << "Refit chi2: " << chi2 << " Ndf: " << ndf << endl;
		cout << "Prob: " << TMath::Prob(chi2,ndf) << endl;

//		traj.printTrajectory(1);

		return chi2;
	}
	catch(const genfit::Exception& e)
	{
		return -1;
	}
}

//TODO document
double MillepedeCaller::MC_GBL_refit(unsigned int n_tracks, double smearing_sigma, unsigned int min_hits)
{
	double chi2, lostweight;
	int ndf;
	vector<vector<TVector3>> tracks(n_tracks);
	for(unsigned int i = 0; i < n_tracks; ++i)
	{
		tracks[i] = MC_gen_track();
	}

	unsigned int fitted = 0;
	for(int i = 0; i < tracks.size(); ++i)
	{
		auto track = tracks[i];
		vector<gbl::GblPoint> hitlist = MC_list_hits(track,fitted,smearing_sigma,min_hits);
		if(hitlist.size() < min_hits)
		{
			continue;
		}
		gbl::GblTrajectory traj(hitlist, false);
		traj.milleOut(*m_gbl_mille_binary);
		traj.fit(chi2, ndf, lostweight);
		cout << "MC chi2: " << chi2 << " Ndf: " << ndf << endl;
		cout << "Prob: " << TMath::Prob(chi2,ndf) << endl;
		++fitted;
	}
	cout << "Fitted " << fitted << " out of " << n_tracks << " tracks." << endl;

	return 0.0;
}


void MillepedeCaller::write_resolution_function(const char* filename,
		const genfit::Track& track,
		std::vector<MufluxSpectrometerHit>* raw_hits) const
{
	/*
	 * File output of distance(track-wire) vs. residuals
	 */
	ofstream resolutionfunction("resolution.ascii", ios::app);
	vector<TVector3> linear_model = linear_model_wo_scatter(track);
	vector<MufluxSpectrometerHit*>* fitted = nullptr;
	vector<genfit::TrackPoint*> points = track.getPointsWithMeasurement();
	size_t n_points = points.size();
	if (raw_hits)
	{
		fitted = new vector<MufluxSpectrometerHit*>(0);
		for (unsigned int i = 0; i < n_points; ++i)
		{
			genfit::TrackPoint* point = points[i];
			int id = point->getRawMeasurement()->getDetId();
			for (unsigned int j = 0; j < raw_hits->size(); ++j)
			{
				MufluxSpectrometerHit raw = raw_hits->at(j);
				if (raw.GetDetectorID() == id)
				{
					fitted->push_back((MufluxSpectrometerHit*) &(raw_hits->at(j)));
				}
			}
		}
	}


	for(size_t i = 0; i < n_points; ++i)
	{
		genfit::TrackPoint* point = points[i];
		//get coords of upper and lower end of hit tube
		genfit::AbsMeasurement* raw_measurement = point->getRawMeasurement();
		int det_id = raw_measurement->getDetId();
		TVectorD raw = raw_measurement->getRawHitCoords();
		TVector3 vbot(raw[0],raw[1],raw[2]);
		TVector3 vtop(raw[3],raw[4],raw[5]);
		double measurement = raw[6]; //rt distance [cm]

		TVector3 closest_approach = calc_shortest_distance(vtop,vbot,linear_model[0],linear_model[1]);

		//write ascii output for testing
		resolutionfunction << closest_approach.Mag() << "\t" << closest_approach.Mag() - measurement;
		if(fitted)
		{
			resolutionfunction << "\t" << fitted->at(i)->GetTimeOverThreshold();
		}
		resolutionfunction << endl;

	}
	resolutionfunction.close();
	if(fitted) delete fitted;
}

//Reimplementation of python function
/**
 * Calculates a vector the represents the shortest distance from the sense wire to the track. This vector is perpendicular both
 * to the track as well as the sense wire.
 *
 * @brief Calculates shortest distance from sense wire to track
 *
 * @author Stefan Bieschke
 * @date Oct. 01, 2019
 * @version 1.1
 *
 * @param wire_top Top position of the sense wire as TVector3 with x,y,z components
 * @param wire_bot Bottom position of the sense wire as TVector3 with x,y,z components
 * @param track_pos Some position on the track. Could be anything but must be on the (straight) track or track segment
 * @param track_mom Momentum vector of the (straight) track or track segment
 * @param PCA_on_wire Coordinates of point of closest approach (PCA) on the wire (defaults to nullptr)
 * @param PCA_on_track Coordinates of point of closest approach (PCA) on the track (defaults to nullptr)
 *
 * @return TVector3 of shortest distance pointing from the wire to the track
 */
TVector3 MillepedeCaller::calc_shortest_distance(const TVector3& wire_top,
		const TVector3& wire_bot, const TVector3& track_pos,
		const TVector3& track_mom, TVector3* PCA_on_wire,
		TVector3* PCA_on_track) const
{
	TVector3 wire_dir = wire_top - wire_bot;

	//construct a helper plane that contains one of the straights and is parallel to the other one
	TVector3 plane_pos = track_pos - wire_bot;
	TVector3 plane_dir_1(track_mom);
	TVector3 plane_dir_2(-1 * wire_dir);

	//Construct components of equation system M * x = c where M is the coefficient matrix, x the solution and c the const_vector below
	TVectorD const_vector(2);
	TMatrixD coeff_matrix(2,2);

	const_vector[0] = -(plane_pos.Dot(track_mom));
	const_vector[1] = -(plane_pos.Dot(wire_dir));

	coeff_matrix[0][0] = plane_dir_1.Dot(track_mom);
	coeff_matrix[0][1] = plane_dir_2.Dot(track_mom);
	coeff_matrix[1][0] = plane_dir_1.Dot(wire_dir);
	coeff_matrix[1][1] = plane_dir_2.Dot(wire_dir);

	TDecompLU solvable_matrix(coeff_matrix);
	TVectorD result(const_vector);
	int rc = solvable_matrix.Solve(result);

	TVector3 PCA_track(track_pos + result[0] * track_mom);
	TVector3 PCA_wire(wire_bot + result[1] * wire_dir);

	//if TVector3 pointers were passed as parameter the actual locations of closest approach are stored there
	if(PCA_on_wire)
	{
		*(PCA_on_wire) = PCA_wire;
	}
	if(PCA_on_track)
	{
		*(PCA_on_track) = PCA_track;
	}

	return TVector3(PCA_track - PCA_wire);
}

/**
 * Calculate the rotation matrix in a way such that a given vector v is the x-axis of the rotated coordinate frame.
 *
 * @brief Rotation of a given vector in lab frame so that it is new x axis
 *
 * @author Stefan Bieschke
 * @date Aug. 09, 2019
 * @version 1.0
 *
 * @param v TVector3 that is meant to be the new x axis of the rotated coordinate frame
 *
 * @return TRotation containing the rotation matrix
 */
TRotation MillepedeCaller::calc_rotation_of_vector(const TVector3& v) const
{
	TRotation rot;
	rot.SetXAxis(v);

	return rot;
}

/**
 * Convert a TRotation to a its 3x3 rotation matrix given as TMatrixD
 *
 * @brief Convert a TRotation to a (3x3) TMatrixD
 *
 * @author Stefan Bieschke
 * @date Sep. 09, 2019
 * @version 1.0
 *
 * @param rot TRotation object for that the rotation matrix is needed
 *
 * @return TMatrixD object with dimensions 3x3
 */
TMatrixD MillepedeCaller::rot_to_matrix(const TRotation& rot) const
{
	TMatrixD result(3,3);
	for(uint8_t i = 0; i < 3; ++i)
	{
		for(uint8_t j = 0; j < 3; ++j)
		{
			result[i][j] = rot[i][j];
		}
	}
	return result;
}

/**
 * Calculate a position and a direction vector for a linear model that models a given track from genfit. While the track
 * from genfit might take scattering in detector material into account, this linear model doesn't consider scatterers.
 * The way this works is to calculate and return a straight track between the first and the last hit of the passed track.
 *
 * @brief Linear model of the passed track
 *
 * @author Stefan Bieschke
 * @date Sep. 09, 2019
 * @version 1.0
 *
 * @param track Track that is meant to be modeled
 *
 * @return std::vector of length two, whos first entry is the position vector, second one being the direction vector
 */
vector<TVector3> MillepedeCaller::linear_model_wo_scatter(const genfit::Track& track) const
{
	vector<TVector3> result(2);

	//get first and last fitted states
	size_t n_hits = track.getNumPointsWithMeasurement();
	genfit::StateOnPlane first_hit = track.getFittedState(0);
	genfit::StateOnPlane last_hit = track.getFittedState(n_hits - 1);

	//position is fitted position on first hit
	TVector3 pos = first_hit.getPos();

	//direction is difference between positions at first and last hits
	TVector3 dir = last_hit.getPos() - first_hit.getPos();
	result[0] = pos;
	result[1] = dir;

	return result;
}

/**
 * Calculate the projection matrix (dimension 2x2) from the global reference frame (x,y,z) to the measurement system (u,v,w) for an
 * individual hit.
 *
 * @brief Calculate projection matrix global to measumrent system
 *
 * @author Stefan Bieschke
 * @date Oct. 07, 2019
 * @version 1.0
 *
 * @param fit_system_base_vectors Matrix with dimensions 2x3 with the base vectors of the local fit system expressed in global system. E.g [(1., 0., 0.),(0.,1.,0.)] for x,y and track in general z direction.
 * @param rotation_global_to_measurement Rotation matrix rotating the global frame to the measurement frame
 *
 * @return TMatrixD (dimension 2x2) containing the projection of the fit system (x,y) on the measurement system
 */
TMatrixD MillepedeCaller::calc_projection_matrix(
		const TMatrixD& fit_system_base_vectors,
		const TMatrixD& rotation_global_to_measurement) const
{
	TMatrixD result(2,2); //projection matrix has dimensions 2x2
	TMatrixD measurement_to_global(rotation_global_to_measurement); //copy rotational matrix
//	measurement_to_global.Invert(); //invert matrix in place
	measurement_to_global.ResizeTo(3,2); //TODO check if this is correct, want to skip column normal to measurement direction
	result.Mult(fit_system_base_vectors,measurement_to_global);
	result.Invert();
	return result;
}


void MillepedeCaller::print_model_parameters(const vector<TVector3>& model) const
{
	double slope_x = model[1].X() / model[1].Z();
	double slope_y = model[1].Y() / model[1].Z();

	cout << "Printing linear track model parameters" << endl;
	cout << "Initial fit position:" << endl;
	model[0].Print();
	cout << "Track direction:" << endl;
	model[1].Print();
	cout << "Model parameters: x0, y0, slope x, slope y" << endl;
	cout << "(" << model[0].X() << ", " << model[0].Y() << ", " << slope_x  << ", " << slope_y <<")" << endl;
}

/**
 * Checks, if the points on a genfit seed track are ordered by arc length. The arc length
 * is the distance from an early point on the track, e.g the vertex or the very first hit.
 * This doesn't matter if the arc length itself is not used but the GblPoints must only
 * be ordered correctly.
 *
 * @brief Checks if points on seed track are ordered by arc length
 *
 * @author Stefan Bieschke
 * @date Nov. 6, 2019
 * @version 1.0
 *
 * @param track genfit seed track
 * @return true if points are correctly ordered, false else
 *
 * @note This is mainly for debugging and maintenance purposes, not recommended to be called in productive code
 */
bool MillepedeCaller::check_ordered_by_arclen(const genfit::Track& track) const
{
	size_t n_points = track.getNumPointsWithMeasurement();
	double z = 0.0;

	for(size_t i = 0; i < n_points; ++i)
	{
		double z_next = track.getFittedState(i).getPos().Z();
		if(z_next > z)
		{
			z = z_next;
		}
		else
		{
			return false;
		}
	}

	return true;
}

void MillepedeCaller::print_seed_hits(const genfit::Track& track) const
{
	size_t n_hits = track.getNumPointsWithMeasurement();
	const vector< genfit::TrackPoint* > points = track.getPoints();

	for(size_t i = 0; i < n_hits; ++i)
	{
		int det_id = points[i]->getRawMeasurement()->getDetId();
		cout << "Hit: " << i << "\t" << "ID: " << det_id << endl;
	}
}


/**
 * Produce a list of GblPoint objects from which a GblTrajectory can be constructed. It is passed a simple MC track, an event id,
 * the sigma of a gaussian smearing of the hits as well as (optionally) the minimum number of hits required. If the track has
 * too few hits, an empty list will be returned.
 *
 * @brief Produce a sorted list of GblPoint objects from which a GblTrajectory can be built
 *
 * @author Stefan Bieschke
 * @date Nov. 22, 2019
 * @version 1.0
 *
 * @param mc_track_model Simple MC track (for example generated with @c MC_gen_track()
 * @param event id Number identifying the track
 * @param smearing_sigma Sigma of a gaussian distribution (mu = 0) with which the hits are smeared
 * @param min_hits minimum number of hits required. If track has less hits, returned list will be empty
 *
 * @return std::vector<gbl::GblPoint> containing all hits, except track has less than @c min_hits hits, then it will be empty
 */
vector<gbl::GblPoint> MillepedeCaller::MC_list_hits(const vector<TVector3>& mc_track_model, int event_id, double smearing_sigma, unsigned int min_hits)
{
	//apply gaussian smearing of measured hit
	normal_distribution<double> gaussian_smear(0,smearing_sigma); //mean 0, sigma 350um in cm

//	vector<pair<int,double>> hits = MC_gen_hits(mc_track_model[0], mc_track_model[1], &(m_modules["T3bX"]));
	vector<pair<int,double>> hits = MC_gen_hits(mc_track_model[0], mc_track_model[1]);
	if(hits.size() < min_hits)
	{
		return {};
	}

	vector<gbl::GblPoint> gbl_hits = {};

	TMatrixD fit_system_base_vectors(2,3);
	fit_system_base_vectors.Zero();
	fit_system_base_vectors[0][0] = 1.0; 	//first row vector for x direction
	fit_system_base_vectors[1][1] = 1.0; 	//second row vector for y direction



	//zero arc length GBL point for first hit
	TMatrixD* unity = new TMatrixD(5,5);
	unity->UnitMatrix();

	TVector3 vbot, vtop;
	MufluxSpectrometer::TubeEndPoints(hits[0].first, vtop, vbot);
	TVector3 PCA_wire;
	TVector3 PCA_track;
	TVector3 closest_approach = calc_shortest_distance(vtop, vbot,mc_track_model[0],mc_track_model[1], &PCA_wire, &PCA_track);

	gbl_hits.push_back(gbl::GblPoint(*unity));
	TRotation rot = calc_rotation_of_vector(closest_approach);
	TMatrixD rot_mat = rot_to_matrix(rot);
	TMatrixD projection_matrix = calc_projection_matrix(fit_system_base_vectors,rot_mat);
	TVectorD rotated_residual(2);
	double smear = gaussian_smear(m_mersenne_twister);
	rotated_residual[0] = closest_approach.Mag() - (hits[0].second + smear);
	rotated_residual[1] = 0;
	TVectorD precision(rotated_residual);
	precision[0] = 1.0 / (smearing_sigma * smearing_sigma);
	gbl_hits.back().addMeasurement(projection_matrix,rotated_residual,precision);

	//add labels and derivatives for first hit
	vector<int> label = labels(MODULE,hits[0].first);
	TMatrixD* globals = calc_global_parameters(PCA_track,mc_track_model);
	gbl_hits.back().addGlobals(label, *globals);
	delete globals;

	for(size_t i = 1; i < hits.size(); ++i)
	{
		TVector3 PCA_track_old = PCA_track;
		MufluxSpectrometer::TubeEndPoints(hits[i].first, vtop, vbot);
		TVector3 closest_approach = calc_shortest_distance(vtop,vbot,mc_track_model[0],mc_track_model[1],&PCA_wire,&PCA_track);
		TMatrixD* jacobian = calc_jacobian(PCA_track_old, PCA_track);
		gbl_hits.push_back(gbl::GblPoint(*jacobian));
		TRotation rot = calc_rotation_of_vector(closest_approach);
		TMatrixD rot_mat = rot_to_matrix(rot);
		TMatrixD projection_matrix = calc_projection_matrix(fit_system_base_vectors,rot_mat);
		TVectorD rotated_residual(2);
		smear = gaussian_smear(m_mersenne_twister);
		rotated_residual[0] = closest_approach.Mag() - (hits[i].second + smear);
		rotated_residual[1] = 0;
		TVectorD precision(rotated_residual);
		precision[0] = 1.0 / (smearing_sigma * smearing_sigma);
		gbl_hits.back().addMeasurement(projection_matrix,rotated_residual,precision);
		delete jacobian;

		//calculate labels and global derivatives for hit
		vector<int> label = labels(MODULE,hits[i].first);
		TMatrixD* globals = calc_global_parameters(PCA_track,mc_track_model);
		gbl_hits.back().addGlobals(label, *globals);
		delete globals;
	}

	return gbl_hits;
}

/**
 * Generates a very simple toy MC track. There is no physics involved, this is just meant to test the GBL refit. The
 * MC tracks return a vector<TVector3> containing two TVector3 objects. The first being a position in front of T1, the second
 * one a direction, the track propagates from there. These are only straight lines, no scattering or any physics involved.
 *
 * @brief Generate a very simple toy MC track
 *
 * @author Stefan Bieschke
 * @date Nov. 22, 2019
 * @version 1.0
 *
 * @return std::vector<TVector3> two element vector with first element position (e.g vertex) and second element direction
 */
vector<TVector3> MillepedeCaller::MC_gen_track()
{
	/*
	 * Track locations uniformly distributed in front of first plane and in front of first plane of T4.
	 *
	 * area for x and y is calculated as follows:
	 * allowed x area is centerpos +- 21 cm for T1, (centerpos T4d - 21cm) to (centerpos T4a + 21cm) for T4
	 * allowed y area is centerpos +- 21 cm for T1, centerpos +- 75 cm for T4
	 *
	 * Positions below calculated with: https://github.com/StBies/FairShip/blob/a7ead53a260a1b796371eb09ec406421c027ad18/charmdet/drifttubeMonitoring.py#L955
     * Module: T1X Centerpos = 1.8225      -1.995  24.99
     * Module: T1U Centerpos = 1.33890796894       0.782207722399  56.20825
     * Module: T2V Centerpos = 1.71329088596       0.663508489501  93.68275
     * Module: T2X Centerpos = 1.03        -2.0875 125.082
     * Module: T3aX Centerpos = 76.9353333333       10.512  583.19
	 * Module: T3bX Centerpos = 25.987      10.289  585.065
     * Module: T3cX Centerpos = -24.539     10.28075        585.1175
     * Module: T3dX Centerpos = -75.0575    10.23975        585.315
     * Module: T4aX Centerpos = 76.74       10.163  748.2675
     * Module: T4bX Centerpos = 26.3425     9.939   748.1975
     * Module: T4cX Centerpos = -24.0075    9.8685  748.27
	 * Module: T4dX Centerpos = -74.48      9.92025 748.4325
	 */
	uniform_real_distribution<double> uniform(0.0, 1.0);
	double offset_x_beginning = -19.1775 + 42.0 * uniform(m_mersenne_twister);
	double offset_y_beginning = -22.995 + 42.0 * uniform(m_mersenne_twister);
	double offset_x_end = -95.48 + 193.22 * uniform(m_mersenne_twister);
	double offset_y_end = -65.0 + 150.0 * uniform(m_mersenne_twister);

	double z_beginning = 17.0;
	double z_end = 739.0;

	TVector3 beginning = TVector3(offset_x_beginning, offset_y_beginning, z_beginning);
	TVector3 end = TVector3(offset_x_end, offset_y_end, z_end);
	TVector3 direction = end - beginning;
	vector<TVector3> result = {beginning, direction};

	return result;
}

/**
 * Generate unsmeared hits for a given MC track. The track is passed with its startpoint and direction to this method.
 * When generating hits the closest distance to the wire is calculated for every wire in the detector to this track. If the
 * distance is smaller 18.15mm, a hit is written. The hits are stored in a vector<pair<int,double>>. The vector is as long as
 * the number of hits that were found. Each element is a pair of int and double, the int being the hit detectorID, the double
 * being the unsmeared hit distance (drift distance). The hits are sorted by arc length (which is the distance from the very
 * first hit).
 *
 * @brief Generate MC hits for a given track
 *
 * @author Stefan Bieschke
 * @date Nov. 22, 2019
 * @version 1.0
 *
 * @param start Starting point (like vertex) of MC track
 * @param direction Direction of the MC track
 *
 * @return std::vector<std::pair<int,double>> list of pairs of int and double, the first one being the detectorID, the second one the rt distance
 */
vector<pair<int,double>> MillepedeCaller::MC_gen_hits(const TVector3& start, const TVector3& direction, const vector<int>* shifted_det_ids)
{
	vector<pair<int,double>> result(0);
	TVector3 wire_end_top;
	TVector3 wire_end_bottom;
	TVector3 wire_to_track;

	//check distance to every tube
	for(int id : m_tube_ids)
	{
		MufluxSpectrometer::TubeEndPoints(id, wire_end_top, wire_end_bottom);
		if(shifted_det_ids)
		{
			TVector3 translation(-0.07, 0.01, 0.08);
//			TVector3 rotation(TMath::Pi() / 180, - 0.4 * TMath::Pi() / 180 , 0.3 * TMath::Pi() / 180);
			if(shifted_det_ids->end() != find(shifted_det_ids->begin(),shifted_det_ids->end(),id))
			{
				wire_end_bottom = wire_end_bottom + translation;
				wire_end_top = wire_end_top + translation;
			}
		}
		wire_to_track = calc_shortest_distance(wire_end_top, wire_end_bottom, start, direction, nullptr, nullptr);
		if(wire_to_track.Mag() < 1.815)
		{
			result.push_back(pair<int,double>(id,wire_to_track.Mag()));
		}
	}

	//sort with lambda comparison
	sort(result.begin(), result.end(), [&](pair<int,double> element1, pair<int,double> element2){
		MufluxSpectrometer::TubeEndPoints(element1.first, wire_end_top, wire_end_bottom);
		double z1 = (wire_end_bottom + ((wire_end_top - wire_end_bottom)* 0.5)).Z();
		MufluxSpectrometer::TubeEndPoints(element2.first, wire_end_top, wire_end_bottom);
		double z2 = (wire_end_bottom + ((wire_end_top - wire_end_bottom)* 0.5)).Z();
		return z1 < z2;
		});

	return result;
}

/**
 * Calculate the jacobi matrix for a linear model from the points of closest approach on the track for two consecutive hits.
 *
 * @brief Calculate jacobian between two consecutive hits
 *
 * @author Stefan Bieschke
 * @date November 22, 2019
 * @version 1.0
 *
 * @param PCA_1 Point of closest approach on the track for the first hit
 * @param PCA_2 Point of closest approach on the track for the second hit
 *
 * @return Pointer to a heap object of TMatrixD type with dimensions 5x5
 *
 * @warning Heap object without auto deletion
 */
TMatrixD* MillepedeCaller::calc_jacobian(const TVector3& PCA_1, const TVector3& PCA_2) const
{
	TMatrixD* jacobian = new TMatrixD(5,5);

	// 1.) init unity matrix
	jacobian->UnitMatrix();

	//2.) enter non-zero partial differentials
	//2.1) get the two points on track where reconstruction happened
	double dz = PCA_2.Z() - PCA_1.Z();

	//2.2) enter dx and dy to jacobian
	(*jacobian)[3][1] = dz;
	(*jacobian)[4][2] = dz;

	return jacobian;
}
