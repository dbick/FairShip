#include "MillepedeCaller.h"
#include <iostream>

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
MillepedeCaller::MillepedeCaller(const char *outFileName, bool asBinary, bool writeZero)
: mille(outFileName, asBinary, writeZero)
{

}

/**
 * Default destructor. Tears down the object
 *
 * @brief Destructor
 */
MillepedeCaller::~MillepedeCaller()
{

}

/**
 * Call the mille function of the Mille object contained in this class. For a documentation of the parameters, see the Millepede documentation at:
 * http://www.desy.de/~kleinwrt/MP2/doc/html/draftman_page.html
 *
 * @brief Mille caller
 *
 * @author Stefan Bieschke
 * @date May 9, 2019
 * @version 1.0
 *
 * @param n_local_derivatives Number of local derivatives.
 * @param local_derivatives Pointer to values of the local derivatives. Array of floats as long as n_local_derivatives states.
 * @param n_global_derivatives Number of global derivatives.
 * @param global_derivatives Pointer to values of the global derivatives. Array of floats as long as n_global_derivatives states.
 * @param label Array of unique labels to identify global parameters. This must be as long as n_global_derivatives states.
 * @param measured_residual Residual of the hit for that mille is called.
 * @param sigma Sigma of this hit.
 */
void MillepedeCaller::call_mille(int n_local_derivatives,
					const float *local_derivatives,
					int n_global_derivatives,
					const float *global_derivatives,
					const int *label,
					float measured_residual,
					float sigma)
{
	mille.mille(n_local_derivatives,local_derivatives,n_global_derivatives,global_derivatives,label,measured_residual,sigma);
}

vector<gbl::GblPoint> MillepedeCaller::list_hits(const genfit::Track& track) const
{
	vector<gbl::GblPoint> result(0);

	size_t n_points = track.getNumPointsWithMeasurement();
	vector<genfit::TrackPoint* > points = track.getPointsWithMeasurement();

	for(size_t i = 0; i < n_points; i++)
	{
		genfit::TrackPoint* point = points[i];
		genfit::AbsMeasurement* raw = point->getRawMeasurement();
		int detId = raw->getDetId();
	}

	//TODO continue implementing

//	    for i in range(n_points):
//	        point = points[i]
//	        raw_measurement = point.getRawMeasurement()
//	        det_id = raw_measurement.getDetId()
//	        rt_dist = raw_measurement.getRawHitCoords()[6] * u.mm #rt distance stored here
//	        # 2.) Parse detector id to module
//	        module_id = parse_det_id(det_id)
//	        module = dtmodules[module_id['module']]
//	        # 3.) Find correct drift tube in module
//	        for j in range(len(module.get_tubes())):
//	            tube = module.get_tubes()[j]
//	            if tube._ID == det_id:
//	                break
//	        tube = module.get_tubes()[j]

	return result;
}



const int* MillepedeCaller::labels() const
{
	return new int[100];
}

TMatrixD* MillepedeCaller::calc_jacobian(const genfit::Track& track, const unsigned int hit_id_1, const unsigned int hit_id_2) const
{
	TMatrixD* jacobian = new TMatrixD(5,5);

	return jacobian;
}

map<double,TMatrixD*> MillepedeCaller::jacobians_with_arclength(const genfit::Track& track) const
{
	map<double,TMatrixD*> result;


	//calculate length of the track between the two hits (in GBL terms arc length)
	TVector3 fitted_pos_1 = track.getFittedState(hit_id_1);
	TVector3 fitted_pos_2 = track.getFittedState(hit_id_2);
	TVector3 between_hits = fitted_pos_2 - fitted_pos_1;
	double distance = between_hits.Mag();

	TMatrixD* jacobian = calc_jacobian(track, hit_id_1, hit_id_2);



	return result;
}
