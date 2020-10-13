/*
 * GBLseedtrack.cpp
 *
 *  Created on: 29.04.2020
 *      Author: bieschke
 *  Modified on: 12.06.2020
 *      Author: Daniel Bick
 */

#include "GBLseedtrack.h"
#include "TDecompLU.h"

using std::vector;
using std::pair;

GBL_seed_track::GBL_seed_track(const genfit::Track& track)
{
	vector<genfit::TrackPoint*> points = track.getPoints();
	size_t n_points = points.size();

	m_hits.resize(n_points);

	for(size_t i = 0; i < n_points; ++i)
	{
		genfit::TrackPoint* point = points[i];
		genfit::AbsMeasurement* raw_measurement = point->getRawMeasurement();
		int det_id = raw_measurement->getDetId();
		double drift_radius = raw_measurement->getRawHitCoords()[6];
		pair<int,double> hit_info = std::make_pair(det_id, drift_radius);

		m_hits[i] = hit_info;
	}

	m_position = track.getFittedState().getPos();
	m_direction = track.getFittedState().getMom();
}


GBL_seed_track::GBL_seed_track(const vector<TVector3>& pos_mom, const vector<pair<int,double>> hits)
{
	m_position = pos_mom[0];
	m_direction = pos_mom[1];
	m_hits.resize(hits.size());
	//TODO copy ctor invokation
	for(size_t i = 0; i < m_hits.size(); ++i)
	{
		m_hits[i] = hits[i];
	}
}


/*
 * Constructor for seed track from pattern reco. Creates seedrack without hits yet, only momentum and direction is provided. Hits should then be added by add_hit(). 
 *
 * @brief Constructor
 *
 * @param position A point on the track
 * @param direction Direction of the track
 *
 * @author Daniel Bick
 * @date 12.06.2020
 * @version 1.0
 * 
 */
GBL_seed_track::GBL_seed_track(TVector3 position, TVector3 direction)
{
	m_position = position;
	m_direction = direction;
}





GBL_seed_track::GBL_seed_track::~GBL_seed_track()
{
	// TODO Auto-generated destructor stub
}

const size_t GBL_seed_track::get_number_of_hits() const
{
	return m_hits.size();
}


const TVector3& GBL_seed_track::get_position() const
{
	return m_position;
}

const TVector3& GBL_seed_track::get_direction() const
{
	return m_direction;
}

const vector<pair<int,double>>& GBL_seed_track::get_hits() const
{
	return m_hits;
}



const vector<int> GBL_seed_track::get_hit_detIDs() const
{
	vector<int> result(m_hits.size());

	for(size_t i = 0; i < result.size(); ++i)
	{
		result[i] = m_hits[i].first;
	}

	return result;
}




/*
 * When the seed track is created with direcion and momentum only, hits can be added from the pattern reco
 *
 * @brief Adding hits from pattern reco
 *
 * @param detectorID The detector ID of the hit
 * @param driftradius The drift radius of the hit
 *
 * @author Daniel Bick
 * @date 12.06.2020
 * @version 1.0
 */
void GBL_seed_track::add_hit(int detectorID, double driftradius)
{
  std::pair<int, double> hit(detectorID, driftradius);
  m_hits.push_back(hit);
}

/*
 * Retreive the point of closes approach of the track to any wire coordinates
 *
 * @brief Get point of closest approach on track to wire
 *
 * @param vtop Top coordinate of the wire
 * @param vtop Bottom coordinate of the wire
 *
 * @author Daniel Bick, derived from Stefans Code in MillepedeCaller
 * @date 17.06.2020
 * @version 1.0
 */
TVector3 GBL_seed_track::PCA_track(TVector3 vtop, TVector3 vbot){
  
  TVector3 wire_dir = vtop - vbot;
  
  //construct a helper plane that contains one of the straights and is parallel to the other one
  TVector3 plane_pos = m_position - vbot;
  TVector3 plane_dir_1(m_direction);
  TVector3 plane_dir_2(-1 * wire_dir);
  
  //Construct components of equation system M * x = c where M is the coefficient matrix, x the solution and c the const_vector below
  TVectorD const_vector(2);
  TMatrixD coeff_matrix(2,2);
  
  const_vector[0] = -(plane_pos.Dot(m_direction));
  const_vector[1] = -(plane_pos.Dot(wire_dir));
  
  coeff_matrix[0][0] = plane_dir_1.Dot(m_direction);
  coeff_matrix[0][1] = plane_dir_2.Dot(m_direction);
  coeff_matrix[1][0] = plane_dir_1.Dot(wire_dir);
  coeff_matrix[1][1] = plane_dir_2.Dot(wire_dir);
  
  TDecompLU solvable_matrix(coeff_matrix);
  TVectorD result(const_vector);
  int rc = solvable_matrix.Solve(result);
  
  return (m_position + result[0] * m_direction);

  //in case PCA_wire is needed: return would be  (vbot + result[1] * wire_dir);
  
}


std::pair<TVector3,TVector3> GBL_seed_track::PCA(TVector3 vtop, TVector3 vbot){
  TVector3 wire_dir = vtop - vbot;
  
  //construct a helper plane that contains one of the straights and is parallel to the other one
  TVector3 plane_pos = m_position - vbot;
  TVector3 plane_dir_1(m_direction);
  TVector3 plane_dir_2(-1 * wire_dir);
  
  //Construct components of equation system M * x = c where M is the coefficient matrix, x the solution and c the const_vector below
  TVectorD const_vector(2);
  TMatrixD coeff_matrix(2,2);
  
  const_vector[0] = -(plane_pos.Dot(m_direction));
  const_vector[1] = -(plane_pos.Dot(wire_dir));
  
  coeff_matrix[0][0] = plane_dir_1.Dot(m_direction);
  coeff_matrix[0][1] = plane_dir_2.Dot(m_direction);
  coeff_matrix[1][0] = plane_dir_1.Dot(wire_dir);
  coeff_matrix[1][1] = plane_dir_2.Dot(wire_dir);
  
  TDecompLU solvable_matrix(coeff_matrix);
  TVectorD result(const_vector);
  int rc = solvable_matrix.Solve(result);
  
  TVector3 track = (m_position + result[0] * m_direction);
  TVector3 wire = (vbot + result[1] * wire_dir);

  std::pair<TVector3,TVector3> PCA_points(track, wire);
  
  return PCA_points;
  
}
