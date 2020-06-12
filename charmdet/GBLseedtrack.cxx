/*
 * GBLseedtrack.cpp
 *
 *  Created on: 29.04.2020
 *      Author: bieschke
 *  Modified on: 12.06.2020
 *      Author: Daniel Bick
 */

#include "GBLseedtrack.h"

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
 * 
 */
void GBL_seed_track::add_hit(int detectorID, double driftradius)
{
  std::pair<int, double> hit(detectorID, driftradius);
  m_hits.push_back(hit);
}
