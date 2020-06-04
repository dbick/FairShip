/*
 * GBLseedtrack.cpp
 *
 *  Created on: 29.04.2020
 *      Author: bieschke
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
