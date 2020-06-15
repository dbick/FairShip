/*
 * GBLseedtrack.h
 *
 *  Created on: 29.04.2020
 *      Author: bieschke
 */

#include "Track.h"
#include "TVector3.h"
#include <vector>
#ifndef CHARMDET_GBLSEEDTRACK_H_
#define CHARMDET_GBLSEEDTRACK_H_

class GBL_seed_track
{
public:
	GBL_seed_track(const genfit::Track& track);
	GBL_seed_track(TVector3 position, TVector3 direction);
	virtual ~GBL_seed_track();

	const size_t get_number_of_hits() const;
	const TVector3& get_position() const;
	const TVector3& get_direction() const;

	const std::vector<std::pair<int,double>>& get_hits() const;

	const std::vector<int> get_hit_detIDs() const;

	void add_hit(int detectorID, double driftradius);

	
private:
	TVector3 m_position;
	TVector3 m_direction;

	std::vector<std::pair<int,double>> m_hits;
};

#endif /* CHARMDET_GBLSEEDTRACK_H_ */
