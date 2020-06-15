#ifndef MUFLUXSPECTROMETERDTPATREC_H
#define MUFLUXSPECTROMETERDTPATREC_H

#include "MufluxSpectrometerRTRelation.h"
#include "MufluxSpectrometerDTSurvey.h"
#include "MufluxSpectrometerHit.h"
#include "MufluxSpectrometer.h"
#include "TTreeReaderArray.h"
#include "GBLseedtrack.h"

struct tangent2d {
  double alpha;
  double p;
  int closehits;
  double avres; //average residual
  bool operator<(tangent2d &t){
        return closehits > t.closehits;
  }
} ;


int GetStationB(int DetectorID);
int GetView(int DetectorID);


Double_t YslopeFromProjections(Double_t alpha1, Double_t beta1, Double_t alpha2, Double_t beta2);
Double_t XslopeFromProjections(Double_t alpha1, Double_t beta1, Double_t alpha2, Double_t beta2);
Double_t Y0FromProjections(Double_t alpha1, Double_t p1, Double_t beta1, Double_t alpha2, Double_t p2, Double_t beta2);
Double_t X0FromProjections(Double_t alpha1, Double_t p1, Double_t beta1, Double_t alpha2, Double_t p2, Double_t beta2);


tangent2d SimplePattern2D(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits, MufluxSpectrometerRTRelation &RTRel, int bin_station, int view);

GBL_seed_track *seedtrack(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits, MufluxSpectrometerRTRelation &RTRel);

#endif
