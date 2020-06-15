#ifndef MUFLUXSPECTROMETERDTTOOLS_H
#define MUFLUXSPECTROMETERDTTOOLS_H

#include "TTreeReaderArray.h"
#include "MufluxSpectrometerHit.h"
#include "MufluxSpectrometerRTRelation.h"
#include "MufluxSpectrometerDTPatRec.h"
#include "TCanvas.h"


void DrawDTEvent(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits);
void DrawDTEvent(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits,MufluxSpectrometerRTRelation &RTRel);
void DrawDTEvent(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits,MufluxSpectrometerRTRelation &RTRel, bool rtrel);


void DrawDTGeom(TCanvas *& disp);
void DrawDTGeomT12(TCanvas *& disp);
void DrawDTGeomT34(TCanvas *& disp);
void DrawDTHits(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits, TCanvas *& disp);
void DrawDTHitsToT(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits,TCanvas *& disp);
void DrawDTHitsRT(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits,MufluxSpectrometerRTRelation &RTRel, TCanvas *& disp);
void DrawDTTangent(tangent2d tangent, TCanvas *& disp);
void DrawDTStereoTangent(tangent2d tangent, TCanvas *& disp, int station);



#endif
