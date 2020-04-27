#ifndef DANIELANA_H
#define DANIELANA_H

#include "TROOT.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH1D.h"
#include "MufluxSpectrometerHit.h"
#include "ScintillatorHit.h"
#include "MufluxSpectrometerRTRelation.h"
#include "MufluxSpectrometerDTTools.h"



void helloDaniel();
void dtAna(TTree *TData);
void dtAnaChain(TTreeReader *t);
TH1D *FilterDTSpectrum(TTreeReader *&t);

void PrintSpectrometerHit(MufluxSpectrometerHit &hit);
void PrintScintillatorHit(ScintillatorHit &hit);

int GetStationB(int DetectorID);
int GetView(int DetectorID);

tangent2d SimplePattern2D(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits, MufluxSpectrometerRTRelation &RTRel, int station, int view);

#endif
