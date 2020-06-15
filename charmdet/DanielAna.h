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
#include "GBLseedtrack.h"


void helloDaniel();
void dtAna(TTree *TData);
void dtAnaChain(TTreeReader *t);
void dtAnaChain(TTreeReader *t, int event);
void dtPatAna(TTreeReader *t);

void dtPatSeed(TTreeReader *t);

TH1D *FilterDTSpectrum(TTreeReader *&t);

void PrintSpectrometerHit(MufluxSpectrometerHit &hit);
void PrintScintillatorHit(ScintillatorHit &hit);

//tangent2d SimplePattern2D(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits, MufluxSpectrometerRTRelation &RTRel, int station, int view);

void checkAlignment();
void drawModule(Int_t module);
void test();
#endif
