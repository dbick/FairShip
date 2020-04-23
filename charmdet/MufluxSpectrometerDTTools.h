#ifndef MUFLUXSPECTROMETERDTTOOLS_H
#define MUFLUXSPECTROMETERDTTOOLS_H

#include "TTreeReaderArray.h"
#include "MufluxSpectrometerHit.h"
#include "MufluxSpectrometerRTRelation.h"

void DrawDTEvent(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits);
void DrawDTEvent(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits,MufluxSpectrometerRTRelation &RTRel);
void DrawDTEvent(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits,MufluxSpectrometerRTRelation &RTRel, bool rtrel);



#endif
