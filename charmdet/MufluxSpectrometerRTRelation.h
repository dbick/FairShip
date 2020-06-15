#ifndef MUFLUXSPECTROMETERRTREL_H
#define MUFLUXSPECTROMETERRTREL_H

#include "TH1D.h"

class MufluxSpectrometerRTRelation : public TH1D {

 public:
  
  MufluxSpectrometerRTRelation();
  MufluxSpectrometerRTRelation(const char *name,const char *title,Int_t nbins,Double_t xlow,Double_t xup);
  MufluxSpectrometerRTRelation(const TH1D& HDriftTimes);

  //MufluxSpectrometerRTRelation(const MufluxSpectrometerRTRelation &point);
  //MufluxSpectrometerRTRelation operator=(const MufluxSpectrometerRTRelation &point);
  
  //virtual ~MufluxSpectrometerRTRelation();
  
  void InitRTRelation(const TH1D &HDriftTimes);
  Double_t GetRadius(Double_t drifttime);

 private:

  ClassDef(MufluxSpectrometerRTRelation, 2);

};



#endif
