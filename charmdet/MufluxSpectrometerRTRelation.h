#ifndef MUFLUXSPECTROMETERRTREL_H
#define MUFLUXSPECTROMETERRTREL_H

#include "TH1D.h"



//class MufluxSpectrometerRTRelation {
class MufluxSpectrometerRTRelation : public TH1D {
 public:
  MufluxSpectrometerRTRelation();

  MufluxSpectrometerRTRelation(const char *name,const char *title,Int_t nbins,Double_t xlow,Double_t xup);

  void InitRTRelation(const TH1D &HDriftTimes);
  Double_t GetRadius(Double_t drifttime);

  /*
  //MufluxSpectrometerRTRelation(const MufluxSpectrometerRTRelation &point);
  //MufluxSpectrometerRTRelation operator=(const MufluxSpectrometerRTRelation &point);

  
  
  virtual ~MufluxSpectrometerRTRelation();
  //MufluxSpectrometerRTRelation() = default;

  void SetRTRelation(TH1D *theHRTRelation);
  void InitRTRelation(TH1D *HDriftTimes);
  Double_t GetRadius(Double_t drifttime);
  //TH1D *RTRelationHistrogram(){return HRTRelation;};
  
  */  
 private:
   
  //TH1D *HRTRelation;
  ClassDef(MufluxSpectrometerRTRelation, 2);
};



#endif
