#include "MufluxSpectrometerRTRelation.h"

MufluxSpectrometerRTRelation::MufluxSpectrometerRTRelation() : TH1D(){}

MufluxSpectrometerRTRelation::MufluxSpectrometerRTRelation(const char *name,const char *title,Int_t nbins,Double_t xlow,Double_t xup) : TH1D(name,title,nbins,xlow,xup) {}

MufluxSpectrometerRTRelation::MufluxSpectrometerRTRelation(const TH1D& HDriftTimes) : TH1D(HDriftTimes){
  Reset();
  SetTitle("RT Relation");
  SetName("RTRel");
  int nbins=GetNbinsX();
  Double_t norm=1.815/HDriftTimes.Integral(1,nbins);
  
  for(int i = 1; i <= nbins; i++){
    SetBinContent(i,norm*HDriftTimes.Integral(1,i));
  }

  
}


void MufluxSpectrometerRTRelation::InitRTRelation(const TH1D &HDriftTimes){
  //const int nbins = HDriftTimes.GetNbinsX();
  //const int low = HDriftTimes.GetBinLowEdge(1);
  //const int high = HDriftTimes.GetBinLowEdge(nbins+1);
  
  int nbins=GetNbinsX();
  
  
  Double_t norm=1.815/HDriftTimes.Integral(1,nbins);

  for(int i = 1; i <= nbins; i++){
    SetBinContent(i,norm*HDriftTimes.Integral(1,i));
  }


  
}

Double_t MufluxSpectrometerRTRelation::GetRadius(Double_t drifttime){
  return GetBinContent(FindBin(drifttime));
}

