#include "MufluxSpectrometerRTRelation.h"
#include <stdio.h>
#include "iostream"


MufluxSpectrometerRTRelation::MufluxSpectrometerRTRelation() : TH1D(){}

MufluxSpectrometerRTRelation::MufluxSpectrometerRTRelation(const char *name,const char *title,Int_t nbins,Double_t xlow,Double_t xup) : TH1D(name,title,nbins,xlow,xup) {}


void MufluxSpectrometerRTRelation::SetRTRelation(TH1D *theHRTRelation){
  //HRTRelation = theHRTRelation;
}


void MufluxSpectrometerRTRelation::InitRTRelation(TH1D *HDriftTimes){
  const int nbins = HDriftTimes->GetNbinsX();
  const int low = HDriftTimes->GetBinLowEdge(1);
  const int high = HDriftTimes->GetBinLowEdge(nbins+1);

    
  //TH1D * pHRTRelation=HRTRelation;
  // TH1D* pHRTRelation = new TH1D("HRTRelation","RT Relation",nbins,low,high);
  
  Double_t norm=1.815/HDriftTimes->Integral(1,nbins);

  for(int i = 1; i <= nbins; i++){
    //pHRTRelation->SetBinContent(i,norm*HDriftTimes->Integral(1,i));
    SetBinContent(i,norm*HDriftTimes->Integral(1,i));
  }


  std::cout << "test " << GetRadius(100) << std::endl;
  
  //HRTRelation=pHRTRelation;
  
}


Double_t MufluxSpectrometerRTRelation::GetRadius(Double_t drifttime){
  //return -1;

  //return HRTRelation->GetBinContent(HRTRelation->FindBin(drifttime));
  
  return GetBinContent(FindBin(drifttime));
}
