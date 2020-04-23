#include "DanielAna.h"
#include <iostream>
#include "MufluxSpectrometerHit.h"
#include "ScintillatorHit.h"
#include "MufluxSpectrometer.h"
#include "MufluxSpectrometerRTRelation.h"

#include "TTreeReaderArray.h"
#include "TVector3.h"
#include "TH1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TLine.h"

void helloDaniel(){
  std::cout << "Hello Daniel" << std::endl;
}


void dtAna(TTree *TData){
  //TData->Draw("Digi_MufluxSpectrometerHits.fdigi");
  TClonesArray    *cDigi_MufluxSpectrometerHits;
  cDigi_MufluxSpectrometerHits = 0;
  TBranch        *b_Digi_MufluxSpectrometerHits; 
  
  
  TData->SetBranchAddress("Digi_MufluxSpectrometerHits", &cDigi_MufluxSpectrometerHits, &b_Digi_MufluxSpectrometerHits);

  TData->GetEntry(0);
  std::cout << b_Digi_MufluxSpectrometerHits->GetEntries() << std::endl; 
  //b_Digi_MufluxSpectrometerHits->Inspect();

  b_Digi_MufluxSpectrometerHits->GetEntry(0);

  TData->Draw("Digi_MufluxSpectrometerHits.fdigi");
  
}


void dtAnaChain(TTreeReader *t){
  
  /*
    TTree* fChain = t->GetTree();
    
    TClonesArray    *cDigi_MufluxSpectrometerHits;
    cDigi_MufluxSpectrometerHits = 0;
    TBranch        *b_Digi_MufluxSpectrometerHits; 
    
    fChain->SetBranchAddress("Digi_MufluxSpectrometerHits", &cDigi_MufluxSpectrometerHits, &b_Digi_MufluxSpectrometerHits);
  */
  
  TTreeReaderArray <MufluxSpectrometerHit> Digi_MufluxSpectrometerHits(*t, "Digi_MufluxSpectrometerHits");
  TTreeReaderArray <MufluxSpectrometerHit> Digi_LateMufluxSpectrometerHits(*t, "Digi_LateMufluxSpectrometerHits");
  TTreeReaderArray <MufluxSpectrometerHit> Digi_Triggers(*t, "Digi_Triggers");
  
  
  TVector3 *vtop = new TVector3();
  TVector3 *vbot = new TVector3();

  int event=0;
  TH1D *HDriftTimes = new TH1D("HDriftTimes","Drift Time Spectrum",3000,-500,2500);
  TH1D *HToT = new TH1D("HToT","Time Over Threshold",600,0,600);

  
  while (t->Next()) {
    // behaves like an iterator
    
    int n =  Digi_MufluxSpectrometerHits.GetSize();
   
    
    int trigger_n =  Digi_Triggers.GetSize();
    


    bool goodevent=true;
    for (int i=0; i<trigger_n; ++i){
      //if(Digi_Triggers[i].GetFlags() !=1)goodevent=false;
      if(event==71832) std::cout << Digi_Triggers[i].GetFlags() << std::endl;
    }

    if(trigger_n!=5)goodevent=false;
    
    if(!goodevent){
      //std::cout << event << " Number of hits: " << n << std::endl;
      //std::cout << event << " Number of trigger hits: " << trigger_n << std::endl;
    }
    
    if(goodevent)
    for (int i = 0; i < n; ++i) {


      
      //std::cout << Digi_MufluxSpectrometerHits[i].GetDetectorID() << " " ;
      MufluxSpectrometerHit* hit = &(Digi_MufluxSpectrometerHits[i]);
      hit->MufluxSpectrometerEndPoints(*vbot,*vtop);
      // std::cout << hit->GetDetectorID() << " " << hit->GetDigi() << " " << hit->GetTimeOverThreshold();
      //std::cout << std::endl << " " << vbot->x() << " " << vbot->y() << " " << vbot->z() << std::endl << " " << vtop->x() << " " << vtop->y() << " " << vtop->z() ;
      //std::cout << std::endl;
      if(hit->GetTimeOverThreshold()>40&&!(hit->GetTimeOverThreshold()>167.1999&&hit->GetTimeOverThreshold()<=167.2)){
	HDriftTimes->Fill(hit->GetDigi());

      }

	HToT->Fill(hit->GetTimeOverThreshold());        
      
      
    }


    //break;
    event++;

  }
  //HDriftTimes->Draw();
  HToT->Draw();

  /*
  const int nbins = HDriftTimes->GetNbinsX();
  const int low = HDriftTimes->GetBinLowEdge(1);
  const int high = HDriftTimes->GetBinLowEdge(nbins+1);
  TH1D *HRTRelation = new TH1D("HRTRelation","RT Relation",nbins,low,high);
  

  Double_t norm=1.815/HDriftTimes->Integral(1,nbins);

  for(int i = 1; i <= nbins; i++){
    HRTRelation->SetBinContent(i,norm*HDriftTimes->Integral(1,i));
  }

  HRTRelation->Draw();
  */
  
  
  //MufluxSpectrometerRTRelation* RTRel;
  //  TH1D *HDriftTimes = new TH1D("HDriftTimes","Drift Time Spectrum",3000,-500,2500);

  //MufluxSpectrometerRTRelation* RTRel = new MufluxSpectrometerRTRelation("RTRel","RT Relation",3000,-500,2500);

 MufluxSpectrometerRTRelation* RTRel = new MufluxSpectrometerRTRelation(*HDriftTimes);
  
 //RTRel->InitRTRelation(*HDriftTimes);
  std::cout << RTRel->GetRadius(100) << std::endl;;
  //RTRel->GetEntries();
    
  
}




TH1D *FilterDTSpectrum(TTreeReader *t){
  
  TTreeReaderArray <MufluxSpectrometerHit> Digi_MufluxSpectrometerHits(*t, "Digi_MufluxSpectrometerHits");
  TTreeReaderArray <MufluxSpectrometerHit> Digi_LateMufluxSpectrometerHits(*t, "Digi_LateMufluxSpectrometerHits");
  TTreeReaderArray <ScintillatorHit> Digi_Triggers(*t, "Digi_Triggers");
  TTreeReaderArray <ScintillatorHit> Digi_MasterTrigger(*t, "Digi_MasterTrigger");
  TTreeReaderArray <ScintillatorHit> Digi_Scintillators(*t, "Digi_Scintillators");
  
  
  int event=0;
  TH1D *HDriftTimes = new TH1D("HDriftTimes","Drift Time Spectrum",3000,-500,2500);
  TH1D *HToT = new TH1D("HToT","Time Over Threshold",600,0,600);
  
  
  while (t->Next()) {
    
    int n =  Digi_MufluxSpectrometerHits.GetSize();
    int late_n =  Digi_LateMufluxSpectrometerHits.GetSize();
    int trigger_n =  Digi_Triggers.GetSize();
    int master_n =  Digi_MasterTrigger.GetSize();
    int scint_n =  Digi_Scintillators.GetSize();
    
    if(event==71832){
    for (int i = 0; i < n; ++i) {
      //std::cout << Digi_MufluxSpectrometerHits[i].GetDetectorID() << " " ;
      MufluxSpectrometerHit* hit = &(Digi_MufluxSpectrometerHits[i]);
      std::cout << "ID: " << hit->GetDetectorID() << " FLAGS " << hit->GetFlags() << " Time " << hit->GetDigi() << " (" << hit->GetTimeOverThreshold() << ") ";
      for(int j=0;j<late_n;j++){
	MufluxSpectrometerHit* latehit = &(Digi_LateMufluxSpectrometerHits[j]);
	if(hit->GetDetectorID()==latehit->GetDetectorID()){
	  std::cout << "   " << latehit->GetDigi() << " (" << latehit->GetTimeOverThreshold() << ") ";
	}
      }
      std::cout << std::endl;
    } 
    for (int i = 0; i < trigger_n; ++i) {
      ScintillatorHit* thit = &(Digi_Triggers[i]);
      std::cout << "Trigger ID: " << thit->GetDetectorID() << " FLAGS " << thit->GetFlags() << " Time " << thit->GetDigi() << " (" << thit->GetTimeOverThreshold() << ") " << std::endl;
    }

    for (int i = 0; i < scint_n; ++i) {
      ScintillatorHit* shit = &(Digi_Scintillators[i]);
      std::cout << "Scint ID: " << shit->GetDetectorID() << " FLAGS " << shit->GetFlags() << " Time " << shit->GetDigi() << " (" << shit->GetTimeOverThreshold() << ") " << std::endl;
    }

    DrawDTEvent(Digi_MufluxSpectrometerHits);   
    }
    //break;
    event++;
  }
  
  
  return HDriftTimes;
}


  
void DrawDTEvent(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits){
 
    
  gStyle->SetLineScalePS(1);
  TCanvas *disp = new TCanvas("disp","Event Display",1600,480);
  disp->Range(0,-120,800,120);

  TCanvas *zoomdisp = new TCanvas("zoomdisp","Event Display Zoom",960,768);
  zoomdisp->Range(0,-72,180,72);
  
  Int_t DetectorID;

  MufluxSpectrometer* mySpectrometer= new MufluxSpectrometer();
  TVector3 *vtop = new TVector3();
  TVector3 *vbot = new TVector3();

  //Draw Detector Geometry

  //stations 1 and 2 have two views with one module each
  for(int station=1;station<=2;station++){
    for(int vnb=0;vnb<2;vnb++){
      for(int pnb=0;pnb<2;pnb++){
	for(int lnb=0;lnb<2;lnb++){
	  for(int wire=1;wire<=12;wire++){
	    DetectorID=station*10000000+vnb*1000000+pnb*100000+lnb*10000+2000+wire;
	    mySpectrometer->TubeEndPoints(DetectorID, *vbot, *vtop);
	    //std::cout << DetectorID << " " << vtop->x() << " " << vtop->y() << " " << vtop->z() << std::endl;
	    double scale=1;
	    if(station+vnb==2)scale=1./cos(60*3.14159/180);//scale stereo views
	    TEllipse *el = new TEllipse((vtop->z()+vbot->z())/2,scale*(vtop->x()+vbot->x())/2,2.,2.);
	    el->Draw();
	    zoomdisp->cd();
	    el->Draw();
	    disp->cd();
	  }
	}
      }
    }
  }


   //stations 3 and 4 have one view with four modules
  for(int station=3;station<=4;station++){
    for(int vnb=0;vnb<1;vnb++){
      for(int pnb=0;pnb<2;pnb++){
	for(int lnb=0;lnb<2;lnb++){
	  for(int wire=1;wire<=48;wire++){
	    DetectorID=station*10000000+vnb*1000000+pnb*100000+lnb*10000+2000+wire;
	    mySpectrometer->TubeEndPoints(DetectorID, *vbot, *vtop);
	    //std::cout << DetectorID << " " << vtop->x() << " " << vtop->y() << " " << vtop->z() << std::endl;
	    TEllipse *el = new TEllipse((vtop->z()+vbot->z())/2,(vtop->x()+vbot->x())/2,2.,2.);
	    el->Draw();
	  }
	}
      }
    }
  }


  //Draw Hits



 

  int n =  Digi_MufluxSpectrometerHits.GetSize();
  for(int i=0;i<n;i++){
    MufluxSpectrometerHit* hit = &(Digi_MufluxSpectrometerHits[i]);
    DetectorID = hit->GetDetectorID();

    int station=DetectorID/10000000;
    int vnb=(DetectorID%10000000)/1000000;
    
    //hit->MufluxSpectrometerEndPoints(*vbot,*vtop);
    //std::cout << "Method 1: " << vtop->x() << " " << vtop->y() << " " << vtop->z() << std::endl;

    mySpectrometer->TubeEndPoints(DetectorID, *vbot, *vtop);

    //std::cout << "Method 2: " << vtop->x() << " " << vtop->y() << " " << vtop->z() << std::endl;
    double scale=1;
    if(station+vnb==2)scale=1./cos(60*3.14159/180);//scale stereo views
    TEllipse *hitel = new TEllipse((vtop->z()+vbot->z())/2,scale*(vtop->x()+vbot->x())/2,2.,2.);
    hitel->SetFillColor(2);
    if(station+vnb==2)hitel->SetFillColor(3+vnb);
    hitel->Draw();
    zoomdisp->cd();
    hitel->Draw();
    disp->cd();

    
    if(station<3){
      TLine *tube = new TLine(
			      -vtop->x()+400,
			      vtop->y(),
			      -vbot->x()+400,
			      vbot->y()
			      );
      tube->SetLineColor(2);
      if(station+vnb==2)tube->SetLineColor(3+vnb);
      tube->Draw();
    }
  }
  
}
