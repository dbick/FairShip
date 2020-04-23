#include "DanielAna.h"
#include <iostream>
#include "MufluxSpectrometerHit.h"
#include "ScintillatorHit.h"
#include "MufluxSpectrometer.h"
#include "MufluxSpectrometerRTRelation.h"
#include "MufluxSpectrometerDTTools.h"

#include "TTreeReaderArray.h"
#include "TVector3.h"
#include "TH1.h"
#include "TCanvas.h"


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
  TTreeReaderArray <MufluxSpectrometerHit> Digi_MufluxSpectrometerHits(*t, "Digi_MufluxSpectrometerHits");
  TTreeReaderArray <MufluxSpectrometerHit> Digi_LateMufluxSpectrometerHits(*t, "Digi_LateMufluxSpectrometerHits");
  TTreeReaderArray <MufluxSpectrometerHit> Digi_Triggers(*t, "Digi_Triggers");
  */

  /*
  
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
  */
  
  
  //MufluxSpectrometerRTRelation* RTRel;
  //  TH1D *HDriftTimes = new TH1D("HDriftTimes","Drift Time Spectrum",3000,-500,2500);

  //MufluxSpectrometerRTRelation* RTRel = new MufluxSpectrometerRTRelation("RTRel","RT Relation",3000,-500,2500);

  TH1D *HDriftTimes = FilterDTSpectrum(t);

  TCanvas *c1 = new TCanvas();
  HDriftTimes->Draw();
  
  MufluxSpectrometerRTRelation *RTRel = new MufluxSpectrometerRTRelation(*HDriftTimes);

  TCanvas *c2 = new TCanvas();
  RTRel->Draw();
  //RTRel->InitRTRelation(*HDriftTimes);
  //RTRel->GetEntries();

  /*
  t->Restart();
  //t->Next();
  //t->SetEntry(71832);
  t->SetEntry(2);
 
  DrawDTEvent(Digi_MufluxSpectrometerHits,*RTRel);
  //DrawDTEvent(Digi_MufluxSpectrometerHits);   
  */
  
}




TH1D *FilterDTSpectrum(TTreeReader *t){
  
  
  TTreeReaderArray <MufluxSpectrometerHit> Digi_MufluxSpectrometerHits(*t, "Digi_MufluxSpectrometerHits");
  TTreeReaderArray <MufluxSpectrometerHit> Digi_LateMufluxSpectrometerHits(*t, "Digi_LateMufluxSpectrometerHits");
  TTreeReaderArray <ScintillatorHit> Digi_Triggers(*t, "Digi_Triggers");
  TTreeReaderArray <ScintillatorHit> Digi_MasterTrigger(*t, "Digi_MasterTrigger");
  TTreeReaderArray <ScintillatorHit> Digi_Scintillators(*t, "Digi_Scintillators");
  
  
  TH1D *HDriftTimes = new TH1D("HDriftTimes","Drift Time Spectrum",3000,-500,2500);
  TH1D *HToT = new TH1D("HToT","Time Over Threshold",600,0,600);

  t->Restart();
  while (t->Next()) {
    
    int n =  Digi_MufluxSpectrometerHits.GetSize();
    int late_n =  Digi_LateMufluxSpectrometerHits.GetSize();
    int trigger_n =  Digi_Triggers.GetSize();
    int master_n =  Digi_MasterTrigger.GetSize();
    int scint_n =  Digi_Scintillators.GetSize();

    //first we check if all triggers were good

    bool triggersGood = true;
    for (int i = 0; i < trigger_n; ++i) {
      ScintillatorHit* thit = &(Digi_Triggers[i]);
      if(thit->GetFlags()!=1){
	triggersGood = false;
	continue;
      }
    }

    //ignore event if trigger is bad
    if(!triggersGood) continue;


    //now we check all hits, including late hits
    
    for (int i = 0; i < n; ++i) {
      MufluxSpectrometerHit* hit = &(Digi_MufluxSpectrometerHits[i]);

      bool hitGood=true;
      
      for(int j=0;j<late_n;j++){
	MufluxSpectrometerHit* latehit = &(Digi_LateMufluxSpectrometerHits[j]);
	if(hit->GetDetectorID()==latehit->GetDetectorID()){
	  if(latehit->GetFlags()!=1){
	    hitGood=false;
	    continue;
	  }
	}

	//we don't need to look at remaining late hits if we have a bad one!
	if(!hitGood)continue;
      }


      if(hitGood&&hit->GetFlags()==1&&hit->GetTimeOverThreshold()>40){
	HDriftTimes->Fill(hit->GetDigi());
	HToT->Fill(hit->GetTimeOverThreshold());
      }
      
    } 
    
  }


  return HDriftTimes;
}


void PrintSpectrometerHit(MufluxSpectrometerHit &hit){
  std::cout << " FLAGS " << hit.GetFlags() << " Time " << hit.GetDigi() << " (" << hit.GetTimeOverThreshold() << ") ";
}
void PrintScintillatorHit(ScintillatorHit &hit){
  std::cout << " FLAGS " << hit.GetFlags() << " Time " << hit.GetDigi() << " (" << hit.GetTimeOverThreshold() << ") " << std::endl;
}
