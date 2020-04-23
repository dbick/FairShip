#include <iostream>

#include "MufluxSpectrometerDTTools.h"
#include "MufluxSpectrometerHit.h"
#include "ScintillatorHit.h"
#include "MufluxSpectrometer.h"
#include "MufluxSpectrometerRTRelation.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TTreeReaderArray.h"
#include "TVector3.h"


void DrawDTEvent(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits){
  MufluxSpectrometerRTRelation *RTRel;
  DrawDTEvent(Digi_MufluxSpectrometerHits,*RTRel, false);
}

void DrawDTEvent(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits, MufluxSpectrometerRTRelation &RTRel){
  DrawDTEvent(Digi_MufluxSpectrometerHits, RTRel, true);
}

void DrawDTEvent(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits, MufluxSpectrometerRTRelation &RTRel, bool rtrel){
  
  
  gStyle->SetLineScalePS(1);
  TCanvas *disp = new TCanvas("disp","Event Display",1600,480);
  disp->Range(0,-120,800,120);
  
  TCanvas *zoomdisp;
  TCanvas *zoomdispback;
  
  if(rtrel){
    zoomdisp = new TCanvas("zoomdisp","Event Display zoom T1/2",960,768);
    zoomdisp->Range(0,-72,180,72);
    zoomdispback = new TCanvas("zoomdispback","Event Display zoom T3/4",960,960);
    zoomdispback->Range(560,-120,800,120);
  }
  disp->cd();
  
  
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
	    if(rtrel){
	      zoomdisp->cd();
	      el->Draw();
	      disp->cd();
	    }
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
	    if(rtrel){
	      zoomdispback->cd();
	      el->Draw();
	      disp->cd();
	    }
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
    double radius=2;
    if(rtrel) radius = RTRel.GetRadius(hit->GetDigi());
    
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

    //draw hits in zoomed display if we have an rt relation
    if(rtrel){
      std::cout << "Drawing radius " << radius << " for time " << hit->GetDigi() << std::endl;
      
      TEllipse *hitelr = new TEllipse((vtop->z()+vbot->z())/2,scale*(vtop->x()+vbot->x())/2,radius*10./9,radius*10./9);
      if(station+vnb==2)hitelr->SetFillColor(3+vnb);
      else hitelr->SetFillColor(2);
      if(station<3) zoomdisp->cd();
      else zoomdispback->cd();
      hitelr->Draw();
      disp->cd();
    }
    
    //draw xy projection
    if(station<3){
      TLine *tube = new TLine(-vtop->x()+400,vtop->y(),->x()+400,vbot->y());
      if(station+vnb==2)tube->SetLineColor(3+vnb);
      else tube->SetLineColor(2);
      tube->Draw();
    }
  }
  
}
