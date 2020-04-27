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
#include "TLatex.h"
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
	    mySpectrometer->TubeEndPoints(DetectorID, *vtop, *vbot);
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

  //label "directions"
  TLatex text;
  text.DrawLatex(275,80, "Jura");
  text.DrawLatex(425,-80, "Sal#grave{e}ve");
    
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
      TLine *tube = new TLine(-vtop->x()+350,vtop->y(),-vbot->x()+350,vbot->y());
      if(station+vnb==2)tube->SetLineColor(3+vnb);
      else tube->SetLineColor(2);
      tube->Draw();
    }
  }
  
}


void DrawDTGeom(TCanvas *& disp){
  
  gStyle->SetLineScalePS(1);
  disp = new TCanvas("disp","Event Display",1600,480);
  disp->Range(0,-120,800,120);
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
	    mySpectrometer->TubeEndPoints(DetectorID, *vtop, *vbot);
	    //std::cout << DetectorID << " " << vtop->x() << " " << vtop->y() << " " << vtop->z() << std::endl;
	    double scale=1;
	    if(station+vnb==2)scale=1./cos(60*3.14159/180);//scale stereo views
	    TEllipse *el = new TEllipse((vtop->z()+vbot->z())/2,scale*(vtop->x()+vbot->x())/2,2.,2.);
	    el->Draw();
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

  //label "directions"
  TLatex text;
  text.DrawLatex(275,80, "Jura");
  text.DrawLatex(425,-80, "Sal#grave{e}ve");

}


void DrawDTGeomT12(TCanvas *& disp){
  
  disp = new TCanvas("zoomdisp12","Event Display zoom T1/2",960,768);
  disp->Range(0,-72,180,72);
  disp->cd();
  
  Double_t DetectorID;
  
  MufluxSpectrometer* mySpectrometer= new MufluxSpectrometer();
  TVector3 *vtop = new TVector3();
  TVector3 *vbot = new TVector3();
  
  for(int station=1;station<=2;station++){
    for(int vnb=0;vnb<2;vnb++){
      for(int pnb=0;pnb<2;pnb++){
	for(int lnb=0;lnb<2;lnb++){
	  for(int wire=1;wire<=12;wire++){
	    DetectorID=station*10000000+vnb*1000000+pnb*100000+lnb*10000+2000+wire;
	    mySpectrometer->TubeEndPoints(DetectorID, *vtop, *vbot);
	    double scale=1;
	    if(station+vnb==2)scale=1./cos(60*3.14159/180);//scale stereo views
	    TEllipse *el = new TEllipse((vtop->z()+vbot->z())/2,scale*(vtop->x()+vbot->x())/2,2.,2.);
	    el->Draw();
	  }
	}
      }
    }
  }
} 


void DrawDTGeomT34(TCanvas *& disp){

  disp = new TCanvas("zoomdisp34","Event Display zoom T3/4",960,960);
  disp->Range(560,-120,800,120);
  disp->cd();

  Double_t DetectorID;
  
  MufluxSpectrometer* mySpectrometer= new MufluxSpectrometer();
  TVector3 *vtop = new TVector3();
  TVector3 *vbot = new TVector3();
  
  for(int station=3;station<=4;station++){
    for(int vnb=0;vnb<1;vnb++){
      for(int pnb=0;pnb<2;pnb++){
	for(int lnb=0;lnb<2;lnb++){
	  for(int wire=1;wire<=48;wire++){
	    DetectorID=station*10000000+vnb*1000000+pnb*100000+lnb*10000+2000+wire;
	    mySpectrometer->TubeEndPoints(DetectorID, *vbot, *vtop);
	    TEllipse *el = new TEllipse((vtop->z()+vbot->z())/2,(vtop->x()+vbot->x())/2,2.,2.);
	    el->Draw();
	  }
	}
      }
    }
  }
}


void DrawDTHits(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits,TCanvas *& disp){
  
  disp->cd();
  
  Int_t DetectorID;
  
  MufluxSpectrometer* mySpectrometer= new MufluxSpectrometer();
  TVector3 *vtop = new TVector3();
  TVector3 *vbot = new TVector3();
    
  //Draw Hits
  
  int n =  Digi_MufluxSpectrometerHits.GetSize();
  for(int i=0;i<n;i++){
    MufluxSpectrometerHit* hit = &(Digi_MufluxSpectrometerHits[i]);
    DetectorID = hit->GetDetectorID();
    double radius=2;
    //if(rtrel) radius = RTRel.GetRadius(hit->GetDigi());
    
    int station=DetectorID/10000000;
    int vnb=(DetectorID%10000000)/1000000;
    
    mySpectrometer->TubeEndPoints(DetectorID, *vbot, *vtop);

    double scale=1;
    if(station+vnb==2)scale=1./cos(60*3.14159/180);//scale stereo views

    TEllipse *hitel = new TEllipse((vtop->z()+vbot->z())/2,scale*(vtop->x()+vbot->x())/2,2.,2.);
    hitel->SetFillColor(2);
    if(station+vnb==2)hitel->SetFillColor(3+vnb);
    hitel->Draw();

    
    //draw xy projection
    if(station<3){
      TLine *tube = new TLine(-vtop->x()+350,vtop->y(),-vbot->x()+350,vbot->y());
      if(station+vnb==2)tube->SetLineColor(3+vnb);
      else tube->SetLineColor(2);
      tube->Draw();
    }
  }
}



void DrawDTHitsRT(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits,MufluxSpectrometerRTRelation &RTRel, TCanvas *& disp){
  
  disp->cd();
  
  Int_t DetectorID;
  
  MufluxSpectrometer* mySpectrometer= new MufluxSpectrometer();
  TVector3 *vtop = new TVector3();
  TVector3 *vbot = new TVector3();
    
  //Draw Hits
  
  int n =  Digi_MufluxSpectrometerHits.GetSize();
  for(int i=0;i<n;i++){
    MufluxSpectrometerHit* hit = &(Digi_MufluxSpectrometerHits[i]);
    DetectorID = hit->GetDetectorID();
    double radius=2;
    radius = RTRel.GetRadius(hit->GetDigi());
    
    int station=DetectorID/10000000;
    int vnb=(DetectorID%10000000)/1000000;
    
    mySpectrometer->TubeEndPoints(DetectorID, *vbot, *vtop);

    double scale=1;
    if(station+vnb==2)scale=1./cos(60*3.14159/180);//scale stereo views

    TEllipse *hitel = new TEllipse((vtop->z()+vbot->z())/2,scale*(vtop->x()+vbot->x())/2,radius*2./1.8,radius*2./1.8);
    hitel->SetFillColor(2);
    if(station+vnb==2)hitel->SetFillColor(3+vnb);
    hitel->Draw();

  }
}




void DrawDTTangent(tangent2d tangent, TCanvas *& disp){
  disp->cd();

  double z1=0;
  double z2=800;
  double x1=(tangent.p-z1*sin(tangent.alpha))/cos(tangent.alpha);
  double x2=(tangent.p-z2*sin(tangent.alpha))/cos(tangent.alpha);

  TLine *track = new TLine(z1,x1,z2,x2);
  track->Draw();
}

void DrawDTStereoTangent(tangent2d tangent, TCanvas *& disp, int station){
  disp->cd();

  double z1=0;
  double z2=800;
  int color=0;
  if(station==1){
    z1=30;
    z2=80;
    color=2;
  }
  if(station==2){
    z1=65;
    z2=115;
    color=3;
  }
  double x1=(tangent.p-z1*sin(tangent.alpha))/cos(tangent.alpha);
  double x2=(tangent.p-z2*sin(tangent.alpha))/cos(tangent.alpha);
  TLine *track = new TLine(z1,x1,z2,x2);
  track->SetLineColor(color);
  track->Draw();
}

