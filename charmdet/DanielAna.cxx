#include "DanielAna.h"
#include <iostream>
#include "math.h"
#include "MufluxSpectrometerHit.h"
#include "ScintillatorHit.h"
#include "MufluxSpectrometer.h"
#include "MufluxSpectrometerRTRelation.h"
#include "MufluxSpectrometerDTTools.h"
#include "MufluxSpectrometerDTSurvey.h"
#include "MufluxSpectrometerDTPatRec.h"
#include "MillepedeCaller.h"
#include "GblTrajectory.h"

#include "TTreeReaderArray.h"
#include "TVector3.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TMarker.h"
#include "TTree.h"
#include "TFile.h"

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
  dtAnaChain(t, 0);
}

void dtAnaChain(TTreeReader *t, int event){

  
  MufluxSpectrometerDTSurvey *surv = new MufluxSpectrometerDTSurvey();
  surv->Init();

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

  //void *temp=t;
  
  TH1D *HDriftTimes = FilterDTSpectrum(t);

  //t=(TTreeReader*)temp;
  
  //TCanvas *c1 = new TCanvas();
  //HDriftTimes->Draw();
  
  MufluxSpectrometerRTRelation *RTRel = new MufluxSpectrometerRTRelation(*HDriftTimes);

  //TCanvas *c2 = new TCanvas();
  //RTRel->Draw();


  //RTRel->InitRTRelation(*HDriftTimes);
  //RTRel->GetEntries();

  //t->Restart();
  TTreeReaderArray <MufluxSpectrometerHit> Digi_MufluxSpectrometerHits(*t, "Digi_MufluxSpectrometerHits");
  
  
  //t->Next();
  //t->SetEntry(71832);
  //t->SetEntry(2);
  t->SetEntry(event);
    
  tangent2d front = SimplePattern2D(Digi_MufluxSpectrometerHits,*RTRel,0b0011,0);

  TCanvas *fdisp;
  DrawDTGeomT12(fdisp);
  DrawDTHitsRT(Digi_MufluxSpectrometerHits,*RTRel,fdisp);
  DrawDTTangent(front,fdisp);

  
  tangent2d back = SimplePattern2D(Digi_MufluxSpectrometerHits,*RTRel,0b1100,0);

  TCanvas *bdisp;
  DrawDTGeomT34(bdisp);
  DrawDTHitsRT(Digi_MufluxSpectrometerHits,*RTRel,bdisp);
  DrawDTTangent(back,bdisp);

  
  tangent2d all = SimplePattern2D(Digi_MufluxSpectrometerHits,*RTRel,0b1111,0);

  TCanvas *disp;
  DrawDTGeom(disp);
  DrawDTHits(Digi_MufluxSpectrometerHits,disp);
  //DrawDTHitsToT(Digi_MufluxSpectrometerHits,disp);
  DrawDTTangent(all,disp);

  DrawDTTangent(front,disp);
  DrawDTTangent(back,disp);

  DrawDTTangent(all,fdisp);
  DrawDTTangent(back,fdisp);
  
 
 
  //SimplePattern2D(Digi_MufluxSpectrometerHits,*RTRel,0b1111,0);

  tangent2d stereo1 = SimplePattern2D(Digi_MufluxSpectrometerHits,*RTRel,0b0011,1);
  tangent2d stereo2 = SimplePattern2D(Digi_MufluxSpectrometerHits,*RTRel,0b0011,2);

  DrawDTStereoTangent(stereo1,fdisp,1);
  DrawDTStereoTangent(stereo2,fdisp,2);
  
    
  std::cout << "Front" << std::endl << front.closehits << " \t" << front.alpha*180./3.14159 << " \t" << front.p << " \t" << front.avres << " \t" << front.avres*front.closehits/(front.closehits-2) <<  std::endl;;

  std::cout << "Back" << std::endl << back.closehits << " \t" << back.alpha*180./3.14159 << " \t" << back.p << " \t" << back.avres << " \t" << back.avres*back.closehits/(back.closehits-2) <<  std::endl;;

    std::cout << "All" << std::endl << all.closehits << " \t" << all.alpha*180./3.14159 << " \t" << all.p << " \t" << all.avres << " \t" << all.avres*all.closehits/(all.closehits-2) <<  std::endl;;

  std::cout << "S1" << std::endl << stereo1.closehits << " \t" << stereo1.alpha*180./3.14159 << " \t" << stereo1.p << " \t" << stereo1.avres << " \t" << stereo1.avres*stereo1.closehits/(stereo1.closehits-2) <<  std::endl;;

  std::cout << "S2" << std::endl << stereo2.closehits << " \t" << stereo2.alpha*180./3.14159 << " \t" << stereo2.p << " \t" << stereo2.avres << " \t" << stereo2.avres*stereo2.closehits/(stereo2.closehits-2) <<  std::endl;;  

  std::cout << all.p/all.alpha << " " << stereo1.p/stereo1.alpha-stereo2.p/stereo2.alpha << std::endl;
  std::cout << "sum tan stereo: " << tan(stereo1.alpha)+tan(stereo2.alpha) << " \t tan vertical: " << tan(all.alpha) << std::endl;

  std::cout << "Alpha in vertical modules " << all.alpha*180./TMath::Pi() << std::endl;
  std::cout << "Expected Alpha from Stereo Moudles  " << atan(-tan(stereo1.alpha)-tan(stereo2.alpha))*180./TMath::Pi() << std::endl;


  /*
  double beta=-60./180*TMath::Pi();

  beta=surv->DTSurveyStereoAngle(11002001);
  
  std::cout << "Stereo 1 -- y0 " << (-cos(beta)*all.p/cos(all.alpha)+stereo1.p/cos(stereo1.alpha))/sin(beta) << "\t\t delta_y " << (-tan(stereo1.alpha)+tan(all.alpha)*cos(beta))/sin(beta) << std::endl;

  //beta=-60./180*TMath::Pi();

  beta=surv->DTSurveyStereoAngle(20002001);
  std::cout << "Stereo 2 -- y0 " << (-cos(beta)*all.p/cos(all.alpha)+stereo2.p/cos(stereo2.alpha))/sin(beta) << "\t\t delta_y " << (-tan(stereo2.alpha)+tan(all.alpha)*cos(beta))/sin(beta) << std::endl;
  */


  Double_t beta1=surv->DTSurveyStereoAngle(11002001);
  Double_t beta2=surv->DTSurveyStereoAngle(20002001);

  Double_t yslope = YslopeFromProjections(stereo1.alpha, beta1, stereo2.alpha, beta2);
  Double_t y0 = Y0FromProjections(stereo1.alpha, stereo1.p, beta1, stereo2.alpha, stereo2.p, beta2);

  Double_t yslope1 = YslopeFromProjections(stereo1.alpha, beta1, front.alpha, 0);
  Double_t y01 = Y0FromProjections(stereo1.alpha, stereo1.p, beta1, front.alpha, front.p, 0);

  Double_t yslope2 = YslopeFromProjections(stereo2.alpha, beta2, front.alpha, 0);
  Double_t y02 = Y0FromProjections(stereo2.alpha, stereo2.p, beta2, front.alpha, front.p, 0);

  std::cout << "Slope:\t" << yslope << " \t" << yslope1 << " \t" << yslope2 << std::endl;
  std::cout << "Y0:\t" << y0 << " \t" << y01 << " \t" << y02 << std::endl;
  
  
}




void dtPatAna(TTreeReader *t){

  
  MufluxSpectrometerDTSurvey *surv = new MufluxSpectrometerDTSurvey();
  surv->Init();
  Double_t beta1=surv->DTSurveyStereoAngle(11002001);
  Double_t beta2=surv->DTSurveyStereoAngle(20002001);
  
  TH1D *HDriftTimes = FilterDTSpectrum(t);

  MufluxSpectrometerRTRelation *RTRel = new MufluxSpectrometerRTRelation(*HDriftTimes);


  TTreeReaderArray <MufluxSpectrometerHit> Digi_MufluxSpectrometerHits(*t, "Digi_MufluxSpectrometerHits");
  
  TFile *fout = new TFile("pattern.root","recreate");
  TTree *TPat = new TTree("TPat","Test for pattern reco");

 
  tangent2d front;
  tangent2d back;
  tangent2d all;

  tangent2d stereo1;
  tangent2d stereo2;

 

  Double_t yslope,yslope1,yslope2,xslope12,xslope;
  Double_t y0,y01,y02,x012,x0;

  Double_t yslope1_front,yslope2_front,xslope_front,xslope_back;
  Double_t y01_front,y02_front,x0_front,x0_back;

  int event=0;
  
  TPat->Branch("falpha",&front.alpha,"falpha/D");
  TPat->Branch("fp",&front.p,"fp/D");
  TPat->Branch("balpha",&back.alpha,"balpha/D");
  TPat->Branch("bp",&back.p,"bp/D");
  TPat->Branch("alpha",&all.alpha,"alpha/D");
  TPat->Branch("p",&all.p,"p/D");
  TPat->Branch("s1alpha",&stereo1.alpha,"s1alpha/D");
  TPat->Branch("s1p",&stereo1.p,"s1p/D");
  TPat->Branch("s2alpha",&stereo2.alpha,"s2alpha/D");
  TPat->Branch("s2p",&stereo2.p,"s2p/D");
  TPat->Branch("event",&event,"event/I");
  TPat->Branch("yslope",&yslope,"yslope/D");
  TPat->Branch("yslope1",&yslope1,"yslope1/D");
  TPat->Branch("yslope2",&yslope2,"yslope2/D");
  TPat->Branch("xslope12",&xslope12,"xslope12/D");
  TPat->Branch("xslope",&xslope,"xslope/D");
  TPat->Branch("y0",&y0,"y0/D");
  TPat->Branch("y01",&y01,"y01/D");
  TPat->Branch("y02",&y02,"y02/D");
  TPat->Branch("x012",&x012,"x012/D");
  TPat->Branch("x0",&x0,"x0/D");

  TPat->Branch("yslope1_front",&yslope1_front,"yslope1_front/D");
  TPat->Branch("yslope2_front",&yslope2_front,"yslope2_front/D");
  TPat->Branch("xslope_front",&xslope_front,"xslope_front/D");
  TPat->Branch("xslope_back",&xslope_back,"xslope_back/D");
  TPat->Branch("y01_front",&y01_front,"y01_front/D");
  TPat->Branch("y02_front",&y02_front,"y02_front/D");
  TPat->Branch("x0_front",&x0_front,"x0_front/D");
  TPat->Branch("x0_back",&x0_back,"x0_back/D");
  
  
  //t->Next();
  //t->SetEntry(71832);
  //t->SetEntry(2);
  //t->SetEntry(event);

  while (t->Next()) {
  
  front = SimplePattern2D(Digi_MufluxSpectrometerHits,*RTRel,0b0011,0);
  back = SimplePattern2D(Digi_MufluxSpectrometerHits,*RTRel,0b1100,0);
  all = SimplePattern2D(Digi_MufluxSpectrometerHits,*RTRel,0b1111,0);

  stereo1 = SimplePattern2D(Digi_MufluxSpectrometerHits,*RTRel,0b0011,1);
  stereo2 = SimplePattern2D(Digi_MufluxSpectrometerHits,*RTRel,0b0011,2);

  if(stereo1.p>=0&&stereo2.p>=0){
    yslope = YslopeFromProjections(stereo1.alpha, beta1, stereo2.alpha, beta2);
    y0 = Y0FromProjections(stereo1.alpha, stereo1.p, beta1, stereo2.alpha, stereo2.p, beta2);
    xslope12 = XslopeFromProjections(stereo1.alpha, beta1, stereo2.alpha, beta2);
    x012 = X0FromProjections(stereo1.alpha, stereo1.p, beta1, stereo2.alpha, stereo2.p, beta2);
  }
  else{
    yslope=0;
    y0=0;
    xslope12=0;
    x012=0;
  }

  if(stereo1.p>=0&&all.p>=0){
    yslope1 = YslopeFromProjections(stereo1.alpha, beta1, all.alpha, 0);
    y01 = Y0FromProjections(stereo1.alpha, stereo1.p, beta1, all.alpha, all.p, 0);
  }
  else{
    yslope1=0;
    y01=0;
  }


  if(stereo2.p>=0&&all.p>=0){
    yslope2 = YslopeFromProjections(stereo2.alpha, beta2, all.alpha, 0);
    y02 = Y0FromProjections(stereo2.alpha, stereo2.p, beta2, all.alpha, all.p, 0);}
  else{
    yslope2=0;
    y02=0;
  }

  
  if(stereo1.p>=0&&front.p>=0){
    yslope1_front = YslopeFromProjections(stereo1.alpha, beta1, front.alpha, 0);
    y01_front = Y0FromProjections(stereo1.alpha, stereo1.p, beta1, front.alpha, front.p, 0);
  }
  else{
    yslope1_front=0;
    y01_front=0;
  }
  
  
  if(stereo2.p>=0&&front.p>=0){
    yslope2_front = YslopeFromProjections(stereo2.alpha, beta2, front.alpha, 0);
    y02_front = Y0FromProjections(stereo2.alpha, stereo2.p, beta2, front.alpha, front.p, 0);}
  else{
    yslope2_front=0;
    y02_front=0;
  }


  if(all.p>=0){
    xslope=-tan(all.alpha);
    x0=all.p/cos(all.alpha);
  }
  else{
    xslope=0;
    x0=0;
  }

  if(front.p>=0){
    xslope_front=-tan(front.alpha);
    x0_front=front.p/cos(front.alpha);
  }
  else{
    xslope_front=0;
    x0_front=0;
  }

  
   if(back.p>=0){
    xslope_back=-tan(back.alpha);
    x0_back=back.p/cos(back.alpha);
  }
  else{
    xslope_back=0;
    x0_back=0;
  }
  
  
  TPat->Fill();
  event++;
  
  //if(event>=10000)break;

  }


  //create quality plots



   TCanvas *c5 = new TCanvas("c5","y0",800,800);

  TPat->Draw("y01:y02>>h3(200,-100,100,200,-100,100)","p>=0&&s1p>=0&&s2p>=0","colz");
  
  TH1 *h3=TPat->GetHistogram();
  h3->SetTitle("Comparison of y_{0} calculated individually for each stereo module");
  h3->GetXaxis()->SetTitle("y_{0} [cm] from T2a and vertical modules");
  h3->GetYaxis()->SetTitle("y_{0} [cm] from T1b and vertical modules");
  h3->GetYaxis()->SetTitleOffset(1.2); 
  

  TCanvas *c5b = new TCanvas("c5b","y0 stereo stereo",800,800);
    
  TPat->Draw("y02:y0>>h3b(200,-100,100,200,-100,100)","p>=0&&s1p>=0&&s2p>=0","colz");
  
  TH1 *h3b=TPat->GetHistogram();
  h3b->SetTitle("Comparison of y_{0} T2a-vertical vs. T1b-T2a");
  h3b->GetXaxis()->SetTitle("y_{0} [cm] from T2a and T1b");
  h3b->GetYaxis()->SetTitle("y_{0} [cm] from T2a and vertical modules");
  h3b->GetYaxis()->SetTitleOffset(1.2);


  TCanvas *c5c = new TCanvas("c5c","y0 stereo stereo",800,800);
  
  
  TPat->Draw("y01:y0>>h3c(200,-100,100,200,-100,100)","p>=0&&s1p>=0&&s2p>=0","colz");
  
  TH1 *h3c=TPat->GetHistogram();
  h3c->SetTitle("Comparison of y_{0} T1b-vertical vs. T1b-T2a");
  h3c->GetXaxis()->SetTitle("y_{0} [cm] from T2a and T1b");
  h3c->GetYaxis()->SetTitle("y_{0} [cm] from T1b and vertical modules");
  h3c->GetYaxis()->SetTitleOffset(1.2);
  
  
  TCanvas *c6 = new TCanvas("c6","#Delta y",800,800);
 
  TPat->Draw("yslope1:yslope2>>h4(100,-.3,.3,100,-.3,.3)","p>=0&&s1p>=0&&s2p>=0","colz");
  TH1 *h4=TPat->GetHistogram();

  h4->SetTitle("Comparison of #Delta y calculated individually for each stereo module");
  h4->GetXaxis()->SetTitle("#Delta y from T2a and vertical modules");
  h4->GetYaxis()->SetTitle("#Delta y from T1b and vertical modules");
  h4->GetYaxis()->SetTitleOffset(1.4);




  
  TCanvas *c6b = new TCanvas("c6b","#Delta y stereo",800,800);

  TPat->Draw("yslope1:yslope>>h4b(100,-.3,.3,100,-.3,.3)","p>=0&&s1p>=0&&s2p>=0","colz");
  TH1 *h4b=TPat->GetHistogram();

  h4b->SetTitle("Comparison of #Delta y calculated from  T1b-vertical vs. T1b-T2a");
  h4b->GetXaxis()->SetTitle("#Delta y from T1b and T2a");
  h4b->GetYaxis()->SetTitle("#Delta y from T1b and vertical modules");
  h4b->GetYaxis()->SetTitleOffset(1.4);

  
  TCanvas *c6c = new TCanvas("c6c","#Delta y stereo",800,800);

  TPat->Draw("yslope2:yslope>>h4c(100,-.3,.3,100,-.3,.3)","p>=0&&s1p>=0&&s2p>=0","colz");
  TH1 *h4c=TPat->GetHistogram();

  h4c->SetTitle("Comparison of #Delta y calculated from  T2a-vertical vs. T1b-T2a");
  h4c->GetXaxis()->SetTitle("#Delta y from T1b and T2a");
  h4c->GetYaxis()->SetTitle("#Delta y from T2a and vertical modules");
  h4c->GetYaxis()->SetTitleOffset(1.4);

  TCanvas *c7 = new TCanvas("c7","#Delta x",800,800);
  
  TPat->Draw("xslope12:-tan(alpha)>>h5(100,-.3,.3,100,-.3,.3)","p>=0&&s1p>=0&&s2p>=0","colz");
  TH1 *h5=TPat->GetHistogram();

  h5->SetTitle("Comparison of #Delta x calculated from vertical modules vs. T1b-T2a");
  h5->GetXaxis()->SetTitle("#Delta x from T1b and T2a");
  h5->GetYaxis()->SetTitle("#Delta x from vertical modules");
  h5->GetYaxis()->SetTitleOffset(1.4);



  TCanvas *c8 = new TCanvas("c8","x0",800,800);
  
  TPat->Draw("x012:p/cos(alpha)>>h6(200,-100,100,200,-100,100)","p>=0&&s1p>=0&&s2p>=0","colz");
  
  TH1 *h6=TPat->GetHistogram();
  h6->SetTitle("Comparison of x_{0} vertical vs. T1b-T2a");
  h6->GetXaxis()->SetTitle("x_{0} [cm] from T2a and T1b");
  h6->GetYaxis()->SetTitle("x_{0} [cm] vertical modules");
  h6->GetYaxis()->SetTitleOffset(1.2);
  
  
  h3->Write();
  h3b->Write();
  h3c->Write();
  h4->Write();
  h4b->Write();
  h4c->Write();
  h5->Write();
  h6->Write();

  
  fout->Write();
  fout->Close();
}



void dtPatSeed(TTreeReader *t){

  
  MufluxSpectrometerDTSurvey *surv = new MufluxSpectrometerDTSurvey();
  surv->Init();
  Double_t beta1=surv->DTSurveyStereoAngle(11002001);
  Double_t beta2=surv->DTSurveyStereoAngle(20002001);
  
  TH1D *HDriftTimes = FilterDTSpectrum(t);

  MufluxSpectrometerRTRelation *RTRel = new MufluxSpectrometerRTRelation(*HDriftTimes);


  MillepedeCaller *mpc= new MillepedeCaller("bla.milleout");
  
  TTreeReaderArray <MufluxSpectrometerHit> Digi_MufluxSpectrometerHits(*t, "Digi_MufluxSpectrometerHits");
  //t->SetEntry(15);
  //{
  while (t->Next()) {
    GBL_seed_track *seed = seedtrack(Digi_MufluxSpectrometerHits,*RTRel);
    if(seed!=nullptr){
      gbl::GblTrajectory trajectory = mpc->perform_GBL_refit(*seed, 5e-2);
      
      if(trajectory.isValid()){
	TVectorD localPar(5);
	TMatrixDSym localCov(5,5);
	
	std::vector<int> hits = seed->get_hit_detIDs();

	TVector3 vtop;
	TVector3 vbot;

	TVector3 pos=seed->get_position();
	TVector3 mom=seed->get_direction();
	

	TVectorT<double> parameters(5);
	TMatrixTSym<double> covariance(5, 5);

	for (unsigned int i = 1; i <= trajectory.getNumPoints(); ++i){

	  Int_t detID=hits[i-1];
	  surv->TubeEndPointsSurvey(detID, vtop, vbot);

	  TVector3 PCA_track=seed->PCA_track(vbot,vtop);
	  
	  trajectory.getResults(i, parameters, covariance);


	  TVector3 fitpos(parameters[3],parameters[4],0);
	  fitpos+=PCA_track;
	  TVector3 fitmom(parameters[1],parameters[2],0);
	  fitmom+=mom;
	  
	  Double_t zfactor=-(fitpos[2])/(fitmom[2]);

	  TVector3 reference=fitpos+zfactor*fitmom;

	  std::cout << "Hit " << i << " \t" << detID << " \t" << reference[0] << " \t" << reference[1] << " \t" << reference[2] << std::endl;	
	  /*
	  std::cout << "Hit: " << i << std::endl;
	  for (unsigned int j = 0; j < parameters.GetNrows(); ++j)
	    {
	      std::cout << "Parameter: " << j << " value: " << parameters[j] << std::endl;
	    }
	  */
	}
		
      }
      
    
      delete seed;
    }
    
  }
}



TH1D *FilterDTSpectrum(TTreeReader *&t){
  
  
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

  //leave t in fresh state
  t->Restart();
  return HDriftTimes;
}


void PrintSpectrometerHit(MufluxSpectrometerHit &hit){
  std::cout << " FLAGS " << hit.GetFlags() << " Time " << hit.GetDigi() << " (" << hit.GetTimeOverThreshold() << ") ";
}
void PrintScintillatorHit(ScintillatorHit &hit){
  std::cout << " FLAGS " << hit.GetFlags() << " Time " << hit.GetDigi() << " (" << hit.GetTimeOverThreshold() << ") " << std::endl;
}


void checkAlignment(){
  
  
  MufluxSpectrometerDTSurvey *surv = new MufluxSpectrometerDTSurvey();
  surv->Init();
  MufluxSpectrometer* mySpectrometer= new MufluxSpectrometer();
  TVector3 *vtop = new TVector3();
  TVector3 *vbot = new TVector3();


  TVector3 *stop = new TVector3();
  TVector3 *sbot = new TVector3();

  bool table=false;
  
  for(int module=0;module<12;module++){
    std::cout << "Module " << module << std::endl;
    int station=4;
    if(module<8) station=3;
    if(module<4) station=2;
    if(module<2) station=1;
    int view=0;
    if(module==1||module==3)view=1;
    for(int pnb=0;pnb<2;pnb++){
      for(int lnb=0;lnb<2;lnb++){
	for(int tube=0;tube<12;tube++){
	  int snb=tube+1;
	  if(module>=4)snb+=(module%4)*12;
	  int DetectorID=station*1e7+view*1e6+pnb*1e5+lnb*1e4+2000+snb;
	  mySpectrometer->TubeEndPoints(DetectorID, *vbot, *vtop);
	  surv->TubeEndPointsSurvey(DetectorID, *stop, *sbot);

	  if(!table){
	    // std::cout << DetectorID << " TOP   " << stop->x()-vtop->x() << "  \t" << stop->y()-vtop->y() << "  \t" << stop->z()-vtop->z() << "  \t\t BOT   " << sbot->x()-vbot->x() << "  \t" << sbot->y()-vbot->y() << "  \t" << sbot->z()-vbot->z() << "\t\t\t\t\t" << vtop->x()-vbot->x() << std::endl;
	    
	  
	  std::cout << DetectorID << "\t TOP \t survey \t FairShip \t Delta \t\t BOT\t survey \t FairShip \t Delta " << std::endl;
	  std::cout << "       x \t\t "  << stop->x() << "   \t " << vtop->x() << "   \t " <<  stop->x()-vtop->x() << "   \t\t "  << sbot->x() << "    \t " << vbot->x() << "    \t " <<  sbot->x()-vbot->x() << std:: endl;
	  std::cout << "       y \t\t "  << stop->y() << "   \t " << vtop->y() << "   \t " <<  stop->y()-vtop->y() << "   \t\t "  << sbot->y() << "    \t " << vbot->y() << "    \t " <<  sbot->y()-vbot->y() << std:: endl;
	  std::cout << "       z \t\t "  << stop->z() << "   \t " << vtop->z() << "   \t " <<  stop->z()-vtop->z() << "   \t\t "  << sbot->z() << "    \t " << vbot->z() << "    \t " <<  sbot->z()-vbot->z() << std:: endl;

	  }
	  
	  
	}
      }
      
      
    }
    
    if(table)std::cout << module << " & " << stop->x()-vtop->x() << " & " << stop->y()-vtop->y() << " & " << stop->z()-vtop->z() << " & " << sbot->x()-vbot->x() << " & " << sbot->y()-vbot->y() << " & " << sbot->z()-vbot->z() << " \\\\" << std::endl;
    else  std::cout << "-------------------------------------------------------------------------------------------------------------------------" << std::endl;
  }
  
}



void drawModule(Int_t module){

  TCanvas *cmod = new TCanvas("cmod","Module",1000,1000);

  TVector3 *vtop=new TVector3();
  TVector3 *vbot=new TVector3();

  Int_t station=module;
  Int_t offset=0;
  if(module>30){
    station=(module/10)*10;
    offset=(module%10)*12;
  }
  
  Int_t DetectorID=station*1e6+2007+offset;

  MufluxSpectrometerDTSurvey *s = new MufluxSpectrometerDTSurvey();
  MufluxSpectrometer *m = new MufluxSpectrometer();

  s->Init();

  s->TubeEndPointsSurvey(DetectorID,*vbot,*vtop);

  Double_t x,y;

  x=-(vbot->x()+vtop->x())/2;
  y=(vbot->y()+vtop->y())/2;
  
  cmod->Range(x-100,y-100,x+100,y+100);
  cmod->cd();

  for(int pnb=0;pnb<2;pnb++){
    for(int lnb=0;lnb<2;lnb++){
      for(int snb=0;snb<12;snb++){
	DetectorID=station*1e6+pnb*1e5+lnb*1e4+2000+offset+snb+1;
	s->TubeEndPointsSurvey(DetectorID,*vbot,*vtop);

	TMarker *t=new TMarker(-vtop->x(),vtop->y(),7);
	t->Draw();
	TMarker *b=new TMarker(-vbot->x(),vbot->y(),7);
	b->Draw();


	m->TubeEndPoints(DetectorID,*vbot,*vtop);

	TMarker *mt=new TMarker(-vtop->x(),vtop->y(),7);
	mt->SetMarkerColor(kRed);
	mt->Draw();

	TMarker *mb=new TMarker(-vbot->x(),vbot->y(),7);
	mb->SetMarkerColor(kRed);
	mb->Draw();
      }
    }
  }

  
}


void test(){
  
  MufluxSpectrometerDTSurvey *surv = new MufluxSpectrometerDTSurvey();
  surv->Init();

  TVector3 *v1=new TVector3();
  TVector3 *v2=new TVector3();
  TubeEndPoints pp=surv->TubeEndPointsSurvey(10012001);
  surv->TubeEndPointsSurvey(10012001,*v1,*v2);

  
  std::cout << pp.top.x() << " " << v1->x() << std::endl;
  

}
  
