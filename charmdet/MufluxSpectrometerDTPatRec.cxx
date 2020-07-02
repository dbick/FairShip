#include "MufluxSpectrometerDTPatRec.h"





int GetStationB(int DetectorID){ //retuns binary station
  return 1<<(DetectorID/10000000-1);
}

int GetView(int DetectorID){
  int station=GetStationB(DetectorID);
  int vnb=(DetectorID%10000000)/1000000;
  if(station==1&&vnb==1)return 1;
  if(station==2&&vnb==0)return 2;
  return 0;
}




Double_t YslopeFromProjections(Double_t alpha1, Double_t stereo1, Double_t alpha2, Double_t stereo2){
  return (tan(alpha2)/cos(stereo2)-tan(alpha1)/cos(stereo1))/(tan(stereo1)-tan(stereo2));
}
Double_t XslopeFromProjections(Double_t alpha1, Double_t stereo1, Double_t alpha2, Double_t stereo2){
  return (tan(alpha2)/sin(stereo2)-tan(alpha1)/sin(stereo1))/(1/tan(stereo1)-1/tan(stereo2));
}
Double_t Y0FromProjections(Double_t alpha1, Double_t p1, Double_t stereo1, Double_t alpha2, Double_t p2, Double_t stereo2){
  return (p1/(cos(stereo1)*cos(alpha1))-p2/(cos(stereo2)*cos(alpha2)))/(tan(stereo1)-tan(stereo2));
}
Double_t X0FromProjections(Double_t alpha1, Double_t p1, Double_t stereo1, Double_t alpha2, Double_t p2, Double_t stereo2){
  return (p1/(sin(stereo1)*cos(alpha1))-p2/(sin(stereo2)*cos(alpha2)))/(1/tan(stereo1)-1/tan(stereo2));
}




tangent2d SimplePattern2D(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits, MufluxSpectrometerRTRelation &RTRel, int bin_station, int view){

  bool b_survey=true;

  MufluxSpectrometer* mySpectrometer= new MufluxSpectrometer();
  MufluxSpectrometerDTSurvey *survey = new MufluxSpectrometerDTSurvey();
  
  TVector3 *vtop = new TVector3();
  TVector3 *vbot = new TVector3();

  Double_t x[2],z[2],r[2];
  MufluxSpectrometerHit* hit;

  std::list<tangent2d> listoftangents;
  
  
  int n =  Digi_MufluxSpectrometerHits.GetSize();
  for(int i=0;i<n;i++){
    hit = &(Digi_MufluxSpectrometerHits[i]);

    //int station=hit->GetDetectorID()/10000000;
    //int vnb=(hit->GetDetectorID()%10000000)/1000000;
    //if(station>2||station+vnb==2)continue; //only T1/2 no stereo

    if(hit->GetTimeOverThreshold()<40)continue;
    
    if(GetView(hit->GetDetectorID())!=view)continue;
    if(GetStationB(hit->GetDetectorID())&bin_station==0)continue;
    
    if(!b_survey) mySpectrometer->TubeEndPoints(hit->GetDetectorID(), *vtop, *vbot);
    else survey->TubeEndPointsSurvey(hit->GetDetectorID(), *vtop, *vbot);

    Double_t stereo_angle=60./180*TMath::Pi();
    
    if(view!=0){
      Double_t angle;
      if(hit->GetDetectorID()/10000000==1){
	if(b_survey) angle = survey->DTSurveyStereoAngle(hit->GetDetectorID());
	else angle=stereo_angle;
	vbot->RotateZ(-angle);
	vtop->RotateZ(-angle);
      }
      if(hit->GetDetectorID()/10000000==2){
	if(b_survey) angle = survey->DTSurveyStereoAngle(hit->GetDetectorID());
	else angle=-stereo_angle;
	vbot->RotateZ(-angle);
	vtop->RotateZ(-angle);
      }
    }
    
    x[0]=(vtop->x()+vbot->x())/2;
    z[0]=(vbot->z()+vbot->z())/2;
    r[0]=RTRel.GetRadius(hit->GetDigi());

    for(int j=i+1;j<n;j++){
      hit = &(Digi_MufluxSpectrometerHits[j]);

      //station=hit->GetDetectorID()/10000000;
      //vnb=(hit->GetDetectorID()%10000000)/1000000;
      //if(station>2||station+vnb==2)continue; //only T1/2 no stereo

      if(GetView(hit->GetDetectorID())!=view)continue;
      if((GetStationB(hit->GetDetectorID())&bin_station)==0)continue;
      
      if(!b_survey) mySpectrometer->TubeEndPoints(hit->GetDetectorID(), *vtop, *vbot);
      else survey->TubeEndPointsSurvey(hit->GetDetectorID(), *vtop, *vbot);

      if(view!=0){
	Double_t angle;
	if(hit->GetDetectorID()/10000000==1){
	  if(b_survey) angle = survey->DTSurveyStereoAngle(hit->GetDetectorID());
	  else angle=stereo_angle;
	  vbot->RotateZ(-angle);
	  vtop->RotateZ(-angle);
	}
	if(hit->GetDetectorID()/10000000==2){
	  if(b_survey) angle = survey->DTSurveyStereoAngle(hit->GetDetectorID());
	  else angle=-stereo_angle;
	  vbot->RotateZ(-angle);
	  vtop->RotateZ(-angle);
	}
      }
      
      
      x[1]=(vtop->x()+vbot->x())/2;
      z[1]=(vbot->z()+vbot->z())/2;
      r[1]=RTRel.GetRadius(hit->GetDigi());

    
      /*
      std::cout << "Comparing two hits" << std::endl;
      std::cout << "\t Hit 0 at " << x[0] << ", " << z[0] << "\t\t" << "radius " << r[0] << std::endl;
      std::cout << "\t Hit 1 at " << x[1] << ", " << z[1] << "\t\t" << "radius " << r[1] << std::endl; 
      */
     

      Double_t deltax=x[0]-x[1];
      Double_t deltaz=z[0]-z[1];
      
      Double_t distance=sqrt(pow(deltax,2)+pow(deltaz,2));
      //if(distance<5)continue; //skip adjacent tubes since high uncertainty due to lever arm 

      //std::cout << "\t Distance is " << distance <<std::endl;
      
      Double_t deltar[2];
      deltar[0]=r[0]-r[1];
      deltar[1]=r[0]+r[1];

      Int_t pm[2];
      pm[0]=1;
      pm[1]=-1;
      
      for(int k=0;k<2;k++){
	for(int l=0;l<2;l++){
	  for(int m=0;m<2;m++){
	    int o=4*k+2*l+m;
	    tangent2d tangent;
	    tangent.alpha=2*atan((deltaz+pm[l]*sqrt(pow(deltaz,2)-pow(deltar[k],2)+pow(deltax,2)))/(pm[m]*deltar[k]+deltax));
	    //m is the "sign" used for r[0] 
	    tangent.p=x[0]*cos(tangent.alpha)+z[0]*sin(tangent.alpha)-pm[m]*r[0];
	    /*
	    std::cout << "\t\t" << alpha_single[o]*180./3.14159 << "   \t" << p_single[o] ;//<< std::endl;
	    p_single[o]=x[1]*cos(alpha_single[o])+z[1]*sin(alpha_single[o])-pm[(m+k)%2]*r[1];
	    std::cout << "   \t" << p_single[o] << std::endl;
	    */
	    if(tangent.p>=0){
	      //std::cout << "\t\t" << tangent.alpha*180./3.14159 << "   \t" << tangent.p << std::endl;

	      tangent.closehits=0;
	      tangent.avres=0;
	      //std::cout << "\t\t\t Distance to other hits: "; 
	      //distance to other hits
	      for(int ii=0;ii<n;ii++){
		hit = &(Digi_MufluxSpectrometerHits[ii]);

		//station=hit->GetDetectorID()/10000000;
		//vnb=(hit->GetDetectorID()%10000000)/1000000;
		//if(station+vnb==2)continue; //no stereo
		if(GetView(hit->GetDetectorID())!=view)continue;
		if((GetStationB(hit->GetDetectorID())&bin_station)==0)continue;

		if(!b_survey) mySpectrometer->TubeEndPoints(hit->GetDetectorID(), *vtop, *vbot);
		else survey->TubeEndPointsSurvey(hit->GetDetectorID(), *vtop, *vbot);
		
		if(view!=0){
		  Double_t angle;
		  if(hit->GetDetectorID()/10000000==1){
		    if(b_survey) angle = survey->DTSurveyStereoAngle(hit->GetDetectorID());
		    else angle=stereo_angle;
		    vbot->RotateZ(-angle);
		    vtop->RotateZ(-angle);
		  }
		  if(hit->GetDetectorID()/10000000==2){
		    if(b_survey) angle = survey->DTSurveyStereoAngle(hit->GetDetectorID());
		    else angle=-stereo_angle;
		    vbot->RotateZ(-angle);
		    vtop->RotateZ(-angle);
		  }
		}

		   
		Double_t xhit=(vtop->x()+vbot->x())/2;
		Double_t zhit=(vtop->z()+vbot->z())/2;

		Double_t track_distance = xhit*cos(tangent.alpha)+zhit*sin(tangent.alpha)-tangent.p;
		Double_t residual = fabs(track_distance)-RTRel.GetRadius(hit->GetDigi());
		if(fabs(residual)<.9&&fabs(track_distance)<2.1){ //2.1 referes to half tube distance... realistic is 1.815, but for pattern reco ok
		  tangent.closehits++;
		  tangent.avres+=fabs(residual);
		}
		//std::cout << "\t" << residual;
	      }
	      tangent.avres/=tangent.closehits;
	      //we need at least 3 hits for a valid tangent
	      if(tangent.closehits>2){
		listoftangents.push_back(tangent);
	      }
	      //std::cout << "\t\t\t" << tangent.closehits << " close hits, average " << tangent.avres << std::endl;
	    }
	  }
	}
      }
    }
  }


    
  vbot->Delete();
  vtop->Delete();
  
  delete mySpectrometer;
  delete survey;
  
  listoftangents.sort();
  
  Int_t npass=0;
  while(listoftangents.size()>0){
    Int_t npass=NTubesPassed2D(listoftangents.front(), view, bin_station);

    if(npass<2*listoftangents.front().closehits) return listoftangents.front();
    else listoftangents.pop_front();;
  }

  
  tangent2d falsetangent;
  falsetangent.p=-1; // real p always positive
  falsetangent.alpha=-10; // real alpha does not exceed 2pi
  falsetangent.closehits=0;
  falsetangent.avres=0;
  return falsetangent;

}


tangent2d SimplePattern2D(std::vector <MufluxSpectrometerHit> List_MufluxSpectrometerHits, MufluxSpectrometerRTRelation &RTRel, int bin_station, int view){

  bool b_survey=true;

  MufluxSpectrometer* mySpectrometer= new MufluxSpectrometer();
  MufluxSpectrometerDTSurvey *survey = new MufluxSpectrometerDTSurvey();
  
  TVector3 *vtop = new TVector3();
  TVector3 *vbot = new TVector3();

  Double_t x[2],z[2],r[2];
  MufluxSpectrometerHit* hit;

  std::list<tangent2d> listoftangents;
  
  
  int n =  List_MufluxSpectrometerHits.size();
  for(int i=0;i<n;i++){
    hit = &(List_MufluxSpectrometerHits[i]);

    //int station=hit->GetDetectorID()/10000000;
    //int vnb=(hit->GetDetectorID()%10000000)/1000000;
    //if(station>2||station+vnb==2)continue; //only T1/2 no stereo

    if(hit->GetTimeOverThreshold()<40)continue;
    
    if(GetView(hit->GetDetectorID())!=view)continue;
    if(GetStationB(hit->GetDetectorID())&bin_station==0)continue;
    
    if(!b_survey) mySpectrometer->TubeEndPoints(hit->GetDetectorID(), *vtop, *vbot);
    else survey->TubeEndPointsSurvey(hit->GetDetectorID(), *vtop, *vbot);

    Double_t stereo_angle=60./180*TMath::Pi();
    
    if(view!=0){
      Double_t angle;
      if(hit->GetDetectorID()/10000000==1){
	if(b_survey) angle = survey->DTSurveyStereoAngle(hit->GetDetectorID());
	else angle=stereo_angle;
	vbot->RotateZ(-angle);
	vtop->RotateZ(-angle);
      }
      if(hit->GetDetectorID()/10000000==2){
	if(b_survey) angle = survey->DTSurveyStereoAngle(hit->GetDetectorID());
	else angle=-stereo_angle;
	vbot->RotateZ(-angle);
	vtop->RotateZ(-angle);
      }
    }
    
    x[0]=(vtop->x()+vbot->x())/2;
    z[0]=(vbot->z()+vbot->z())/2;
    r[0]=RTRel.GetRadius(hit->GetDigi());

    for(int j=i+1;j<n;j++){
      hit = &(List_MufluxSpectrometerHits[j]);

      //station=hit->GetDetectorID()/10000000;
      //vnb=(hit->GetDetectorID()%10000000)/1000000;
      //if(station>2||station+vnb==2)continue; //only T1/2 no stereo

      if(GetView(hit->GetDetectorID())!=view)continue;
      if((GetStationB(hit->GetDetectorID())&bin_station)==0)continue;
      
      if(!b_survey) mySpectrometer->TubeEndPoints(hit->GetDetectorID(), *vtop, *vbot);
      else survey->TubeEndPointsSurvey(hit->GetDetectorID(), *vtop, *vbot);

      if(view!=0){
	Double_t angle;
	if(hit->GetDetectorID()/10000000==1){
	  if(b_survey) angle = survey->DTSurveyStereoAngle(hit->GetDetectorID());
	  else angle=stereo_angle;
	  vbot->RotateZ(-angle);
	  vtop->RotateZ(-angle);
	}
	if(hit->GetDetectorID()/10000000==2){
	  if(b_survey) angle = survey->DTSurveyStereoAngle(hit->GetDetectorID());
	  else angle=-stereo_angle;
	  vbot->RotateZ(-angle);
	  vtop->RotateZ(-angle);
	}
      }
      
      
      x[1]=(vtop->x()+vbot->x())/2;
      z[1]=(vbot->z()+vbot->z())/2;
      r[1]=RTRel.GetRadius(hit->GetDigi());

    
      /*
      std::cout << "Comparing two hits" << std::endl;
      std::cout << "\t Hit 0 at " << x[0] << ", " << z[0] << "\t\t" << "radius " << r[0] << std::endl;
      std::cout << "\t Hit 1 at " << x[1] << ", " << z[1] << "\t\t" << "radius " << r[1] << std::endl; 
      */
     

      Double_t deltax=x[0]-x[1];
      Double_t deltaz=z[0]-z[1];
      
      Double_t distance=sqrt(pow(deltax,2)+pow(deltaz,2));
      //if(distance<5)continue; //skip adjacent tubes since high uncertainty due to lever arm 

      //std::cout << "\t Distance is " << distance <<std::endl;
      
      Double_t deltar[2];
      deltar[0]=r[0]-r[1];
      deltar[1]=r[0]+r[1];

      Int_t pm[2];
      pm[0]=1;
      pm[1]=-1;
      
      for(int k=0;k<2;k++){
	for(int l=0;l<2;l++){
	  for(int m=0;m<2;m++){
	    int o=4*k+2*l+m;
	    tangent2d tangent;
	    tangent.alpha=2*atan((deltaz+pm[l]*sqrt(pow(deltaz,2)-pow(deltar[k],2)+pow(deltax,2)))/(pm[m]*deltar[k]+deltax));
	    //m is the "sign" used for r[0] 
	    tangent.p=x[0]*cos(tangent.alpha)+z[0]*sin(tangent.alpha)-pm[m]*r[0];
	    /*
	    std::cout << "\t\t" << alpha_single[o]*180./3.14159 << "   \t" << p_single[o] ;//<< std::endl;
	    p_single[o]=x[1]*cos(alpha_single[o])+z[1]*sin(alpha_single[o])-pm[(m+k)%2]*r[1];
	    std::cout << "   \t" << p_single[o] << std::endl;
	    */
	    if(tangent.p>=0){
	      //std::cout << "\t\t" << tangent.alpha*180./3.14159 << "   \t" << tangent.p << std::endl;

	      tangent.closehits=0;
	      tangent.avres=0;
	      //std::cout << "\t\t\t Distance to other hits: "; 
	      //distance to other hits
	      for(int ii=0;ii<n;ii++){
		hit = &(List_MufluxSpectrometerHits[ii]);

		//station=hit->GetDetectorID()/10000000;
		//vnb=(hit->GetDetectorID()%10000000)/1000000;
		//if(station+vnb==2)continue; //no stereo
		if(GetView(hit->GetDetectorID())!=view)continue;
		if((GetStationB(hit->GetDetectorID())&bin_station)==0)continue;

		if(!b_survey) mySpectrometer->TubeEndPoints(hit->GetDetectorID(), *vtop, *vbot);
		else survey->TubeEndPointsSurvey(hit->GetDetectorID(), *vtop, *vbot);
		
		if(view!=0){
		  Double_t angle;
		  if(hit->GetDetectorID()/10000000==1){
		    if(b_survey) angle = survey->DTSurveyStereoAngle(hit->GetDetectorID());
		    else angle=stereo_angle;
		    vbot->RotateZ(-angle);
		    vtop->RotateZ(-angle);
		  }
		  if(hit->GetDetectorID()/10000000==2){
		    if(b_survey) angle = survey->DTSurveyStereoAngle(hit->GetDetectorID());
		    else angle=-stereo_angle;
		    vbot->RotateZ(-angle);
		    vtop->RotateZ(-angle);
		  }
		}

		   
		Double_t xhit=(vtop->x()+vbot->x())/2;
		Double_t zhit=(vtop->z()+vbot->z())/2;

		Double_t track_distance = xhit*cos(tangent.alpha)+zhit*sin(tangent.alpha)-tangent.p;
		Double_t residual = fabs(track_distance)-RTRel.GetRadius(hit->GetDigi());
		if(fabs(residual)<.9&&fabs(track_distance)<2.1){ //2.1 referes to half tube distance... realistic is 1.815, but for pattern reco ok
		  tangent.closehits++;
		  tangent.avres+=fabs(residual);
		}
		//std::cout << "\t" << residual;
	      }
	      tangent.avres/=tangent.closehits;
	      //we need at least 3 hits for a valid tangent
	      if(tangent.closehits>2){
		listoftangents.push_back(tangent);
	      }
	      //std::cout << "\t\t\t" << tangent.closehits << " close hits, average " << tangent.avres << std::endl;
	    }
	  }
	}
      }
    }
  }


    
  vbot->Delete();
  vtop->Delete();
  
  delete mySpectrometer;
  delete survey;
  
  listoftangents.sort();
  
  Int_t npass=0;
  while(listoftangents.size()>0){
    Int_t npass=NTubesPassed2D(listoftangents.front(), view, bin_station);

    if(npass<2*listoftangents.front().closehits) return listoftangents.front();
    else listoftangents.pop_front();;
  }

  
  tangent2d falsetangent;
  falsetangent.p=-1; // real p always positive
  falsetangent.alpha=-10; // real alpha does not exceed 2pi
  falsetangent.closehits=0;
  falsetangent.avres=0;
  return falsetangent;

}



Int_t NTubesPassed2D(tangent2d tangent, Int_t view, Int_t bin_station){
  Int_t npass=0;
  
  MufluxSpectrometerDTSurvey *survey = new MufluxSpectrometerDTSurvey();

  Double_t angle=0;
  if(view==1)angle=survey->DTSurveyStereoAngle(11002001);
  else if(view==2)angle=survey->DTSurveyStereoAngle(20002001);
    
  for (auto const& tube : survey->TubeList()) {
    
    if(GetView(tube)!=view)continue;
    if((GetStationB(tube)&bin_station)==0)continue;
    
    //TVector3 tubeaty0=survey->AtY(tube,0);
    TVector3 ttop;
    TVector3 tbot;
    survey->TubeEndPointsSurvey(tube, ttop, tbot);

    if(view!=0){
      tbot.RotateZ(-angle);
      ttop.RotateZ(-angle);
    }
    
    Double_t xwire=(ttop.x()+tbot.x())/2;
    Double_t zwire=(ttop.z()+tbot.z())/2;
    
    double dist= fabs(xwire*cos(tangent.alpha)+zwire*sin(tangent.alpha)-tangent.p);
    if(dist<1.815)npass++;
  }

  delete survey;
  
  return npass;
}



GBL_seed_track *seedtrack(TTreeReaderArray <MufluxSpectrometerHit> &Digi_MufluxSpectrometerHits, MufluxSpectrometerRTRelation &RTRel){
  
  //front = SimplePattern2D(Digi_MufluxSpectrometerHits,*RTRel,0b0011,0);
  //back = SimplePattern2D(Digi_MufluxSpectrometerHits,*RTRel,0b1100,0);
  tangent2d all = SimplePattern2D(Digi_MufluxSpectrometerHits,RTRel,0b1111,0);
  if(all.p<0)return nullptr;
  
  tangent2d stereo1 = SimplePattern2D(Digi_MufluxSpectrometerHits,RTRel,0b0011,1);
  if(stereo1.p<0)return nullptr;
  tangent2d stereo2 = SimplePattern2D(Digi_MufluxSpectrometerHits,RTRel,0b0011,2);
  if(stereo2.p<0)return nullptr;

  MufluxSpectrometerDTSurvey *surv = new MufluxSpectrometerDTSurvey();
  
  Double_t beta1=surv->DTSurveyStereoAngle(11002001);
  Double_t beta2=surv->DTSurveyStereoAngle(20002001);
  
  Double_t xslope=-tan(all.alpha);
  Double_t x0=all.p/cos(all.alpha);
  
  Double_t yslope = YslopeFromProjections(stereo1.alpha, beta1, stereo2.alpha, beta2);
  Double_t y0 = Y0FromProjections(stereo1.alpha, stereo1.p, beta1, stereo2.alpha, stereo2.p, beta2);
  
  
  TVector3 position,direction;
  position.SetXYZ(x0,y0,0);
  direction.SetXYZ(xslope,yslope,1);
 
  GBL_seed_track *seed=new GBL_seed_track(position, direction);

  TVector3 *vtop=new TVector3();
  TVector3 *vbot=new TVector3();

  MufluxSpectrometerHit* hit;
  
  for(int ii=0;ii<Digi_MufluxSpectrometerHits.GetSize();ii++){
    hit = &(Digi_MufluxSpectrometerHits[ii]);

    surv->TubeEndPointsSurvey(hit->GetDetectorID(), *vtop, *vbot);

    tangent2d tangent=all;
    int view=GetView(hit->GetDetectorID());
    if(view!=0){
      Double_t angle = surv->DTSurveyStereoAngle(hit->GetDetectorID());
      vbot->RotateZ(-angle);
      vtop->RotateZ(-angle);
      if(view==1)tangent=stereo1;
      else tangent=stereo2;
    }
    
    
    Double_t xhit=(vtop->x()+vbot->x())/2;
    Double_t zhit=(vtop->z()+vbot->z())/2;
    
    Double_t track_distance = xhit*cos(tangent.alpha)+zhit*sin(tangent.alpha)-tangent.p;
    Double_t residual = fabs(track_distance)-RTRel.GetRadius(hit->GetDigi());
    if(fabs(residual)<.9&&fabs(track_distance)<2.385){ //2.1 referes to half tube distance... realistic is 1.815, but for pattern reco ok, 4.2-1.815 means that git is not in next tube
      seed->add_hit(hit->GetDetectorID(),RTRel.GetRadius(hit->GetDigi()));
    }

  }
  


  delete vtop;
  delete vbot;
  delete surv;
  return seed;
  
}



GBL_seed_track *seedtrack(std::vector <MufluxSpectrometerHit> Digi_MufluxSpectrometerHits, MufluxSpectrometerRTRelation &RTRel){

  std::vector<MufluxSpectrometerHit> HitsT12X;
  std::vector<MufluxSpectrometerHit> HitsT1U;
  std::vector<MufluxSpectrometerHit> HitsT2V;
  std::vector<MufluxSpectrometerHit> HitsT34X;
  std::vector<MufluxSpectrometerHit> HitsT1234X;
  
  //MufluxSpectrometerHit* hit;
  for(int i=0;i<Digi_MufluxSpectrometerHits.size();i++){

    Int_t ModID=Digi_MufluxSpectrometerHits[i].GetDetectorID()/1e6;
    
    if(ModID>=30){
      HitsT1234X.push_back(Digi_MufluxSpectrometerHits[i]);
      HitsT34X.push_back(Digi_MufluxSpectrometerHits[i]);
    }
    else if(ModID==10||ModID==21){
      HitsT1234X.push_back(Digi_MufluxSpectrometerHits[i]);
      HitsT12X.push_back(Digi_MufluxSpectrometerHits[i]);
    }
    else if(ModID==11){
      HitsT1U.push_back(Digi_MufluxSpectrometerHits[i]);
    }
    else if(ModID==20){
      HitsT2V.push_back(Digi_MufluxSpectrometerHits[i]);
    }	
    
  }
  
  
  
  //front = SimplePattern2D(Digi_MufluxSpectrometerHits,*RTRel,0b0011,0);
  //back = SimplePattern2D(Digi_MufluxSpectrometerHits,*RTRel,0b1100,0);
  //tangent2d all = SimplePattern2D(Digi_MufluxSpectrometerHits,RTRel,0b1111,0);
  tangent2d all = SimplePattern2D(HitsT1234X,RTRel,0b1111,0);
  if(all.p<0)return nullptr;
  
  //tangent2d stereo1 = SimplePattern2D(Digi_MufluxSpectrometerHits,RTRel,0b0011,1);
  tangent2d stereo1 = SimplePattern2D(HitsT1U,RTRel,0b0011,1);
  if(stereo1.p<0)return nullptr;
  //tangent2d stereo2 = SimplePattern2D(Digi_MufluxSpectrometerHits,RTRel,0b0011,2);
  tangent2d stereo2 = SimplePattern2D(HitsT2V,RTRel,0b0011,2);
  if(stereo2.p<0)return nullptr;

  MufluxSpectrometerDTSurvey *surv = new MufluxSpectrometerDTSurvey();
  
  Double_t beta1=surv->DTSurveyStereoAngle(11002001);
  Double_t beta2=surv->DTSurveyStereoAngle(20002001);
  
  Double_t xslope=-tan(all.alpha);
  Double_t x0=all.p/cos(all.alpha);
  
  Double_t yslope = YslopeFromProjections(stereo1.alpha, beta1, stereo2.alpha, beta2);
  Double_t y0 = Y0FromProjections(stereo1.alpha, stereo1.p, beta1, stereo2.alpha, stereo2.p, beta2);
  
  
  TVector3 position,direction;
  position.SetXYZ(x0,y0,0);
  direction.SetXYZ(xslope,yslope,1);
 
  GBL_seed_track *seed=new GBL_seed_track(position, direction);

  TVector3 *vtop=new TVector3();
  TVector3 *vbot=new TVector3();

  MufluxSpectrometerHit* hit;
  
  for(int ii=0;ii<Digi_MufluxSpectrometerHits.size();ii++){
    hit = &(Digi_MufluxSpectrometerHits[ii]);
    
    surv->TubeEndPointsSurvey(hit->GetDetectorID(), *vtop, *vbot);

    tangent2d tangent=all;
    int view=GetView(hit->GetDetectorID());
    if(view!=0){
      Double_t angle = surv->DTSurveyStereoAngle(hit->GetDetectorID());
      vbot->RotateZ(-angle);
      vtop->RotateZ(-angle);
      if(view==1)tangent=stereo1;
      else tangent=stereo2;
    }
    
    
    Double_t xhit=(vtop->x()+vbot->x())/2;
    Double_t zhit=(vtop->z()+vbot->z())/2;
    
    Double_t track_distance = xhit*cos(tangent.alpha)+zhit*sin(tangent.alpha)-tangent.p;
    Double_t residual = fabs(track_distance)-RTRel.GetRadius(hit->GetDigi());
    if(fabs(residual)<.9&&fabs(track_distance)<2.385){ //2.1 referes to half tube distance... realistic is 1.815, but for pattern reco ok, 4.2-1.815 means that git is not in next tube
      seed->add_hit(hit->GetDetectorID(),RTRel.GetRadius(hit->GetDigi()));
    }

  }
  


  delete vtop;
  delete vbot;
  delete surv;
  return seed;
  
}

