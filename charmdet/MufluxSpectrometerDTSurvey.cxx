#include "MufluxSpectrometerDTSurvey.h"
#include "ShipUnit.h"
#include "iostream"

/**
 * Constructor. 
 *
 * @brief Constructor
 *
 * @author Daniel Bick
 * @date May 20, 2020
 * @version 1.0
 *
 */
MufluxSpectrometerDTSurvey::MufluxSpectrometerDTSurvey(){}


/**
 * Initializes the std map for the tube end points based on the survey
 *
 * @brief Initializes tube end points
 *
 * @author Daniel Bick
 * @date May 20, 2020
 * @version 1.0
 *
 */
void MufluxSpectrometerDTSurvey::Init(){
  
  TubeEndPoints endpoints;

  for(int module=0;module<12;module++){
    int station=4;
    if(module<8) station=3;
    if(module<4) station=2;
    if(module<2) station=1;
    int view=0;
    if(module==1||module==3)view=1;
    int BaseDetectorID=station*1e7+view*1e6+1;
    if(station>2)BaseDetectorID+=(module%4)*12;
    DTSurveyPoints spoint=DTSurveyMuflux(BaseDetectorID)+DTSurveyDistanceToRefTube(BaseDetectorID);

    /*
    std::cout << BaseDetectorID << std::endl;
    std::cout << "\t" << spoint.top_jur.x() << " \t" << spoint.bot_jur.x() << " \t" << spoint.top_sal.x() << " \t"<< spoint.bot_sal.x() << std::endl;
    std::cout << "\t" << spoint.top_jur.y() << " \t" << spoint.bot_jur.y() << " \t" << spoint.top_sal.y() << " \t"<< spoint.bot_sal.y() << std::endl;
    std::cout << "\t" << spoint.top_jur.z() << " \t" << spoint.bot_jur.z() << " \t" << spoint.top_sal.z() << " \t"<< spoint.bot_sal.z() << std::endl;
    */

    TubeEndPoints reftube;
    
    reftube.top.SetXYZ((spoint.top_jur.x()+spoint.top_sal.x())/2,
			(spoint.top_jur.y()+spoint.top_sal.y())/2,
			(spoint.top_jur.z()+spoint.top_sal.z())/2);
    reftube.bot.SetXYZ((spoint.bot_jur.x()+spoint.bot_sal.x())/2,
			(spoint.bot_jur.y()+spoint.bot_sal.y())/2,
			(spoint.bot_jur.z()+spoint.bot_sal.z())/2);


    
    for(int pnb=0;pnb<2;pnb++){
      for(int lnb=0;lnb<2;lnb++){
	for(int tube=0;tube<12;tube++){

	  int snb=tube+1;
	  if(module>=4)snb+=(module%4)*12;

	  int DetectorID=station*1e7+view*1e6+pnb*1e5+lnb*1e4+2000+snb;

	  TVector3 staggering=DTStaggering(DetectorID);

	  TubeEndPoints endpoint;
	  endpoint.top=reftube.top+staggering;
	  endpoint.bot=reftube.bot+staggering;

	  DTSurveyEndPointMap[DetectorID]=endpoint;

	  /*
	  std::cout << DetectorID << std::endl;
	  std::cout << "\t" << endpoint.top.x() << "  \t" << endpoint.bot.x() << std::endl;
	  std::cout << "\t" << endpoint.top.y() << "  \t" << endpoint.bot.y() << std::endl;
	  std::cout << "\t" << endpoint.top.z() << "  \t" << endpoint.bot.z() << std::endl;
	  */
	}
      }
    }
  }
  DTSurveyIsInitialized=true;
}



/**
 * Test function to print the coordinates of the tube end points
 *
 * @brief Print coordinates of the tube end points
 *
 * @param DetectorID ID of the tube to test
 *
 * @author Daniel Bick
 * @date May 20, 2020
 * @version 1.0
 *
 */
void MufluxSpectrometerDTSurvey::Test(Int_t DetectorID){
  /*
  std::cout << DTSurveyEndPointMap.find(DetectorID)->second.top.x() << std::endl;
  std::cout << DTSurveyEndPointMap.find(DetectorID)->second.top.y() << std::endl;
  std::cout << DTSurveyEndPointMap.find(DetectorID)->second.top.z() << std::endl;
  std::cout << DTSurveyEndPointMap.find(DetectorID)->second.bot.x() << std::endl;
  std::cout << DTSurveyEndPointMap.find(DetectorID)->second.bot.y() << std::endl;
  std::cout << DTSurveyEndPointMap.find(DetectorID)->second.bot.z() << std::endl;
  */
  

  double xt=DTSurveyEndPointMap[DetectorID].top.x();
  double xb=DTSurveyEndPointMap[DetectorID].bot.x();
  double yt=DTSurveyEndPointMap[DetectorID].top.y();
  double yb=DTSurveyEndPointMap[DetectorID].bot.y();
  double zt=DTSurveyEndPointMap[DetectorID].top.z();
  double zb=DTSurveyEndPointMap[DetectorID].bot.z();

  TVector3 center=AtY(DetectorID,0);
  
  std::cout << DetectorID << std::endl;
  std::cout << "   \t TOP      \t BOT          \t y=0" << std::endl;
  std::cout << "  x\t" << xt << "       \t"  << xb << "       \t" << center.x() << std::endl;
  std::cout << "  y\t" << yt << "       \t"  << yb << "       \t0" << std::endl;
  std::cout << "  z\t" << zt << "       \t"  << zb << "       \t" << center.z() << std::endl;
}


/**
 * Returns a unique number for each physical module
 *
 * @brief Returns a unique number for each physical module
 *
 * @param DetectorID ID of any tube in that module
 *
 * @author Daniel Bick
 * @date May 20, 2020
 * @version 1.0
 *
 */
Int_t MufluxSpectrometerDTSurvey::Module(Int_t DetectorID){
  return DetectorID/1e6+(DetectorID%100-1)/12; //unique for physical modules
}

/**
 * Returns DTSurveyPoints struct with the original survey points of a mdoule
 *
 * @brief survey points of a module
 *
 * @param DetectorID ID of any tube in that module
 *
 * @author Daniel Bick
 * @date May 20, 2020
 * @version 1.0
 *
 */
DTSurveyPoints  MufluxSpectrometerDTSurvey::DTSurveyMuflux(Int_t DetectorID){
  DTSurveyPoints SurveyPoints;

  Int_t module=Module(DetectorID);
 
  if(module==10){
    
    SurveyPoints.top_jur.SetX((0.2443*1000)*ShipUnit::mm);
    SurveyPoints.top_jur.SetY((0.7102*1000)*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ((9.0527*1000)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX((.2436*1000)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY((-.6503*1000)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ((9.0578*1000)*ShipUnit::mm);
    
    SurveyPoints.top_sal.SetX((-.2078*1000)*ShipUnit::mm);
    SurveyPoints.top_sal.SetY((.7092*1000)*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ((9.0502*1000)*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX((-.2075*1000)*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY((-.6495*1000)*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ((9.0564*1000)*ShipUnit::mm);
    
  }
  
  else if(module==11){
    
    SurveyPoints.top_jur.SetX(-504.544*ShipUnit::mm);
    SurveyPoints.top_jur.SetY(564.599*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ(9537.42*ShipUnit::mm);
    SurveyPoints.top_sal.SetX(-729.384*ShipUnit::mm);
    SurveyPoints.top_sal.SetY(172.827*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ(9539*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX(444.63*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY(-499.492*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ(9541.94*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX(669.285*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY(-107.278*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ(9543.47*ShipUnit::mm);

  }//module11

  else if(module==20){

    SurveyPoints.top_jur.SetX(758.838*ShipUnit::mm);
    SurveyPoints.top_jur.SetY(173.029*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ(9742.35*ShipUnit::mm);
    SurveyPoints.top_sal.SetX(532.668*ShipUnit::mm);
    SurveyPoints.top_sal.SetY(565.716*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ(9743.27*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX(-637.445*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY(-110.358*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ(9738.21*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX(-412.303*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY(-501.888*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ(9740.98*ShipUnit::mm);

  }//module==20

  else if(module==21){
    
    SurveyPoints.top_jur.SetX((0.2385*1000)*ShipUnit::mm);
    SurveyPoints.top_jur.SetY((0.7078*1000)*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ((10.2314*1000)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX((.2361*1000)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY((-.6495*1000)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ((10.2302*1000)*ShipUnit::mm);
  
    SurveyPoints.top_sal.SetX((-.2121*1000)*ShipUnit::mm);
    SurveyPoints.top_sal.SetY((.7087*1000)*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ((10.2276*1000)*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX((-.2157*1000)*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY((-.6488*1000)*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ((10.2285*1000)*ShipUnit::mm);
  
  }//module==21

  else if(module==33){
  
    SurveyPoints.top_jur.SetX(925.397*ShipUnit::mm);
    SurveyPoints.top_jur.SetY(893.085*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ(14573.2*ShipUnit::mm);
    SurveyPoints.top_sal.SetX(589.417*ShipUnit::mm);
    SurveyPoints.top_sal.SetY(891.359*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ(14574.1*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX(592.618*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY(-682.233*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ(14570.3*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX(928.585*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY(-681.743*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ(14571.1*ShipUnit::mm);

  }

  else if(module==32){
    
    SurveyPoints.top_jur.SetX(421.229*ShipUnit::mm);
    SurveyPoints.top_jur.SetY(890.753*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ(14574.6*ShipUnit::mm);
    SurveyPoints.top_sal.SetX(85.282*ShipUnit::mm);
    SurveyPoints.top_sal.SetY(890.549*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ(14575*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX(88.457*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY(-685.388*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ(14568.6*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX(424.517*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY(-684.353*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ(14569.9*ShipUnit::mm);

  }

  else if(module==31){

    SurveyPoints.top_jur.SetX(-83.89*ShipUnit::mm);
    SurveyPoints.top_jur.SetY(889.943*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ(14575.5*ShipUnit::mm);
    SurveyPoints.top_sal.SetX(-419.692*ShipUnit::mm);
    SurveyPoints.top_sal.SetY(888.852*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ(14576.9*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX(-417.093*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY(-683.992*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ(14569.3*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX(-81.229*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY(-683.565*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ(14568.5*ShipUnit::mm);
 
  }

  else if(module==30){

    SurveyPoints.top_jur.SetX(-589.556*ShipUnit::mm);
    SurveyPoints.top_jur.SetY(890.769*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ(14578.1*ShipUnit::mm);
    SurveyPoints.top_sal.SetX(-925.574*ShipUnit::mm);
    SurveyPoints.top_sal.SetY(889.619*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ(14581.2*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX(-921.584*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY(-684.471*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ(14569.2*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX(-585.804*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY(-686.33*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ(14569.6*ShipUnit::mm);
  
  }

  else if(module==43){
  
    SurveyPoints.top_jur.SetX(920.702*ShipUnit::mm);
    SurveyPoints.top_jur.SetY(889.806*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ(16545*ShipUnit::mm);
    SurveyPoints.top_sal.SetX(584.439*ShipUnit::mm);
    SurveyPoints.top_sal.SetY(888.298*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ(16545.7*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX(582.355*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY(-686.731*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ(16541.9*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX(918.349*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY(-684.858*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ(16543.6*ShipUnit::mm);

  }

  else if(module==42){
   
    SurveyPoints.top_jur.SetX(416.735*ShipUnit::mm);
    SurveyPoints.top_jur.SetY(887.258*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ(16546.1*ShipUnit::mm);
    SurveyPoints.top_sal.SetX(80.336*ShipUnit::mm);
    SurveyPoints.top_sal.SetY(886.136*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ(16547.5*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX(78.512*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY(-688.315*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ(16538.9*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX(414.365*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY(-687.523*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ(16540.9*ShipUnit::mm);
  
  }


  else if(module==41){
    
    SurveyPoints.top_jur.SetX(-87.554*ShipUnit::mm);
    SurveyPoints.top_jur.SetY(886.099*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ(16548*ShipUnit::mm);
    SurveyPoints.top_sal.SetX(-423.404*ShipUnit::mm);
    SurveyPoints.top_sal.SetY(886.185*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ(16549.8*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX(-424.753*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY(-688.486*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ(16539.6*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX(-88.796*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY(-689.051*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ(16538.9*ShipUnit::mm);

  }

  else if(module==40){
  
    SurveyPoints.top_jur.SetX(-591.889*ShipUnit::mm);
    SurveyPoints.top_jur.SetY(886.582*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ(16550.8*ShipUnit::mm);
    SurveyPoints.top_sal.SetX(-928.063*ShipUnit::mm);
    SurveyPoints.top_sal.SetY(886.854*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ(16551.8*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX(-929.54*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY(-687.767*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ(16540.2*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX(-593.653*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY(-688.852*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ(16540*ShipUnit::mm);

  }


  else{

    std::cout << "WARNING: MufluxSpectrometerDTSurvey::DTSurveyMuflux could not convert detector id " << DetectorID << " to valid module number" << std::endl;
  }

  //z correction for survey reference system
  TVector3 zcorrect;
  zcorrect.SetXYZ(0,0,-8.89609*ShipUnit::m);
  SurveyPoints.top_jur+=zcorrect;
  SurveyPoints.top_sal+=zcorrect;
  SurveyPoints.bot_jur+=zcorrect;
  SurveyPoints.bot_sal+=zcorrect;
  
  return SurveyPoints;
}


/**
 * A struct DTSurveyPoints is returned with the distance of the survey points in the module to a refrence tube (tube 1, layer 0, plane 0 for non-roated modules, tube 12 layer plane 1 otherwise)
 *
 * @brief Hard coded distance of reference tube to survey points for all modules
 *
 * @param DetectorID ID of any tube in the module
 *
 * @author Daniel Bick
 * @date May 20, 2020
 * @version 1.0
 *
 */
DTSurveyPoints MufluxSpectrometerDTSurvey::DTSurveyDistanceToRefTube(Int_t DetectorID){
  DTSurveyPoints distance; //distance to lower right tube (upper left for rotated modules)

  Int_t module=Module(DetectorID);
  double mm=ShipUnit::mm;


    
    if(module==10){
      
      distance.top_jur.SetXYZ(-11*42*mm-10*mm,-170*mm-8.35*mm,30.3*mm);
      distance.top_sal.SetXYZ(-10*mm,-170*mm-8.35*mm,30.3*mm);
      
      distance.bot_jur.SetXYZ(-11*42*mm-10*mm,70*mm+8.65*mm,30.3*mm);
      distance.bot_sal.SetXYZ(-10*mm,70*mm+8.65*mm,30.3*mm);
    
    }

    else if(module==11){
      
      distance.top_jur.SetXYZ(-11*42*mm-10*mm,-170*mm,-143.7*mm);
      distance.top_sal.SetXYZ(-10*mm,-170*mm,-143.7*mm);
      
      distance.bot_jur.SetXYZ(-11*42*mm-10*mm,70*mm,-143.7*mm);
      distance.bot_sal.SetXYZ(-10*mm,70*mm,-143.7*mm);

      double stereo=DTSurveyStereoAngle(DetectorID);//-1.051;
      
      distance.top_jur.RotateZ(-stereo);
      distance.top_sal.RotateZ(-stereo);
      
      distance.bot_jur.RotateZ(-stereo);
      distance.bot_sal.RotateZ(-stereo);
      
    }

    else if(module==20){ //rotated

      distance.top_sal.SetXYZ(11*42*mm+10*mm,-170*mm,143.7*mm);
      distance.top_jur.SetXYZ(10*mm,-170*mm,143.7*mm);
      
      distance.bot_sal.SetXYZ(11*42*mm+10*mm,70*mm,143.7*mm);
      distance.bot_jur.SetXYZ(10*mm,70*mm,143.7*mm);

      double stereo=DTSurveyStereoAngle(DetectorID);
      
      distance.top_jur.RotateZ(-stereo);
      distance.top_sal.RotateZ(-stereo);
      
      distance.bot_jur.RotateZ(-stereo);
      distance.bot_sal.RotateZ(-stereo);
      
    }

    else if(module==21){ //rotated

      distance.top_sal.SetXYZ(11*42*mm+10*mm,-170*mm-7.3*mm,-30.3*mm);
      distance.top_jur.SetXYZ(10*mm,-170*mm-7.3*mm,-30.3*mm);
      
      distance.bot_sal.SetXYZ(11*42*mm+10*mm,70*mm+8.45*mm,-30.3*mm);
      distance.bot_jur.SetXYZ(10*mm,70*mm+8.45*mm,-30.3*mm);

    }

    else if(module==30||module==31||module==32||module==33){

      distance.top_sal.SetXYZ(-63*mm,12.5*mm,70*mm+43.3*mm);
      distance.top_jur.SetXYZ(-399*mm,12.5*mm,70*mm+43.3*mm);

      distance.bot_sal.SetXYZ(-63*mm,-12.5*mm,70*mm+43.3*mm);
      distance.bot_jur.SetXYZ(-399*mm,-12.5*mm,70*mm+43.3*mm);
      
    }
  
    else if(module==40||module==41||module==42||module==43){
      
      distance.top_sal.SetXYZ(-52*mm,12.5*mm,-70*mm-43.3*mm-113.4*mm);
      distance.top_jur.SetXYZ(-388*mm,12.5*mm,-70*mm-43.3*mm-113.4*mm);
      
      distance.bot_sal.SetXYZ(-52*mm,-12.5*mm,-70*mm-43.3*mm-113.4*mm);
      distance.bot_jur.SetXYZ(-388*mm,-12.5*mm,-70*mm-43.3*mm-113.4*mm);
      
    }
  
    else {
      std::cout << "WARNING: MufluxSpectrometerDTSurvey::DTSurveyDistanceToRefTube could not convert detector id " << DetectorID << " to valid module number" << std::endl;
    }
    
    return distance;
}



/**
 * Here the staggering of the module is implemented. A vector is returned with the distance of the tube and the refrence tube (tube 1, layer 0, plane 0 for non-roated modules, tube 12 layer plane 1 otherwise)
 *
 * @brief Hard coded distance of reference tube to givent tube in module
 *
 * @param DetectorID ID of the tube in the module
 *
 * @author Daniel Bick
 * @date May 20, 2020
 * @version 1.0
 *
 */
 TVector3 MufluxSpectrometerDTSurvey::DTStaggering(Int_t DetectorID){
  TVector3 staggering;
  Int_t pnb=(DetectorID%1000000)/100000;
  Int_t lnb=(DetectorID%100000)/10000;
  Int_t tube=((DetectorID-1)%100)%12;

  Int_t station=DetectorID/1000000;

  if(station==20||station==21){//rotated
    pnb=(pnb+1)&1;
    lnb=(lnb+1)&1;
    tube=11-tube;
  }
  
  double x,z;

  z=pnb*77*ShipUnit::mm+lnb*36.4*ShipUnit::mm;
  x=tube*42*ShipUnit::mm+lnb*21*ShipUnit::mm+pnb*10*ShipUnit::mm;

  if((pnb+lnb)==2)x-=42*ShipUnit::mm;

  if(station==20||station==21) staggering.SetXYZ(-x,0,-z);
  else staggering.SetXYZ(x,0,z);

  if(station==11){
    double stereo=DTSurveyStereoAngle(DetectorID);
    staggering.RotateZ(-stereo);
  }

  if(station==20){
    double stereo=DTSurveyStereoAngle(DetectorID);
    staggering.RotateZ(-stereo);
  }
  
  return staggering;
}


/**
 * Gives you the position of the wire top and bottom calculated from survey
 *
 * @brief returns position of the wire
 *
 * @param DetectorID ID of the tube in the module
 * @param vtop position of the top wire end
 * @param vbot position of the bottom wire end
 *
 * @author Daniel Bick
 * @date May 20, 2020
 * @version 1.0
 *
 */
void MufluxSpectrometerDTSurvey::TubeEndPointsSurvey(Int_t DetectorID, TVector3 &vtop, TVector3 &vbot){
  vtop=DTSurveyEndPointMap[DetectorID].top;
  vbot=DTSurveyEndPointMap[DetectorID].bot;
}

/**
 * Calculates the rotation of the module around z from the two long sides (distance of the survey points)
 *
 * @brief returns rotation of module around z
 *
 * @param DetectorID ID of any tube in the module
 *
 * @author Daniel Bick
 * @date May 20, 2020
 * @version 1.0
 *
 */
Double_t MufluxSpectrometerDTSurvey::DTSurveyStereoAngle(Int_t DetectorID){
  
  DTSurveyPoints spoint=DTSurveyMuflux(DetectorID);
  
  Double_t angle_jur=atan((spoint.top_jur.x()-spoint.bot_jur.x())/(spoint.top_jur.y()-spoint.bot_jur.y()));
  Double_t angle_sal=atan((spoint.top_sal.x()-spoint.bot_sal.x())/(spoint.top_sal.y()-spoint.bot_sal.y()));

  //std::cout << angle_jur << " " << angle_sal << std::endl;

  return (angle_jur+angle_sal)/2;
  
}

/**
 * Calculates tube length based on survey data
 *
 * @brief Calculates tube length based on survey data
 *
 * @param DetectorID ID of the tube 
 *
 * @author Daniel Bick
 * @date May 20, 2020
 * @version 1.0
 *
 */
Double_t MufluxSpectrometerDTSurvey::DTSurveyTubeLength(Int_t DetectorID){
  if(!DTSurveyIsInitialized) std::cout << "Warning: MufluxSpectrometerDTSurvey is not initialized!" << std::endl;
  TVector3 vtop=DTSurveyEndPointMap[DetectorID].top;
  TVector3 vbot=DTSurveyEndPointMap[DetectorID].bot;
  TVector3 wire=vtop-vbot;
  return wire.Mag();
}

/**
 * Returns point on wire at given y
 *
 * @brief Returns point on wire at given y
 *
 * @param DetectorID ID of the tube 
 * @param y the y coordinate of the point
 *
 * @author Daniel Bick
 * @date May 20, 2020
 * @version 1.0
 *
 */
TVector3 MufluxSpectrometerDTSurvey::AtY(Int_t DetectorID, Double_t y){
  
  Double_t xt=DTSurveyEndPointMap[DetectorID].top.x();
  Double_t xb=DTSurveyEndPointMap[DetectorID].bot.x();
  Double_t yt=DTSurveyEndPointMap[DetectorID].top.y();
  Double_t yb=DTSurveyEndPointMap[DetectorID].bot.y();
  Double_t zt=DTSurveyEndPointMap[DetectorID].top.z();
  Double_t zb=DTSurveyEndPointMap[DetectorID].bot.z();

  TVector3 wirepos;

  wirepos.SetXYZ(xb+(xt-xb)*(y-yb)/(yt-yb),y,zb+(zt-zb)*(y-yb)/(yt-yb));

  return wirepos;
}

