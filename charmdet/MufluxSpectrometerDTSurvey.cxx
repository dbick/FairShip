#include "MufluxSpectrometerDTSurvey.h"
#include "ShipUnit.h"
#include "iostream"
MufluxSpectrometerDTSurvey::MufluxSpectrometerDTSurvey(){}


void MufluxSpectrometerDTSurvey::Init(){
  
  TubeEndPoints enden;
  
  enden.top.SetXYZ(1,2,3);
  enden.bot.SetXYZ(4,5,6);
  
  EndPointMap.insert(std::pair<int,TubeEndPoints>(10001000,enden));
  enden.top.SetXYZ(11,21,31);
  enden.bot.SetXYZ(41,51,61);

  EndPointMap[10001001]=enden;
}


void MufluxSpectrometerDTSurvey::Test(){
  std::cout << EndPointMap.find(10001000)->second.top.x() << std::endl;
  std::cout << EndPointMap.find(10001000)->second.top.y() << std::endl;
  std::cout << EndPointMap.find(10001000)->second.top.z() << std::endl;
  std::cout << EndPointMap.find(10001000)->second.bot.x() << std::endl;
  std::cout << EndPointMap.find(10001000)->second.bot.y() << std::endl;
  std::cout << EndPointMap.find(10001000)->second.bot.z() << std::endl;
  std::cout << EndPointMap[10001001].top.x() << std::endl;
  std::cout << EndPointMap[10001001].top.y() << std::endl;
  std::cout << EndPointMap[10001001].top.z() << std::endl;
  std::cout << EndPointMap[10001001].bot.x() << std::endl;
  std::cout << EndPointMap[10001001].bot.y() << std::endl;
  std::cout << EndPointMap[10001001].bot.z() << std::endl;
}

Int_t MufluxSpectrometerDTSurvey::Module(Int_t DetectorID){
  return DetectorID/1e6+(DetectorID%1000)/12; //unique for physical modules
}

DTSurveyPoints  MufluxSpectrometerDTSurvey::DTSurveyMuflux(Int_t DetectorID){
  DTSurveyPoints SurveyPoints;

  Int_t module=Module(DetectorID);
 
  if(module==10){
    
    SurveyPoints.top_jur.SetX((0.2443*1000)*ShipUnit::mm);
    SurveyPoints.top_jur.SetY((0.7102*1000-8.35-150-20)*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ((9.0527*1000)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX((.2436*1000)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY((-.6503*1000+8.65+50+20)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ((9.0578*1000)*ShipUnit::mm);
    
    SurveyPoints.top_sal.SetX((-.2078*1000)*ShipUnit::mm);
    SurveyPoints.top_sal.SetY((.7092*1000-8.35-150-20)*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ((9.0502*1000)*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX((-.2075*1000)*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY((-.6495*1000+8.65+50+20)*ShipUnit::mm);
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

    std::cout << "WARNING: MufluxSpectrometerDTSurvey::DTSurveyMuflux bolts not correctly implemented for stereo module survey data" << std::endl;
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

    std::cout << "WARNING: MufluxSpectrometerDTSurvey::DTSurveyMuflux bolts not correctly implemented for stereo module survey data" << std::endl;
  }//module==20

  else if(module==21){
    
    SurveyPoints.top_jur.SetX((0.2385*1000)*ShipUnit::mm);
    SurveyPoints.top_jur.SetY((0.7078*1000-7.3-150-20)*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ((10.2314*1000)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX((.2361*1000)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY((-.6495*1000+8.45+50+20)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ((10.2302*1000)*ShipUnit::mm);
  
    SurveyPoints.top_sal.SetX((-.2121*1000)*ShipUnit::mm);
    SurveyPoints.top_sal.SetY((.7087*1000-7.3-150-20)*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ((10.2276*1000)*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX((-.2157*1000)*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY((-.6488*1000+8.45+50+20)*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ((10.2285*1000)*ShipUnit::mm);
  
  }//module==21

  else if(module==30){
  
    SurveyPoints.top_jur.SetX(925.397*ShipUnit::mm);
    SurveyPoints.top_jur.SetY(893.085*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ((14573.2+70)*ShipUnit::mm);
    SurveyPoints.top_sal.SetX(589.417*ShipUnit::mm);
    SurveyPoints.top_sal.SetY(891.359*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ((14574.1+70)*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX(592.618*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY(-682.233*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ((14570.3+70)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX(928.585*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY(-681.743*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ((14571.1+70)*ShipUnit::mm);

  }

  else if(module==31){
    
    SurveyPoints.top_jur.SetX(421.229*ShipUnit::mm);
    SurveyPoints.top_jur.SetY(890.753*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ((14574.6+70)*ShipUnit::mm);
    SurveyPoints.top_sal.SetX(85.282*ShipUnit::mm);
    SurveyPoints.top_sal.SetY(890.549*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ((14575+70)*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX(88.457*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY(-685.388*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ((14568.6+70)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX(424.517*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY(-684.353*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ((14569.9+70)*ShipUnit::mm);

  }

  else if(module==32){

    SurveyPoints.top_jur.SetX(-83.89*ShipUnit::mm);
    SurveyPoints.top_jur.SetY(889.943*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ(14575.5+70*ShipUnit::mm);
    SurveyPoints.top_sal.SetX(-419.692*ShipUnit::mm);
    SurveyPoints.top_sal.SetY(888.852*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ(14576.9+70*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX(-417.093*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY(-683.992*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ(14569.3+70*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX(-81.229*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY(-683.565*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ(14568.5+70*ShipUnit::mm);
 
  }

  else if(module==33){

    SurveyPoints.top_jur.SetX(-589.556*ShipUnit::mm);
    SurveyPoints.top_jur.SetY(890.769*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ((14578.1+70)*ShipUnit::mm);
    SurveyPoints.top_sal.SetX(-925.574*ShipUnit::mm);
    SurveyPoints.top_sal.SetY(889.619*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ((14581.2+70)*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX(-921.584*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY(-684.471*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ((14569.2+70)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX(-585.804*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY(-686.33*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ((14569.6+70)*ShipUnit::mm);
  
  }

  else if(module==40){
  
    SurveyPoints.top_jur.SetX(920.702*ShipUnit::mm);
    SurveyPoints.top_jur.SetY(889.806*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ((16545-70)*ShipUnit::mm);
    SurveyPoints.top_sal.SetX(584.439*ShipUnit::mm);
    SurveyPoints.top_sal.SetY(888.298*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ((16545.7-70)*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX(582.355*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY(-686.731*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ((16541.9-70)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX(918.349*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY(-684.858*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ((16543.6-70)*ShipUnit::mm);

  }

  else if(module==41){
   
    SurveyPoints.top_jur.SetX(416.735*ShipUnit::mm);
    SurveyPoints.top_jur.SetY(887.258*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ((16546.1-70)*ShipUnit::mm);
    SurveyPoints.top_sal.SetX(80.336*ShipUnit::mm);
    SurveyPoints.top_sal.SetY(886.136*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ((16547.5-70)*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX(78.512*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY(-688.315*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ((16538.9-70)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX(414.365*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY(-687.523*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ((16540.9-70)*ShipUnit::mm);
  
  }


  else if(module==42){
    
    SurveyPoints.top_jur.SetX(-87.554*ShipUnit::mm);
    SurveyPoints.top_jur.SetY(886.099*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ((16548-70)*ShipUnit::mm);
    SurveyPoints.top_sal.SetX(-423.404*ShipUnit::mm);
    SurveyPoints.top_sal.SetY(886.185*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ((16549.8-70)*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX(-424.753*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY(-688.486*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ((16539.6-70)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX(-88.796*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY(-689.051*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ((16538.9-70)*ShipUnit::mm);

  }

  else if(module==43){
  
    SurveyPoints.top_jur.SetX(-591.889*ShipUnit::mm);
    SurveyPoints.top_jur.SetY(886.582*ShipUnit::mm);
    SurveyPoints.top_jur.SetZ((16550.8-70)*ShipUnit::mm);
    SurveyPoints.top_sal.SetX(-928.063*ShipUnit::mm);
    SurveyPoints.top_sal.SetY(886.854*ShipUnit::mm);
    SurveyPoints.top_sal.SetZ((16551.8-70)*ShipUnit::mm);
    SurveyPoints.bot_sal.SetX(-929.54*ShipUnit::mm);
    SurveyPoints.bot_sal.SetY(-687.767*ShipUnit::mm);
    SurveyPoints.bot_sal.SetZ((16540.2-70)*ShipUnit::mm);
    SurveyPoints.bot_jur.SetX(-593.653*ShipUnit::mm);
    SurveyPoints.bot_jur.SetY(-688.852*ShipUnit::mm);
    SurveyPoints.bot_jur.SetZ((16540-70)*ShipUnit::mm);

  }


  else{

    std::cout << "WARNING: MufluxSpectrometerDTSurvey::DTSurveyMuflux could not convert detector id " << DetectorID << " to valid module number" << std::endl;
  }
  
  return SurveyPoints;
}


DTSurveyPoints MufluxSpectrometerDTSurvey::DTSurveyDistanceToRefTube(Int_t DetectorID){
  DTSurveyPoints distance; //distance to lower right tube (upper left for rotated modules)

  Int_t module=Module(DetectorID);
  double mm=ShipUnit::mm;


    
    if(module==10){
      
      distance.top_jur.SetXYZ(-11*42*mm,-170*mm-8.35*mm,30.3*mm);
      distance.top_sal.SetXYZ(-10,-170*mm-8.35*mm,30.3*mm);
      
      distance.bot_jur.SetXYZ(-11*42*mm,70*mm+8.65*mm,30.3*mm);
      distance.bot_sal.SetXYZ(-10*mm,70*mm+8.65*mm,30.3*mm);
    
    }

    else if(module==11){

      distance.top_jur.SetXYZ(-11*42*mm,-170*mm,-143.7*mm);
      distance.top_sal.SetXYZ(-10*mm,-170*mm,-143.7*mm);
      
      distance.bot_jur.SetXYZ(-11*42*mm,70*mm,-143.7*mm);
      distance.bot_sal.SetXYZ(-10*mm,70*mm,-143.7*mm);

    }

    else if(module==20){ //rotated

      distance.top_sal.SetXYZ(11*42*mm,-170*mm,143.7*mm);
      distance.top_jur.SetXYZ(10*mm,-170*mm,143.7*mm);
      
      distance.bot_sal.SetXYZ(11*42*mm,70*mm,143.7*mm);
      distance.bot_jur.SetXYZ(10*mm,70*mm,143.7*mm);

    }

    else if(module==21){ //rotated

      distance.top_sal.SetXYZ(11*42*mm,-170*mm-7.3*mm,-30.3*mm);
      distance.top_jur.SetXYZ(10,-170*mm-7.3*mm,-30.3*mm);
      
      distance.bot_sal.SetXYZ(11*42*mm,70*mm+8.45*mm,-30.3*mm);
      distance.bot_jur.SetXYZ(10*mm,70*mm+8.45*mm,-30.3*mm);

    }

    else if(module==30||module==31||module==32||module==33){

      distance.top_sal.SetXYZ(-63*mm,0,70*mm+43.3*mm);
      distance.top_jur.SetXYZ(-399*mm,0,70*mm+43.3*mm);

      distance.bot_sal.SetXYZ(-63*mm,0,70*mm+43.3*mm);
      distance.bot_jur.SetXYZ(-399*mm,0,70*mm+43.3*mm);
      
    }
  
    else if(module==40||module==41||module==42||module==43){
      
      distance.top_sal.SetXYZ(-52*mm,0,-70*mm-43.3*mm-113.4*mm);
      distance.top_jur.SetXYZ(-388*mm,0,-70*mm-43.3*mm-113.4*mm);
      
      distance.bot_sal.SetXYZ(-52*mm,0,-70*mm-43.3*mm-113.4*mm);
      distance.bot_jur.SetXYZ(-388*mm,0,-70*mm-43.3*mm-113.4*mm);
      
    }
  
    else {
      std::cout << "WARNING: MufluxSpectrometerDTSurvey::DTSurveyDistanceToRefTube could not convert detector id " << DetectorID << " to valid module number" << std::endl;
    }
    
    return distance;
}
