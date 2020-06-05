#ifndef MUFLUXSPECTROMETERDTSURVEY_H
#define MUFLUXSPECTROMETERDTSURVEY_H

#include <map>
#include <TVector3.h>

/**
 * Struct to store the four survey points (vectors) of the modules
 *
 * @author Daniel Bick
 * @date May. 20, 2020
 */
struct DTSurveyPoints{
  TVector3 top_jur;
  TVector3 top_sal;
  TVector3 bot_jur;
  TVector3 bot_sal;

  inline DTSurveyPoints operator +(DTSurveyPoints point){
    return{
      top_jur+point.top_jur,
      top_sal+point.top_sal,
      bot_jur+point.bot_jur,
      bot_sal+point.bot_sal
    };
  };
  inline DTSurveyPoints operator -(DTSurveyPoints point){
    return{
      top_jur-point.top_jur,
      top_sal-point.top_sal,
      bot_jur-point.bot_jur,
      bot_sal-point.bot_sal
    };
  };
};


/**
 * Struct for the wire end position to be stored in a std map
 *
 * @author Daniel Bick
 * @date May. 20, 2020
 */
typedef struct{
  TVector3 top;
  TVector3 bot;
}TubeEndPoints;


static std::map<int,TubeEndPoints> DTSurveyEndPointMap; 
static bool DTSurveyIsInitialized;

/**
 * A class to query the hard coded syrvey positions of the drift tubes
 *
 * @author Daniel Bick
 * @date May. 20, 2020
 */
class MufluxSpectrometerDTSurvey{
  
 public:
  
  MufluxSpectrometerDTSurvey();
  virtual ~MufluxSpectrometerDTSurvey();
  void Init();
  void Test(Int_t DetectorID);
  Int_t Module(Int_t DetectorID);
  DTSurveyPoints DTSurveyMuflux(Int_t DetectorID);
  DTSurveyPoints DTSurveyDistanceToRefTube(Int_t DetectorID);
  TVector3 DTStaggering(Int_t DetectorID);
  void TubeEndPointsSurvey(Int_t DetectorID, TVector3 &vtop, TVector3 &vbot);
  TubeEndPoints TubeEndPointsSurvey(Int_t DetectorID);
  Double_t DTSurveyStereoAngle(Int_t DetectorID);
  Double_t DTSurveyTubeLength(Int_t DetectorID);
  TVector3 AtY(Int_t DetectorID, Double_t y);
 private:

  ClassDef(MufluxSpectrometerDTSurvey, 2);
  
};


#endif
