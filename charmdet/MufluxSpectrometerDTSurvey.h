#ifndef MUFLUXSPECTROMETERDTSURVEY_H
#define MUFLUXSPECTROMETERDTSURVEY_H

#include <map>
#include <TVector3.h>

typedef struct{
  TVector3 top_jur;
  TVector3 top_sal;
  TVector3 bot_jur;
  TVector3 bot_sal;
}DTSurveyPoints;


typedef struct{
  TVector3 top;
  TVector3 bot;
}TubeEndPoints;

static std::map<int,TubeEndPoints> EndPointMap; 

class MufluxSpectrometerDTSurvey{
  
 public:
  
  MufluxSpectrometerDTSurvey();
  ~MufluxSpectrometerDTSurvey();
  void Init();
  void Test();
  Int_t Module(Int_t DetectorID);
  DTSurveyPoints DTSurveyMuflux(Int_t DetectorID);
  DTSurveyPoints DTSurveyDistanceToRefTube(Int_t DetectorID);

 private:

  ClassDef(MufluxSpectrometerDTSurvey, 2);
  
};


#endif
