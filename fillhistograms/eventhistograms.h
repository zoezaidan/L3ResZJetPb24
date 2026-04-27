#ifndef __eventhistograms_h__
#define __eventhistograms_h__

#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

class eventhistograms {

public:
  Float_t isMC;

  // Event quality histograms
  TH1D *event_vz;
  TH1D *event_hiBin;
  // TODO: rho

  TH2D *event_pthatwsgenweight;

  TDirectory *dir;

  eventhistograms(TDirectory *dir, bool ismc);

  //  ~eventhistograms();
  //  this->blabla

  void Write();

private:
};

#endif
