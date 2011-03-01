#ifndef __EXVMEASLAYER__
#define __EXVMEASLAYER__
//*************************************************************************
//* ===================
//*  EXVMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by TVTrackHit.
//* (Requires)
//* (Provides)
//*     class EXVMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*************************************************************************
//
#include "TVector3.h"
#include "kaltest/TKalMatrix.h"
#include "kaltest/TCylinder.h"
#include "kaltest/TVMeasLayer.h"
#include "kaltest/TAttDrawable.h"
#include "kaltest/KalTrackDim.h"
#include "TString.h"

class TVTrackHit;
class TNode;

namespace EVENT{
  class TrackerHit ;
} 


class EXVMeasLayer : public TVMeasLayer, public TAttDrawable {
public:
  static Bool_t kActive;
  static Bool_t kDummy;

  // Ctors and Dtor

  EXVMeasLayer(TMaterial &min,
	       TMaterial &mout,
	       Bool_t     type = EXVMeasLayer::kActive,
	       int layerID = -1 , 
	       const Char_t    *name = "MeasL");
  virtual ~EXVMeasLayer();
  
  virtual void ProcessHit(const TVector3  &xx,
			  TObjArray &hits, EVENT::TrackerHit* hit=0) = 0;
  
  virtual TVTrackHit* createHit( EVENT::TrackerHit* hit=0 ) = 0;

  inline int getLayerID(){ return _layerID ; } 
  
  inline TString GetMLName () const { return fName;    }
  
  inline TNode  *GetNodePtr() const { return fNodePtr; }
  inline void    SetNodePtr(TNode *nodep) { fNodePtr = nodep; }
  
  
private:
  TString  fName;      // layer name
  TNode   *fNodePtr;   // node pointer
  
  int _layerID ;


  //ClassDef(EXVMeasLayer,1) 	// Sample measurement layer class
};

#endif
