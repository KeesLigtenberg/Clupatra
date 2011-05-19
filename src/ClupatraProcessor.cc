#include "ClupatraProcessor.h"

#include <time.h>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <math.h>
#include <cmath>

//---- MarlinUtil 
#include "NNClusters_clupa.h"
#include "ClusterShapes.h"
#include "MarlinCED.h"

//---- LCIO ---
#include "IMPL/LCCollectionVec.h"
#include "IMPL/TrackImpl.h"
#include "IMPL/TrackerHitImpl.h"
#include "EVENT/SimTrackerHit.h"
#include "IMPL/LCFlagImpl.h"
#include "UTIL/Operators.h"
#include "UTIL/LCTOOLS.h"


//---- GEAR ----
#include "marlin/Global.h"
#include "gear/GEAR.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"
#include "gear/BField.h"

#include "LCIterator.h"

#include "KalTest.h"

using namespace lcio ;
using namespace marlin ;


typedef GenericCluster<TrackerHit> HitCluster ;
typedef GenericHit<TrackerHit>     Hit ;
typedef GenericHitVec<TrackerHit>  HitVec ;


// copy_if (missing from STL)
template <class In, class Out, class Pred> Out copy_if(In first, In last, Out res, Pred p){
  
  while( first != last ){

    if( p( *first) ){

      *res++ = first ;
      ++first ;
    }
  }
  return res ;
}


// delete helper
template<class P>  void delete_ptr(P* p) { delete p;}

/** helper class that maps array to gear::Vector3D */
struct VecFromArray{
  gear::Vector3D _v ;
  VecFromArray( const double* v) : _v( v[0] , v[1] , v[2] ) {}
  VecFromArray( const float* v) : _v( v[0] , v[1] , v[2] ) {}
  const gear::Vector3D& v(){ return _v ; }
} ;

/** helper function to restrict the range of the azimuthal angle to ]-pi,pi]*/
inline double toBaseRange( double phi){
  while( phi <= -M_PI ){  phi += 2. * M_PI ; }
  while( phi >   M_PI ){  phi -= 2. * M_PI ; }
  return phi ;
}

/** helper class to compute the chisquared of two points in rho and z coordinate */
class Chi2_RPhi_Z{
  double _sigsr, _sigsz ;
public :
  Chi2_RPhi_Z(double sigr, double sigz) : _sigsr( sigr * sigr ) , _sigsz( sigz * sigz ){}
  double operator()( const gear::Vector3D& v0, const gear::Vector3D& v1) {

    //    return (v0 - v1 ).r() ;

    //double dRPhi = v0.rho() * v0.phi() - v1.rho() * v1.phi() ;

    double dPhi = std::abs(  v0.phi() - v1.phi() )  ;
    if( dPhi > M_PI )
      dPhi = 2.* M_PI - dPhi ;

    double dRPhi =  dPhi *  v0.rho() ; 

    double dZ = v0.z() - v1.z() ;

    return  dRPhi * dRPhi / _sigsr + dZ * dZ / _sigsz  ;
  }
};

/** helper class to compute the chisquared of two points in rho and z coordinate */
struct Chi2_RPhi_Z_Hit{
  double operator()( const TrackerHit* h, const gear::Vector3D& v1) {


    gear::Vector3D v0( h->getPosition()[0] ,  h->getPosition()[1] ,  h->getPosition()[2] ) ;

    double sigsr =  sqrt( h->getCovMatrix()[0] + h->getCovMatrix()[2] ) ;
    double sigsz =  h->getCovMatrix()[5] ;
    // double sigsr =  0.01 ; 
    // double sigsz =  0.1 ;
    

    double dPhi = std::abs(  v0.phi() - v1.phi() )  ;
    if( dPhi > M_PI )
      dPhi = 2.* M_PI - dPhi ;

    double dRPhi =  dPhi *  v0.rho() ; 

    double dZ = v0.z() - v1.z() ;

    return  dRPhi * dRPhi / sigsr + dZ * dZ / sigsz  ;
  }
};

// helper class to assign additional parameters to TrackerHits
struct HitInfoStruct{
  HitInfoStruct() :layerID(-1), usedInTrack(false) {}
  int layerID ;
  bool usedInTrack ;
  double chi2Residual ;
} ;
struct HitInfo : LCOwnedExtension<HitInfo, HitInfoStruct> {} ;


//------------------------------------------------------
// function to extract position for Kaltest:
TVector3 hitPosition( Hit* h)  { 
  return TVector3( h->first->getPosition()[0],   
  		   h->first->getPosition()[1],
  		   h->first->getPosition()[2]  ) ; 
}   

// function to extract layerID from generic Hit:
int hitLayerID( const Hit* h, int offset=0) { return  h->first->ext<HitInfo>()->layerID + offset  ; } 

// functor for layer ID
class HitLayerID{
  int _off ;
  HitLayerID(){}
public:
  HitLayerID( int off) : _off(off) {}
  int operator()(const Hit* h){ return hitLayerID( h, _off) ; } 
} ;

struct LCIOTrackerHit{ EVENT::TrackerHit* operator()( Hit* h) { return h->first ; }   } ;


//---------------------------------------------------
// helper for sorting cluster wrt layerID
template <bool SortDirection>
struct LayerSort{
  bool operator()( const Hit* l, const Hit* r) {
    return hitLayerID( l ) < hitLayerID( r ) ; 
  }
} ;
template<>
struct LayerSort<KalTest::OrderIncoming>{
  bool operator()( const Hit* l, const Hit* r) {
    return hitLayerID( r ) < hitLayerID( l ) ; 
  }
} ;

//------ ordering of KalTracks 
struct KalTrackLengthSort {
  bool operator()( const KalTrack* t0, const KalTrack* t1) {
    return ( t0->getNHits() >= t1->getNHits() );
  }
};



//------------------------------
//helpers for z ordering of hits
struct TrackerHitCast{
  TrackerHit* operator()(LCObject* o) { return (TrackerHit*) o ; }
};

struct ZSort {
  bool operator()( const TrackerHit* l, const TrackerHit* r) {
    return ( l->getPosition()[2] < r->getPosition()[2] );
  }
};



void printZ(TrackerHit* h) { 
  std::cout << h->getPosition()[2] << ", " ;
  if(!( h->id() % 30 )) std::cout << std::endl ;
}



//-------------------------------
template <class T>
void delete_elements(T* t) { delete t ; }

//-------------------------------

//-------------------------------------------------------------------------
template <bool HitOrder, bool FitOrder, bool PropagateIP=false>

struct KalTestFitter{

  KalTest* _kt ; 
  
  KalTestFitter(KalTest* k) : _kt( k ) {}
  
  KalTrack* operator() (HitCluster* clu) {  
    
    static HitLayerID tpcLayerID( _kt->indexOfFirstLayer( KalTest::DetID::TPC )  )  ;
    
    clu->sort( LayerSort<HitOrder>() ) ;
    
    
    // need to reverse the order for incomming track segments (curlers)
    // assume particle comes from IP
    Hit* hf = clu->front() ;
    Hit* hb = clu->back() ;

    bool reverse_order = ( ( HitOrder ==  KalTest::OrderOutgoing ) ?    
			   ( std::abs( hf->first->getPosition()[2] ) > std::abs( hb->first->getPosition()[2]) + 3. )   :   
			   ( std::abs( hf->first->getPosition()[2] ) < std::abs( hb->first->getPosition()[2]) + 3. )   ) ;
    
    // reverse_order = false ;


    KalTrack* trk = _kt->createKalTrack() ;

    trk->setCluster<HitCluster>( clu ) ;
    

    if( PropagateIP  && HitOrder == KalTest::OrderOutgoing ) {
      
      trk->addIPHit() ;
    }  
    
    // // ----- debug ----
    // std::set<int> layers ;
    // for( HitCluster::iterator it=clu->begin() ; it != clu->end() ; ++it){
    //   if( layers.find( tpcLayerID( *it ) ) != layers.end()  )
    // 	std::cout << " +++++++++++++++++++ duplicate layerID in addHits : " <<  tpcLayerID( *it ) << std::endl ;
    //   layers.insert( tpcLayerID( *it ) ) ;
    // }
    // // ---- end debug ----------
    
    if( reverse_order )
      trk->addHits( clu->rbegin() , clu->rend() , hitPosition, tpcLayerID , LCIOTrackerHit() ) ; 
    else
      trk->addHits( clu->begin() , clu->end() , hitPosition, tpcLayerID , LCIOTrackerHit() ) ; 
    

    if( PropagateIP  && HitOrder == KalTest::OrderIncoming ) {
      
      trk->addIPHit() ;
    }  

    trk->fitTrack( FitOrder  ) ;
    
    return trk;
  }
};


struct KalTrack2LCIO{
  TrackImpl* operator() (KalTrack* trk) {  
    TrackImpl* lTrk = new TrackImpl ;
    trk->toLCIOTrack( lTrk  ) ;
    return lTrk ;
  }
};

//-------------------------------------------------------------------------
template <class T>
class RCut {
public:
  RCut( double rcut ) : _rcut( rcut ) {}  
  
  // bool operator() (T* hit) {  // DEBUG ....
  //   return   std::abs( hit->getPosition()[2] ) > 2000. ;
  bool operator() (T* hit) {  
    return  ( (std::sqrt( hit->getPosition()[0]*hit->getPosition()[0] +
			  hit->getPosition()[1]*hit->getPosition()[1] )   > _rcut )   ||
	      ( std::abs( hit->getPosition()[2] ) > (500. + _rcut ) )
	      ); 
  }
protected:
  RCut() {} ;
  double _rcut ;
} ;

template <class T>
class RCutInverse {
public:
  RCutInverse( double rcut ) : _rcut( rcut ) {}  
  
  bool operator() (T* hit) {  
    return (  ( std::sqrt( hit->getPosition()[0]*hit->getPosition()[0] +
			   hit->getPosition()[1]*hit->getPosition()[1] )   <= _rcut )   &&
	      (  std::abs( hit->getPosition()[2] ) <= (500. + _rcut ) )
	      ) ;

  }
protected:
  RCutInverse() {} ;
  double _rcut ;
} ;

// template <class T>
// class RCut {
// public:
//   RCut( double rcut ) : _rcut( rcut ) {}  
  
//   // bool operator() (T* hit) {  // DEBUG ....
//   //   return   std::abs( hit->getPosition()[2] ) > 2000. ;
//   bool operator() (T* hit) {  
//     return   std::sqrt( hit->getPosition()[0]*hit->getPosition()[0] +
// 			hit->getPosition()[1]*hit->getPosition()[1] )   > _rcut ; 
//   }
// protected:
//   RCut() {} ;
//   double _rcut ;
// } ;

// template <class T>
// class RCutInverse {
// public:
//   RCutInverse( double rcut ) : _rcut( rcut ) {}  
  
//   bool operator() (T* hit) {  
//     return   std::sqrt( hit->getPosition()[0]*hit->getPosition()[0] +
// 			hit->getPosition()[1]*hit->getPosition()[1] )   <= _rcut ; 
//   }
// protected:
//   RCutInverse() {} ;
//   double _rcut ;
// } ;

//---------------------------------------------------------------------------------

/** Predicate class for identifying clusters with duplicate pad rows - returns true
 * if the fraction of duplicate hits is larger than 'fraction'.
 */
struct DuplicatePadRows{

  unsigned _N ;
  float _f ; 
  DuplicatePadRows(unsigned nPadRows, float fraction) : _N( nPadRows), _f( fraction )  {}

  bool operator()(const HitCluster* cl) const {
 
    // check for duplicate layer numbers
    std::vector<int> hLayer( _N )  ; 
    typedef HitCluster::const_iterator IT ;

    unsigned nHit = 0 ;
    for(IT it=cl->begin() ; it != cl->end() ; ++it ) {
      TrackerHit* th = (*it)->first ;
      ++ hLayer[ th->ext<HitInfo>()->layerID ]   ;
      ++ nHit ;
    } 
    unsigned nDuplicate = 0 ;
    for(unsigned i=0 ; i < _N ; ++i ) {
      if( hLayer[i] > 1 ) 
     	nDuplicate += hLayer[i] ;
    }
    return double(nDuplicate)/nHit > _f ;
  }
};
//TODO: create a faster predicate for no duplicate pad rows ....

//---------------------------------------------------------------------------------

/** Predicate class for 'distance' of NN clustering.
 */
//template <class HitClass, typename PosType > 
class HitDistance{
  typedef TrackerHit HitClass ;
  typedef double PosType ;
public:

  /** Required typedef for cluster algorithm 
   */
  typedef HitClass hit_type ;

  /** C'tor takes merge distance */
  HitDistance(float dCut) : _dCutSquared( dCut*dCut ) , _dCut(dCut)  {} 


  /** Merge condition: true if distance  is less than dCut given in the C'tor.*/ 
  inline bool mergeHits( GenericHit<HitClass>* h0, GenericHit<HitClass>* h1){
    
    if( std::abs( h0->Index0 - h1->Index0 ) > 1 ) return false ;
    
    //     int l0 =  h0->first->ext<HitInfo>()->layerID ;
    //     int l1 =  h1->first->ext<HitInfo>()->layerID ;

    //     //------- don't merge hits from same layer !
    //     if( l0 == l1 )
    //       return false ;

    if( h0->first->ext<HitInfo>()->layerID == h1->first->ext<HitInfo>()->layerID )
      return false ;

    const PosType* pos0 =  h0->first->getPosition() ;
    const PosType* pos1 =  h1->first->getPosition() ;
    
    return 
      ( pos0[0] - pos1[0] ) * ( pos0[0] - pos1[0] ) +
      ( pos0[1] - pos1[1] ) * ( pos0[1] - pos1[1] ) +
      ( pos0[2] - pos1[2] ) * ( pos0[2] - pos1[2] )   
      < _dCutSquared ;
  }
  
protected:
  HitDistance() ;
  float _dCutSquared ;
  float _dCut ;
} ;

class HitDistance_2{
  typedef TrackerHit HitClass ;
  typedef double PosType ;
public:

  /** Required typedef for cluster algorithm 
   */
  typedef HitClass hit_type ;

  /** C'tor takes merge distance */
  HitDistance_2(float dCut) : _dCutSquared( dCut*dCut ) , _dCut(dCut)  {} 


  /** Merge condition: true if distance  is less than dCut given in the C'tor.*/ 
  inline bool mergeHits( GenericHit<HitClass>* h0, GenericHit<HitClass>* h1){
    
    //if( std::abs( h0->Index0 - h1->Index0 ) > 1 ) return false ;

    int l0 =  h0->first->ext<HitInfo>()->layerID ;
    int l1 =  h1->first->ext<HitInfo>()->layerID ;


    //------- don't merge hits from same layer !
    if( l0 == l1 )
      return false ;


    const PosType* pos0 =  h0->first->getPosition() ;
    const PosType* pos1 =  h1->first->getPosition() ;
    
    return inRange<-2,2>(  l0 - l1 )  &&  
      ( pos0[0] - pos1[0] ) * ( pos0[0] - pos1[0] ) +
      ( pos0[1] - pos1[1] ) * ( pos0[1] - pos1[1] ) +
      ( pos0[2] - pos1[2] ) * ( pos0[2] - pos1[2] )   
      < _dCutSquared ;
  }
  
protected:
  HitDistance_2() ;
  float _dCutSquared ;
  float _dCut ;
} ;


template <class T>
struct LCIOTrack{
  
  lcio::Track* operator() (GenericCluster<T>* c) {  
    
    TrackImpl* trk = new TrackImpl ;
    
    double e = 0.0 ;
    int nHit = 0 ;
    for( typename GenericCluster<T>::iterator hi = c->begin(); hi != c->end() ; hi++) {
      
      trk->addHit(  (*hi)->first ) ;
      e += (*hi)->first->getEDep() ;
      nHit++ ;
    }

   
    trk->setdEdx( e/nHit ) ;
    trk->subdetectorHitNumbers().push_back( 1 ) ;  // workaround for bug in lcio::operator<<( Tracks ) - used for picking ....
 
    // FIXME - these are no meaningfull tracks - just a test for clustering tracker hits
    return trk ;
  }

} ;

// helper for creating lcio header for short printout
template <class T> 
const std::string & myheader(){return header(*(T*)(0)); }


void printTrackShort(const LCObject* o){
  
  const Track* trk = dynamic_cast<const Track*> (o) ; 
  
  if( o == 0 ) {
    
    streamlog_out( ERROR ) << "  printTrackShort : dynamic_cast<Track*> failed for LCObject : " << o << std::endl ;
    return  ;
  }
  
  streamlog_out( MESSAGE ) << myheader<Track>()  
			   << lcshort( trk ) << std::endl  ;
  
  
  double r0 = 1. / trk->getOmega() ;
  double d0 = trk->getD0() ;
  double p0 = trk->getPhi() ;
  
  double x0 = ( r0 - d0 ) * sin( p0 ) ;
  double y0 = ( d0 - r0 ) * cos( p0 ) ;
  
  streamlog_out( MESSAGE ) << " circle: r = " << r0 << ", xc = " << x0 << " , yc = " << y0 << std::endl ;
    
}

void printTrackerHit(const LCObject* o){
  
  TrackerHit* trk = const_cast<TrackerHit*> ( dynamic_cast<const TrackerHit*> (o) ) ; 
  
  if( o == 0 ) {
    
    streamlog_out( ERROR ) << "  printTrackerHit : dynamic_cast<TrackerHit*> failed for LCObject : " << o << std::endl ;
    return  ;
  }
  
  streamlog_out( MESSAGE ) << *trk << std::endl 
			   << " err: rPhi" <<  sqrt( trk->getCovMatrix()[0] + trk->getCovMatrix()[2] ) 
			   << " z :  " <<   trk->getCovMatrix()[5] << std::endl 
			   << " chi2 residual to best matching track : " << trk->ext<HitInfo>()->chi2Residual << std::endl ;

    
}


/** Predicate class for track merging with NN clustering.
 */
class TrackStateDistance{
  typedef Track HitClass ;

public:
  
  /** Required typedef for cluster algorithm 
   */
  typedef HitClass hit_type ;
  
  /** C'tor takes merge distance */
  TrackStateDistance(float dCut) : _dCutSquared( dCut*dCut ) , _dCut(dCut)  {} 


  /** Merge condition: true if distance  is less than dCut given in the C'tor.*/ 
  inline bool mergeHits( GenericHit<HitClass>* h0, GenericHit<HitClass>* h1){
    
    if( std::abs( h0->Index0 - h1->Index0 ) > 1 ) return false ;
    
    Track* trk0 = h0->first ;
    Track* trk1 = h1->first ;

    //------- dont' merge complete tracks: ------------------------
    unsigned nh0 = trk0->getTrackerHits().size() ;
    unsigned nh1 = trk1->getTrackerHits().size() ;
    if( nh0 > 220 ) return false ;
    if( nh1 > 220 ) return false ;
    //------------------------------------------------------


    KalTrack* ktrk0 = h0->first->ext<KalTrackLink>() ; 
    KalTrack* ktrk1 = h1->first->ext<KalTrackLink>() ; 

    //---- sanity check on radii-------------------
    double r0 = 1. / trk0->getOmega() ;
    double r1 = 1. / trk1->getOmega() ;

    if( r0 < 300. || r1 < 300. )
      return false ;

    if( std::abs( r0 - r1 ) / std::abs( r0 + r1 )  > 0.02 )  // relative difference larger than 1%
      return false ;
    //---------------------------------------------


    double chi2  = KalTrack::chi2( *ktrk0 , *ktrk1 ) ; 	   

    return chi2  < _dCut ;
      
  }
  
protected:
  TrackStateDistance() ;
  float _dCutSquared ;
  float _dCut ;
} ;



/** helper class for merging track segments, based on circle (and tan lambda) */
class TrackCircleDistance{
  typedef Track HitClass ;

public:
  
  /** Required typedef for cluster algorithm
   */
  typedef HitClass hit_type ;
  
  
  /** C'tor takes merge distance */
  TrackCircleDistance(float dCut) : _dCutSquared( dCut*dCut ) , _dCut(dCut){}
  
  /** Merge condition: ... */
  inline bool mergeHits( GenericHit<HitClass>* h0, GenericHit<HitClass>* h1){
    
    static const double DRMAX = 0.1 ; // make parameter
    static const double DTANLMAX = 0.2 ; //   " 

    if( std::abs( h0->Index0 - h1->Index0 ) > 1 ) return false ;
    
    Track* trk0 = h0->first ;
    Track* trk1 = h1->first ;

    // //------- dont' merge complete tracks: ------------------------
    // unsigned nh0 = trk0->getTrackerHits().size() ;
    // unsigned nh1 = trk1->getTrackerHits().size() ;
    // if( nh0 > 220 ) return false ;
    // if( nh1 > 220 ) return false ;
    // //------------------------------------------------------

    // KalTrack* ktrk0 = h0->first->ext<KalTrackLink>() ; 
    // KalTrack* ktrk1 = h1->first->ext<KalTrackLink>() ; 

    double tl0 = trk0->getTanLambda() ;
    double tl1 = trk1->getTanLambda() ;

    double dtl = 2. * ( tl0 - tl1 ) / ( tl0 + tl1 ) ;
    dtl *= dtl ;
    if(  dtl >  DTANLMAX * DTANLMAX ) 
      return false ;

    double r0 = 1. / trk0->getOmega()  ;
    double r1 = 1. / trk1->getOmega()  ;


    double d0 = trk0->getD0() ;
    double d1 = trk1->getD0() ;

    double p0 = trk0->getPhi() ;
    double p1 = trk1->getPhi() ;

    double x0 = ( r0 - d0 ) * sin( p0 ) ;
    double x1 = ( r1 - d1 ) * sin( p1 ) ;

    double y0 = ( d0 - r0 ) * cos( p0 ) ;
    double y1 = ( d1 - r1 ) * cos( p1 ) ;
    
    double dr = 2. * std::abs( ( r0 -r1 ) / (r0 + r1 ) ) ;

    double distMS = sqrt ( ( x0 - x1 ) * ( x0 - x1 ) + ( y0 - y1 ) * ( y0 - y1 )  ) ;
    
    
    return ( dr < DRMAX && distMS < _dCut*r0 ) ;

  }
  
protected:
  float _dCutSquared ;
  float _dCut ;
} ; 


ClupatraProcessor aClupatraProcessor ;


ClupatraProcessor::ClupatraProcessor() : Processor("ClupatraProcessor") {
  
  // modify processor description
  _description = "ClupatraProcessor : simple nearest neighbour clustering" ;
  
  
  StringVec colDefault ;
  colDefault.push_back("AllTPCTrackerHits" ) ;

  registerInputCollections( LCIO::TRACKERHIT,
			    "HitCollections" , 
			    "Name of the input collections"  ,
			    _colNames ,
			    colDefault ) ;
  
  registerOutputCollection( LCIO::TRACK,
			    "OutputCollection" , 
			    "Name of the output collections"  ,
			    _outColName ,
			    std::string("CluTracks" ) ) ;
  
  
  registerProcessorParameter( "DistanceCut" , 
			      "Cut for distance between hits in mm"  ,
			      _distCut ,
			      (float) 40.0 ) ;
  
  registerProcessorParameter( "MinimumClusterSize" , 
			      "minimum number of hits per cluster"  ,
			      _minCluSize ,
			      (int) 3) ;
  

  registerProcessorParameter( "DuplicatePadRowFraction" , 
			      "allowed fraction of hits in same pad row per track"  ,
			      _duplicatePadRowFraction,
			      (float) 0.01 ) ;

  registerProcessorParameter( "RCut" , 
 			      "Cut for r_min in mm"  ,
 			      _rCut ,
 			      (float) 0.0 ) ;
  
}


void ClupatraProcessor::init() { 

  // usually a good idea to
  printParameters() ;


  _kalTest = new KalTest( *marlin::Global::GEAR ) ;

  _kalTest->setOption( KalTest::CFG::ownsHits , true ) ;

  _kalTest->init() ;

  _nRun = 0 ;
  _nEvt = 0 ;
}

void ClupatraProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void ClupatraProcessor::processEvent( LCEvent * evt ) { 

  clock_t start =  clock() ; 

  GenericHitVec<TrackerHit> h ;
  GenericHitVec<TrackerHit> hSmallR ; 
  
  GenericClusterVec<TrackerHit> cluList ;
  
  RCut<TrackerHit> rCut( _rCut ) ;
  RCutInverse<TrackerHit> rCutInverse( _rCut ) ;
  
  ZIndex<TrackerHit,200> zIndex( -2750. , 2750. ) ; 
  
  //  NNDistance< TrackerHit, double> dist( _distCut )  ;
  HitDistance dist0( _distCut ) ;
  HitDistance dist( 20. ) ;
  //  HitDistance_2 dist_2( 20. ) ;
  
  LCIOTrack<TrackerHit> converter ;
  
  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;
  unsigned nPadRows = padLayout.getNRows() ;

  // create a vector of generic hits from the collection applying a cut on r_min
  for( StringVec::iterator it = _colNames.begin() ; it !=  _colNames.end() ; it++ ){  
    
    LCCollectionVec* col =  dynamic_cast<LCCollectionVec*> (evt->getCollection( *it )  ); 
    
    
    //--- assign the layer number to the TrackerHits
    
    int nHit = col->getNumberOfElements() ;
    for(int i=0 ; i < nHit ; ++i ) {
      
      TrackerHitImpl* th = (TrackerHitImpl*) col->getElementAt(i) ;
      gear::Vector3D v( th->getPosition()[0],th->getPosition()[1], 0 ) ; 
      int padIndex = padLayout.getNearestPad( v.rho() , v.phi() ) ;
      
      th->ext<HitInfo>() = new HitInfoStruct ;

      th->ext<HitInfo>()->layerID = padLayout.getRowNumber( padIndex ) ;
      

      //      std::cout << " layer ID " <<  th->ext<HitInfo>()->layerID << " ..... " << th->getType() << std::endl ;

      //       //--- for fixed sized rows this would also work...
      //       float rMin = padLayout.getPlaneExtent()[0] ;
      //       float rMax = padLayout.getPlaneExtent()[1] ;
      //       float nRow  = padLayout.getNRows() ;
      //       int lCheck =  ( v.rho() - rMin ) / ((rMax - rMin ) /nRow ) ;

      //       streamlog_out( DEBUG ) << " layerID : " << th->ext<HitInfo>()->layerID 
      // 			     << " r: " << v.rho() 
      // 			     << " lCheck : " << lCheck 
      // 			     << " phi : " << v.phi()
      // 			     << " rMin : " << rMin 
      // 			     << " rMax : " << rMax 
      // 			     << std::endl ;

    } //-------------------- end assign layernumber ---------
    
    //addToGenericHitVec( h , col , rCut , zIndex ) ;
    std::list< TrackerHit*> hitList ;
    TrackerHitCast cast ;
    ZSort zsort ;
    std::transform(  col->begin(), col->end(), std::back_inserter( hitList ), cast ) ;

    hitList.sort( zsort ) ;
    //    std::for_each( hitList.begin() , hitList.end() , printZ ) ;

    addToGenericHitVec( h, hitList.begin() , hitList.end() , rCut ,  zIndex ) ;

    // create a vector with the hits at smaller R
    addToGenericHitVec( hSmallR, hitList.begin() , hitList.end() , rCutInverse ,  zIndex ) ;
  }  
  
  // cluster the sorted hits  ( if |diff(z_index)|>1 the loop is stopped)
  cluster_sorted( h.begin() , h.end() , std::back_inserter( cluList )  , &dist0 , _minCluSize ) ;
  //cluster( h.begin() , h.end() , std::back_inserter( cluList )  , &dist , _minCluSize ) ;
  
  streamlog_out( DEBUG ) << "   ***** clusters: " << cluList.size() << std::endl ; 

  LCCollectionVec* allClu = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform(cluList.begin(), cluList.end(), std::back_inserter( *allClu ) , converter ) ;
  evt->addCollection( allClu , "AllTrackClusters" ) ;



  // find 'odd' clusters that have duplicate hits in pad rows
  GenericClusterVec<TrackerHit> ocs ;

  split_list( cluList, std::back_inserter(ocs),  DuplicatePadRows( nPadRows, _duplicatePadRowFraction  ) ) ;


  LCCollectionVec* oddCol = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform( ocs.begin(), ocs.end(), std::back_inserter( *oddCol ) , converter ) ;
  evt->addCollection( oddCol , "OddClu_1" ) ;


  streamlog_out( DEBUG ) << "   ***** clusters: " << cluList.size() 
			 << "   ****** oddClusters " << ocs.size() 
			 << std::endl ; 



  //-------------------- split up cluster with duplicate rows 

  GenericClusterVec<TrackerHit> sclu ; // new split clusters

  std::vector< GenericHit<TrackerHit>* > oddHits ;
  oddHits.reserve( h.size() ) ;

  typedef GenericClusterVec<TrackerHit>::iterator GCVI ;


  //========================== first iteration ================================================
  for( GCVI it = ocs.begin() ; it != ocs.end() ; ++it ){
    (*it)->takeHits( std::back_inserter( oddHits )  ) ;
    delete (*it) ;
  }
  ocs.clear() ;

  int _nRowForSplitting = 10 ; //FIXME:  make proc param
  // reset the hits index to row ranges for reclustering
  unsigned nOddHits = oddHits.size() ;
  for(unsigned i=0 ; i< nOddHits ; ++i){
    int layer =  oddHits[i]->first->ext<HitInfo>()->layerID  ;
    oddHits[i]->Index0 =   2 * int( layer / _nRowForSplitting ) ;
  }

  //----- recluster in pad row ranges
  cluster( oddHits.begin(), oddHits.end() , std::back_inserter( sclu ), &dist , _minCluSize ) ;

  LCCollectionVec* oddCol2 = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform( sclu.begin(), sclu.end(), std::back_inserter( *oddCol2 ) , converter ) ;
  evt->addCollection( oddCol2 , "OddClu_2" ) ;


  streamlog_out( DEBUG ) << "   ****** oddClusters fixed" << sclu.size() 
			 << std::endl ; 
 
  //--------- remove pad row range clusters where merge occured 
  split_list( sclu, std::back_inserter(ocs), DuplicatePadRows( nPadRows, _duplicatePadRowFraction  ) ) ;


  LCCollectionVec* oddCol3 = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform( ocs.begin(), ocs.end(), std::back_inserter( *oddCol3 ) , converter ) ;
  evt->addCollection( oddCol3 , "OddClu_3" ) ;


  for( GCVI it = ocs.begin() ; it != ocs.end() ; ++it ){
    (*it)->takeHits( std::back_inserter( oddHits )  ) ;
    delete (*it) ;
  }
  ocs.clear() ;


  //   //========================== second iteration in shifted pad row ranges ================================================


  oddHits.clear() ;
  for( GCVI it = sclu.begin() ; it != sclu.end() ; ++it ){
    (*it)->takeHits( std::back_inserter( oddHits )  ) ;
    delete (*it) ;
  }
  sclu.clear() ;

  //  int _nRowForSplitting = 10 ; //FIXME:  make proc param
  // reset the hits index to row ranges for reclustering
  nOddHits = oddHits.size() ;

  streamlog_out( DEBUG ) << "   left over odd hits for second iteration of pad row range clustering " << nOddHits << std::endl ;

  for(unsigned i=0 ; i< nOddHits ; ++i){
    int layer =  oddHits[i]->first->ext<HitInfo>()->layerID  ;
    oddHits[i]->Index0 =  2 * int( 0.5 +  ( (float) layer / (float) _nRowForSplitting ) ) ;
  }
  
  //----- recluster in pad row ranges
  cluster( oddHits.begin(), oddHits.end() , std::back_inserter( sclu ), &dist , _minCluSize ) ;

  LCCollectionVec* oddCol2_1 = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform( sclu.begin(), sclu.end(), std::back_inserter( *oddCol2_1 ) , converter ) ;
  evt->addCollection( oddCol2_1 , "OddClu_2_1" ) ;


  streamlog_out( DEBUG ) << "   ****** oddClusters fixed" << sclu.size() 
			 << std::endl ; 
 
  //--------- remove pad row range clusters where merge occured 
  split_list( sclu, std::back_inserter(ocs), DuplicatePadRows( nPadRows, _duplicatePadRowFraction  ) ) ;


  LCCollectionVec* oddCol3_1 = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform( ocs.begin(), ocs.end(), std::back_inserter( *oddCol3_1 ) , converter ) ;
  evt->addCollection( oddCol3_1 , "OddClu_3_1" ) ;

  //----------------end  split up cluster with duplicate rows 
  
  for( GCVI it = ocs.begin() ; it != ocs.end() ; ++it ){
    (*it)->takeHits( std::back_inserter( oddHits )  ) ;
    delete (*it) ;
  }
  ocs.clear() ;
  

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  // --- recluster the good clusters w/ all pad rows

  oddHits.clear() ;
  for( GCVI it = sclu.begin() ; it != sclu.end() ; ++it ){
    (*it)->takeHits( std::back_inserter( oddHits )  ) ;
    delete (*it) ;
  }
  sclu.clear() ;

  //   reset the index for 'good' hits coordinate again...
  nOddHits = oddHits.size() ;
  for(unsigned i=0 ; i< nOddHits ; ++i){
    oddHits[i]->Index0 = zIndex ( oddHits[i]->first ) ;
  }

  cluster( oddHits.begin(), oddHits.end() , std::back_inserter( sclu ), &dist , _minCluSize ) ;

  LCCollectionVec* oddCol4 = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform( sclu.begin(), sclu.end(), std::back_inserter( *oddCol4 ) , converter ) ;
  evt->addCollection( oddCol4 , "OddClu_4" ) ;

  // --- end recluster the good clusters w/ all pad rows

  // merge the good clusters to final list
  cluList.merge( sclu ) ;
  
  LCCollectionVec* cluCol = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform( cluList.begin(), cluList.end(), std::back_inserter( *cluCol ) , converter ) ;
  evt->addCollection( cluCol , "CluTrackSegments" ) ;


  //DEBUG ..... check if there are really no duplicate pad rows ...
  ocs.clear() ; 
  split_list( cluList, std::back_inserter(ocs), DuplicatePadRows( nPadRows, _duplicatePadRowFraction  ) ) ;
  LCCollectionVec* dupCol = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform( ocs.begin(), ocs.end(), std::back_inserter( *dupCol ) , converter ) ;
  evt->addCollection( dupCol , "DuplicatePadRowCluster" ) ;

  streamlog_out( DEBUG ) << "   DuplicatePadRowCluster.sizE() : " << dupCol->getNumberOfElements() << std::endl ;

  //================================================================================
   
  // create vector with left over hits
  std::vector< Hit* > leftOverHits ;
  leftOverHits.reserve(  h.size() ) ;

  typedef HitVec::const_iterator GHVI ;

  for( GHVI it = h.begin(); it != h.end() ; ++it ){

    if ( (*it)->second == 0 ) leftOverHits.push_back( *it ) ;
  }

  // add all hits that failed the rcut 
  std::copy( hSmallR.begin() , hSmallR.end() , std::back_inserter( leftOverHits )  ) ;


  //  GenericClusterVec<TrackerHit> mergedClusters ; // new split clusters


  //*********************************************************
  //   run KalTest on track segments (clusters)
  //*********************************************************

  streamlog_out( DEBUG ) <<  "************* fitted segments and KalTest tracks : **********************************" 
			 << std::endl ;


  std::list< KalTrack* > ktracks ;

  //  KalTestFitter<KalTest::OrderIncoming, KalTest::FitForward > fitter( _kalTest ) ;
  KalTestFitter<KalTest::OrderOutgoing, KalTest::FitBackward > fitter( _kalTest ) ;
    
  std::transform( cluList.begin(), cluList.end(), std::back_inserter( ktracks ) , fitter ) ;
  
  std::for_each( ktracks.begin(), ktracks.end(), std::mem_fun( &KalTrack::findXingPoints ) ) ;
  
  
  LCCollectionVec* trksegs = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform( ktracks.begin(), ktracks.end(), std::back_inserter( *trksegs ) , KalTrack2LCIO() ) ;
  evt->addCollection( trksegs , "KalTrackSegments" ) ;

  

  //===========  merge track segments based on xing points ==================================================



  //=========== assign left over hits ... ==================================================================
  
  static const bool use_best_track = false ;

  if( use_best_track ) {

    streamlog_out( DEBUG ) << "  ------ assign left over hits - best matching track for every hit ..."  << std::endl ;

    Chi2_RPhi_Z ch2rz( 0.1 , 1. ) ; // fixme - need proper errors ....
    Chi2_RPhi_Z_Hit ch2rzh ;


    HitLayerID  tpcLayerID( _kalTest->indexOfFirstLayer( KalTest::DetID::TPC )  ) ;
    
    for( GHVI ih = leftOverHits.begin() ; ih != leftOverHits.end() ; ++ih ){
      
      Hit* hit = *ih ;
      VecFromArray hPos(  hit->first->getPosition() ) ;
      
      double ch2Min = 999999999999999. ;
      KalTrack* bestTrk = 0 ;
      
      for( std::list< KalTrack* >::iterator it = ktracks.begin() ; it != ktracks.end() ; ++it ){
	
	const gear::Vector3D* kPos = (*it)->getXingPointForLayer( tpcLayerID( hit ) ) ;
	
	// double rh  =  hPos.v().rho() ;
	// double rk  =  kPos->rho() ;
	// if( std::abs( rh - rk ) > 0.1 ) {
	// 	streamlog_out( WARNING ) << " --- different radii for hit and crossing point : " <<  tpcLayerID( hit ) << ": " << rh << " - " << rk 
	// 				 <<  *kPos  << std::endl ;
	// } 
	
	if( kPos != 0 ){
	  
	  //	  double ch2 = ch2rz( hPos.v() , *kPos )  ;
	  double ch2 = ch2rzh( hit->first , *kPos )  ;
	  
	  if( ch2 < ch2Min ){
	    
	    ch2Min = ch2 ;
	    bestTrk = *it ;
	  }
	  
	}
	
	// else {
	// 	streamlog_out( MESSAGE ) << " --- no crossing point found for layer : " <<  tpcLayerID( hit ) << ": " << hPos.v() << std::endl ;
	// }
	
      }
      if( bestTrk ) {
	
	const gear::Vector3D* kPos = bestTrk->getXingPointForLayer( tpcLayerID( hit ) ) ;
	
	// double rh  =  hPos.v().rho() ;
	// double rk  =  kPos->rho() ;
	// if( std::abs( rh - rk ) > 0.1 ) {
	// 	streamlog_out( WARNING ) << "  different radii for hit and crossing point : " << rh << " - " << rk << std::endl ;
	// } 
	
	//      if( std::abs( hPos.v().rho() - kPos->rho() ) < 0.5 &&   std::abs( hPos.v().z() - kPos->z() ) < 5. ) {
	
	if(  (  hPos.v() - *kPos ).r()  < 3. ) {   // check for bad outliers... FIXME: need proper criterion here .....
	  
	  
	  HitCluster* clu = bestTrk->getCluster< HitCluster >() ;
	  
	  streamlog_out( DEBUG ) << " ---- assigning left over hit : " << hPos.v() << " <-> " << *kPos  
				 <<   " dist: " <<  (  hPos.v() - *kPos ).r()  << std::endl ;
	  
	  clu->addHit( hit ) ;
	}	
	else 
	  streamlog_out( DEBUG ) << " ---- NOT assigning left over hit : " << hPos.v() << " <-> " << *kPos << std::endl ;
      }
      else
	streamlog_out( DEBUG ) << " ---- NO best track found ??? ---- " << std::endl ;
      
    }
    

    //        ==========================================================================================
  } else { // ================== use best matching hit for every track segment =========================
    //        ==========================================================================================
    


    streamlog_out( DEBUG1 ) << "  ------ assign left over hits - best matching hit for every track ..."  << std::endl ;
    
    HitLayerID  tpcLayerID( _kalTest->indexOfFirstLayer( KalTest::DetID::TPC )  ) ;
    

    //------------- create vector of left over hits per layer
    typedef std::list<Hit*> HitList ;
    typedef std::vector< HitList > HitListVector ;
    HitListVector hitsInLayer( _kalTest->maxLayerIndex() ) ;
    
    
    for( GHVI ih = leftOverHits.begin() ; ih != leftOverHits.end() ; ++ih ) {
      
      Hit* hit = *ih ;
      //      std::cout << " ++++++  layerId: " << tpcLayerID( hit ) << " max layer index : " <<  _kalTest->maxLayerIndex() << std::endl  ;
      hitsInLayer[ tpcLayerID( hit ) ].push_back( hit )  ;
    }
    //-----------------------------
    
    std::map< HitCluster* , KalTrack* > clu2trkMap ;

    const bool use_segment_hits = false ; //true ;
    
    if( use_segment_hits  ){
      
      // store first and last hit of every segment in map with leftover hits in this layer
      
      for( GenericClusterVec<TrackerHit>::iterator icv = cluList.begin() ; icv != cluList.end() ; ++ icv ) {
	
	Hit* h0 = (*icv)->front() ;
	Hit* h1 = (*icv)->back() ;
	
	hitsInLayer[ tpcLayerID( h0 ) ].push_back( h0 )  ;
	hitsInLayer[ tpcLayerID( h1 ) ].push_back( h1 )  ;
      }
      
      // sort the tracks wrt. lenghts (#hits)
      ktracks.sort( KalTrackLengthSort() ) ;

      // store assoaciation between cluster and track 
      for( std::list< KalTrack* >::iterator it = ktracks.begin() ; it != ktracks.end() ; ++it ){
	HitCluster* c = (*it)->getCluster< HitCluster >() ;
	clu2trkMap[ c ] = *it ;
      }	   
    }
    //-------------------------------
    

    //    Chi2_RPhi_Z ch2rz( 0.1 , 1. ) ; // fixme - need proper errors 
    Chi2_RPhi_Z_Hit  ch2rzh ;
    

    for( std::list< KalTrack* >::iterator it = ktracks.begin() ; it != ktracks.end() ; ++it ){
      
      KalTrack* theTrack = *it ;
      if( theTrack == 0 ) 
	continue ;
      
      
      // ----- define chi2 cut    ~15 for 1 GeV pt 
      double chi2Cut = 100000. / ( std::log(1.) - std::log( std::abs(theTrack->getOmega()) ) ) ;


      streamlog_out( DEBUG3 ) << " ------- searching for leftover hits for track : " << theTrack 
			      << "   chi2 cut : " << chi2Cut  << " -  omega : " << theTrack->getOmega() <<  std::endl ;
      
      int xpLayer = 0 ;
      
      // const PointList& xptList = theTrack->getXingPoints() ;
      // for(PointList::const_iterator itXP = xptList.begin() ; itXP != xptList.end() ; ++itXP , xpLayer++ ) {
      // 	const gear::Vector3D* kPos =  *itXP ;
      
      PointList& xpVec = theTrack->getXingPoints() ;
      for( unsigned ixp=0 ; ixp < xpVec.size() ; ++ixp, xpLayer++  ) {
      	const gear::Vector3D* kPos =  xpVec[ixp]  ;
	
	if( kPos == 0 ) {   // we don't have a xing point
	  continue ;
	}
	
       	double ch2Min = 10e99 ;
	Hit* bestHit = 0 ;
	
	HitList& hLL = hitsInLayer.at( xpLayer ) ;
	
	for( HitList::const_iterator ih = hLL.begin() ; ih != hLL.end() ; ++ih ){
	  
	  Hit* hit = *ih ;
	  
	  //VecFromArray hPos(  hit->first->getPosition() ) ;
	  //double ch2 = ch2rz( hPos.v() , *kPos )  ;
	  double ch2 = ch2rzh( hit->first , *kPos )  ;

	  if( ch2 < ch2Min ){
	    
	    ch2Min = ch2 ;
	    bestHit = hit ;
	  }
	}
	
	if( bestHit != 0 ) {
	  
	  VecFromArray hPos(  bestHit->first->getPosition() ) ;
	  
	  //	  if( ch2Min  <  6. ) { // Sum( pdf(ch2,ndf==2) )_0^6 ~ 95% )
	  //	  if( ch2Min  <  20. ) { // Sum( pdf(ch2,ndf==2) )_0^20 ~ 99.x% ?? ) // FIXME: need steering parameter and optimize value
	  
	  
	  bestHit->first->ext<HitInfo>()->chi2Residual = ch2Min ;


	  if( ch2Min  < chi2Cut ) { 
	    
	    streamlog_out( DEBUG1 ) <<   " ---- assigning left over hit : " << hPos.v() << " <-> " << *kPos
				    <<   " dist: " <<  (  hPos.v() - *kPos ).r()
				    <<   " chi2: " <<  ch2Min 
				    <<   "  hit errors :  rphi=" <<  sqrt( bestHit->first->getCovMatrix()[0] + bestHit->first->getCovMatrix()[2] ) 
				    <<	 "  z= " <<  sqrt( bestHit->first->getCovMatrix()[5] )
				    << std::endl ;
	    
	    
	    if( bestHit->second != 0 ) { //--------------------------------------------------------------------------------
	      
	      // hit is already part of a track segment 
	      
	      
	      HitCluster* c = bestHit->second  ;
	      KalTrack* trk = clu2trkMap[ c ] ;
	      

	      if( trk == theTrack ) {
		streamlog_out( ERROR ) << " =======================  found best matching hit from track itself: " << *bestHit->first 
				       <<     std::endl  
				       <<  "      track has  " << trk->getNHits()  << " hits " << std::endl ;

		for( unsigned ii=0 ; ii < xpVec.size() ; ++ii) {
		  if( xpVec[ii] ) 
		    streamlog_out( ERROR ) << "  xing pt : "  << ii << " - " << *xpVec[ii]  ;
		}
		
		
		for( HitCluster::iterator its = c->begin(); its != c->end() ; ++its ){
		  Hit* hit = *its ;
		  VecFromArray hPos(  hit->first->getPosition() ) ;
		  streamlog_out( ERROR ) << "  hit  : layer: "  <<   tpcLayerID( hit )   << " - " << hPos.v()  ;
		}
		

	      } else {

		
		streamlog_out( DEBUG3 ) << " +++++++++ found best hit already part of a track segment !!!!!! " 
					<< " trk : " << trk  << " #hits: " << trk->getNHits() 
					<< " cluster " << c << c->size() 
					<< std::endl ;   
		
		
		unsigned goodHits(0), allHits(0) ;
		
		double chi2Max = 10. ; // fixme parameter
		
		for( HitCluster::iterator its = c->begin(); its != c->end() ; ++its ){
		  
		  ++allHits ;
		  
		  Hit* hit = *its ;
		  VecFromArray hPos(  hit->first->getPosition() ) ;
		  
		  const gear::Vector3D* kPos = theTrack->getXingPointForLayer( tpcLayerID( hit ) ) ;
		  
		  if( kPos != 0 ) {
		    
		    //double chi2 = ch2rz( hPos.v() , *kPos )  ;
		    double chi2 = ch2rzh( hit->first , *kPos )  ;

		    streamlog_out( DEBUG3 ) << " +++++++++ chi2 : " << chi2 << hPos.v() 
					    << " +++++++++                  " << *kPos 
					    << " +++++++++  hit id " << std::hex << hit->first->id() << std::dec 
					    << std::endl ;
		    
		    if( chi2 < chi2Max ){
		      
		      ++goodHits ;
		    }
		  }
		}
		
		double goodFraction = double( goodHits ) / double(  allHits ) ;
		
		streamlog_out( DEBUG3 ) << " +++++++++ fraction of matching hits : " << goodFraction 
					<< std::endl ;   
		
		
		// ---------------------  merge the track segements -------------------------
		
		if( goodFraction > 0.5  ) { // fixme: what is reasonable here - make parameter ...
		  
		  
		  for( HitCluster::iterator its = c->begin(); its != c->end() ; ++its ){

		    delete  xpVec[  tpcLayerID( *its ) ] ; // erase crossing points for these hit
		    xpVec[  tpcLayerID( *its ) ]  = 0 ;   
		  }
		  HitCluster* clu = theTrack->getCluster< HitCluster >() ;
		  
		  // merge the cluster into the larger one and delete it - remove the hits from the hitsinlayer vector first
		  
		  Hit* h0 = c->front() ;
		  Hit* h1 = c->back() ;
		  
		  hitsInLayer[ tpcLayerID( h0 ) ].remove( h0 )  ;
		  hitsInLayer[ tpcLayerID( h1 ) ].remove( h1 )  ;
		  
		  clu->mergeClusters( c ) ;
		  
		  cluList.remove( c  ) ;
		  
		  
		  streamlog_out( DEBUG3) << " ************ found matching segment, merged all hits: delete cluster : " << c 
					 << " and track : " << trk << std::endl ;
		  
		  delete c ;
		  
		  ktracks.remove( trk ) ;
		  
		} //-------------------------------------------------------------

	      }

		
	    }  else  {  //--------------------------------------------------------------------------------
	      
	      hLL.remove(  bestHit ) ;
	      
	      HitCluster* clu = theTrack->getCluster< HitCluster >() ;
	      
	      streamlog_out( DEBUG3) << "    ************ found matching hit, add to  cluster : " << clu  << std::endl ;
	      
	      clu->addHit( bestHit ) ;
	    }


	  }
	} 
          // else {
	  //	  streamlog_out( DEBUG1 ) << "????????????????????? no best Hit found xing pnt  : chi2  " << *xpVec[ixp]  << " : " << ch2Min << std::endl ;
	  //	}

      }
    }
    
  }

   //================================================================================================================ 
  


  //std::transform( ktracks.begin(), ktracks.end(), std::back_inserter( *kaltracks ) , KalTrack2LCIO() ) ;

  std::list< KalTrack* > newKTracks ;

  //KalTestFitter<KalTest::OrderIncoming, KalTest::FitForward, KalTest::PropagateToIP > ipFitter( _kalTest ) ;

  //  KalTestFitter<KalTest::OrderOutgoing, KalTest::FitBackward, KalTest::PropagateToIP > ipFitter( _kalTest ) ;

  //FIXME: DEBUG - non ip fitter
  KalTestFitter<KalTest::OrderOutgoing, KalTest::FitBackward > ipFitter( _kalTest ) ;
  
  std::transform( cluList.begin(), cluList.end(), std::back_inserter( newKTracks ) , ipFitter  ) ;


  LCCollectionVec* kaltracks = new LCCollectionVec( LCIO::TRACK ) ;
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  kaltracks->setFlag( trkFlag.getFlag()  ) ;
  
  std::transform( newKTracks.begin(), newKTracks.end(), std::back_inserter( *kaltracks ) , KalTrack2LCIO() ) ;
  //  std::transform( cluList.begin(), cluList.end(), std::back_inserter( *kaltracks ) , converter ) ;

  evt->addCollection( kaltracks , _outColName ) ;
  
 

 //================================================================================================================ 
  //   merge track segments based on track parameters and errors ...
  //
  static const int merge_track_segments = true ;

  if( merge_track_segments ) {

    GenericHitVec<Track> trkVec ;
    GenericClusterVec<Track> trkCluVec ;
    LCCollectionVec* mergedTracks = new LCCollectionVec( LCIO::TRACK ) ;
  
    addToGenericHitVec( trkVec , kaltracks ,  AllwaysTrue()  ) ;

    //    TrackStateDistance trkMerge( 50. ) ;
    TrackCircleDistance trkMerge( 0.1 ) ; 

    cluster( trkVec.begin() , trkVec.end() , std::back_inserter( trkCluVec ), &trkMerge  , 2 ) ;


    streamlog_out( DEBUG4 ) << " ===== merged tracks - # cluster: " << trkCluVec.size()   << "  ============================== " << std::endl ;
    
    for( GenericClusterVec<Track>::iterator it= trkCluVec.begin() ; it != trkCluVec.end() ; ++it) {
      
      streamlog_out( DEBUG2 ) <<  myheader<Track>() << std::endl ;
      
      GenericCluster<Track>* trkClu = *it ;
      
      std::list<Track*> mergedTrk ;
      for( GenericCluster<Track>::iterator itC = trkClu->begin() ; itC != trkClu->end() ; ++ itC ){
	
	streamlog_out( DEBUG2 ) << lcshort(  (*itC)->first ) << std::endl ; 
	
	mergedTrk.push_back( (*itC)->first ) ; 
      }

      TrackImpl* trk = new TrackImpl ;
      Track* bestTrk = 0 ;
      double chi2Min = 99999999999999999. ;
      for( std::list<Track*>::iterator itML = mergedTrk.begin() ; itML != mergedTrk.end() ; ++ itML ){
	
	const TrackerHitVec& hV = (*itML)->getTrackerHits() ;

	for(unsigned i=0 ; i < hV.size() ; ++i){

	  trk->addHit( hV[i] ) ;
	  double chi2ndf = (*itML)->getChi2() / (*itML)->getNdf() ;

	  if( chi2ndf < chi2Min ){
	    bestTrk = (*itML) ;
	    chi2Min = chi2ndf ;
	  }
	}
      }
      if( bestTrk != 0 ){ 

	trk->setD0( bestTrk->getD0() ) ;
	trk->setOmega( bestTrk->getOmega() ) ;
	trk->setPhi( bestTrk->getPhi() ) ;
	trk->setZ0( bestTrk->getZ0() ) ;
	trk->setTanLambda( bestTrk->getTanLambda() ) ;
	trk->setCovMatrix( bestTrk->getCovMatrix()  ) ;
	// ...
	
      }
      else{
	streamlog_out( ERROR ) << "   no best track found for merged tracks ... !? " << std::endl ; 
      }
      mergedTracks->addElement( trk )  ;

    }

    // add all tracks that have not been merged :
    for( GenericHitVec<Track>::iterator it = trkVec.begin(); it != trkVec.end() ;++it){

      if( (*it)->second == 0 ){

	mergedTracks->addElement(  new TrackImpl( *dynamic_cast<TrackImpl*>( (*it)->first ) ) ) ;
      }
    }


    evt->addCollection( mergedTracks , "MergedKalTracks" ) ;
  }


  //------ register some debugging print funtctions for picking in CED :

  CEDPickingHandler::getInstance().registerFunction( LCIO::TRACKERHIT , &printTrackerHit ) ; 
  CEDPickingHandler::getInstance().registerFunction( LCIO::TRACK , &printTrackShort ) ; 
  //======================================================================================================

 
  //========== cleanup KalTracks ========
  std::for_each( ktracks.begin() , ktracks.end() , delete_ptr<KalTrack> ) ;

  //FIXME: memory leak - need for debugging....

  //  std::for_each( newKTracks.begin() , newKTracks.end() , delete_ptr<KalTrack> ) ;
  //=====================================



  //========  create collections of used and unused TPC hits ===========================================

  LCCollectionVec* usedHits = new LCCollectionVec( LCIO::TRACKERHIT ) ;   ;
  LCCollectionVec* unUsedHits = new LCCollectionVec( LCIO::TRACKERHIT ) ;   ;
  usedHits->setSubset() ;
  unUsedHits->setSubset() ;
  usedHits->reserve( h.size() ) ;
  unUsedHits->reserve( h.size() ) ;
  //  typedef GenericHitVec<TrackerHit>::iterator GHVI ;
  for( GHVI it = h.begin(); it != h.end() ;++it){
    if( (*it)->second != 0 ){
      usedHits->push_back( (*it)->first ) ;
    } else {
      unUsedHits->push_back( (*it)->first ) ;          
    }
  }
  for( GHVI it = hSmallR.begin(); it != hSmallR.end() ;++it){
    if( (*it)->second != 0 ){
      usedHits->push_back( (*it)->first ) ;
    } else {
      unUsedHits->push_back( (*it)->first ) ;          
    }
  }
  evt->addCollection( usedHits ,   "UsedTPCCluTrackerHits" ) ;
  evt->addCollection( unUsedHits , "UnUsedTPCCluTrackerHits" ) ;
  
  //========================================================================================================
  
  _nEvt ++ ;

  clock_t end = clock () ; 
  
  streamlog_out( DEBUG )  << "---  clustering time: " 
 			  <<  double( end - start ) / double(CLOCKS_PER_SEC) << std::endl  ;
  
}


/*************************************************************************************************/
void ClupatraProcessor::check( LCEvent * evt ) { 
  /*************************************************************************************************/

  std::string colName( "MergedKalTracks"  ) ;
  // std::string colName(  _outColName ) ;


  bool checkForDuplicatePadRows =  false ; //true ;
  bool checkForMCTruth =  true ;

  bool checkForSplitTracks =  true ;   // WARNING: DEBUG only - this requires the kaltracks to not be deleted in processEvent !!!!!!!!! 


  streamlog_out( MESSAGE ) <<  " check called.... " << std::endl ;

  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& pL = gearTPC.getPadLayout() ;


  //====================================================================================
  // check for duplicate padRows 
  //====================================================================================

  if( checkForDuplicatePadRows ) {

    LCCollectionVec* oddCol = new LCCollectionVec( LCIO::TRACK ) ;
    oddCol->setSubset( true ) ;
    // try iterator class ...

    LCIterator<Track> trIt( evt, colName ) ;
    while( Track* tr = trIt.next()  ){

      
      // check for duplicate layer numbers
      std::vector<int> hitsInLayer( pL.getNRows() ) ; 
      const TrackerHitVec& thv = tr->getTrackerHits() ;
      typedef TrackerHitVec::const_iterator THI ;
      for(THI it = thv.begin() ; it  != thv.end() ; ++it ) {
	TrackerHit* th = *it ;
	++ hitsInLayer.at( th->ext<HitInfo>()->layerID )   ;
      } 
      unsigned nHit = thv.size() ;
      unsigned nDouble = 0 ;
      for(unsigned i=0 ; i < hitsInLayer.size() ; ++i ) {
	if( hitsInLayer[i] > 1 ){
	  ++nDouble ;
	  streamlog_out( DEBUG4 ) << " &&&&&&&&&&&&&&&&&&&&&&&&&& duplicate hit in layer : " << i << std::endl ;
	}
      }
      if( double(nDouble) / nHit > _duplicatePadRowFraction ){
	//if( nDouble  > 0){
	streamlog_out( DEBUG4 ) << " oddTrackCluster found with "<< 100. * double(nDouble) / nHit 
				<< "% of double hits " << std::endl ;
	oddCol->addElement( tr ) ;
      }
    }
    evt->addCollection( oddCol , "OddCluTracks" ) ;
  }
  //====================================================================================
  // check Monte Carlo Truth via SimTrackerHits 
  //====================================================================================

  if( checkForMCTruth ) {
 

    LCCollectionVec* oddCol = new LCCollectionVec( LCIO::TRACK ) ;
    oddCol->setSubset( true ) ;

    LCCollectionVec* splitCol = new LCCollectionVec( LCIO::TRACK ) ;
    splitCol->setSubset( true ) ;
    
    typedef std::map<Track* , unsigned > TRKMAP ; 
    
    typedef std::map< MCParticle* , TRKMAP > MCPTRKMAP ;
    MCPTRKMAP mcpTrkMap ;
    
    typedef std::map< MCParticle* , unsigned > MCPMAP ;
    MCPMAP hitMap ;
    
    
    // if( streamlog_level( DEBUG4) )
    //   LCTOOLS::printTracks( evt->getCollection("KalTestTracks") ) ;


    LCIterator<Track> trIt( evt, colName  ) ;  
    //    "KalTestTracks" ) ;
    //    LCIterator<Track> trIt( evt, _outColName ) ;
    //    LCIterator<Track> trIt( evt, "TPCTracks" ) ;

    while( Track* tr = trIt.next()  ){
      
      MCPMAP mcpMap ;

      const TrackerHitVec& thv = tr->getTrackerHits() ;
      typedef TrackerHitVec::const_iterator THI ;

      // get relation between mcparticles and tracks
      for(THI it = thv.begin() ; it  != thv.end() ; ++it ) {

	TrackerHit* th = *it ;
	// FIXME:
	// we know that the digitizer puts the sim hit into the raw hit pointer
	// but of course the proper way is to go through the LCRelation ...
	SimTrackerHit* sh = (SimTrackerHit*) th->getRawHits()[0] ;
	MCParticle* mcp = sh->getMCParticle() ;

	
	hitMap[ mcp ] ++ ;   // count all hits from this mcp
	
	mcpMap[ mcp ]++ ;    // count hits from this mcp for this track
	
	mcpTrkMap[ mcp ][ tr ]++ ;  // map between mcp, tracks and hits
	
      } 

      // check for tracks with hits from several mcparticles
      unsigned nHit = thv.size() ;
      unsigned maxHit = 0 ; 
      for( MCPMAP::iterator it= mcpMap.begin() ;
	   it != mcpMap.end() ; ++it ){
	if( it->second  > maxHit ){
	  maxHit = it->second ;
	}
      }

      if( double(maxHit) / nHit < 0.99 ){ // What is acceptable here ???
	//if( nDouble  > 0){
	streamlog_out( MESSAGE ) << " oddTrackCluster found with only "
				 << 100.*double(maxHit)/nHit 
				 << "% of hits  form one MCParticle " << std::endl ;
	oddCol->addElement( tr ) ;
      }
    }
    evt->addCollection( oddCol , "OddMCPTracks" ) ;
    
    
    if( checkForSplitTracks ) {
      
      streamlog_out( DEBUG ) << " checking for split tracks - mcptrkmap size : " <<  mcpTrkMap.size() << std::endl ;
      
      // check for split tracks 
      for( MCPTRKMAP::iterator it0 = mcpTrkMap.begin() ; it0 != mcpTrkMap.end() ; ++it0){
	
	streamlog_out( DEBUG ) << " checking for split tracks - map size : " <<  it0->second.size() << std::endl ;
	
	
	if( it0->second.size() > 1 ) {
	  
	  
	  typedef std::list< EVENT::Track* > TL ;
	  TL trkList ;
	  
	  for( TRKMAP::iterator it1 = it0->second.begin() ; it1 != it0->second.end() ; ++it1){
	    
	    double totalHits = hitMap[ it0->first ]  ; // total hits for this track 
	    
	    double thisMCPHits = it1->second ;     //  hits from this mcp
	    
	    double ratio =  thisMCPHits / totalHits  ;
	    
	    streamlog_out( DEBUG ) << " checking for split tracks - ratio : " 
				   << thisMCPHits << " / " << totalHits << " = " << ratio << std::endl ;
	    
	    if( ratio > 0.03 && ratio < 0.95 ){
	      // split track
	      
	      splitCol->addElement( it1->first ) ; 
	      
	      trkList.push_back( it1->first ) ;
	    } 
	  }
	  // chi2 between split track segments :
	  // for( TRKMAP::iterator ist0 = it0->second.begin() ; ist0 != it0->second.end() ; ++ist0){
	    
	  //   KalTrack* sptrk0 = ist0->first->ext<KalTrackLink>() ; 
	    
	  //   TRKMAP::iterator ist0_pp = ist0 ;
	  //   ++ist0_pp ;

	  //   for( TRKMAP::iterator ist1 = ist0_pp ; ist1 != it0->second.end() ; ++ist1){
	  
	  //     KalTrack* sptrk1 = ist1->first->ext<KalTrackLink>() ; 
	      
	  //     double chi2 =  KalTrack::chi2( *sptrk0 ,  *sptrk1 ) ;
	      
	  //     streamlog_out( DEBUG4 ) << " *********************  chi2 between split tracks : "  << chi2 << std::endl 
	  // 			      << myheader< Track >() << std::endl 
	  // 			      << lcshort( ist0->first )  << std::endl 
	  // 			      << lcshort( ist1->first )	 << std::endl ; 
	      
	  //   }
	  // }



	  streamlog_out( DEBUG2 ) << " ------------------------------------------------------ " << std::endl ;
	  
	  for( TL::iterator it0 = trkList.begin() ; it0 != trkList.end() ; ++it0 ){
	    
	    
	    //	    KalTrack* trk0 = (*it0)->ext<KalTrackLink>() ; 
	    
	    HelixClass hel ;
	    hel.Initialize_Canonical( (*it0)->getPhi(),
				      (*it0)->getD0(),
				      (*it0)->getZ0(),
				      (*it0)->getOmega(),
				      (*it0)->getTanLambda(),
				      3.50 ) ;
	    
	    streamlog_out( DEBUG1 ) << hel.getXC() << "\t"
				    << hel.getYC() << "\t"
				    << hel.getRadius() << "\t" 
				    << hel.getTanLambda() << std::endl ; 
	    
	    
	    // streamlog_out( DEBUG1 ) << (*it0)->getPhi() << "\t"
	    // 			  << (*it0)->getD0()  << "\t"
	    // 			  << (*it0)->getOmega()  << "\t"
	    // 			  << (*it0)->getZ0()  << "\t"
	    // 			  << (*it0)->getTanLambda()  << "\t"
	    // 			  << std::endl ;
	    
	    //	    streamlog_out( DEBUG1 ) << " trk0 : " << *trk0 << std::endl ;
	    
	    // TL::iterator its = it0 ;
	    // ++its ;
	    
	    // for( TL::iterator it1 =  its ; it1 != trkList.end() ; ++it1 ){
	      
	    //   KalTrack* trk1 = (*it1)->ext<KalTrackLink>() ; 
	      
	    //   streamlog_out( DEBUG1 ) << "    - trk0 : " << *trk0 << std::endl ;
	    //   streamlog_out( DEBUG1 ) << "    - trk1 : " << *trk1 << std::endl ;
	      
	    //   double chi2 =  KalTrack::chi2( *trk0 ,  *trk1 ) ;
	      
	    //   streamlog_out( DEBUG1 ) << " +++++++++++++++++  chi2 between split tracks : " 
	    // 			      << trk0 << " - " << trk1 << " : " << chi2 << std::endl ; 
	      
	      
	    // }
	  }
	  
	}
      }
      evt->addCollection( splitCol , "SplitTracks" ) ;
    }

  }
  //====================================================================================

}


void ClupatraProcessor::end(){ 
  
  //   std::cout << "ClupatraProcessor::end()  " << name() 
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;
  
  delete _kalTest ;
}


//====================================================================================================
