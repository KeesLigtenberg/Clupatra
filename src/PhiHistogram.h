/*
 * PhiHistogram.h
 *
 *  Created on: Jan 16, 2017
 *      Author: cligtenb
 */

#ifndef SRC_PHIHISTOGRAM_H_
#define SRC_PHIHISTOGRAM_H_

#include <memory>

#include "TH1I.h"
#include "TMath.h"

namespace clupatra_new {

	class PhiHistogram {
	public:
		PhiHistogram(int nbins);

		void Fill(Hit* h);

		//maximum functions of TH1 are buffered, so there is no need to worry about the number of calls to maximum
		double getMaximum(int minHits=1);
		int getMaximumSize();
		//return the cluster in the maximum bin and hits 'binrange' to the left and to the right of the maximum, e.g. 1 returns 3 bins!
		Clusterer::cluster_type getMaximimumCluster(int binrange) { return getMaximimumCluster(1,binrange); } ;
		Clusterer::cluster_type getMaximimumCluster( int minHits, int binrange);

		static 	TrackStateImpl MakeTrackStateFromPhi( double phi);

	private:
		//works for CLHEP, DDSurface and GEAR vectors
		template<class T>
		void fillHistogram(const T& pos);
		TH1I _hist;

		void fillCluster(Hit* h);
		std::vector< std::unique_ptr<Clusterer::cluster_type> > _cluVector;

	};

	PhiHistogram::PhiHistogram(int nbins) :
			_hist("phiHist", "histogram for PhiHistogram", nbins, -TMath::Pi() , TMath::Pi() ),
			_cluVector(nbins+1)
			{

	}

	void PhiHistogram::Fill(Hit* h) {
		fillHistogram(h->first->pos);
		fillCluster(h);
	}

	template<class T>
	void PhiHistogram::fillHistogram(const T& pos) {
		_hist.Fill(pos.phi());
	}

	double PhiHistogram::getMaximum(int minHits) {
		int maxbin=_hist.GetMaximumBin();
		if(_hist.GetBinContent(maxbin)>minHits) {
			return _hist.GetBinCenter(maxbin);
		} else {
			return 0;
		}
	}

	Clusterer::cluster_type PhiHistogram::getMaximimumCluster(int minHits, int binrange) {
		int maxbin=_hist.GetMaximumBin();
		if(_hist.GetBinContent(maxbin)>minHits) {
			Clusterer::cluster_type cluster=*_cluVector[maxbin];
		    _cluVector[maxbin]->clear();
		    _hist.SetBinContent(maxbin, 0);

		    //loop over neighbouring bins, but not the bin itself
			for(int i=-binrange; i<=binrange; i++) {
			    if(!i) continue;

			    //put bins in the axis range
			    int bin=i+maxbin;
			    bin %= _hist.GetNbinsX();
			    if(bin<0) bin+=_hist.GetNbinsX();

			    //check if cluster is not empty
			   if(!_hist.GetBinContent(bin)) continue;

			    cluster.mergeClusters( _cluVector[bin].get() );
			    _cluVector[bin]->clear();
			    _hist.SetBinContent(bin, 0);
			}
		    return cluster;
		} else {
			//return empty cluster
			return Clusterer::cluster_type();
		}
	}

inline TrackStateImpl PhiHistogram::MakeTrackStateFromPhi(double phi) {
		const float originReference[3] = {0,0,0}; //reference point is origin
		const float d0error=3000, //set a large d0 error
				phierror=0.1,
				unknown=1E16; //extremely large error for unknowns
		return TrackStateImpl(EVENT::TrackState::AtOther, /*location*/
				0., /*d0*/
				phi, /*phi0*/
				1E-6, /*omega*/
				0., /*z0*/
				-.5, /*tanLambda*/
				{ d0error,
				  unknown, phierror,
				  unknown, unknown, unknown,
				  unknown, unknown, unknown, unknown,
				  unknown, unknown, unknown, unknown, unknown,}, /*covMatrix*/
				originReference ); /*reference*/
}

int PhiHistogram::getMaximumSize() {
	int maxbin=_hist.GetMaximumBin();
	return _hist.GetBinContent(maxbin);
}

	void PhiHistogram::fillCluster(Hit* h) {
		int bin=_hist.GetXaxis()->FindBin(h->first->pos.phi() );
		if(bin<=0 or bin>=(_hist.GetNbinsX()+1) ) return; //overflow or underflow
		if(_cluVector[bin]) { //already a cluster at this bin
			streamlog_out(DEBUG0)<<"there are "<<_hist.GetBinContent(bin)<<" hits in the histogram bin"<<std::endl;
			streamlog_out(DEBUG0)<<"cluster in this bin has "<<_hist.GetBinContent(bin)<<" elements"<<std::endl;
			_cluVector[bin]->addElement(h);
		} else { //no cluster yet
			_cluVector[bin]=std::unique_ptr<Clusterer::cluster_type>(new Clusterer::cluster_type(h) );
		}
	}

}

#endif /* SRC_PHIHISTOGRAM_H_ */
