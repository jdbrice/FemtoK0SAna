#ifndef FEMTOK0S_ANA_H
#define FEMTOK0S_ANA_H

#include "TreeAnalyzer.h"

// FemtoDstFormat
#include "FemtoDstFormat/BranchReader.h"
#include "FemtoDstFormat/TClonesArrayReader.h"
#include "FemtoDstFormat/FemtoEvent.h"
#include "FemtoDstFormat/FemtoTrack.h"
#include "FemtoDstFormat/FemtoTrackHelix.h"
#include "FemtoDstFormat/FemtoMtdPidTraits.h"
#include "FemtoDstFormat/FemtoBTofPidTraits.h"
#include "FemtoDstFormat/FemtoTrackProxy.h"

#include "vendor/loguru.h"

#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/StDcaGeometry.h"
#include "StarClassLibrary/PhysicalConstants.h"
#include "StarClassLibrary/StThreeVectorD.hh"

#include <vector>
#include <memory>

const float BFIELD = -4.984500;
const float MASS_PI = 0.13957; // in GeV

struct FullTrack {
	FemtoTrack *t = nullptr;
	FemtoTrackHelix *h = nullptr;
	FemtoMtdPidTraits * m = nullptr;
	FemtoBTofPidTraits *b = nullptr;
	StPhysicalHelixD helix;
	void calcHelix(){
		StThreeVectorF gMom, gOrigin;
		gMom.set( this->h->mPar[0], this->h->mPar[1], this->h->mPar[2] );
		gOrigin.set( this->h->mPar[3], this->h->mPar[4], this->h->mPar[5] );
		StPhysicalHelixD gHelix(gMom, gOrigin, BFIELD * kilogauss, static_cast<float>(this->t->charge()) );
		this->helix = gHelix;
	}
};

class FemtoK0SAna : public TreeAnalyzer {
protected:
	FemtoEvent *_event;
	FemtoTrackProxy _proxy;
	FemtoTrackProxy _proxy2;

	BranchReader<FemtoEvent> _rEvent;
	TClonesArrayReader<FemtoTrack> _rTracks;
	TClonesArrayReader<FemtoTrackHelix> _rHelices;
	TClonesArrayReader<FemtoMtdPidTraits> _rMtdPid;
	TClonesArrayReader<FemtoBTofPidTraits> _rBTofPid;
	

	vector<shared_ptr<FullTrack> > pip, pim;

public:
	FemtoK0SAna() {}
	~FemtoK0SAna() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		_rEvent.setup( chain, "Event" );
		_rTracks.setup( chain, "Tracks" );
		_rHelices.setup( chain, "Helices" );
		_rMtdPid.setup( chain, "MtdPidTraits" );
		_rBTofPid.setup( chain, "BTofPidTraits" );

	}

protected:

	void analyzePairs( FemtoEvent * _event ){

		StThreeVectorF pVtx;
		pVtx.set( _event->mPrimaryVertex_mX1, _event->mPrimaryVertex_mX2, _event->mPrimaryVertex_mX3 );
		TLorentzVector lvp, lvn, lv;
		for ( auto pos : pip ){
			for ( auto neg : pim ){

				size_t nMTD = 0;
				if ( nullptr != pos->m )
					nMTD++;
				if ( nullptr != neg->m )
					nMTD++;
				if ( nMTD < 1 ) continue;
				  
				pair<double, double> pathLengths = pos->helix.pathLengths( neg->helix );
				StThreeVectorD pNegPosAtDca = pos->helix.at(pathLengths.first);
				StThreeVectorD pPosPosAtDca = neg->helix.at(pathLengths.second);
				
				StThreeVectorD dcaVector = pNegPosAtDca-pPosPosAtDca;
				StThreeVectorD secVtxPos = pPosPosAtDca + ( dcaVector * 0.5 );
				StThreeVectorD decLenVec = secVtxPos - pVtx;

				book->fill( "decayLength", decLenVec.mag() );
				book->fill( "dcaVec", dcaVector.mag() );
				
				if (decLenVec.mag() < 2.7) continue; // cut on path length
				if (dcaVector.mag() > 1.5) continue; // cut on track mutual dca

				StThreeVectorD pNegMomAtDca = pos->helix.momentumAt(pathLengths.first, BFIELD * kilogauss );
				StThreeVectorD pPosMomAtDca = neg->helix.momentumAt(pathLengths.second, BFIELD * kilogauss );
				StThreeVectorD K0sMomAtDCA = pNegMomAtDca + pPosMomAtDca;
				lvp.SetPtEtaPhiM( pNegMomAtDca.perp(), pNegMomAtDca.pseudoRapidity(), pNegMomAtDca.phi(), MASS_PI );
				lvn.SetPtEtaPhiM( pPosMomAtDca.perp(), pPosMomAtDca.pseudoRapidity(), pPosMomAtDca.phi(), MASS_PI );
				lv = lvp+lvn;

				double pointingAngle = K0sMomAtDCA.angle(decLenVec);
				book->fill( "pointingAngle", pointingAngle );

				double K0pT = lv.Pt();
				if (TMath::Abs(pointingAngle) > 0.106+0.056-0.1123*K0pT+0.025*K0pT*K0pT) continue;
				

				

				if ( 0 == nMTD  )
					book->fill( "pt_mass", lv.M(), lv.Pt() );
				else if ( 1 == nMTD  )
					book->fill( "pt_mass_smtd", lv.M(), lv.Pt() );
				else if ( 2 == nMTD  )
					book->fill( "pt_mass_bmtd", lv.M(), lv.Pt() );


			} // loop neg
		}// loop pos

	}


	virtual void analyzeEvent(){
		_event = _rEvent.get();

		StThreeVectorF pVtx;
		pVtx.set( _event->mPrimaryVertex_mX1, _event->mPrimaryVertex_mX2, _event->mPrimaryVertex_mX3 );
		size_t nTracks = _rTracks.N();
		LOG_F( 9, "mult=%lu", nTracks );

		pip.clear();
		pim.clear();


		for (size_t i = 0; i < nTracks; i++ ){

			auto ft = shared_ptr<FullTrack>( new FullTrack() );

			ft->t = _rTracks.get( i );
			ft->h = _rHelices.get( ft->t->mHelixIndex );
			ft->calcHelix();

			if ( ft->t->mBTofPidTraitsIndex >= 0 )
				ft->b = _rBTofPid.get( ft->t->mBTofPidTraitsIndex );
			if ( ft->t->mMtdPidTraitsIndex >= 0 )
				ft->m = _rMtdPid.get( ft->t->mMtdPidTraitsIndex );

			if ( ft->t->charge() > 0 )
				pip.push_back( ft );
			else 
				pim.push_back( ft );
			
			continue;

			LOG_F( INFO, "dca = %f", ft->h->mDCA  );

			StDcaGeometry dcaGeom;
			Float_t errMat[12] = {0};
			dcaGeom.set(ft->h->mPar, errMat);
			// StPhysicalHelixD gHelix = dcaGeom.helix();
			StThreeVectorF gMom, gOrigin;
			gMom.set( ft->h->mPar[0], ft->h->mPar[1], ft->h->mPar[2] );
			gOrigin.set( ft->h->mPar[3], ft->h->mPar[4], ft->h->mPar[5] );
			StPhysicalHelixD gHelix(gMom, gOrigin, BFIELD * kilogauss, static_cast<float>(ft->t->charge()) );

			
			float gP = gHelix.momentumAt(gHelix.pathLength(pVtx), BFIELD * kilogauss ).mag();
			LOG_F( INFO, "gP = %f, pP = %f", gP, ft->t->momentum().Mag() );

			gHelix.moveOrigin( gHelix.pathLength(pVtx) );

			StThreeVectorF origin = gHelix.origin();
			StThreeVectorF HMom   = gHelix.momentum( BFIELD * kilogauss); // taken from PicoEvent

			StThreeVectorF diff   = origin - pVtx;
			Float_t gdca = diff.mag();
			LOG_F( INFO, "DCA calc = %f", gdca );


		} // loop on Tracks

		if ( pip.size() < 1 || pim.size() < 1 )
			return;
		
		analyzePairs( _event );


	} // analyze Event



};

#endif
