#ifndef FEMTO_TRACK_H
#define FEMTO_TRACK_H

#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "FemtoDstBranch.h"

class FemtoTrack : public FemtoDstBranch
{
public:
	virtual const char* classname() const { return "FemtoTrack"; }
	FemtoTrack(){
		reset();
	}
	~FemtoTrack(){
	}

	void reset(){
		mPt                 = 0.0;
		mEta                = 0.0;
		mPhi                = 0.0;
		mId                 = 0.0;
		mDedx               = 0.0;
		mNHitsFit           = 0.0;
		mNHitsMax           = 0.0;
		mNHitsDedx          = 0.0;
		mNSigmaPion         = 0.0;
		mNSigmaKaon         = 0.0;
		mNSigmaProton       = 0.0;
		mNSigmaElectron     = 0.0;
		mDCA                = 0.0;
		
		// pidTraits
		mBTofPidTraitsIndex = -1;
		mMtdPidTraitsIndex  = -1;
		mEmcPidTraitsIndex  = -1;
		mHelixIndex         = -1;
		mMcIndex            = -1;
	}

	void copy( FemtoTrack * that ) {
		this->mPt                 = that->mPt;
		this->mEta                = that->mEta;
		this->mPhi                = that->mPhi;
		this->mId                 = that->mId;
		this->mDedx               = that->mDedx;
		this->mNHitsFit           = that->mNHitsFit;
		this->mNHitsMax           = that->mNHitsMax;
		this->mNHitsDedx          = that->mNHitsDedx;
		this->mNSigmaPion         = that->mNSigmaPion;
		this->mNSigmaKaon         = that->mNSigmaKaon;
		this->mNSigmaProton       = that->mNSigmaProton;
		this->mNSigmaElectron     = that->mNSigmaElectron;
		this->mDCA                = that->mDCA;

		this->mBTofPidTraitsIndex = that->mBTofPidTraitsIndex;
		this->mMtdPidTraitsIndex  = that->mMtdPidTraitsIndex;
		this->mEmcPidTraitsIndex  = that->mEmcPidTraitsIndex;
		this->mHelixIndex         = that->mHelixIndex;
		this->mMcIndex            = that->mMcIndex;
	}

	int charge() { return mNHitsFit > 0 ? 1 : -1; }
	float nSigmaPion() { return mNSigmaPion / 100.0; }
	void nSigmaPion( float tmp ){ this->mNSigmaPion = tmp * 100; }
	float nSigmaKaon() { return mNSigmaKaon / 100.0; }
	void nSigmaKaon( float tmp ){ this->mNSigmaKaon = tmp * 100; }
	float nSigmaProton() { return mNSigmaProton / 100.0; }
	void nSigmaProton( float tmp ){ this->mNSigmaProton = tmp * 100; }
	float nSigmaElectron() { return mNSigmaElectron / 100.0; }
	void nSigmaElectron( float tmp ){ this->mNSigmaElectron = tmp * 100; }

	float gDCA() { return mDCA / 1000.0; }
	void gDCA( float dca ) 
	{ 
		if (dca < 0 ) dca = 50;		// dont default to zero, something is wrong 
		this->mDCA = dca * 1000; 
	}

	TVector3 momentum() { TVector3 v; v.SetPtEtaPhi( mPt, mEta, mPhi ); return v; } 
	void momentum( float pt, float eta, float phi ) { 
		this->mPt  = pt;
		this->mEta = eta;
		this->mPhi = phi;
	}
	TLorentzVector lv( double m ) { 
		TLorentzVector lv;
		lv.SetPtEtaPhiM( mPt, mEta, mPhi, m );
		return lv;
	}

	float dEdx() { return mDedx / 1000.0; }
	void dEdx( float tmp ) { this->mDedx = tmp * 1000; }
	
	template <typename T>
	T v( std::string name ){
		if ( name == "mPt" ) return mPt;
		if ( name == "mEta" ) return mEta;
		if ( name == "mPhi" ) return mPhi;
		if ( name == "mId" ) return mId;
		if ( name == "mDedx" ) return mDedx;
		if ( name == "mNHitsFit" ) return mNHitsFit;
		if ( name == "mNHitsMax" ) return mNHitsMax;
		if ( name == "mNHitsDedx" ) return mNHitsDedx;
		if ( name == "mNSigmaPion" ) return nSigmaPion();
		if ( name == "mNSigmaKaon" ) return nSigmaKaon();
		if ( name == "mNSigmaProton" ) return nSigmaProton();
		if ( name == "mNSigmaElectron" ) return nSigmaElectron();
		if ( name == "mDCA" ) return gDCA();

		if ( name == "mBTofPidTraitsIndex" ) return mBTofPidTraitsIndex;
		if ( name == "mMtdPidTraitsIndex" ) return mMtdPidTraitsIndex;
		if ( name == "mEmcPidTraitsIndex" ) return mEmcPidTraitsIndex;
		if ( name == "mHelixIndex" ) return mHelixIndex;
		if ( name == "mMcIndex" ) return mMcIndex;
		return 0;
	}
	

	Float_t 	mPt;				// primary track px
	Float_t 	mEta;				// primary track py
	Float_t 	mPhi;				// primary track pz
	UShort_t 	mId;				// track Id
	UShort_t 	mDedx;				// dEdx*1000
	Char_t 		mNHitsFit;			// q*nHitsFit - TPC
	Char_t 		mNHitsMax;			// nHitsMax - TPC
	UChar_t 	mNHitsDedx;			// nHitsDedx - TPC
	Short_t 	mNSigmaPion;		// nsigmaPi * 100
	Short_t 	mNSigmaKaon;		// nsigmaK * 100
	Short_t 	mNSigmaProton;		// nsigmaP * 100
	Short_t 	mNSigmaElectron;	// nsigmaE * 100
	UShort_t 	mDCA;				// global track gDCA * 1000
	
	// pidTraits
	Short_t 	mBTofPidTraitsIndex; // index of the BTOF pidTratis in the event
	Short_t 	mMtdPidTraitsIndex;  // index of the MTD pidTratis in the event
	Short_t 	mEmcPidTraitsIndex;  // index of the EMC pidTratis in the event
	Short_t 	mHelixIndex;		 // index of associated Track Helix
	Short_t 	mMcIndex;			 // index of associated Track Helix
	
	ClassDef( FemtoTrack, 3 )
};

#endif
