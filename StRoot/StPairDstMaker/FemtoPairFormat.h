#ifndef FEMTO_PAIR_H
#define FEMTO_PAIR_H

#include "TObject.h"

class FemtoDstBranch : public TObject
{
public:
    // virtual const char* branchName() const = 0;
    virtual const char* classname() const = 0;

    virtual ~FemtoDstBranch() {};

    ClassDef( FemtoDstBranch, 1 )
};

class FemtoPair
{
public:
    virtual const char* classname() const { return "FemtoPair"; }
    FemtoPair() {

    }
    ~FemtoPair() {}

    void reset(){
        mRunID                = -999;
        mVertexZ              = -999;
        mDeltaVertexZ         = -999;
        mGRefMult             = 0;
        mZDCEast              = 0;
        mZDCWest              = 0;
        mChargeSum            = -3;
        mPt                   = 0.0;
        mEta                  = 0.0;
        mPhi                  = 0.0;
        mMass                 = 0.0;
        mRapidity             = 0.0;
        d1_mYLocal            = 0.0;
        d1_mZLocal            = 0.0;
        d1_mTof               = 0.0;
        d1_mMatchFlag         = 0.0;
        d1_mLength            = 0.0;
        d1_mPt                = 0.0;
        d1_mEta               = 0.0;
        d1_mPhi               = 0.0;
        d1_mgPt               = 0.0;
        d1_mgEta              = 0.0;
        d1_mgPhi              = 0.0;
        d1_mId                = 0.0;
        d1_mNHitsFit          = 0.0;
        d1_mNHitsMax          = 0.0;
        d1_mNHitsDedx         = 0.0;
        d1_mNSigmaPion        = 0.0;
        d1_mNSigmaKaon        = 0.0;
        d1_mNSigmaProton      = 0.0;
        d1_mNSigmaElectron    = 0.0;
        d1_mDCA               = 0.0;
        d1_mPid               = 0.0;
        d2_mYLocal            = 0.0;
        d2_mZLocal            = 0.0;
        d2_mTof               = 0.0;
        d2_mMatchFlag         = 0.0;
        d2_mLength            = 0.0;
        d2_mPt                = 0.0;
        d2_mEta               = 0.0;
        d2_mPhi               = 0.0;
        d2_mId                = 0.0;
        d2_mNHitsFit          = 0.0;
        d2_mNHitsMax          = 0.0;
        d2_mNHitsDedx         = 0.0;
        d2_mNSigmaPion        = 0.0;
        d2_mNSigmaKaon        = 0.0;
        d2_mNSigmaProton      = 0.0;
        d2_mNSigmaElectron    = 0.0;
        d2_mDCA               = 0.0;
        d2_mPid               = 0.0;
    }
    void copy( FemtoPair * that ){
        this->mRunID                = that->mRunID;
        this->mVertexZ              = that->mVertexZ;
        this->mZDCEast              = that->mZDCEast;
        this->mZDCWest              = that->mZDCWest;
        this->mDeltaVertexZ         = that->mDeltaVertexZ;
        this->mGRefMult             = that->mGRefMult;
        this->mChargeSum            = that->mChargeSum;
        this->mPt                   = that->mPt;
        this->mEta                  = that->mEta;
        this->mPhi                  = that->mPhi;
        this->mMass                 = that->mMass;
        this->mRapidity             = that->mRapidity;
        this->d1_mYLocal            = that->d1_mYLocal;
        this->d1_mZLocal            = that->d1_mZLocal;
        this->d1_mTof               = that->d1_mTof;
        this->d1_mMatchFlag         = that->d1_mMatchFlag;
        this->d1_mLength            = that->d1_mLength;
        this->d1_mPt                = that->d1_mPt;
        this->d1_mEta               = that->d1_mEta;
        this->d1_mPhi               = that->d1_mPhi;
        this->d1_mgPt               = that->d1_mgPt;
        this->d1_mgEta              = that->d1_mgEta;
        this->d1_mgPhi              = that->d1_mgPhi;
        this->d1_mId                = that->d1_mId;
        this->d1_mNHitsFit          = that->d1_mNHitsFit;
        this->d1_mNHitsMax          = that->d1_mNHitsMax;
        this->d1_mNHitsDedx         = that->d1_mNHitsDedx;
        this->d1_mNSigmaPion        = that->d1_mNSigmaPion;
        this->d1_mNSigmaKaon        = that->d1_mNSigmaKaon;
        this->d1_mNSigmaProton      = that->d1_mNSigmaProton;
        this->d1_mNSigmaElectron    = that->d1_mNSigmaElectron;
        this->d1_mDCA               = that->d1_mDCA;
        this->d1_mPid               = that->d1_mPid;
        this->d2_mYLocal            = that->d2_mYLocal;
        this->d2_mZLocal            = that->d2_mZLocal;
        this->d2_mTof               = that->d2_mTof;
        this->d2_mMatchFlag         = that->d2_mMatchFlag;
        this->d2_mLength            = that->d2_mLength;
        this->d2_mPt                = that->d2_mPt;
        this->d2_mEta               = that->d2_mEta;
        this->d2_mPhi               = that->d2_mPhi;
        this->d2_mgPt               = that->d2_mgPt;
        this->d2_mgEta              = that->d2_mgEta;
        this->d2_mgPhi              = that->d2_mgPhi;
        this->d2_mId                = that->d2_mId;
        this->d2_mNHitsFit          = that->d2_mNHitsFit;
        this->d2_mNHitsMax          = that->d2_mNHitsMax;
        this->d2_mNHitsDedx         = that->d2_mNHitsDedx;
        this->d2_mNSigmaPion        = that->d2_mNSigmaPion;
        this->d2_mNSigmaKaon        = that->d2_mNSigmaKaon;
        this->d2_mNSigmaProton      = that->d2_mNSigmaProton;
        this->d2_mNSigmaElectron    = that->d2_mNSigmaElectron;
        this->d2_mDCA               = that->d2_mDCA;
        this->d2_mPid               = that->d2_mPid;
    }

    // Event
    Int_t       mRunID;                 // Run number
    Float_t     mVertexZ;               // Event VertexZ
    Float_t     mDeltaVertexZ;          // delta vz
    UShort_t    mGRefMult;              // global RefMult
    UShort_t    mZDCEast;               // ZDC East
    UShort_t    mZDCWest;               // ZDC West


    // Pair
    Char_t    mChargeSum;
    Float_t   mPt;
    Float_t   mEta;
    Float_t   mPhi;
    Float_t   mMass;
    Float_t   mRapidity;

    // BTOF Stuff
    Float_t   d1_mYLocal;
    Float_t   d1_mZLocal;
    Float_t   d1_mTof;
    Char_t    d1_mMatchFlag;
    Float_t   d1_mLength;

    // Track Stuff
    Float_t     d1_mPt;                 // primary track px
    Float_t     d1_mEta;                // primary track py
    Float_t     d1_mPhi;                // primary track pz
    UShort_t    d1_mId;                 // track Id
    Char_t      d1_mNHitsFit;           // q*nHitsFit - TPC
    Char_t      d1_mNHitsMax;           // nHitsMax - TPC
    UChar_t     d1_mNHitsDedx;          // nHitsDedx - TPC
    Float_t     d1_mNSigmaPion;         // nsigmaPi * 100
    Float_t     d1_mNSigmaKaon;
    Float_t     d1_mNSigmaProton;
    Float_t     d1_mNSigmaElectron;     // nsigmaElectron
    Float_t     d1_mDCA;
    Float_t     d1_mPid;

    // BTOF Stuff
    Float_t   d2_mYLocal;
    Float_t   d2_mZLocal;
    Float_t   d2_mTof;
    Char_t    d2_mMatchFlag;
    Float_t   d2_mLength;

    // Track Stuff
    Float_t     d2_mPt;                 // primary track px
    Float_t     d2_mEta;                // primary track py
    Float_t     d2_mPhi;                // primary track pz
    UShort_t    d2_mId;                 // track Id
    Char_t      d2_mNHitsFit;           // q*nHitsFit - TPC
    Char_t      d2_mNHitsMax;           // nHitsMax - TPC
    UChar_t     d2_mNHitsDedx;          // nHitsDedx - TPC
    Float_t     d2_mNSigmaPion;         // nsigmaPi * 100
    Float_t     d2_mNSigmaKaon;
    Float_t     d2_mNSigmaProton;
    Float_t     d2_mNSigmaElectron;     // nsigmaElectron
    Float_t     d2_mDCA;
    Float_t     d2_mPid;

    Float_t     d1_mgPt;                    // global track px
    Float_t     d1_mgEta;               // global track py
    Float_t     d1_mgPhi;               // global track pz

    Float_t     d2_mgPt;                    // global track px
    Float_t     d2_mgEta;               // global track py
    Float_t     d2_mgPhi;               // global track pz

    ClassDef( FemtoPair, 3 )
};


#endif
