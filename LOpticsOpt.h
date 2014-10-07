///////////////////////////////////////////////////////////////////////////////
//
// LOpticsOpt
//
// HRS optics matrix optimization class
// Based on THaVDC
//
// Units used:
//   For X, Y, and Z coordinates of track    -  meters
//   For Theta and Phi angles of track       -  tan(angle)
//   For Momentums, Masses                   -  GeV, GeV/c^2
//
// Author: Jin Huang <jinhuang@jlab.org>
//         Chao Gu <cg2ja@jlab.org>
//
// Modification:
//   Jun 25, 2010 Updated for APEX optics calibration
//   Aug 01, 2013 Updated for G2P optics calibration (Chao Gu)
//
///////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_LOpticsOpt
#define ROOT_LOpticsOpt

#include <vector>

#include "TMath.h"
#include "TString.h"
#include "TRotation.h"

#include "THaTrackingDetector.h"

///////////////////////////////////////////////////////////////////////////////
// Debug Definitions
// place this section below any other head files
///////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG_LEVEL
#undef DEBUG_LEVEL
#endif
//     DEBUG_LEVEL;
//     =0      or not define: no debug, full speed
//     >=1     enable debug extra warning (suggested setting)
//     >=2     above + enable debug assert
//     >=3     above + enable debug extra info
//     >=4     above + massive info (in a for or while)
//     >=5     Decode dump
#define DEBUG_LEVEL 3
#include "DebugDef.h"

///////////////////////////////////////////////////////////////////////////////
class TClonesArray;
class TList;
class TTree;
class TVector3;

class THaEvData;
class THaMatrixElement;

class LOpticsOpt : public THaTrackingDetector {
public:
    LOpticsOpt(const char* name = "Optimizer", const char* description = "Optimizer for HRS Optics", THaApparatus* a = NULL);
    virtual ~LOpticsOpt();

    ///////////////////////////////////////////////////////////////////////////
    // Database input/output
    ///////////////////////////////////////////////////////////////////////////
    TString OldComments;
    Int_t LoadDataBase(TString DataBaseName); // Database file -> Memory
    Int_t SaveDataBase(TString DataBaseName); // Memory -> Database file

    virtual void Print(const Option_t* opt) const;

    Int_t Matrix2Array(Double_t Array[], Bool_t FreeParaFlag[] = NULL, Int_t UseFPOff = 0); // fCurrentMatrixElems -> Array
    Int_t Matrix2Array(Double_t Array[], const std::vector<THaMatrixElement>& Matrix, Bool_t FreeParaFlag[] = NULL, Int_t UseFPOff = 0);

    Int_t Array2Matrix(const Double_t Array[], Int_t UseFPOff = 0); // Array -> fCurrentMatrixElems
    Int_t Array2Matrix(const Double_t Array[], std::vector<THaMatrixElement>& Matrix, Int_t UseFPOff = 0);

    ///////////////////////////////////////////////////////////////////////////
    // Data storage
    ///////////////////////////////////////////////////////////////////////////
    void DCS2FCS(const Double_t* det, Double_t* rot);

    enum {
        kMaxNEventData = 50, kMaxNRawData = 2000000, kPreCalPow = 10, kMaxDataGroup = 12 * 12 * 7 * 7, kMaxDataPerGroup = 10000
    };

    Int_t LoadRawData(TString DataFileName, Int_t NLoad = kMaxNRawData, Int_t MaxDataPerGroup = kMaxDataPerGroup); // load data to Rawdata[]

    typedef struct {
        Double_t Data[kMaxNEventData]; // [CommonIdx]
        Double_t powers[kPreCalPow][5]; // {(x), th, y, ph, abs(th) }
    } EventData;
    EventData fRawData[kMaxNRawData]; // [fNRawData]
    Int_t fNRawData;

    enum CommonIdx {
        kCutID = 0, // cut ID in order of tree2ascii cut file
        kDetX = 1, // L.tr.d_x
        kDetTh = 2, // L.tr.d_th
        kDetY = 3, // L.tr.d_y
        kDetPh = 4, // L.tr.d_ph
        kBeamE = 5, // Beam energy
        kBeamX = 6, // urb.x or rb.x
        kBeamY = 7, // urb.y or rb.y
        kSimX = 8,
        kSimTh = 9,
        kSimY = 10,
        kSimPh = 11,
        kSimOrTh = 12, // use this to calculate kSimDp
        kSimOrPh = 13
    };

    enum ExtraSieveIdx {
        kRotX = 14, // L.tr.r_x
        kRotTh, // L.tr.r_th
        kRotY, // L.tr.r_y
        kRotPh, // L.tr.r_ph
        kRealTh, // real target th from survey (or sim)
        kRealPh, // real target ph from survey (or sim)
        kRealX, // real target x from survey, beam (or sim)
        kRealY, // real target y from survey, beam (or sim)
        kRealThMatrix, // expected target th before extended target corrections
        kRealPhMatrix, // expected target ph before extended target corrections
        kCalcTh, // calculated th from matrix
        kCalcPh // calculated ph from matrix
    };

    enum ExtraVertexIdx {
        kFoilID = 30,
        kRealReactZ, // expected ReactZ
        kCalcY, // calculated y_tg
        kCalcReactZ // calculated ReactZ
    };

    enum ExtraDpIdx {
        kKineID = 30, // Delta Scan Kinematics ID
        kCentralp, // Central Momentum
        kRadiLossDp, // Radiation Loss for this event in unit of dp
        kScatterAngle, // Scattering Angle
        kDpKinOffsets, // = dp - dp_kin, based on sieve hole survey
        kRealDpKin, // expected dp_kin, before radiation correction
        kRealDpKinMatrix, // real dp kin before extended target corrections
        kCalcDpKinMatrix, // calculated dp kin before extended target corrections
        kCalcDpKin, // calculated dp_kin, before radiation correction
        kRealDpKinExcitations // first index of expected dp_kins for all excitation states
    };

    ///////////////////////////////////////////////////////////////////////////
    // Optimization related Commands
    ///////////////////////////////////////////////////////////////////////////
    const TVector3 GetSieveHoleTCS(Int_t Col, Int_t Row);

    void PrepareSieve(void);
    Double_t SumSquareDTh(Int_t UseFPOff = 0);
    Double_t SumSquareDPhi(Int_t UseFPOff = 0);
    TList* CheckSieve(Int_t PlotType = 0);

    void PrepareVertex(void);
    Double_t SumSquareDY(Int_t UseFPOff = 0);
    TList* CheckY(Int_t PlotType = 0);

    Double_t fArbitaryDpKinShift[100]; // compensate bias due to dp event selections, array of [KineID]
    void PrepareDp(void);
    Double_t SumSquareDp(Int_t UseFPOff = 0);
    TList* CheckDp(Int_t PlotType = 0);

    TRotation fTCSInHCS; // transformations vector from TCS to HCS
    TVector3 fPointingOffset; // Optical point in lab coordinate system

    Double_t ScatMom(Double_t DM, Double_t Ma, Double_t P0, Double_t Theta);

    ///////////////////////////////////////////////////////////////////////////
    // declarations for target vertex reconstruction
    ///////////////////////////////////////////////////////////////////////////

    enum {
        kPORDER = 7
    };

    friend class THaMatrixElement;
    std::vector<THaMatrixElement>* fCurrentMatrixElems;
    // initial matrix elements
    std::vector<THaMatrixElement> fTMatrixElems;
    std::vector<THaMatrixElement> fDMatrixElems;
    std::vector<THaMatrixElement> fPMatrixElems;
    std::vector<THaMatrixElement> fPTAMatrixElems; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fYMatrixElems;
    std::vector<THaMatrixElement> fYTAMatrixElems; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fFPMatrixElems; // Focus Plane offsets
    std::vector<THaMatrixElement> fLMatrixElems; // Path-length corrections (meters)

    void CalcMatrix(const Double_t x, std::vector<THaMatrixElement> &matrix);
    Double_t CalcTargetVar(const std::vector<THaMatrixElement> &matrix, const Double_t powers[][5]);

    ///////////////////////////////////////////////////////////////////////////
    // Inherited from THaTrackingDetector
    ///////////////////////////////////////////////////////////////////////////
    virtual Int_t Decode(const THaEvData&);
    virtual Int_t CoarseTrack(TClonesArray&);
    virtual Int_t FineTrack(TClonesArray&);
    virtual EStatus Init(const TDatime&);

private:
    ClassDef(LOpticsOpt, 0) // HRS Optics Optimizer
};

///////////////////////////////////////////////////////////////////////////////
// class for storing matrix element data
///////////////////////////////////////////////////////////////////////////////

class THaMatrixElement {
public:

    THaMatrixElement();

    Bool_t match(const THaMatrixElement& rhs) const;

    Bool_t iszero; // whether the element is zero
    std::vector<Int_t> pw; // exponents of matrix element
    //   e.g. D100 = { 1, 0, 0 }
    Int_t order;
    Double_t v; // its computed value
    std::vector<Double_t> poly; // the associated polynomial

    void SkimPoly(); //reduce order to highest non-zero poly

    Int_t OptOrder; //order optimize to
};

inline Double_t LOpticsOpt::ScatMom(Double_t DM, Double_t Ma, Double_t P0, Double_t Theta)
{
    // Calc Scattered Electron Momentum
    // Assuming Electron with P0, Scattering on a fix target of Mass Ma Assuming
    // recoil system is a resonance of DM and scattering angle is Theta.
    return (-DM * DM - 2 * Ma * DM + 2 * Ma * P0) / (2 * (Ma + P0 - P0 * TMath::Cos(Theta)));
}

#endif
