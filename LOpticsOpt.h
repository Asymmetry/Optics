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
//
// Modification:
//   Jun 25, 2010 Updated for APEX optics calibration
//   Aug 01, 2013 Updated for G2P optics calibration (Chao Gu)
//
///////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_LOpticsOpt
#define ROOT_LOpticsOpt

#include <vector>

#include "TRotation.h"
#include "TMath.h"

#include "THaTrackingDetector.h"
#include "THaString.h"

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
class TCanvas;
class THaTrack;
class TClonesArray;
class TTree;
class TVector3;

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

    UInt_t Matrix2Array(Double_t Array[], Bool_t FreeParaFlag[] = NULL, Int_t UseFPOff = 0) // fCurrentMatrixElems -> Array
    {
        assert(fCurrentMatrixElems);
        return Matrix2Array(Array, (*fCurrentMatrixElems), FreeParaFlag, UseFPOff);
    }
    UInt_t Matrix2Array(Double_t Array[], const std::vector<THaMatrixElement> &Matrix, Bool_t FreeParaFlag[] = NULL, Int_t UseFPOff = 0);

    UInt_t Array2Matrix(const Double_t Array[], Int_t UseFPOff = 0) // Array -> fCurrentMatrixElems
    {
        assert(fCurrentMatrixElems);
        return Array2Matrix(Array, (*fCurrentMatrixElems), UseFPOff);
    }
    UInt_t Array2Matrix(const Double_t Array[], std::vector<THaMatrixElement> &Matrix, Int_t UseFPOff = 0);

    ///////////////////////////////////////////////////////////////////////////
    // Data storage
    ///////////////////////////////////////////////////////////////////////////
    void DCS2FCS(const Double_t* det, Double_t* rot);
    UInt_t LoadRawData(TString DataFileName, UInt_t NLoad = MaxNRawData, UInt_t MaxDataPerGroup = (UInt_t) - 1); // load data to Rawdata[]

    enum {
        MaxNEventData = 50, MaxNRawData = 2000000, kNUM_PRECOMP_POW = 10, kMaxDataGroup = 180 * 5 * 5
    };

    typedef struct {
        Double_t Data[MaxNEventData]; // [CommonIdx]
        Double_t powers[kNUM_PRECOMP_POW][5]; // {(x), th, y, ph, abs(th) }
    } EventData;
    EventData fRawData[MaxNRawData]; // [fNRawData]
    UInt_t fNRawData;
    UInt_t fNCalibData; // for dp calib only

    enum CommonIdx {
        kCutID = 0, // cut ID in order of tree2ascii cut file
        kX = 1, // L.tr.d_x
        kTh = 2, // L.tr.d_th
        kY = 3, // L.tr.d_y
        kPh = 4, // L.tr.d_ph
        kurb_e = 5, // Beam energy
        kBeamX = 6, // urb.x or rb.x
        kBeamY = 7, // urb.y or rb.y
        kSimX = 8,
        kSimTh = 9,
        kSimY = 10,
        kSimPh = 11,
        kOrinTh = 12,
        kOrinPh = 13
    };

    enum ExtraSieveIdx {
        kRotX = 30, // L.tr.r_x
        kRotTh, // L.tr.r_th
        kRotY, // L.tr.r_y
        kRotPh, // L.tr.r_ph
        kRealTh, // real target th from survey
        kRealPh, // real target ph from survey
        kRealTgX, // real target x from survey, beam
        kRealTgY, // real target y from survey, beam
        kRealThMatrix, // expected target th before extended target corrections
        kRealPhMatrix,
        kCalcTh, // calculated th from matrix
        kCalcPh, // calculated ph from matrix
        kSieveX,
        kSieveY,
        kSieveZ,
        kBeamZ
    };

    enum ExtraDpIdx {
        kExtraDataFlag = 14, //Whether this event is for optimization; 0=used for optimization, 1=for plotting only
        kKineID, //Delta Scan Kinematics ID
        kCentralp, //Central Momentum
        kRadiLossDp, //Radiation Loss for this event in unit of dp
        kScatterAngle, //Scattering Angle
        kDpKinOffsets, //=dp-dp_kin, based on sieve hole survey
        kRealDpKin, //expected dp_kin, before radiation correction
        kRealDpKinMatrix, //real dp kin before extended target corrections
        kCalcDpKinMatrix, //calculated dp kin before extended target corrections
        kCalcDpKin, //calculated dp_kin, before radiation correction
        kRealDpKinExcitations/*Do not append more index*/ //first index of expected dp_kins for all excitation states
    };

    ///////////////////////////////////////////////////////////////////////////
    // Optimization related Commands
    ///////////////////////////////////////////////////////////////////////////
    const TVector3 GetSieveHoleTCS(UInt_t Col, UInt_t Row);

    void PrepareSieve(void);
    void PrepareSieveWithField(void);
    Double_t SumSquareDTh(Int_t UseFPOff = 0);
    Double_t SumSquareDPhi(Int_t UseFPOff = 0);
    TCanvas* CheckSieve(Int_t PlotKine = 0, UInt_t PlotFoilID = 0);

    Double_t fArbitaryDpKinShift[100]; // compensate bias due to dp event selections, array of [KineID]
    void PrepareDp(void);
    void PrepareDpWithField(void);
    Double_t SumSquareDp(Bool_t IncludeExtraData = kFALSE);
    TCanvas* CheckDp(void);
    TCanvas* CheckDpGlobal(void);

    TRotation fTCSInHCS; // transformations vector from TCS to HCS
    TVector3 fPointingOffset; // Optical point in lab coordinate system

    inline Double_t ScatMom(Double_t DM, Double_t Ma, Double_t P0, Double_t Theta)
    {
        // Calc Scattered Electron Momentum
        // Assuming Electron with P0, Scattering on a fix target of Mass Ma Assuming
        // recoil system is a resonance of DM and scattering angle is Theta.
        return (-DM * DM - 2 * Ma * DM + 2 * Ma * P0) / (2 * (Ma + P0 - P0 * TMath::Cos(Theta)));
    }

    ///////////////////////////////////////////////////////////////////////////
    // declarations for target vertex reconstruction
    ///////////////////////////////////////////////////////////////////////////

    enum ECoordTypes {
        kTransport, kRotatingTransport
    };

    enum EFPMatrixElemTags {
        T000 = 0, Y000, P000
    };

    enum {
        kPORDER = 7
    };

    friend class THaMatrixElement;
    std::vector<THaMatrixElement> * fCurrentMatrixElems;
    // initial matrix elements
    std::vector<THaMatrixElement> fTMatrixElems;
    std::vector<THaMatrixElement> fDMatrixElems;
    std::vector<THaMatrixElement> fPMatrixElems;
    std::vector<THaMatrixElement> fPTAMatrixElems; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fYMatrixElems;
    std::vector<THaMatrixElement> fYTAMatrixElems; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fFPMatrixElems; // Focus Plane offsets
    std::vector<THaMatrixElement> fLMatrixElems; // Path-length corrections (meters)

    void CalcMatrix(const double x, std::vector<THaMatrixElement> &matrix);
    // Double_t DoPoly(const int n, const std::vector<double> &a, const double x);
    // Double_t PolyInv(const double x1, const double x2, const double xacc, const double y, const int norder, const std::vector<double> &a);
    Double_t CalcTargetVar(const std::vector<THaMatrixElement> &matrix, const double powers[][5]);

    ///////////////////////////////////////////////////////////////////////////
    // Inherited from THaTrackingDetector
    ///////////////////////////////////////////////////////////////////////////

    virtual Int_t Decode(const THaEvData&)
    {
        return 0;
    }

    virtual Int_t CoarseTrack(TClonesArray&)
    {
        return 0;
    }

    virtual Int_t FineTrack(TClonesArray&)
    {
        return 0;
    }

    virtual EStatus Init(const TDatime&)
    {
        return fStatus = kOK;
    };

private:
    ClassDef(LOpticsOpt, 0) // HRS Optics Optimizer
};

///////////////////////////////////////////////////////////////////////////////
// class for storing matrix element data
///////////////////////////////////////////////////////////////////////////////

class THaMatrixElement {
public:

    THaMatrixElement() : iszero(true), pw(3), order(0), v(0), poly(LOpticsOpt::kPORDER), OptOrder(0)
    {
    }
    bool match(const THaMatrixElement& rhs) const;

    bool iszero; // whether the element is zero
    std::vector<int> pw; // exponents of matrix element
    //   e.g. D100 = { 1, 0, 0 }
    int order;
    double v; // its computed value
    std::vector<double> poly; // the associated polynomial

    void SkimPoly(); //reduce order to highest non-zero poly

    UInt_t OptOrder; //order optimize to
};
#endif
