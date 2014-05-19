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

#include <cstdio>
#include <cstdlib>
#include <map>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TDatime.h"
#include "TGraphErrors.h"
#include "TClonesArray.h"
#include "TList.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TLatex.h"
#include "TVector3.h"
#include "TLine.h"
#include "TArrow.h"
#include "TMath.h"
#include "TString.h"

#include "THaGlobals.h"
#include "THaEvData.h"
#include "THaDetMap.h"
#include "THaTrack.h"
#include "THaScintillator.h"
#include "THaSpectrometer.h"
#include "VarDef.h"

#include "LOpticsOpt.h"

#ifdef WITH_DEBUG
#include <iostream>
#endif

using namespace std;
using THaString::Split;

///////////////////////////////////////////////////////////////////////////////
// Input Sections
///////////////////////////////////////////////////////////////////////////////

#include "InputG2PL.h"

///////////////////////////////////////////////////////////////////////////////
// Constructors
///////////////////////////////////////////////////////////////////////////////

LOpticsOpt::LOpticsOpt(const char* name, const char* description, THaApparatus* apparatus) :
THaTrackingDetector(name, description, apparatus)
{
    fPrefix = new char[1000];
    sprintf(fPrefix, "%s", Prefix);

    fCurrentMatrixElems = NULL;

    TVector3 TCSX(0, -1, 0);
    TVector3 TCSZ(TMath::Sin(HRSAngle), 0, TMath::Cos(HRSAngle));
    TVector3 TCSY = TCSZ.Cross(TCSX);
    fTCSInHCS.RotateAxes(TCSX, TCSY, TCSZ);

    fPointingOffset.SetXYZ(-MissPointZ * TMath::Sin(HRSAngle) * TMath::Cos(HRSAngle), MissPointY, MissPointZ * TMath::Sin(HRSAngle) * TMath::Sin(HRSAngle));
    DEBUG_INFO("LOpticsOpt", "Read in configuration " + InputID);
    DEBUG_INFO("LOpticsOpt", "HRS @ %f Degree, PointingOffset = (%f,%f,%f), SievePos = (%f,%f,%f)", HRSAngle / TMath::Pi()*180, fPointingOffset.X(), fPointingOffset.Y(), fPointingOffset.Z(), SieveOffX, SieveOffY, ZPos);

    fNRawData = 0;

    for (UInt_t i = 0; i < 100; i++)
        fArbitaryDpKinShift[i] = 0;
}

LOpticsOpt::~LOpticsOpt()
{
    // Destructor.
}

///////////////////////////////////////////////////////////////////////////////
// Database I/O
///////////////////////////////////////////////////////////////////////////////

Int_t LOpticsOpt::LoadDataBase(TString DataBaseName)
{
    static const char* const here = "LoadDataBase";
    OldComments = "";

    FILE* file = fopen(DataBaseName, "r");
    if (!file) {
        Error("LoadDataBase", "%s can not be opened", DataBaseName.Data());
        assert(0); //
        return kFileError;
    } else DEBUG_INFO("LoadDataBase", "Parsing Database %s", DataBaseName.Data());

    const int LEN = 200;
    char buff[LEN];

    // Look for the section [<prefix>.global] in the file, e.g. [ R.global ]
    TString tag(fPrefix);
    Ssiz_t tagpos = tag.Index(".");
    if (tagpos != kNPOS)
        tag = tag(0, tagpos + 1);
    else
        tag.Append(".");
    tag.Prepend("[");
    tag.Append("global]");
    TString tag2(tag);
    tag.ToLower();

    bool found = false;
    while (!found && fgets(buff, LEN, file) != NULL) {
        // read in comments
        TString line = buff;
        if (line.BeginsWith("#")) {
            OldComments += line;
            // OldComments += "\n";
        }

        line = ::Compress(buff); // strip blanks
        if (line.EndsWith("\n")) line.Chop();

        line.ToLower();
        if (tag == line) found = true;
    }
    if (!found) {
        Error(Here(here), "Database entry %s not found!", tag2.Data());
        fclose(file);
        assert(0); //
        return kInitError;
    }

    // We found the section, now read the data
    fgets(buff, LEN, file); // Skip constant line
    fgets(buff, LEN, file); // Skip comment line

    fTMatrixElems.clear();
    fDMatrixElems.clear();
    fPMatrixElems.clear();
    fPTAMatrixElems.clear();
    fYMatrixElems.clear();
    fYTAMatrixElems.clear();
    fLMatrixElems.clear();

    fFPMatrixElems.clear();
    fFPMatrixElems.resize(3);

    typedef vector<string>::size_type vsiz_t;
    map<string, vsiz_t> power;
    power["t"] = 3; // transport to focal-plane tensors
    power["y"] = 3;
    power["p"] = 3;
    power["D"] = 3; // focal-plane to target tensors
    power["T"] = 3;
    power["Y"] = 3;
    power["YTA"] = 4;
    power["P"] = 3;
    power["PTA"] = 4;
    power["L"] = 4; // pathlength from z=0 (target) to focal plane (meters)
    power["XF"] = 5; // forward: target to focal-plane (I think)
    power["TF"] = 5;
    power["PF"] = 5;
    power["YF"] = 5;

    map<string, vector<THaMatrixElement>*> matrix_map;
    matrix_map["t"] = &fFPMatrixElems;
    matrix_map["y"] = &fFPMatrixElems;
    matrix_map["p"] = &fFPMatrixElems;
    matrix_map["D"] = &fDMatrixElems;
    matrix_map["T"] = &fTMatrixElems;
    matrix_map["Y"] = &fYMatrixElems;
    matrix_map["YTA"] = &fYTAMatrixElems;
    matrix_map["P"] = &fPMatrixElems;
    matrix_map["PTA"] = &fPTAMatrixElems;
    matrix_map["L"] = &fLMatrixElems;

    map<string, int> fp_map;
    fp_map["t"] = 0;
    fp_map["y"] = 1;
    fp_map["p"] = 2;

    // Read in as many of the matrix elements as there are.
    // Read in line-by-line, so as to be able to handle tensors of
    // different orders.
    while (fgets(buff, LEN, file)) {
        string line(buff);
        // Erase trailing newline
        if (line.size() > 0 && line[line.size() - 1] == '\n') {
            buff[line.size() - 1] = 0;
            line.erase(line.size() - 1, 1);
        }
        // Split the line into whitespace-separated fields
        vector<string> line_spl = Split(line);

        // Stop if the line does not start with a string referring to
        // a known type of matrix element. In particular, this will
        // stop on a subsequent timestamp or configuration tag starting with "["
        if (line_spl.empty()) continue; //ignore empty lines
        const char* w = line_spl[0].c_str();
        vsiz_t npow = power[w];
        if (npow == 0) break;

#if DEBUG_LEVEL>=4
        cout << "Matrix Line = ";
        for (Ssiz_t i = 1; (UInt_t) i < (UInt_t) line_spl.size(); i++) {
            cout << i << "(" << line_spl[i].c_str() << "), ";
        }
        cout << endl;
#endif

        // Looks like a good line, go parse it.
        THaMatrixElement ME;
        ME.pw.resize(npow);
        ME.iszero = true;
        ME.order = 0;
        vsiz_t pos;
        for (pos = 1; pos <= npow && pos < line_spl.size(); pos++) {
            ME.pw[pos - 1] = atoi(line_spl[pos].c_str());
        }
        vsiz_t p_cnt;
        for (p_cnt = 0; pos < line_spl.size() && p_cnt < kPORDER && pos <= npow + kPORDER; pos++, p_cnt++) {
            ME.poly[p_cnt] = atof(line_spl[pos].c_str());
            if (ME.poly[p_cnt] != 0.0) {
                ME.iszero = false;
                ME.order = p_cnt + 1;
            }
        }
        if (p_cnt < 1) {
            Error(Here(here), "Could not read in Matrix Element %s%d%d%d!", w, ME.pw[0], ME.pw[1], ME.pw[2]);
            Error(Here(here), "Line looks like: %s", line.c_str());
            fclose(file);
            return kInitError;
        }

        // order optimize to
        ME.OptOrder = atoi(line_spl[line_spl.size() - 1].c_str());

        // Don't bother with all-zero matrix elements
        if (ME.iszero) continue;

        // Add this matrix element to the appropriate array
        vector<THaMatrixElement> *mat = matrix_map[w];
        if (mat) {
            // Special checks for focal plane matrix elements
            if (mat == &fFPMatrixElems) {
                if (ME.pw[0] == 0 && ME.pw[1] == 0 && ME.pw[2] == 0) {
                    THaMatrixElement& m = (*mat)[fp_map[w]];
                    if (m.order > 0) {
                        Warning(Here(here), "Duplicate definition of focal plane matrix element: %s. Using first definition.", buff);
                    } else
                        m = ME;
                } else
                    Warning(Here(here), "Bad coefficients of focal plane matrix element %s", buff);
            } else {
                // All other matrix elements are just appended to the respective array
                // but ensure that they are defined only once!
                bool match = false;
                for (vector<THaMatrixElement>::iterator it = mat->begin(); it != mat->end() && !(match = it->match(ME)); it++) {
                }
                if (match) {
                    Warning(Here(here), "Duplicate definition of matrix element: %s. Using first definition.", buff);
                } else
                    mat->push_back(ME);
            }
        } else if (fDebug > 0)
            Warning(Here(here), "Not storing matrix for: %s !", w);
    }

    CalcMatrix(1., fLMatrixElems); // tensor without explicit polynomial in x_fp

    fIsInit = true;
    fclose(file);
    return kOK;
}

Int_t LOpticsOpt::SaveDataBase(TString DataBaseName)
{
    // Output database in memory to new database file
    // WARNING: Hard coded text included

    DEBUG_INFO("SaveDataBase", "Saving to %s", DataBaseName.Data());

    typedef vector<THaMatrixElement>::size_type vsiz_t;

    FILE* file = fopen(DataBaseName, "w");
    if (!file) {
        Info("SaveDataBase", "Error Open in %s", DataBaseName.Data());
        return kFileError;
    }
    TDatime dt;
    fprintf(file, "# -------------------------------------------------------------");
    fprintf(file, "\n");
    fprintf(file, "# Optimized by Chao Gu @ %s", dt.AsString());
    fprintf(file, "\n");
    fprintf(file, "# Saved to %s", DataBaseName.Data());
    fprintf(file, "\n");
    fprintf(file, OldComments);

    fprintf(file, DatabaseHeader);

    DEBUG_INFO("SaveDataBase", "Transport Matrix:  D-terms");
    for (vsiz_t i = 0; i < fDMatrixElems.size(); i++) {
        fprintf(file, "D ");
        const THaMatrixElement& m = fDMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            fprintf(file, "%d ", m.pw[j]);
        }
        int j;
        for (j = 0; j < m.order; j++) {
            fprintf(file, " %13.6e", m.poly[j]);
        }
        for (; j < kPORDER; j++) {
            fprintf(file, " %13.6e", 0.0);
        }
        fprintf(file, "  %d", m.OptOrder);
        fprintf(file, "\n");
    }

    DEBUG_INFO("SaveDataBase", "Transport Matrix:  T-terms");
    for (vsiz_t i = 0; i < fTMatrixElems.size(); i++) {
        fprintf(file, "T ");
        const THaMatrixElement& m = fTMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            fprintf(file, "%d ", m.pw[j]);
        }
        int j;
        for (j = 0; j < m.order; j++) {
            fprintf(file, " %13.6e", m.poly[j]);
        }
        for (; j < kPORDER; j++) {
            fprintf(file, " %13.6e", 0.0);
        }
        fprintf(file, "  %d", m.OptOrder);
        fprintf(file, "\n");
    }

    DEBUG_INFO("SaveDataBase", "Transport Matrix:  P-terms");
    for (vsiz_t i = 0; i < fPMatrixElems.size(); i++) {
        fprintf(file, "P ");
        const THaMatrixElement& m = fPMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            fprintf(file, "%d ", m.pw[j]);
        }
        int j;
        for (j = 0; j < m.order; j++) {
            fprintf(file, " %13.6e", m.poly[j]);
        }
        for (; j < kPORDER; j++) {
            fprintf(file, " %13.6e", 0.0);
        }
        fprintf(file, "  %d", m.OptOrder);
        fprintf(file, "\n");
    }

    DEBUG_INFO("SaveDataBase", "Transport Matrix:  Y-terms");
    for (vsiz_t i = 0; i < fYMatrixElems.size(); i++) {
        fprintf(file, "Y ");
        const THaMatrixElement& m = fYMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            fprintf(file, "%d ", m.pw[j]);
        }
        int j;
        for (j = 0; j < m.order; j++) {
            fprintf(file, " %13.6e", m.poly[j]);
        }
        for (; j < kPORDER; j++) {
            fprintf(file, " %13.6e", 0.0);
        }
        fprintf(file, "  %d", m.OptOrder);
        fprintf(file, "\n");
    }

    fprintf(file, DatabaseFooter);

    fclose(file);

    return kOK;
}

void LOpticsOpt::Print(const Option_t* opt) const
{
    // Print current matrix

    THaTrackingDetector::Print(opt);
    typedef vector<THaMatrixElement>::size_type vsiz_t;

    // Print out the optics matrices, to verify they make sense
    //    printf("Matrix FP (t000, y000, p000)\n");
    //    for (vsiz_t i = 0; i < fFPMatrixElems.size(); i++) {
    //        const THaMatrixElement& m = fFPMatrixElems[i];
    //        for (vsiz_t j = 0; j < m.pw.size(); j++) {
    //            printf("  %2d", m.pw[j]);
    //        }
    //        for (int j = 0; j < m.order; j++) {
    //            printf("  %g", m.poly[j]);
    //        }
    //        printf(" : Opt -> %d", m.OptOrder);
    //        printf("\n");
    //    }

    printf("LOpticsOpt::Print: Transport Matrix:  D-terms\n");
    for (vsiz_t i = 0; i < fDMatrixElems.size(); i++) {
        const THaMatrixElement& m = fDMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            printf("  %2d", m.pw[j]);
        }
        for (int j = 0; j < m.order; j++) {
            printf("\t%g", m.poly[j]);
        }
        printf(" : Opt -> %d", m.OptOrder);
        if ((UInt_t) m.order != m.OptOrder) printf(" != Matrix Order !!");
        printf("\n");
    }

    printf("LOpticsOpt::Print: Transport Matrix:  T-terms\n");
    for (vsiz_t i = 0; i < fTMatrixElems.size(); i++) {
        const THaMatrixElement& m = fTMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            printf("  %2d", m.pw[j]);
        }
        for (int j = 0; j < m.order; j++) {
            printf("\t%g", m.poly[j]);
        }
        printf(" : Opt -> %d", m.OptOrder);
        if ((UInt_t) m.order != m.OptOrder) printf(" != Matrix Order !!");
        printf("\n");
    }

    printf("LOpticsOpt::Print: Transport Matrix:  Y-terms\n");
    for (vsiz_t i = 0; i < fYMatrixElems.size(); i++) {
        const THaMatrixElement& m = fYMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            printf("  %2d", m.pw[j]);
        }
        for (int j = 0; j < m.order; j++) {
            printf("\t%g", m.poly[j]);
        }
        printf(" : Opt -> %d", m.OptOrder);
        if ((UInt_t) m.order != m.OptOrder) printf(" != Matrix Order !!");
        printf("\n");
    }

    //    printf("Transport Matrix:  YTA-terms (abs(theta))\n");
    //    for (vsiz_t i = 0; i < fYTAMatrixElems.size(); i++) {
    //        const THaMatrixElement& m = fYTAMatrixElems[i];
    //        for (vsiz_t j = 0; j < m.pw.size(); j++) {
    //            printf("  %2d", m.pw[j]);
    //        }
    //        for (int j = 0; j < m.order; j++) {
    //            printf("\t%g", m.poly[j]);
    //        }
    //        printf(" : Opt -> %d", m.OptOrder);
    //        printf("\n");
    //    }

    printf("LOpticsOpt::Print: Transport Matrix:  P-terms\n");
    for (vsiz_t i = 0; i < fPMatrixElems.size(); i++) {
        const THaMatrixElement& m = fPMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            printf("  %2d", m.pw[j]);
        }
        for (int j = 0; j < m.order; j++) {
            printf("\t%g", m.poly[j]);
        }
        printf(" : Opt -> %d", m.OptOrder);
        if ((UInt_t) m.order != m.OptOrder) printf(" != Matrix Order !!");
        printf("\n");
    }

    //    printf("Transport Matrix:  PTA-terms\n");
    //    for (vsiz_t i = 0; i < fPTAMatrixElems.size(); i++) {
    //        const THaMatrixElement& m = fPTAMatrixElems[i];
    //        for (vsiz_t j = 0; j < m.pw.size(); j++) {
    //            printf("  %2d", m.pw[j]);
    //        }
    //        for (int j = 0; j < m.order; j++) {
    //            printf("\t%g", m.poly[j]);
    //        }
    //        printf(" : Opt -> %d", m.OptOrder);
    //        printf("\n");
    //    }
    //
    //    printf("Matrix L\n");
    //    for (vsiz_t i = 0; i < fLMatrixElems.size(); i++) {
    //        const THaMatrixElement& m = fLMatrixElems[i];
    //        for (vsiz_t j = 0; j < m.pw.size(); j++) {
    //            printf("  %2d", m.pw[j]);
    //        }
    //        for (int j = 0; j < m.order; j++) {
    //            printf("\t%g", m.poly[j]);
    //        }
    //        printf(" : Opt -> %d", m.OptOrder);
    //        printf("\n");
    //    }

    printf("fArbitaryDpKinShift[%d] = {", NKine);
    for (UInt_t KineID = 0; KineID < NKine; KineID++)
        printf("%f  ", fArbitaryDpKinShift[KineID]);
    printf("}\n");

    return;
}

UInt_t LOpticsOpt::Matrix2Array(Double_t Array[], const std::vector<THaMatrixElement> &Matrix, Bool_t FreeParaFlag[])
{
    // Matrix -> Array

    typedef vector<THaMatrixElement>::size_type vsiz_t;

    UInt_t idx = 0;

    for (vsiz_t i = 0; i < Matrix.size(); i++) {
        const THaMatrixElement& m = Matrix[i];
        UInt_t j;
        for (j = 0; (int) j < m.order; j++) {
            if (FreeParaFlag) FreeParaFlag[idx] = j < m.OptOrder ? kTRUE : kFALSE;
            Array[idx++] = m.poly[j];
        }
        for (; j < kPORDER; j++) {
            if (FreeParaFlag) FreeParaFlag[idx] = j < m.OptOrder ? kTRUE : kFALSE;
            Array[idx++] = 0;
        }
    }

    DEBUG_INFO("Matrix2Array", "Fill Size = %d", idx);

    return idx;
}

UInt_t LOpticsOpt::Array2Matrix(const Double_t Array[], std::vector<THaMatrixElement> &Matrix)
{
    // Array -> fCurrentMatrixElems

    typedef vector<THaMatrixElement>::size_type vsiz_t;

    UInt_t idx = 0;
    Double_t maxele = 0, sumele = 0;
    vsiz_t max_i = -1;
    int max_j = -1;

    for (vsiz_t i = 0; i < Matrix.size(); i++) {
        THaMatrixElement& m = Matrix[i];
        int j;
        m.order = kPORDER;
        for (j = 0; j < m.order; j++) {
            if (TMath::Abs(Array[idx]) > maxele) {
                max_i = i;
                max_j = j;
            }
            sumele += TMath::Abs(Array[idx]);
            m.poly[j] = Array[idx];
            idx++;
        }
        m.SkimPoly();
    }

    assert(max_i < Matrix.size());
    assert(max_j >= 0);
    THaMatrixElement& m = Matrix[max_i];
    DEBUG_INFO("Array2Matrix", "Load Size = %d, max ele = (%d %d %d %d) = %f, average = %f", idx, max_j, m.pw[0], m.pw[1], m.pw[2], m.poly[max_j], sumele / idx);

    return idx;
}

///////////////////////////////////////////////////////////////////////////////
// Data storage
///////////////////////////////////////////////////////////////////////////////

UInt_t LOpticsOpt::LoadRawData(TString DataFileName, UInt_t NLoad, UInt_t MaxDataPerGroup)
{
    // Load "f51" ascii data file to Rawdata[]

    DEBUG_INFO("LoadRawData", "Loading %s", DataFileName.Data());

    if (BeamShiftX != 0)
        DEBUG_WARNING("LoadRawData", "Shift Beam X = %f", BeamShiftX);

    UInt_t datagrpcnt[kMaxDataGroup] = {0};

    FILE* file = fopen(DataFileName, "r");
    if (!file) return kFileError;

    UInt_t NRead = 0;
    const int LEN = 2000;
    char buff[LEN];

    Double_t NDataRead = 0;
    int NLineRead = 0;

    while (fgets(buff, LEN, file)) {
        NLineRead++;

        if (NLineRead % 100000 == 0) DEBUG_INFO("LoadRawData", "%d/%d Entries Loaded", NRead, NLineRead);

        assert(NRead < MaxNRawData); //too much data if fails

        if (NRead >= NLoad) break;

        Double_t * eventdata = fRawData[NRead].Data;

        string line(buff);
        // Erase trailing newline
        if (line.size() > 0 && line[line.size() - 1] == '\n') {
            buff[line.size() - 1] = 0;
            line.erase(line.size() - 1, 1);
        }
        // Split the line into whitespace-separated fields
        vector<string> line_spl = Split(line);

        assert(line_spl.size() <= MaxNEventData); // array size check
        for (UInt_t idx = 0; idx < line_spl.size(); idx++)
            eventdata[idx] = atof(line_spl[idx].c_str());

        // WARNING : shift beam x
        if (BeamShiftX != 0)
            eventdata[kBeamX] += BeamShiftX;

        // determine whether to save this data
        UInt_t cutid = (UInt_t) eventdata[kCutID];
        assert(cutid < kMaxDataGroup); // too many cuts
        UInt_t & grpcnt = datagrpcnt[cutid];
        grpcnt++;
        if (grpcnt > MaxDataPerGroup) {
            DEBUG_MASSINFO("LoadRawData", "ignore data %d from cutid %d (%d ev total)", NLineRead, cutid, grpcnt);
            continue;
        }

        NDataRead += line_spl.size();

        Double_t(*powers)[5] = fRawData[NRead].powers;
        Double_t x_fp = eventdata[kX];
        Double_t th_fp = eventdata[kTh];
        Double_t y_fp = eventdata[kY];
        Double_t ph_fp = eventdata[kPhi];

        // calculate the powers we need
        for (int i = 0; i < kNUM_PRECOMP_POW; i++) {
            powers[i][0] = pow(x_fp, i);
            powers[i][1] = pow(th_fp, i);
            powers[i][2] = pow(y_fp, i);
            powers[i][3] = pow(ph_fp, i);
            powers[i][4] = pow(TMath::Abs(th_fp), i);
        }

        NRead++;
    }

    fclose(file);
    fNRawData = NRead;
    fNCalibData = NRead; // fNCalibData shall be updated later if only part of data read in is for calibration use

    UInt_t goodstatcut = 0, actcutcnt = 0;
    for (int i = 0; i < kMaxDataGroup; i++) {
        if (datagrpcnt[i] > 0) actcutcnt++;

        if (datagrpcnt[i] > MaxDataPerGroup) goodstatcut++;
    }

    DEBUG_INFO("LoadRawData", "Event Limit/Cut = %d, %d / %d ev read, %d / %d cut have enough ev", MaxDataPerGroup, NRead, NLineRead, goodstatcut, actcutcnt);
    DEBUG_INFO("LoadRawData", "%d events x %f record/event read from %s", fNRawData, NDataRead / fNRawData, DataFileName.Data());

    return NRead;
}

///////////////////////////////////////////////////////////////////////////////
// Optimization related Commands
///////////////////////////////////////////////////////////////////////////////

const TVector3 LOpticsOpt::GetSieveHoleTCS(UInt_t Col, UInt_t Row)
{
    assert(Col < NSieveCol);
    assert(Row < NSieveRow);
    TVector3 SieveHoleTCS(SieveOffX + SieveXbyRow[Row], SieveOffY + SieveYbyCol[Col], ZPos);
    return SieveHoleTCS;
}

void LOpticsOpt::PrepareSieve(void)
{
    // Calculate kRealTh, kRealPhi

    // DEBUG_INFO("PrepareSieve","Entry Point");

    Double_t exttargcorr_th = 0, rms_exttargcorr_th = 0;
    Double_t exttargcorr_ph = 0, rms_exttargcorr_ph = 0;

    for (UInt_t idx = 0; idx < fNRawData; idx++) {
        EventData &eventdata = fRawData[idx];

        UInt_t res = (UInt_t) eventdata.Data[kCutID];
        // const UInt_t KineID = res / (NSieveRow * NSieveCol * NFoils); //starting 0!
        res = res % (NSieveRow * NSieveCol * NFoils);
        const UInt_t FoilID = res / (NSieveRow * NSieveCol); //starting 0!
        res = res % (NSieveRow * NSieveCol);
        const UInt_t Col = res / (NSieveRow); //starting 0!
        const UInt_t Row = res % (NSieveRow); //starting 0!

        assert(FoilID < NFoils); //check array index size

        const TVector3 SieveHoleTCS = GetSieveHoleTCS(Col, Row);
        eventdata.Data[kSieveX] = SieveHoleTCS.X();
        eventdata.Data[kSieveY] = SieveHoleTCS.Y();
        eventdata.Data[kSieveZ] = SieveHoleTCS.Z();

        const TVector3 BeamSpotHCS(eventdata.Data[kBeamX], eventdata.Data[kBeamY], targetfoils[FoilID]);
        eventdata.Data[kBeamZ] = targetfoils[FoilID];

        const TVector3 BeamSpotTCS = fTCSInHCS.Inverse()*(BeamSpotHCS - fPointingOffset);
        const TVector3 MomDirectionTCS = SieveHoleTCS - BeamSpotTCS;

        eventdata.Data[kRealTh] = MomDirectionTCS.X() / MomDirectionTCS.Z();
        eventdata.Data[kRealPhi] = MomDirectionTCS.Y() / MomDirectionTCS.Z();

        const Double_t x_tg = BeamSpotTCS.X() - BeamSpotTCS.Z() * eventdata.Data[kRealTh];
        const Double_t y_tg = BeamSpotTCS.Y() - BeamSpotTCS.Z() * eventdata.Data[kRealPhi];
        eventdata.Data[kRealTgX] = x_tg;
        eventdata.Data[kRealTgY] = y_tg;

        // Expected th ph before ext. target correction
        // fDeltaTh = fThetaCorr * x_tg;
        // Double_t theta = trkifo->GetTheta() + fDeltaTh;
        eventdata.Data[kRealThMatrix] = eventdata.Data[kRealTh] - x_tg * ExtTarCor_ThetaCorr;
        eventdata.Data[kRealPhMatrix] = eventdata.Data[kRealPhi] - y_tg * ExtTarCor_PhiCorr;

        exttargcorr_th += x_tg * ExtTarCor_ThetaCorr;
        rms_exttargcorr_th += x_tg * ExtTarCor_ThetaCorr * x_tg * ExtTarCor_ThetaCorr;
        exttargcorr_ph += x_tg * ExtTarCor_PhiCorr;
        rms_exttargcorr_ph += x_tg * ExtTarCor_PhiCorr * y_tg * ExtTarCor_PhiCorr;

        DEBUG_MASSINFO("PrepareSieve", "Real_Th_Matrix = %f,\t Real_Phi = %f", eventdata.Data[kRealThMatrix], eventdata.Data[kRealPhi]);
    }

    DEBUG_INFO("PrepareSieve", "Average Extended Target Correction: th = %f,\t rms_th = %f", exttargcorr_th / fNRawData, TMath::Sqrt(rms_exttargcorr_th / fNRawData));
    DEBUG_INFO("PrepareSieve", "Average Extended Target Correction: ph = %f,\t rms_ph = %f", exttargcorr_ph / fNRawData, TMath::Sqrt(rms_exttargcorr_ph / fNRawData));

    // Make sure kCalcTh, kCalcPh is filled
    SumSquareDTh();
    SumSquareDPhi();

    DEBUG_INFO("PrepareSieve", "Done!");
}

void LOpticsOpt::PrepareSieveWithField(void)
{
    // Calculate kRealTh, kRealPhi

    // DEBUG_INFO("PrepareSieve","Entry Point");

    Double_t exttargcorr_th = 0, rms_exttargcorr_th = 0;
    Double_t exttargcorr_ph = 0, rms_exttargcorr_ph = 0;

    for (UInt_t idx = 0; idx < fNRawData; idx++) {
        EventData &eventdata = fRawData[idx];

        UInt_t res = (UInt_t) eventdata.Data[kCutID];
        // const UInt_t KineID = res / (NSieveRow * NSieveCol * NFoils); //starting 0!
        res = res % (NSieveRow * NSieveCol * NFoils);
        const UInt_t FoilID = res / (NSieveRow * NSieveCol); //starting 0!
        res = res % (NSieveRow * NSieveCol);
        const UInt_t Col = res / (NSieveRow); //starting 0!
        const UInt_t Row = res % (NSieveRow); //starting 0!

        assert(FoilID < NFoils); //check array index size

        const TVector3 SieveHoleTCS = GetSieveHoleTCS(Col, Row);
        eventdata.Data[kSieveX] = SieveHoleTCS.X();
        eventdata.Data[kSieveY] = SieveHoleTCS.Y();
        eventdata.Data[kSieveZ] = SieveHoleTCS.Z();

        eventdata.Data[kBeamZ] = targetfoils[FoilID];

        eventdata.Data[kRealTh] = eventdata.Data[kSimTh];
        eventdata.Data[kRealPhi] = eventdata.Data[kSimPh];

        const Double_t x_tg = eventdata.Data[kSimX];
        const Double_t y_tg = eventdata.Data[kSimY];
        eventdata.Data[kRealTgX] = x_tg;
        eventdata.Data[kRealTgY] = y_tg;

        // Expected th ph before ext. target correction
        // fDeltaTh = fThetaCorr * x_tg;
        // Double_t theta = trkifo->GetTheta() + fDeltaTh;
        eventdata.Data[kRealThMatrix] = eventdata.Data[kRealTh] - x_tg * ExtTarCor_ThetaCorr;
        eventdata.Data[kRealPhMatrix] = eventdata.Data[kRealPhi] - y_tg * ExtTarCor_PhiCorr;

        exttargcorr_th += x_tg * ExtTarCor_ThetaCorr;
        rms_exttargcorr_th += x_tg * ExtTarCor_ThetaCorr * x_tg * ExtTarCor_ThetaCorr;
        exttargcorr_ph += x_tg * ExtTarCor_PhiCorr;
        rms_exttargcorr_ph += x_tg * ExtTarCor_PhiCorr * y_tg * ExtTarCor_PhiCorr;

        DEBUG_MASSINFO("PrepareSieve", "Real_Th_Matrix = %f,\t Real_Phi = %f", eventdata.Data[kRealThMatrix], eventdata.Data[kRealPhi]);
    }

    DEBUG_INFO("PrepareSieve", "Average Extended Target Correction: th = %f,\t rmsth = %f", exttargcorr_th / fNRawData, TMath::Sqrt(rms_exttargcorr_th / fNRawData));

    // Make sure kCalcTh, kCalcPh is filled
    SumSquareDTh();
    SumSquareDPhi();

    DEBUG_INFO("PrepareSieve", "Done!");
}

Double_t LOpticsOpt::SumSquareDTh()
{
    // return square sum of diff between calculated tg_th and expected tg_th

    Double_t dth = 0; //Difference
    Double_t rmsth = 0; //mean square

    static UInt_t NCall = 0;
    NCall++;

    Double_t theta;

    for (UInt_t idx = 0; idx < fNRawData; idx++) {
        EventData &eventdata = fRawData[idx];

        Double_t x_fp = eventdata.Data[kX];
        const Double_t(*powers)[5] = eventdata.powers;

        // calculate the matrices we need
        // CalcMatrix(x_fp, fDMatrixElems);
        CalcMatrix(x_fp, fTMatrixElems);
        // CalcMatrix(x_fp, fYMatrixElems);
        // CalcMatrix(x_fp, fYTAMatrixElems);
        // CalcMatrix(x_fp, fPMatrixElems);
        // CalcMatrix(x_fp, fPTAMatrixElems);

        // calculate the coordinates at the target
        theta = CalcTargetVar(fTMatrixElems, powers);

        dth += theta - eventdata.Data[kRealThMatrix];
        rmsth += (theta - eventdata.Data[kRealThMatrix])*(theta - eventdata.Data[kRealThMatrix]);

        DEBUG_MASSINFO("SumSquareDTh", "D_Th = %f = \t%f - \t%f", theta - eventdata.Data[kRealThMatrix], theta, eventdata.Data[kRealThMatrix]);

        // save the results
        eventdata.Data[kCalcTh] = theta;
    }

    DEBUG_INFO("SumSquareDTh", "#%d : dth = %f,\t rmsth = %f", NCall, dth / fNRawData, TMath::Sqrt(rmsth / fNRawData));
    //printf("SumSquareDTh: #%d : dth = %f,\t rmsth = %f\n", NCall, dth / fNRawData, TMath::Sqrt(rmsth / fNRawData));

    return rmsth;
}

Double_t LOpticsOpt::SumSquareDPhi()
{
    // return square sum of diff between calculated tg_ph and expected tg_ph

    Double_t dphi = 0; //Difference
    Double_t rmsphi = 0; //mean square

    static UInt_t NCall = 0;
    NCall++;

    Double_t phi;

    for (UInt_t idx = 0; idx < fNRawData; idx++) {
        EventData &eventdata = fRawData[idx];

        Double_t x_fp = eventdata.Data[kX];
        const Double_t(*powers)[5] = eventdata.powers;

        // calculate the matrices we need
        // CalcMatrix(x_fp, fDMatrixElems);
        // CalcMatrix(x_fp, fTMatrixElems);
        // CalcMatrix(x_fp, fYMatrixElems);
        // CalcMatrix(x_fp, fYTAMatrixElems);
        CalcMatrix(x_fp, fPMatrixElems);
        CalcMatrix(x_fp, fPTAMatrixElems);

        // calculate the coordinates at the target
        phi = CalcTargetVar(fPMatrixElems, powers) + CalcTargetVar(fPTAMatrixElems, powers);

        dphi += phi - eventdata.Data[kRealPhMatrix];
        rmsphi += (phi - eventdata.Data[kRealPhMatrix])*(phi - eventdata.Data[kRealPhMatrix]);

        DEBUG_MASSINFO("SumSquareDPhi", "D_Phi = %f = \t%f - \t%f", phi - eventdata.Data[kRealPhMatrix], phi, eventdata.Data[kRealPhMatrix]);

        //save the results
        eventdata.Data[kCalcPh] = phi;
    }

    DEBUG_INFO("SumSquareDPhi", "#%d : dphi = %f,\t rmsphi = %f", NCall, dphi / fNRawData, TMath::Sqrt(rmsphi / fNRawData));
    //printf("SumSquareDPhi: #%d : dphi = %f,\t rmsphi = %f\n", NCall, dphi / fNRawData, TMath::Sqrt(rmsphi / fNRawData));

    return rmsphi;
}

TCanvas * LOpticsOpt::CheckSieve(Int_t PlotKine, UInt_t PlotFoilID)
{
    // Visualize Sieve Plane
    DEBUG_INFO("CheckSieve", "Entry Point");

    const UInt_t nplot = (PlotKine != 0) ? NKine : 1;

    TH2D * HSievePlane[NKine] = {0};
    Double_t x_lim[NKine] = {0};
    Double_t y_lim[NKine] = {0};

    for (UInt_t idx = 0; idx < NKine; idx++) {
        x_lim[idx] = 1.3 * TMath::Max(TMath::Abs(SieveYbyCol[0]), TMath::Abs(SieveYbyCol[NSieveCol - 1]));
        y_lim[idx] = 1.5 * TMath::Max(TMath::Abs(SieveXbyRow[0]), TMath::Abs(SieveXbyRow[NSieveRow - 1]));

        HSievePlane[idx] = new TH2D(Form("Sieve_Kine%d", idx), Form("Sieve Plane Proj. (tg_X vs tg_Y) for Data set #%d", idx), 500, -x_lim[idx], x_lim[idx], 500, -y_lim[idx], y_lim[idx]);

        HSievePlane[idx]->SetXTitle("Sieve H [m]");
        HSievePlane[idx]->SetYTitle("Sieve V [m]");
        assert(HSievePlane[idx]); // assure memory allocation
    }

    Double_t dX = 0, dY = 0;

    enum {
        kEventID, kRealSieveX, kRealSieveY, kCalcSieveX, kCalcSieveY
    };

    Double_t SieveEventID[NKine][NSieveCol][NSieveRow][5] = {
        {
            {
                {0}
            }
        }
    };

    for (UInt_t idx = 0; idx < fNRawData; idx++) {
        const EventData &eventdata = fRawData[idx];

        UInt_t res = (UInt_t) eventdata.Data[kCutID];
        UInt_t KineID = res / (NSieveRow * NSieveCol * NFoils); //starting 0!
        res = res % (NSieveRow * NSieveCol * NFoils);
        const UInt_t FoilID = res / (NSieveRow * NSieveCol); //starting 0!
        res = res % (NSieveRow * NSieveCol);
        const UInt_t Col = res / (NSieveRow); //starting 0!
        const UInt_t Row = res % (NSieveRow); //starting 0!

        assert(KineID < NKine); //array index check

        if (FoilID != PlotFoilID) continue;
        if (PlotKine == 0) KineID = 0;

        const TVector3 SieveHoleTCS = GetSieveHoleTCS(Col, Row);

        Double_t ProjectionX = eventdata.Data[kRealTgX] + (eventdata.Data[kCalcTh] + eventdata.Data[kRealTgX] * ExtTarCor_ThetaCorr) * (SieveHoleTCS.Z());
        Double_t ProjectionY = eventdata.Data[kRealTgY] + (eventdata.Data[kCalcPh] + eventdata.Data[kRealTgY] * ExtTarCor_PhiCorr) * (SieveHoleTCS.Z());

        HSievePlane[KineID]->Fill(ProjectionY, ProjectionX);

        dX += ProjectionX - SieveHoleTCS.X();
        dY += ProjectionY - SieveHoleTCS.Y();

        SieveEventID[KineID][Col][Row][kEventID] = idx;
        SieveEventID[KineID][Col][Row][kRealSieveX] = SieveHoleTCS.X();
        SieveEventID[KineID][Col][Row][kRealSieveY] = SieveHoleTCS.Y();
        SieveEventID[KineID][Col][Row][kCalcSieveX] = ProjectionX;
        SieveEventID[KineID][Col][Row][kCalcSieveY] = ProjectionY;
    }

    DEBUG_INFO("CheckSieve", "Average : D_X = %f,\t D_Y = %f", dX / fNRawData, dY / fNRawData);

    TCanvas * c1 = new TCanvas("SieveCheck", "SieveCheck", 1800, 1100);

    if (nplot <= 1) {
        c1 = new TCanvas("SieveCheck", "SieveCheck", 800, 1100);
        c1->Divide(1, 1);
    } else if (nplot <= 3) {
        c1 = new TCanvas("SieveCheck", "SieveCheck", 1800, 1100);
        c1->Divide(3, 1);
    } else if (nplot <= 6) {
        c1 = new TCanvas("SieveCheck", "SieveCheck", 1800, 1100);
        c1->Divide(3, 2);
    } else {
        c1 = new TCanvas("SieveCheck", "SieveCheck", 1800, 1100);
        c1->Divide(4, 2);
    }

    for (UInt_t idx = 0; idx < nplot; idx++) {
        UInt_t KineID = idx;

        c1->cd(idx + 1);

        assert(HSievePlane[idx]); //pointer check
        HSievePlane[idx]->Draw("COLZ");

        // Draw Sieve
        const Double_t plotwidth = 0.004;
        for (UInt_t Row = 0; Row < NSieveRow; Row++) {
            for (UInt_t Col = 0; Col < NSieveCol; Col++) {
                const Double_t posx = SieveOffY + SieveYbyCol[Col];
                const Double_t posy = SieveOffX + SieveXbyRow[Row];

                TLine *lh = new TLine(posx - plotwidth, posy, posx + plotwidth, posy);
                TLine *lv = new TLine(posx, posy - plotwidth, posx, posy + plotwidth);
                lh -> SetLineColor(kBlack);
                lv -> SetLineColor(kBlack);

                lh -> Draw();
                lv -> Draw();
            }
        }

        // Draw arrows
        for (UInt_t Col = 0; Col < NSieveCol; Col++) {
            for (UInt_t Row = 0; Row < NSieveRow; Row++) {
                if (SieveEventID[KineID][Col][Row][kEventID] > 0) {
                    assert(SieveEventID[KineID][Col][Row][kEventID] < fNRawData); //array index bondary check
                    TArrow * ar2 = new TArrow(SieveEventID[KineID][Col][Row][kCalcSieveY], SieveEventID[KineID][Col][Row][kCalcSieveX], SieveEventID[KineID][Col][Row][kRealSieveY], SieveEventID[KineID][Col][Row][kRealSieveX], 0.008, "|>");
                    ar2->SetAngle(40);
                    ar2->SetLineColor(kMagenta);
                    ar2->SetFillColor(kMagenta);

                    const Double_t ignorelimit = 0.005;
                    if ((ar2->GetX1() - ar2->GetX2())*(ar2->GetX1() - ar2->GetX2())+(ar2->GetY1() - ar2->GetY2())*(ar2->GetY1() - ar2->GetY2()) > ignorelimit * ignorelimit)
                        ar2->Draw();
                }
            }
        }
    }

    return c1;
}

void LOpticsOpt::PrepareDp(void)
{
    // calate expected dp_kin, dp_kin offsets ....
    // Fill up fRawData[].Data[] kKineID thr kRealDpKin

    // print Central Momentums
    printf("HRSCentralMom[%d] (GeV) = {", NKine);
    for (UInt_t KineID = 0; KineID < NKine; KineID++)
        printf("%f  ", HRSCentralMom[KineID]);
    printf("}\n");

    // print radiation loss numbers
    printf("RadiationLossByFoil[%d] (MeV) = {", NFoils);
    for (UInt_t FoilID = 0; FoilID < NFoils; FoilID++)
        printf("%f  ", RadiationLossByFoil[FoilID]*1000);
    printf("}\n");

    Double_t scatang = 0;
    Double_t dpkinoff = 0, dpkinoff_rms = 0;
    fNCalibData = 0;
    Double_t exttargcorr_dp = 0, rms_exttargcorr_dp = 0;

    for (UInt_t idx = 0; idx < fNRawData; idx++) {
        DEBUG_MASSINFO("PrepareDp", "=========== Event %d ===========", idx);

        EventData &eventdata = fRawData[idx];

        // decoding kCutID
        UInt_t res = (UInt_t) eventdata.Data[kCutID];

        const UInt_t ExtraDataFlag = res / (NSieveRow * NSieveCol * NFoils * NKine);
        res = res % (NSieveRow * NSieveCol * NFoils * NKine);
        const UInt_t KineID = res / (NSieveRow * NSieveCol * NFoils); //starting 0!
        res = res % (NSieveRow * NSieveCol * NFoils);
        const UInt_t FoilID = res / (NSieveRow * NSieveCol); //starting 0!
        res = res % (NSieveRow * NSieveCol);
        const UInt_t Col = res / (NSieveRow); //starting 0!
        const UInt_t Row = res % (NSieveRow); //starting 0!
        assert(ExtraDataFlag < 2); //Check flag range. beyond 2 is not used
        assert(KineID < NKine); //check array index size
        assert(FoilID < NFoils); //check array index size

        DEBUG_MASSINFO("PrepareDp", "%d => KineID=%d,\tFoilID=%d,\tCol=%d,\tRow=%d", (UInt_t) eventdata.Data[kCutID], KineID, FoilID, Col, Row);

        // write some variables
        eventdata.Data[kExtraDataFlag] = ExtraDataFlag;
        if (!ExtraDataFlag) fNCalibData++;
        eventdata.Data[kKineID] = KineID;
        eventdata.Data[kCentralp] = HRSCentralMom[KineID];

        const TVector3 SieveHoleTCS = GetSieveHoleTCS(Col, Row);

        TVector3 BeamSpotHCS(eventdata.Data[kBeamX], eventdata.Data[kBeamY], targetfoils[FoilID]);

        // Calculate Real Scattering Angles by Sieve Holes cuts
        TVector3 BeamSpotTCS = fTCSInHCS.Inverse()*(BeamSpotHCS - fPointingOffset);

        TVector3 MomDirectionTCS = SieveHoleTCS - BeamSpotTCS;

        eventdata.Data[kRealTh] = MomDirectionTCS.X() / MomDirectionTCS.Z();
        eventdata.Data[kRealPhi] = MomDirectionTCS.Y() / MomDirectionTCS.Z();

        const Double_t x_tg = BeamSpotTCS.X() - BeamSpotTCS.Z() * eventdata.Data[kRealTh];
        eventdata.Data[kRealTgX] = x_tg;
        eventdata.Data[kRealThMatrix] = eventdata.Data[kRealTh] - x_tg * ExtTarCor_ThetaCorr;

        DEBUG_MASSINFO("PrepareDp", "RealTh = %f,\t RealPh = %f", eventdata.Data[kRealThMatrix], eventdata.Data[kRealPhi]);
        DEBUG_MASSINFO("PrepareDp", "SieveHoleY = %f,\t Mom_Y = %f,\t Mom_Z = %f", SieveHoleTCS.Y(), MomDirectionTCS.Y(), MomDirectionTCS.Z());

        TVector3 MomDirectionHCS = fTCSInHCS*MomDirectionTCS;
        TVector3 BeamDirection(0, 0, 1);
        const Double_t ScatteringAngle = BeamDirection.Angle(MomDirectionHCS);
        eventdata.Data[kScatterAngle] = ScatteringAngle;
        scatang += ScatteringAngle;

        // calculate difference between dp_kin and dp
        // dp_kin + kDpKinOffsets = dp
        const Double_t DM = ExcitationEnergy[KineID];
        const Double_t Ma = GroundNuclearMass;
        const Double_t P0 = eventdata.Data[kurb_e];
        const Double_t DpKinOffsets = (ScatMom(DM, Ma, P0, ScatteringAngle) - ScatMom(DM, Ma, P0, TMath::Abs(HRSAngle))) / eventdata.Data[kCentralp];
        eventdata.Data[kDpKinOffsets] = DpKinOffsets;

        dpkinoff += DpKinOffsets;
        dpkinoff_rms += DpKinOffsets*DpKinOffsets;

        // calculate kRealDpKin, should be same for same kine settings
        eventdata.Data[kRadiLossDp] = RadiationLossByFoil[FoilID] / eventdata.Data[kCentralp];
        eventdata.Data[kRealDpKin] = ScatMom(DM, Ma, P0, TMath::Abs(HRSAngle)) / eventdata.Data[kCentralp] - 1 - eventdata.Data[kRadiLossDp];

        DEBUG_MASSINFO("PrepareDp", "ScatterMom = %f,\t Centralp = %f,\t radloss = %f,\t ebeam = %f", ScatMom(DM, Ma, P0, TMath::Abs(HRSAngle)), eventdata.Data[kCentralp], eventdata.Data[kRadiLossDp], P0);

        eventdata.Data[kRealDpKinMatrix] = eventdata.Data[kRealDpKin] - x_tg / ExtTarCor_DeltaCorr;

        exttargcorr_dp += x_tg / ExtTarCor_DeltaCorr;
        rms_exttargcorr_dp += (x_tg / ExtTarCor_DeltaCorr)*(x_tg / ExtTarCor_DeltaCorr);

        // calculate expected dp_kin for all other excitation states
        for (UInt_t ExcitID = 0; ExcitID < NExcitationStates; ExcitID++) {
            assert(kRealDpKinExcitations + ExcitID < kRealTh); //check array index size
            eventdata.Data[kRealDpKinExcitations + ExcitID] = ScatMom(ExcitationEnergyList[ExcitID], Ma, P0, TMath::Abs(HRSAngle)) / eventdata.Data[kCentralp] - 1;
        }

        DEBUG_MASSINFO("PrepareDp", "ScatterAngle = %f,\t DpKinOffsets = %f", ScatteringAngle / TMath::Pi()*180, DpKinOffsets);
    }

    DEBUG_INFO("PrepareDp", "%d out of %d data is for calibration", fNCalibData, fNRawData);
    DEBUG_INFO("PrepareDp", "Average : ScatteringAngle = %f", scatang / fNRawData / TMath::Pi()*180);
    DEBUG_INFO("PrepareDp", "Average DpKinOffsets = %f, RMS DpKinOffsets = %f", dpkinoff / fNRawData, TMath::Sqrt(dpkinoff_rms / fNRawData));
    DEBUG_INFO("PrepareDp", "Average Extended Target Corretion: dp = %f,\t rms_dp = %f", exttargcorr_dp / fNRawData, TMath::Sqrt(rms_exttargcorr_dp / fNRawData));

    // make sure kCalcTh, kCalcPh is filled, although not necessary
    SumSquareDTh();
    SumSquareDPhi();
    SumSquareDp();

    DEBUG_INFO("PrepareDp", "Done!");
}

void LOpticsOpt::PrepareDpWithField(void)
{
    // calate expected dp_kin, dp_kin offsets ....
    // Fill up fRawData[].Data[] kKineID thr kRealDpKin

    // print Central Momentums
    printf("HRSCentralMom[%d] (GeV) = {", NKine);
    for (UInt_t KineID = 0; KineID < NKine; KineID++)
        printf("%f  ", HRSCentralMom[KineID]);
    printf("}\n");

    // print radiation loss numbers
    printf("RadiationLossByFoil[%d] (MeV) = {", NFoils);
    for (UInt_t FoilID = 0; FoilID < NFoils; FoilID++)
        printf("%f  ", RadiationLossByFoil[FoilID]*1000);
    printf("}\n");

    Double_t scatang = 0;
    Double_t dpkinoff = 0, dpkinoff_rms = 0;
    fNCalibData = 0;
    Double_t exttargcorr_dp = 0, rms_exttargcorr_dp = 0;

    for (UInt_t idx = 0; idx < fNRawData; idx++) {
        DEBUG_MASSINFO("PrepareDp", "=========== Event %d ===========", idx);

        EventData &eventdata = fRawData[idx];

        // decoding kCutID
        UInt_t res = (UInt_t) eventdata.Data[kCutID];

        const UInt_t ExtraDataFlag = res / (NSieveRow * NSieveCol * NFoils * NKine);
        res = res % (NSieveRow * NSieveCol * NFoils * NKine);
        const UInt_t KineID = res / (NSieveRow * NSieveCol * NFoils); //starting 0!
        res = res % (NSieveRow * NSieveCol * NFoils);
        const UInt_t FoilID = res / (NSieveRow * NSieveCol); //starting 0!
        res = res % (NSieveRow * NSieveCol);
        assert(ExtraDataFlag < 2); //Check flag range. beyond 2 is not used
        assert(KineID < NKine); //check array index size
        assert(FoilID < NFoils); //check array index size

        DEBUG_MASSINFO("PrepareDp", "%d => KineID=%d,\tFoilID=%d", (UInt_t) eventdata.Data[kCutID], KineID, FoilID);

        // write some variables
        eventdata.Data[kExtraDataFlag] = ExtraDataFlag;
        if (!ExtraDataFlag) fNCalibData++;
        eventdata.Data[kKineID] = KineID;
        eventdata.Data[kCentralp] = HRSCentralMom[KineID];

        eventdata.Data[kRealTh] = eventdata.Data[kOrinTh];
        eventdata.Data[kRealPhi] = eventdata.Data[kOrinPh];

        const Double_t x_tg = eventdata.Data[kSimX];
        eventdata.Data[kRealTgX] = x_tg;
        eventdata.Data[kRealThMatrix] = eventdata.Data[kRealTh] - x_tg * ExtTarCor_ThetaCorr;

        DEBUG_MASSINFO("PrepareDp", "RealTh = %f,\t RealPh = %f", eventdata.Data[kRealThMatrix], eventdata.Data[kRealPhi]);

        TVector3 MomDirectionTCS(eventdata.Data[kRealTh], eventdata.Data[kRealPhi], 1);

        TVector3 MomDirectionHCS = fTCSInHCS*MomDirectionTCS;
        TVector3 BeamDirection(0, 0, 1);
        const Double_t ScatteringAngle = BeamDirection.Angle(MomDirectionHCS);
        eventdata.Data[kScatterAngle] = ScatteringAngle;
        scatang += ScatteringAngle;

        // calculate difference between dp_kin and dp
        // dp_kin + kDpKinOffsets = dp
        const Double_t DM = ExcitationEnergy[KineID];
        const Double_t Ma = GroundNuclearMass;
        const Double_t P0 = eventdata.Data[kurb_e];
        const Double_t DpKinOffsets = (ScatMom(DM, Ma, P0, ScatteringAngle) - ScatMom(DM, Ma, P0, TMath::Abs(HRSAngle))) / eventdata.Data[kCentralp];
        eventdata.Data[kDpKinOffsets] = DpKinOffsets;

        dpkinoff += DpKinOffsets;
        dpkinoff_rms += DpKinOffsets*DpKinOffsets;

        // calculate kRealDpKin, should be same for same kine settings
        eventdata.Data[kRadiLossDp] = RadiationLossByFoil[FoilID] / eventdata.Data[kCentralp];
        eventdata.Data[kRealDpKin] = ScatMom(DM, Ma, P0, TMath::Abs(HRSAngle)) / eventdata.Data[kCentralp] - 1 - eventdata.Data[kRadiLossDp];

        DEBUG_MASSINFO("PrepareDp", "ScatterMom = %f,\t Centralp = %f,\t radloss = %f,\t ebeam = %f", ScatMom(DM, Ma, P0, TMath::Abs(HRSAngle)), eventdata.Data[kCentralp], eventdata.Data[kRadiLossDp], P0);

        eventdata.Data[kRealDpKinMatrix] = eventdata.Data[kRealDpKin] - x_tg / ExtTarCor_DeltaCorr;

        exttargcorr_dp += x_tg / ExtTarCor_DeltaCorr;
        rms_exttargcorr_dp += (x_tg / ExtTarCor_DeltaCorr)*(x_tg / ExtTarCor_DeltaCorr);

        // calculate expected dp_kin for all other excitation states
        for (UInt_t ExcitID = 0; ExcitID < NExcitationStates; ExcitID++) {
            assert(kRealDpKinExcitations + ExcitID < kRealTh); //check array index size
            eventdata.Data[kRealDpKinExcitations + ExcitID] = ScatMom(ExcitationEnergyList[ExcitID], Ma, P0, TMath::Abs(HRSAngle)) / eventdata.Data[kCentralp] - 1;
        }

        DEBUG_MASSINFO("PrepareDp", "ScatterAngle = %f,\t DpKinOffsets = %f", ScatteringAngle / TMath::Pi()*180, DpKinOffsets);
    }

    DEBUG_INFO("PrepareDp", "%d out of %d data is for calibration", fNCalibData, fNRawData);
    DEBUG_INFO("PrepareDp", "Average : ScatteringAngle = %f", scatang / fNRawData / TMath::Pi()*180);
    DEBUG_INFO("PrepareDp", "Average DpKinOffsets = %f, RMS DpKinOffsets = %f", dpkinoff / fNRawData, TMath::Sqrt(dpkinoff_rms / fNRawData));
    DEBUG_INFO("PrepareDp", "Average Extended Target Corretion: dp = %f,\t rms_dp = %f", exttargcorr_dp / fNRawData, TMath::Sqrt(rms_exttargcorr_dp / fNRawData));

    // make sure kCalcTh, kCalcPh is filled, although not necessary
    SumSquareDTh();
    SumSquareDPhi();
    SumSquareDp();

    DEBUG_INFO("PrepareDp", "Done!");
}

TCanvas * LOpticsOpt::CheckDp()
{
    // Visualize 1D hitogram of dp_kin

    DEBUG_INFO("CheckDp", "Entry Point");

    // calculate Data[kCalcDpKin] for all events
    SumSquareDp(kTRUE);

    const Double_t DpRange = .05;
    const UInt_t NDpRange = 5000;

    TH1D * hDpKinCalib[NKine];
    TH1D * hDpKinAll[NKine];
    Double_t RealDpKin[NKine] = {0};
    Double_t AverCalcDpKin[NKine] = {0};
    UInt_t NEvntDpKin[NKine] = {0};
    Double_t RealDpKinAllExcit[NExcitationStates][NKine] = {
        {0}
    };
    Double_t NewArbitaryDpKinShift[NKine];

    for (UInt_t KineID = 0; KineID < NKine; KineID++) {
        hDpKinCalib[KineID] = new TH1D(Form("hDpKinCalib%d", KineID), Form("Dp_Kin for Delta Scan Kine. #%d (Selected Exct. State)", KineID), NDpRange, -DpRange, DpRange);
        hDpKinAll[KineID] = new TH1D(Form("hDpKinAll%d", KineID), Form("Dp_Kin for Delta Scan Kine. #%d (All Data)", KineID), NDpRange, -DpRange, DpRange);

        assert(hDpKinCalib[KineID]); //pointer check
        assert(hDpKinAll[KineID]); //pointer check
    }

    for (UInt_t idx = 0; idx < fNRawData; idx++) {
        const EventData &eventdata = fRawData[idx];
        const UInt_t ExtraDataFlag = (UInt_t) (eventdata.Data[kExtraDataFlag]);
        assert(ExtraDataFlag == 0 || ExtraDataFlag == 1); //flag definition consistency check
        const UInt_t KineID = (UInt_t) (eventdata.Data[kKineID]);
        assert(KineID < NKine);

        if (!ExtraDataFlag) {
            hDpKinCalib[KineID]->Fill(eventdata.Data[kCalcDpKin] + eventdata.Data[kRadiLossDp]);
            AverCalcDpKin[KineID] += eventdata.Data[kCalcDpKin] + eventdata.Data[kRadiLossDp];
            NEvntDpKin[KineID]++;
        }
        hDpKinAll[KineID]->Fill(eventdata.Data[kCalcDpKin] + eventdata.Data[kRadiLossDp]);

        RealDpKin[KineID] = eventdata.Data[kRealDpKin] + eventdata.Data[kRadiLossDp];

        for (UInt_t ExcitID = 0; ExcitID < NExcitationStates; ExcitID++) {
            assert(kRealDpKinExcitations + ExcitID < kRealTh); //index check
            RealDpKinAllExcit[ExcitID][KineID] = eventdata.Data[kRealDpKinExcitations + ExcitID];
        }
    }

    TCanvas * c1 = new TCanvas("CheckDp", "Check Dp Kin Reconstruction", 1800, 900);

    if (NKine <= 3)
        c1->Divide(3, 1);
    else if (NKine <= 6)
        c1->Divide(3, 2);
    else
        c1->Divide(3, 3);

    UInt_t idx = 1;

    for (UInt_t KineID = 0; KineID < NKine; KineID++) {
        c1->cd(idx++);
        gPad -> SetLogy();

        AverCalcDpKin[KineID] /= NEvntDpKin[KineID];
        DEBUG_MASSINFO("CheckDp", "AverCalcDpKin[%d] = %f", KineID, AverCalcDpKin[KineID]);

        // Histograms
        hDpKinCalib[KineID]->SetLineColor(4);
        hDpKinCalib[KineID]->SetFillColor(4);
        hDpKinCalib[KineID]->SetFillStyle(3008);

        hDpKinAll[KineID]->SetLineColor(1);
        // hDpKinAll[KineID]->SetFillColor(1);

        const Double_t dpRange = 0.01;
        hDpKinCalib[KineID]-> SetAxisRange(RealDpKin[KineID] - dpRange, RealDpKin[KineID] + dpRange);
        hDpKinAll[KineID]->SetAxisRange(RealDpKin[KineID] - dpRange, RealDpKin[KineID] + dpRange);

        hDpKinCalib[KineID]->SetXTitle("radiation corrected dp_kin (angular independent dp)");
        hDpKinAll[KineID]->SetXTitle("radiation corrected dp_kin (angular independent dp)");

        hDpKinAll[KineID]->Draw();
        hDpKinCalib[KineID]->Draw("SAME");

        // expectation lines
        const Double_t MaxPlot = 20000;
        for (UInt_t ExcitID = 0; ExcitID < NExcitationStates; ExcitID++) {
            const Double_t x = RealDpKinAllExcit[ExcitID][KineID];
            TLine *l = new TLine(x, 0, x, +MaxPlot);
            l->SetLineColor(3);
            l->SetLineWidth(2);
            l->Draw();
        }
        TLine *l = new TLine(RealDpKin[KineID], 0, RealDpKin[KineID], +MaxPlot);
        l->SetLineColor(6);
        l->SetLineWidth(2);
        l->Draw();

        // Fits
        const Double_t DefResolution = 1e-4;
        const Double_t FitRangeMultiply = 5;

        TString FitFunc = Form("DpPeak%d", KineID);
        TF1 *f = new TF1(FitFunc, "gaus+[3]+[4]*x", AverCalcDpKin[KineID] - DefResolution*FitRangeMultiply, AverCalcDpKin[KineID] + DefResolution * FitRangeMultiply);
        f->SetParameter(1, AverCalcDpKin[KineID]);
        f->SetParameter(2, DefResolution);
        hDpKinAll[KineID] -> Fit(FitFunc, "RN0");
        // Info("CheckDp","Fit for delta scan #%d peak:",KineID);
        // f->Print();
        f->SetLineColor(2);
        f->Draw("SAME");

        TLatex *t = new TLatex(f->GetParameter(1) + DefResolution, f->GetParameter(0) + f->GetParameter(3) + f->GetParameter(4) * f->GetParameter(1), Form("\\Delta \\pm \\sigma = (%2.1f \\pm %2.1f) \\times 10^{-4}", 10000 * (f->GetParameter(1) - RealDpKin[KineID]), 10000 * f->GetParameter(2)));
        t->SetTextSize(0.05);
        t->SetTextAlign(12);
        t->SetTextColor(2);
        t->Draw();

        NewArbitaryDpKinShift[KineID] = f->GetParameter(1) - RealDpKin[KineID] + fArbitaryDpKinShift[KineID];
    }

    Info("CheckDp", "New set of arbitary dp shifts:");
    for (UInt_t KineID = 0; KineID < NKine; KineID++)
        printf("opt->fArbitaryDpKinShift[%d] = %e;\n", KineID, NewArbitaryDpKinShift[KineID]);

    return c1;
}

TCanvas * LOpticsOpt::CheckDpGlobal()
{
    // Visualize 1D hitogram of dp_kin

    DEBUG_INFO("CheckDp", "Entry Point");

    SumSquareDp(kTRUE);

    const Double_t DpRange = .05;
    const UInt_t NDpRange = 5000;

    TH1D * hDpKinCalib[NKine];
    TH1D * hDpKinAll[NKine];
    Double_t RealDpKin[NKine];
    Double_t AverCalcDpKin[NKine] = {0};
    UInt_t NEvntDpKin[NKine] = {0};
    Double_t RealDpKinAllExcit[NExcitationStates][NKine] = {
        {0}
    };

    for (UInt_t KineID = 0; KineID < NKine; KineID++) {
        hDpKinCalib[KineID] = new TH1D(Form("hDpKinCalib%d", KineID), "Dp_Kin for Delta Scan", NDpRange, -DpRange, DpRange);
        hDpKinAll[KineID] = new TH1D(Form("hDpKinAll%d", KineID), Form("Dp_Kin for Delta Scan Kine. #%d", KineID), NDpRange, -DpRange, DpRange);

        assert(hDpKinCalib[KineID]); //pointer check
        assert(hDpKinAll[KineID]); //pointer check
    }

    for (UInt_t idx = 0; idx < fNRawData; idx++) {
        const EventData &eventdata = fRawData[idx];
        const UInt_t ExtraDataFlag = (UInt_t) (eventdata.Data[kExtraDataFlag]);
        assert(ExtraDataFlag == 0 || ExtraDataFlag == 1); //flag definition consistency check
        const UInt_t KineID = (UInt_t) (eventdata.Data[kKineID]);
        assert(KineID < NKine);

        if (!ExtraDataFlag) {
            hDpKinCalib[KineID]->Fill(eventdata.Data[kCalcDpKin] + eventdata.Data[kRadiLossDp]);
            AverCalcDpKin[KineID] += eventdata.Data[kCalcDpKin] + eventdata.Data[kRadiLossDp];
            NEvntDpKin[KineID]++;
        }
        hDpKinAll[KineID]->Fill(eventdata.Data[kCalcDpKin] + eventdata.Data[kRadiLossDp]);

        RealDpKin[KineID] = eventdata.Data[kRealDpKin] + eventdata.Data[kRadiLossDp];

        for (UInt_t ExcitID = 0; ExcitID < NExcitationStates; ExcitID++) {
            assert(kRealDpKinExcitations + ExcitID < kRealTh); //index check
            RealDpKinAllExcit[ExcitID][KineID] = eventdata.Data[kRealDpKinExcitations + ExcitID];
        }
    }

    TCanvas * c1 = new TCanvas("CheckDp", "Check Dp Kin Reconstruction", 1800, 900);

    for (UInt_t KineID = 0; KineID < NKine; KineID++) {
        c1->cd();
        // gPad -> SetLogy();

        AverCalcDpKin[KineID] /= NEvntDpKin[KineID];
        DEBUG_MASSINFO("CheckDp", "AverCalcDpKin[%d] = %f", KineID, AverCalcDpKin[KineID]);

        // Histograms
        hDpKinCalib[KineID]->SetLineColor(4);
        hDpKinCalib[KineID]->SetFillColor(4);
        hDpKinCalib[KineID]->SetFillStyle(3008);

        hDpKinAll[KineID]->SetLineColor(1);
        // hDpKinAll[KineID]->SetFillColor(1);

        if (KineID == 0) {
            const Double_t dpRange = 0.04;
            hDpKinCalib[KineID]->SetAxisRange(-dpRange, +dpRange);
            hDpKinAll[KineID]->SetAxisRange(RealDpKin[KineID] - dpRange, RealDpKin[KineID] + dpRange);

            hDpKinCalib[KineID]->SetXTitle("dp_kin (energy loss corrected)");

            hDpKinCalib[KineID]->Draw();
        } else {
            hDpKinCalib[KineID]->Draw("same");
        }

        // expectation lines
        const Double_t MaxPlot = 20000;
        for (UInt_t ExcitID = 0; ExcitID < NExcitationStates; ExcitID++) {
            const Double_t x = RealDpKinAllExcit[ExcitID][KineID];
            TLine *l = new TLine(x, 0, x, +MaxPlot);
            l->SetLineColor(3);
            l->SetLineWidth(2);
            l->Draw();
        }
        TLine *l = new TLine(RealDpKin[KineID], 0, RealDpKin[KineID], +MaxPlot);
        l->SetLineColor(6);
        l->SetLineWidth(2);
        l->Draw();
    }

    return c1;
}

Double_t LOpticsOpt::SumSquareDp(Bool_t IncludeExtraData)
{
    // return square sum of diff between calculated dp_kin and expected dp_kin

    Double_t d_dp = 0; //Difference
    Double_t rms_dp = 0; //mean square

    static UInt_t NCall = 0;
    NCall++;

    UInt_t NCalibData = 0;

    if (IncludeExtraData) {
        Warning("SumSquareDp", "Data Beyond selected excitation state is included in this calculation");
    }
    for (UInt_t idx = 0; idx < fNRawData; idx++) {
        Double_t dp, dp_kin;

        EventData &eventdata = fRawData[idx];

        //jump through data beyond selected excitation states
        if (eventdata.Data[kExtraDataFlag] > 0 && !IncludeExtraData) continue;
        NCalibData++;

        Double_t x_fp = eventdata.Data[kX];
        const Double_t(*powers)[5] = eventdata.powers;

        // calculate the matrices we need
        CalcMatrix(x_fp, fDMatrixElems);
        // CalcMatrix(x_fp, fTMatrixElems);
        // CalcMatrix(x_fp, fYMatrixElems);
        // CalcMatrix(x_fp, fYTAMatrixElems);
        // CalcMatrix(x_fp, fPMatrixElems);
        // CalcMatrix(x_fp, fPTAMatrixElems);

        // calculate momentum
        dp = CalcTargetVar(fDMatrixElems, powers);
        dp_kin = dp - eventdata.Data[kDpKinOffsets];

        const UInt_t KineID = (UInt_t) (eventdata.Data[kKineID]);
        assert(KineID < NKine); //check array index size
        const Double_t ArbitaryDpKinShift = fArbitaryDpKinShift[KineID];

        d_dp += dp_kin - eventdata.Data[kRealDpKinMatrix] + ArbitaryDpKinShift;
        rms_dp += (dp_kin - eventdata.Data[kRealDpKinMatrix] + ArbitaryDpKinShift)*(dp_kin - eventdata.Data[kRealDpKinMatrix] + ArbitaryDpKinShift);

        DEBUG_MASSINFO("SumSquareDp", "d_dp = %f = \t%f - \t%f", dp_kin - eventdata.Data[kRealDpKinMatrix], dp_kin, eventdata.Data[kRealDpKinMatrix]);

        // save the results
        eventdata.Data[kCalcDpKinMatrix] = dp_kin;
        eventdata.Data[kCalcDpKin] = dp_kin + eventdata.Data[kRealTgX] / ExtTarCor_DeltaCorr;
    }

    if (!IncludeExtraData)
        assert(fNCalibData == NCalibData); // check number of event for calibration

    DEBUG_INFO("SumSquareDp", "#%d : d_dp = %f, rms_dp = %f", NCall, d_dp / NCalibData, TMath::Sqrt(rms_dp / NCalibData));
    //printf("SumSquareDp: #%d : d_dp = %f, rms_dp = %f\n", NCall, d_dp / NCalibData, TMath::Sqrt(rms_dp / NCalibData));

    return rms_dp;
}

///////////////////////////////////////////////////////////////////////////////
// declarations for target vertex reconstruction
///////////////////////////////////////////////////////////////////////////////

void LOpticsOpt::CalcMatrix(const Double_t x, vector<THaMatrixElement>& matrix)
{
    // calculates the values of the matrix elements for a given location
    // by evaluating a polynomial in x of order it->order with
    // coefficients given by it->poly

    for (vector<THaMatrixElement>::iterator it = matrix.begin();
            it != matrix.end(); it++) {
        it->v = 0.0;

        if (it->order > 0) {
            for (int i = it->order - 1; i >= 1; i--)
                it->v = x * (it->v + it->poly[i]);
            it->v += it->poly[0];
        }
    }
}

Double_t LOpticsOpt::CalcTargetVar(const vector<THaMatrixElement>& matrix, const Double_t powers[][5])
{
    // calculates the value of a variable at the target
    // the x-dependence is already in the matrix, so only 1-3 (or np) used
    Double_t retval = 0.0;
    Double_t v = 0;
    for (vector<THaMatrixElement>::const_iterator it = matrix.begin();
            it != matrix.end(); it++)
        if (it->v != 0.0) {
            v = it->v;
            unsigned int np = it->pw.size(); // generalize for extra matrix elems.
            for (unsigned int i = 0; i < np; i++)
                v *= powers[it->pw[i]][i + 1];
            retval += v;
            //      retval += it->v * powers[it->pw[0]][1]
            //	              * powers[it->pw[1]][2]
            //	              * powers[it->pw[2]][3];
        }

    return retval;
}

///////////////////////////////////////////////////////////////////////////////
// class for storing matrix element data
///////////////////////////////////////////////////////////////////////////////

bool THaMatrixElement::match(const THaMatrixElement& rhs) const
{
    // Compare coefficients of this matrix element to another

    if (pw.size() != rhs.pw.size())
        return false;
    for (vector<int>::size_type i = 0; i < pw.size(); i++) {
        if (pw[i] != rhs.pw[i])
            return false;
    }
    return true;
}

void THaMatrixElement::SkimPoly()
{
    // reduce order to highest non-zero poly

    if (iszero) return;

    while (!poly[order - 1] && order > 0) {
        poly.pop_back();
        order = order - 1;
    }

    if (order == 0) iszero = kTRUE;
}

ClassImp(LOpticsOpt);
