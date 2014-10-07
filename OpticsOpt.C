///////////////////////////////////////////////////////////////////////////////
//
// OpticsOpt
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

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <map>
#include <vector>

#include "TROOT.h"
#include "TError.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TDatime.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLine.h"
#include "TList.h"
#include "TMath.h"
#include "TRotation.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

#include "THaAnalysisObject.h"
#include "THaEvData.h"
#include "THaGlobals.h"
#include "THaTrack.h"
#include "THaScintillator.h"
#include "THaSpectrometer.h"
#include "THaString.h"
#include "THaTrackingDetector.h"
#include "VarDef.h"

#include "OpticsOpt.h"

#ifdef WITH_DEBUG
#include <iostream>
#endif

using namespace std;
using THaString::Split;

///////////////////////////////////////////////////////////////////////////////
// Input Sections
///////////////////////////////////////////////////////////////////////////////

#include "InputG2P.h"

///////////////////////////////////////////////////////////////////////////////
// Constructors
///////////////////////////////////////////////////////////////////////////////

OpticsOpt::OpticsOpt(const Char_t* name, const Char_t* description, THaApparatus* apparatus) : THaTrackingDetector(name, description, apparatus)
{
    fPrefix = new Char_t[1000];
    sprintf(fPrefix, "%s", Prefix);

    fCurrentMatrixElems = NULL;

    TVector3 TCSX(0, -1, 0);
    TVector3 TCSZ(TMath::Sin(HRSAngle), 0, TMath::Cos(HRSAngle));
    TVector3 TCSY = TCSZ.Cross(TCSX);
    fTCSInHCS.RotateAxes(TCSX, TCSY, TCSZ);

    fPointingOffset.SetXYZ(-MissPointZ * TMath::Sin(HRSAngle) * TMath::Cos(HRSAngle), MissPointY, MissPointZ * TMath::Sin(HRSAngle) * TMath::Sin(HRSAngle));
    DEBUG_INFO("OpticsOpt", "Read in configuration " + InputID);
    DEBUG_INFO("OpticsOpt", "HRS @ %f Degree, PointingOffset = (%f,%f,%f), SievePos = (%f,%f,%f)", HRSAngle / TMath::Pi()*180, fPointingOffset.X(), fPointingOffset.Y(), fPointingOffset.Z(), SieveOffX, SieveOffY, ZPos);

    fNRawData = 0;

    for (Int_t i = 0; i < 100; i++) fArbitaryDpKinShift[i] = 0;
}

OpticsOpt::~OpticsOpt()
{
    // Nothing to do
}

///////////////////////////////////////////////////////////////////////////////
// Database I/O
///////////////////////////////////////////////////////////////////////////////

Int_t OpticsOpt::LoadDataBase(TString DataBaseName)
{
    static const Char_t * const here = "LoadDataBase";
    OldComments = "";

    FILE* file = fopen(DataBaseName, "r");
    if (!file) {
        Error("LoadDataBase", "%s can not be opened", DataBaseName.Data());
        assert(0); //
        return kFileError;
    } else DEBUG_INFO("LoadDataBase", "Parsing Database %s", DataBaseName.Data());

    const Int_t LEN = 200;
    Char_t buff[LEN];

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

    Bool_t found = false;
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

    map<string, Int_t> fp_map;
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
        const Char_t* w = line_spl[0].c_str();
        vsiz_t npow = power[w];
        if (npow == 0) break;

#if DEBUG_LEVEL>=4
        cout << "Matrix Line = ";
        for (Ssiz_t i = 1; i < line_spl.size(); i++) {
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
        vector<THaMatrixElement>* mat = matrix_map[w];
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
                Bool_t match = false;
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

Int_t OpticsOpt::SaveDataBase(TString DataBaseName)
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
    fprintf(file, "# Optimized by %s @ %s", User.Data(), dt.AsString());
    fprintf(file, "\n");
    fprintf(file, "# Saved to %s", DataBaseName.Data());
    fprintf(file, "\n");
    fprintf(file, OldComments);

    fprintf(file, DatabaseHeader);

    DEBUG_INFO("SaveDataBase", "Transport Matrix:  FP-terms");
    for (vsiz_t i = 0; i < fFPMatrixElems.size(); i++) {
        switch (i) {
        case 0:
            fprintf(file, "t ");
            break;
        case 1:
            fprintf(file, "y ");
            break;
        case 2:
            fprintf(file, "p ");
            break;
        default:
            break;
        }
        const THaMatrixElement& m = fFPMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            fprintf(file, "%d ", m.pw[j]);
        }
        Int_t j;
        for (j = 0; j < m.order; j++) {
            fprintf(file, " %13.6e", m.poly[j]);
        }
        for (; j < kPORDER; j++) {
            fprintf(file, " %13.6e", 0.0);
        }
        fprintf(file, "  %d", m.OptOrder);
        fprintf(file, "\n");
    }

    DEBUG_INFO("SaveDataBase", "Transport Matrix:  D-terms");
    for (vsiz_t i = 0; i < fDMatrixElems.size(); i++) {
        fprintf(file, "D ");
        const THaMatrixElement& m = fDMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            fprintf(file, "%d ", m.pw[j]);
        }
        Int_t j;
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
        Int_t j;
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
        Int_t j;
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
        Int_t j;
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

void OpticsOpt::Print(const Option_t* opt) const
{
    // Print current matrix

    THaTrackingDetector::Print(opt);
    typedef vector<THaMatrixElement>::size_type vsiz_t;

    // Print out the optics matrices, to verify they make sense
    printf("OpticsOpt::Print: Matrix FP (t000, y000, p000)\n");
    for (vsiz_t i = 0; i < fFPMatrixElems.size(); i++) {
        const THaMatrixElement& m = fFPMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            printf("  %2d", m.pw[j]);
        }
        for (Int_t j = 0; j < m.order; j++) {
            printf("  %g", m.poly[j]);
        }
        printf(" : Opt -> %d", m.OptOrder);
        printf("\n");
    }

    printf("OpticsOpt::Print: Transport Matrix:  D-terms\n");
    for (vsiz_t i = 0; i < fDMatrixElems.size(); i++) {
        const THaMatrixElement& m = fDMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            printf("  %2d", m.pw[j]);
        }
        for (Int_t j = 0; j < m.order; j++) {
            printf("\t%g", m.poly[j]);
        }
        printf(" : Opt -> %d", m.OptOrder);
        if (m.order != m.OptOrder) printf(" != Matrix Order !!");
        printf("\n");
    }

    printf("OpticsOpt::Print: Transport Matrix:  T-terms\n");
    for (vsiz_t i = 0; i < fTMatrixElems.size(); i++) {
        const THaMatrixElement& m = fTMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            printf("  %2d", m.pw[j]);
        }
        for (Int_t j = 0; j < m.order; j++) {
            printf("\t%g", m.poly[j]);
        }
        printf(" : Opt -> %d", m.OptOrder);
        if (m.order != m.OptOrder) printf(" != Matrix Order !!");
        printf("\n");
    }

    printf("OpticsOpt::Print: Transport Matrix:  Y-terms\n");
    for (vsiz_t i = 0; i < fYMatrixElems.size(); i++) {
        const THaMatrixElement& m = fYMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            printf("  %2d", m.pw[j]);
        }
        for (Int_t j = 0; j < m.order; j++) {
            printf("\t%g", m.poly[j]);
        }
        printf(" : Opt -> %d", m.OptOrder);
        if (m.order != m.OptOrder) printf(" != Matrix Order !!");
        printf("\n");
    }

    //    printf("Transport Matrix:  YTA-terms (abs(theta))\n");
    //    for (vsiz_t i = 0; i < fYTAMatrixElems.size(); i++) {
    //        const THaMatrixElement& m = fYTAMatrixElems[i];
    //        for (vsiz_t j = 0; j < m.pw.size(); j++) {
    //            printf("  %2d", m.pw[j]);
    //        }
    //        for (Int_t j = 0; j < m.order; j++) {
    //            printf("\t%g", m.poly[j]);
    //        }
    //        printf(" : Opt -> %d", m.OptOrder);
    //        printf("\n");
    //    }

    printf("OpticsOpt::Print: Transport Matrix:  P-terms\n");
    for (vsiz_t i = 0; i < fPMatrixElems.size(); i++) {
        const THaMatrixElement& m = fPMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            printf("  %2d", m.pw[j]);
        }
        for (Int_t j = 0; j < m.order; j++) {
            printf("\t%g", m.poly[j]);
        }
        printf(" : Opt -> %d", m.OptOrder);
        if (m.order != m.OptOrder) printf(" != Matrix Order !!");
        printf("\n");
    }

    //    printf("Transport Matrix:  PTA-terms\n");
    //    for (vsiz_t i = 0; i < fPTAMatrixElems.size(); i++) {
    //        const THaMatrixElement& m = fPTAMatrixElems[i];
    //        for (vsiz_t j = 0; j < m.pw.size(); j++) {
    //            printf("  %2d", m.pw[j]);
    //        }
    //        for (Int_t j = 0; j < m.order; j++) {
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
    //        for (Int_t j = 0; j < m.order; j++) {
    //            printf("\t%g", m.poly[j]);
    //        }
    //        printf(" : Opt -> %d", m.OptOrder);
    //        printf("\n");
    //    }

    printf("fArbitaryDpKinShift[%d] = {", NKines);
    for (Int_t KineID = 0; KineID < NKines; KineID++)
        printf("%f  ", fArbitaryDpKinShift[KineID]);
    printf("}\n");

    return;
}

Int_t OpticsOpt::Matrix2Array(Double_t Array[], Bool_t FreeParaFlag[], Int_t UseFPOff)
{
    assert(fCurrentMatrixElems);
    return Matrix2Array(Array, (*fCurrentMatrixElems), FreeParaFlag, UseFPOff);
}

Int_t OpticsOpt::Matrix2Array(Double_t Array[], const std::vector<THaMatrixElement> &Matrix, Bool_t FreeParaFlag[], Int_t UseFPOff)
{
    // Matrix -> Array

    typedef vector<THaMatrixElement>::size_type vsiz_t;

    Int_t idx = 0;

    for (vsiz_t i = 0; i < Matrix.size(); i++) {
        const THaMatrixElement& m = Matrix[i];
        Int_t j;
        for (j = 0; (Int_t) j < m.order; j++) {
            if (FreeParaFlag) FreeParaFlag[idx] = j < m.OptOrder ? kTRUE : kFALSE;
            Array[idx++] = m.poly[j];
        }
        for (; j < kPORDER; j++) {
            if (FreeParaFlag) FreeParaFlag[idx] = j < m.OptOrder ? kTRUE : kFALSE;
            Array[idx++] = 0;
        }
    }

    if (UseFPOff != 0) {
        for (vsiz_t i = 0; i < fFPMatrixElems.size(); i++) {
            const THaMatrixElement& m = fFPMatrixElems[i];
            Int_t j;
            for (j = 0; (Int_t) j < m.order; j++) {
                if (FreeParaFlag) FreeParaFlag[idx] = j < m.OptOrder ? kTRUE : kFALSE;
                Array[idx++] = m.poly[j];
            }
            for (; j < kPORDER; j++) {
                if (FreeParaFlag) FreeParaFlag[idx] = j < m.OptOrder ? kTRUE : kFALSE;
                Array[idx++] = 0;
            }
        }
    }

    DEBUG_INFO("Matrix2Array", "Fill Size = %d", idx);

    return idx;
}

Int_t OpticsOpt::Array2Matrix(const Double_t Array[], Int_t UseFPOff)
{
    assert(fCurrentMatrixElems);
    return Array2Matrix(Array, (*fCurrentMatrixElems), UseFPOff);
}

Int_t OpticsOpt::Array2Matrix(const Double_t Array[], std::vector<THaMatrixElement> &Matrix, Int_t UseFPOff)
{
    // Array -> fCurrentMatrixElems

    typedef vector<THaMatrixElement>::size_type vsiz_t;

    Int_t idx = 0;

    for (vsiz_t i = 0; i < Matrix.size(); i++) {
        THaMatrixElement& m = Matrix[i];
        Int_t j;
        m.order = kPORDER;
        for (j = 0; j < m.order; j++) {
            m.poly[j] = Array[idx];
            idx++;
        }
        m.SkimPoly();
    }

    if (UseFPOff != 0) {
        for (vsiz_t i = 0; i < fFPMatrixElems.size(); i++) {
            THaMatrixElement& m = fFPMatrixElems[i];
            Int_t j;
            m.order = kPORDER;
            for (j = 0; j < m.order; j++) {
                m.poly[j] = Array[idx];
                idx++;
            }
            m.SkimPoly();
        }
    }

    DEBUG_INFO("Array2Matrix", "Load Size = %d", idx);

    return idx;
}

///////////////////////////////////////////////////////////////////////////////
// Data storage
///////////////////////////////////////////////////////////////////////////////

void OpticsOpt::DCS2FCS(const Double_t* det, Double_t* rot)
{
    Double_t xdet = det[0];
    Double_t tdet = det[1];
    Double_t ydet = det[2];
    Double_t pdet = det[3];

    Double_t tanrho0 = fFPMatrixElems[0].poly[0];
    Double_t cosrho0 = 1.0 / sqrt(1.0 + tanrho0 * tanrho0);

    Double_t ttra = (tdet + tanrho0) / (1.0 - tdet * tanrho0);
    Double_t ptra = pdet / (cosrho0 * (1.0 - tdet * tanrho0));
    Double_t xtra = xdet * cosrho0 * (1 + ttra * tanrho0);
    Double_t ytra = ydet + tanrho0 * cosrho0 * ptra * xdet;

    Double_t x = xtra;

    CalcMatrix(x, fFPMatrixElems);

    Double_t tanrho = fFPMatrixElems[0].v;
    Double_t cosrho = 1.0 / sqrt(1.0 + tanrho * tanrho);

    Double_t y = ytra - fFPMatrixElems[1].v;
    Double_t t = (tdet + tanrho) / (1.0 - tdet * tanrho);
    Double_t p = (pdet - fFPMatrixElems[2].v) / ((1.0 - tdet * tanrho) * cosrho);

    rot[0] = x;
    rot[1] = t;
    rot[2] = y;
    rot[3] = p;
}

Int_t OpticsOpt::LoadRawData(TString DataFileName, Int_t NLoad, Int_t MaxDataPerGroup)
{
    // Load "f51" ascii data file to Rawdata[]

    DEBUG_INFO("LoadRawData", "Loading %s", DataFileName.Data());

    if (BeamShiftX != 0)
        DEBUG_WARNING("LoadRawData", "Shift Beam X = %f", BeamShiftX);

    Int_t datagrpcnt[kMaxDataGroup];
    memset(datagrpcnt, 0, sizeof (datagrpcnt));

    FILE* file = fopen(DataFileName, "r");
    if (!file) return kFileError;

    Int_t NRead = 0;
    const Int_t LEN = 2000;
    Char_t buff[LEN];

    Double_t NDataRead = 0;
    Int_t NLineRead = 0;

    while (fgets(buff, LEN, file)) {
        NLineRead++;

        if (NLineRead % 100000 == 0) DEBUG_INFO("LoadRawData", "%d/%d Entries Loaded", NRead, NLineRead);

        assert(NRead < kMaxNRawData); //too much data if fails

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

        assert(line_spl.size() <= kMaxNEventData); // array size check
        for (vector<string>::size_type i = 0; i < line_spl.size(); i++)
            eventdata[i] = atof(line_spl[i].c_str());

        // WARNING : shift beam x
        if (BeamShiftX != 0)
            eventdata[kBeamX] += BeamShiftX;

        // determine whether to save this data
        Int_t cutid = eventdata[kCutID];
        assert(cutid < kMaxDataGroup); // too many cuts
        Int_t& grpcnt = datagrpcnt[cutid];
        grpcnt++;
        if (grpcnt > MaxDataPerGroup) {
            DEBUG_MASSINFO("LoadRawData", "ignore data %d from cutid %d (%d ev total)", NLineRead, cutid, grpcnt);
            continue;
        }

        NDataRead += line_spl.size();

        Double_t(*powers)[5] = fRawData[NRead].powers;

        Double_t det[4] = {eventdata[kDetX], eventdata[kDetTh], eventdata[kDetY], eventdata[kDetPh]};
        Double_t rot[4] = {0.0, 0.0, 0.0, 0.0};
        DCS2FCS(det, rot);
        eventdata[kRotX] = rot[0];
        eventdata[kRotTh] = rot[1];
        eventdata[kRotY] = rot[2];
        eventdata[kRotPh] = rot[3];
        Double_t x_fp = eventdata[kRotX];
        Double_t th_fp = eventdata[kRotTh];
        Double_t y_fp = eventdata[kRotY];
        Double_t ph_fp = eventdata[kRotPh];

        // calculate the powers we need
        for (Int_t i = 0; i < kPreCalPow; i++) {
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

    Int_t goodstatcut = 0, actcutcnt = 0;
    for (Int_t i = 0; i < kMaxDataGroup; i++) {
        if (datagrpcnt[i] > 0) actcutcnt++;
        if (datagrpcnt[i] >= MaxDataPerGroup) goodstatcut++;
    }

    DEBUG_INFO("LoadRawData", "Event Limit/Cut = %d, %d / %d ev read, %d / %d cut have enough ev", MaxDataPerGroup, NRead, NLineRead, goodstatcut, actcutcnt);
    DEBUG_INFO("LoadRawData", "%d events x %f record/event read from %s", fNRawData, NDataRead / fNRawData, DataFileName.Data());

    return NRead;
}

///////////////////////////////////////////////////////////////////////////////
// Optimization related Commands
///////////////////////////////////////////////////////////////////////////////

const TVector3 OpticsOpt::GetSieveHoleTCS(Int_t Col, Int_t Row)
{
    assert(Col < NSieveCol);
    assert(Row < NSieveRow);
    TVector3 SieveHoleTCS(SieveOffX + SieveXbyRow[Row], SieveOffY + SieveYbyCol[Col], ZPos);
    return SieveHoleTCS;
}

void OpticsOpt::PrepareSieve(void)
{
    // Calculate kRealTh, kRealPh

    Double_t exttargcorr_th = 0, rms_exttargcorr_th = 0;
    Double_t exttargcorr_ph = 0, rms_exttargcorr_ph = 0;

    for (Int_t i = 0; i < fNRawData; i++) {
        EventData& eventdata = fRawData[i];

        Int_t res = eventdata.Data[kCutID];
        // const Int_t KineID = res / (NSieveRow * NSieveCol * NFoils); //starting 0!
        res = res % (NSieveRow * NSieveCol * NFoils);
        const Int_t FoilID = res / (NSieveRow * NSieveCol); //starting 0!
        res = res % (NSieveRow * NSieveCol);
        const Int_t Col = res / (NSieveRow); //starting 0!
        const Int_t Row = res % (NSieveRow); //starting 0!

        assert(FoilID < NFoils); //check array index size

        const TVector3 SieveHoleTCS = GetSieveHoleTCS(Col, Row);

        const TVector3 BeamSpotHCS(eventdata.Data[kBeamX], eventdata.Data[kBeamY], targetfoils[FoilID]);
        const TVector3 BeamSpotTCS = fTCSInHCS.Inverse()*(BeamSpotHCS - fPointingOffset);

        const TVector3 MomDirectionTCS = SieveHoleTCS - BeamSpotTCS;

        Double_t x_tg, y_tg;

        if (TargetField) {
            eventdata.Data[kRealTh] = eventdata.Data[kSimTh];
            eventdata.Data[kRealPh] = eventdata.Data[kSimPh];
            x_tg = eventdata.Data[kSimX];
            y_tg = eventdata.Data[kSimY];
        } else {
            eventdata.Data[kRealTh] = MomDirectionTCS.X() / MomDirectionTCS.Z();
            eventdata.Data[kRealPh] = MomDirectionTCS.Y() / MomDirectionTCS.Z();
            x_tg = BeamSpotTCS.X() - BeamSpotTCS.Z() * eventdata.Data[kRealTh];
            y_tg = BeamSpotTCS.Y() - BeamSpotTCS.Z() * eventdata.Data[kRealPh];
        }

        eventdata.Data[kRealX] = x_tg;
        eventdata.Data[kRealY] = y_tg;

        // Expected th ph before ext. target correction
        // fDeltaTh = fThetaCorr * x_tg;
        // Double_t theta = trkifo->GetTheta() + fDeltaTh;
        eventdata.Data[kRealThMatrix] = eventdata.Data[kRealTh] - x_tg * ExtTarCor_ThetaCorr;
        eventdata.Data[kRealPhMatrix] = eventdata.Data[kRealPh] - y_tg * ExtTarCor_PhiCorr;

        exttargcorr_th += x_tg * ExtTarCor_ThetaCorr;
        rms_exttargcorr_th += x_tg * ExtTarCor_ThetaCorr * x_tg * ExtTarCor_ThetaCorr;
        exttargcorr_ph += y_tg * ExtTarCor_PhiCorr;
        rms_exttargcorr_ph += y_tg * ExtTarCor_PhiCorr * y_tg * ExtTarCor_PhiCorr;

        DEBUG_MASSINFO("PrepareSieve", "Reference Angle: th = %f,\t phi = %f", eventdata.Data[kRealThMatrix], eventdata.Data[kRealPh]);
    }

    DEBUG_INFO("PrepareSieve", "Average Extended Target Correction: th = %f,\t rmsth = %f", exttargcorr_th / fNRawData, TMath::Sqrt(rms_exttargcorr_th / fNRawData));
    DEBUG_INFO("PrepareSieve", "Average Extended Target Correction: phi = %f,\t rmsphi = %f", exttargcorr_ph / fNRawData, TMath::Sqrt(rms_exttargcorr_ph / fNRawData));

    // Make sure kCalcTh, kCalcPh is filled
    SumSquareDTh();
    SumSquareDPhi();

    DEBUG_INFO("PrepareSieve", "Done!");
}

Double_t OpticsOpt::SumSquareDTh(Int_t UseFPOff)
{
    // return square sum of diff between calculated tg_th and expected tg_th

    Double_t dth = 0; //Difference
    Double_t rmsth = 0; //mean square

    static Int_t NCall = 0;
    NCall++;

    Double_t theta;

    for (Int_t i = 0; i < fNRawData; i++) {
        EventData& eventdata = fRawData[i];

        Double_t(*powers)[5] = eventdata.powers;

        if (UseFPOff != 0) {
            Double_t det[4] = {eventdata.Data[kDetX], eventdata.Data[kDetTh], eventdata.Data[kDetY], eventdata.Data[kDetPh]};
            Double_t rot[4] = {0.0, 0.0, 0.0, 0.0};
            DCS2FCS(det, rot);
            eventdata.Data[kRotX] = rot[0];
            eventdata.Data[kRotTh] = rot[1];
            eventdata.Data[kRotY] = rot[2];
            eventdata.Data[kRotPh] = rot[3];

            // calculate the powers we need
            for (Int_t i = 0; i < kPreCalPow; i++) {
                powers[i][0] = pow(rot[0], i);
                powers[i][1] = pow(rot[1], i);
                powers[i][2] = pow(rot[2], i);
                powers[i][3] = pow(rot[3], i);
                powers[i][4] = pow(TMath::Abs(rot[1]), i);
            }
        }

        Double_t x_fp = eventdata.Data[kRotX];

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

        DEBUG_MASSINFO("SumSquareDTh", "dth = %f = \t%f - \t%f", theta - eventdata.Data[kRealThMatrix], theta, eventdata.Data[kRealThMatrix]);

        // save the results
        eventdata.Data[kCalcTh] = theta;
    }

    DEBUG_INFO("SumSquareDTh", "#%d : dth = %f,\t rmsth = %f", NCall, dth / fNRawData, TMath::Sqrt(rmsth / fNRawData));

    return rmsth;
}

Double_t OpticsOpt::SumSquareDPhi(Int_t UseFPOff)
{
    // return square sum of diff between calculated tg_ph and expected tg_ph

    Double_t dphi = 0; //Difference
    Double_t rmsphi = 0; //mean square

    static Int_t NCall = 0;
    NCall++;

    Double_t phi;

    for (Int_t i = 0; i < fNRawData; i++) {
        EventData& eventdata = fRawData[i];

        Double_t(*powers)[5] = eventdata.powers;

        if (UseFPOff != 0) {
            Double_t det[4] = {eventdata.Data[kDetX], eventdata.Data[kDetTh], eventdata.Data[kDetY], eventdata.Data[kDetPh]};
            Double_t rot[4] = {0.0, 0.0, 0.0, 0.0};
            DCS2FCS(det, rot);
            eventdata.Data[kRotX] = rot[0];
            eventdata.Data[kRotTh] = rot[1];
            eventdata.Data[kRotY] = rot[2];
            eventdata.Data[kRotPh] = rot[3];

            // calculate the powers we need
            for (Int_t i = 0; i < kPreCalPow; i++) {
                powers[i][0] = pow(rot[0], i);
                powers[i][1] = pow(rot[1], i);
                powers[i][2] = pow(rot[2], i);
                powers[i][3] = pow(rot[3], i);
                powers[i][4] = pow(TMath::Abs(rot[1]), i);
            }
        }

        Double_t x_fp = eventdata.Data[kRotX];

        // calculate the matrices we need
        // CalcMatrix(x_fp, fDMatrixElems);
        // CalcMatrix(x_fp, fTMatrixElems);
        // CalcMatrix(x_fp, fYMatrixElems);
        // CalcMatrix(x_fp, fYTAMatrixElems);
        CalcMatrix(x_fp, fPMatrixElems);
        //CalcMatrix(x_fp, fPTAMatrixElems);

        // calculate the coordinates at the target
        //phi = CalcTargetVar(fPMatrixElems, powers) + CalcTargetVar(fPTAMatrixElems, powers);
        phi = CalcTargetVar(fPMatrixElems, powers);

        dphi += phi - eventdata.Data[kRealPhMatrix];
        rmsphi += (phi - eventdata.Data[kRealPhMatrix])*(phi - eventdata.Data[kRealPhMatrix]);

        DEBUG_MASSINFO("SumSquareDPhi", "dphi = %f = \t%f - \t%f", phi - eventdata.Data[kRealPhMatrix], phi, eventdata.Data[kRealPhMatrix]);

        //save the results
        eventdata.Data[kCalcPh] = phi;
    }

    DEBUG_INFO("SumSquareDPhi", "#%d : dphi = %f,\t rmsphi = %f", NCall, dphi / fNRawData, TMath::Sqrt(rmsphi / fNRawData));

    return rmsphi;
}

TList* OpticsOpt::CheckSieve(Int_t PlotType)
{
    // Visualize Sieve Plane

    const Int_t NPlots = (PlotType != 0) ? NFoils * NKines : 1;
    assert(NKines <= 12);

    TH2D * HSievePlane[NPlots];
    const Double_t x_lim = 1.3 * TMath::Max(TMath::Abs(SieveYbyCol[0]), TMath::Abs(SieveYbyCol[NSieveCol - 1]));
    const Double_t y_lim = 1.5 * TMath::Max(TMath::Abs(SieveXbyRow[0]), TMath::Abs(SieveXbyRow[NSieveRow - 1]));

    for (Int_t i = 0; i < NPlots; i++) {
        HSievePlane[i] = new TH2D(Form("Dataset%d", i), Form("Sieve Plane Proj. (tg_X vs tg_Y) for Kine #%d", i / NFoils), 400, -x_lim, x_lim, 400, -y_lim, y_lim);
        HSievePlane[i]->SetXTitle("Sieve H [m]");
        HSievePlane[i]->SetYTitle("Sieve V [m]");
        assert(HSievePlane[i]); // assure memory allocation
    }

    Double_t dX = 0, dY = 0;

    enum {
        kEvNum, kRealSieveX, kRealSieveY, kCalcSieveX, kCalcSieveY
    };

    Double_t SieveEventID[NPlots][NSieveCol][NSieveRow][5];
    memset(SieveEventID, 0, sizeof (SieveEventID));

    for (Int_t i = 0; i < fNRawData; i++) {
        const EventData& eventdata = fRawData[i];

        Int_t res = eventdata.Data[kCutID];
        Int_t PlotID = res / (NSieveRow * NSieveCol); //starting 0!
        res = res % (NSieveRow * NSieveCol);
        const Int_t Col = res / (NSieveRow); //starting 0!
        const Int_t Row = res % (NSieveRow); //starting 0!

        assert(PlotID < NKines * NFoils); //array index check
        if (NPlots == 1) PlotID = 0;

        const TVector3 SieveHoleTCS = GetSieveHoleTCS(Col, Row);

        Double_t ProjectionX = eventdata.Data[kRealX] + (eventdata.Data[kCalcTh] + eventdata.Data[kRealX] * ExtTarCor_ThetaCorr) * (SieveHoleTCS.Z());
        Double_t ProjectionY = eventdata.Data[kRealY] + (eventdata.Data[kCalcPh] + eventdata.Data[kRealY] * ExtTarCor_PhiCorr) * (SieveHoleTCS.Z());

        HSievePlane[PlotID]->Fill(ProjectionY, ProjectionX);

        dX += ProjectionX - SieveHoleTCS.X();
        dY += ProjectionY - SieveHoleTCS.Y();

        SieveEventID[PlotID][Col][Row][kEvNum]++;
        SieveEventID[PlotID][Col][Row][kRealSieveX] += SieveHoleTCS.X();
        SieveEventID[PlotID][Col][Row][kRealSieveY] += SieveHoleTCS.Y();
        SieveEventID[PlotID][Col][Row][kCalcSieveX] += ProjectionX;
        SieveEventID[PlotID][Col][Row][kCalcSieveY] += ProjectionY;
    }

    DEBUG_INFO("CheckSieve", "Average : dX = %f,\t dY = %f", dX / fNRawData, dY / fNRawData);

    TCanvas * c[NFoils];
    Double_t w = 720, h = 960;
    Int_t splitw = 1, splith = 1;
    if (NPlots == 1) {
        c[0] = new TCanvas("SieveCheck", "SieveCheck", w, h);
    } else {
        if (NKines <= 4) {
            w = 1440;
            h = 960;
            splitw = 4;
            splith = 1;
        } else if (NKines <= 8) {
            w = 1440;
            h = 960;
            splitw = 4;
            splith = 2;
        } else if (NKines <= 12) {
            w = 1440;
            h = 640;
            splitw = 6;
            splith = 2;
        }
        for (Int_t i = 0; i < NFoils; i++) {
            c[i] = new TCanvas(Form("SieveCheck%d", i), Form("SieveCheck%d", i), w, h);
            c[i]->Divide(splitw, splith);
        }
    }

    for (Int_t i = 0; i < NPlots; i++) {
        const Int_t FoilID = i % NFoils;
        const Int_t KineID = i / NFoils;
        if (NPlots == 1) c[0]->cd();
        else c[FoilID]->cd(KineID + 1);

        assert(HSievePlane[i]); //pointer check
        HSievePlane[i]->Draw("COLZ");

        // Draw Sieve
        const Double_t plotwidth = 0.004;
        for (Int_t Row = 0; Row < NSieveRow; Row++) {
            for (Int_t Col = 0; Col < NSieveCol; Col++) {
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
        for (Int_t Col = 0; Col < NSieveCol; Col++) {
            for (Int_t Row = 0; Row < NSieveRow; Row++) {
                SieveEventID[i][Col][Row][kRealSieveX] /= SieveEventID[i][Col][Row][kEvNum];
                SieveEventID[i][Col][Row][kRealSieveY] /= SieveEventID[i][Col][Row][kEvNum];
                SieveEventID[i][Col][Row][kCalcSieveX] /= SieveEventID[i][Col][Row][kEvNum];
                SieveEventID[i][Col][Row][kCalcSieveY] /= SieveEventID[i][Col][Row][kEvNum];
                if (SieveEventID[i][Col][Row][kEvNum] > 5) {
                    TArrow * ar2 = new TArrow(SieveEventID[i][Col][Row][kCalcSieveY], SieveEventID[i][Col][Row][kCalcSieveX], SieveEventID[i][Col][Row][kRealSieveY], SieveEventID[i][Col][Row][kRealSieveX], 0.008, "|>");
                    ar2->SetAngle(40);
                    ar2->SetLineColor(kMagenta);
                    ar2->SetFillColor(kMagenta);

                    const Double_t ignorelimit = 0.005;
                    if ((ar2->GetX1() - ar2->GetX2())*(ar2->GetX1() - ar2->GetX2())+(ar2->GetY1() - ar2->GetY2())*(ar2->GetY1() - ar2->GetY2()) > ignorelimit * ignorelimit) ar2->Draw();
                }
            }
        }
    }

    TList* l = new TList();
    if (NPlots == 1) l->Add(c[0]);
    else for (Int_t i = 0; i < NFoils; i++) l->Add(c[i]);

    return l;
}

void OpticsOpt::PrepareVertex(void)
{
    // calculate kRealY

    for (Int_t i = 0; i < fNRawData; i++) {
        EventData& eventdata = fRawData[i];

        Int_t res = eventdata.Data[kCutID];
        // const Int_t KineID = res / (NSieveRow * NSieveCol * NFoils); //starting 0!
        res = res % (NSieveRow * NSieveCol * NFoils);
        const Int_t FoilID = res / (NSieveRow * NSieveCol); //starting 0!
        res = res % (NSieveRow * NSieveCol);
        const Int_t Col = res / (NSieveRow); //starting 0!
        const Int_t Row = res % (NSieveRow); //starting 0!

        assert(FoilID < NFoils); //check array index size

        eventdata.Data[kFoilID] = FoilID;

        const TVector3 SieveHoleTCS = GetSieveHoleTCS(Col, Row);

        const TVector3 BeamSpotHCS(eventdata.Data[kBeamX], eventdata.Data[kBeamY], targetfoils[FoilID]);
        const TVector3 BeamSpotTCS = fTCSInHCS.Inverse()*(BeamSpotHCS - fPointingOffset);

        const TVector3 MomDirectionTCS = SieveHoleTCS - BeamSpotTCS;

        Double_t y_tg;

        if (TargetField) {
            eventdata.Data[kRealPh] = eventdata.Data[kSimPh];
            y_tg = eventdata.Data[kSimY];
        } else {
            eventdata.Data[kRealPh] = MomDirectionTCS.Y() / MomDirectionTCS.Z();
            y_tg = BeamSpotTCS.Y() - BeamSpotTCS.Z() * eventdata.Data[kRealPh];
        }

        eventdata.Data[kRealY] = y_tg;
        eventdata.Data[kRealReactZ] = targetfoils[FoilID];

        DEBUG_MASSINFO("PrepareVertex", "Reference: y = %f,\t z = %f", eventdata.Data[kRealY], eventdata.Data[kRealReactZ]);
    }

    // make sure kCalcTh, kCalcPh is filled
    SumSquareDY();

    DEBUG_INFO("PrepareVertex", "Done!");
}

Double_t OpticsOpt::SumSquareDY(Int_t UseFPOff)
{
    // return square sum of diff between calculated tg_y and expected tg_y

    Double_t dy = 0; //Difference
    Double_t rmsy = 0; //mean square

    static Int_t NCall = 0;
    NCall++;

    Double_t y;

    for (Int_t i = 0; i < fNRawData; i++) {
        EventData &eventdata = fRawData[i];

        Double_t(*powers)[5] = eventdata.powers;

        if (UseFPOff != 0) {
            Double_t det[4] = {eventdata.Data[kDetX], eventdata.Data[kDetTh], eventdata.Data[kDetY], eventdata.Data[kDetPh]};
            Double_t rot[4] = {0.0, 0.0, 0.0, 0.0};
            DCS2FCS(det, rot);
            eventdata.Data[kRotX] = rot[0];
            eventdata.Data[kRotTh] = rot[1];
            eventdata.Data[kRotY] = rot[2];
            eventdata.Data[kRotPh] = rot[3];

            // calculate the powers we need
            for (Int_t i = 0; i < kPreCalPow; i++) {
                powers[i][0] = pow(rot[0], i);
                powers[i][1] = pow(rot[1], i);
                powers[i][2] = pow(rot[2], i);
                powers[i][3] = pow(rot[3], i);
                powers[i][4] = pow(TMath::Abs(rot[1]), i);
            }
        }

        Double_t x_fp = eventdata.Data[kRotX];

        // calculate the matrices we need
        // CalcMatrix(x_fp, fDMatrixElems);
        // CalcMatrix(x_fp, fTMatrixElems);
        CalcMatrix(x_fp, fYMatrixElems);
        //CalcMatrix(x_fp, fYTAMatrixElems);
        // CalcMatrix(x_fp, fPMatrixElems);
        // CalcMatrix(x_fp, fPTAMatrixElems);

        // calculate the coordinates at the target
        //y = CalcTargetVar(fYMatrixElems, powers) + CalcTargetVar(fYTAMatrixElems, powers);
        y = CalcTargetVar(fYMatrixElems, powers);

        const Int_t FoilID = eventdata.Data[kFoilID];
        assert(FoilID < NFoils);

        dy += y - eventdata.Data[kRealY];
        rmsy += (y - eventdata.Data[kRealY])*(y - eventdata.Data[kRealY]);

        DEBUG_MASSINFO("SumSquareDY", "dy = %f = \t%f - \t%f", y - eventdata.Data[kRealY], y, eventdata.Data[kRealY]);

        // save the results
        eventdata.Data[kCalcY] = y;
    }

    DEBUG_INFO("SumSquareDY", "#%d : dy = %f,\t rmsy = %f", NCall, dy / fNRawData, TMath::Sqrt(rmsy / fNRawData));

    return rmsy;
}

TList* OpticsOpt::CheckY(Int_t PlotType)
{
    // Visualize Y spectrum

    const Int_t NPlots = (PlotType != 0) ? NKines : 1;
    assert(NKines <= 12);

    TH1D * HTgY[NPlots];
    TH1D * HTgYReal[NPlots];
    TH1D * HTgYDiff[NPlots];

    Double_t zmin = 1e38, zmax = -1e38;
    for (Int_t i = 0; i < NFoils; i++) {
        if (targetfoils[i] < zmin) zmin = targetfoils[i];
        if (targetfoils[i] > zmax) zmax = targetfoils[i];
    }
    Double_t ylimu = -zmin * TMath::Sin(HRSAngle) + 10e-3;
    Double_t ylimd = -zmax * TMath::Sin(HRSAngle) - 10e-3;

    for (Int_t i = 0; i < NPlots; i++) {
        HTgY[i] = new TH1D(Form("DataSet%d", i), Form("Target Y for Kine #%d", i), 400, ylimd, ylimu);
        HTgYReal[i] = new TH1D(Form("DataSet%d", i), Form("Target Y for Kine #%d", i), 400, ylimd, ylimu);
        HTgY[i]->SetXTitle("Target Y [m]");
        HTgYReal[i]->SetXTitle("Target Y [m]");

        HTgYDiff[i] = new TH1D(Form("DataSet%d", i), Form("Diff of Target Y for Kine #%d", i), 200, -10e-3, 10e-3);
        HTgYDiff[i]->SetXTitle("Cal Y - Real Y [m]");
        assert(HTgY[i]); // assure memory allocation
        assert(HTgYDiff[i]);
    }

    Double_t dy = 0;
    Double_t rmsy = 0;

    for (Int_t i = 0; i < fNRawData; i++) {
        EventData& eventdata = fRawData[i];

        Int_t res = eventdata.Data[kCutID];
        Int_t PlotID = res / (NSieveRow * NSieveCol * NFoils); //starting 0!

        assert(PlotID < NKines); //array index check
        if (NPlots == 1) PlotID = 0;

        HTgY[PlotID]->Fill(eventdata.Data[kCalcY]);
        HTgYReal[PlotID]->Fill(eventdata.Data[kRealY]);
        HTgYDiff[PlotID]->Fill(eventdata.Data[kCalcY] - eventdata.Data[kRealY]);

        dy += eventdata.Data[kCalcY] - eventdata.Data[kRealY];
        rmsy += (eventdata.Data[kCalcY] - eventdata.Data[kRealY])*(eventdata.Data[kCalcY] - eventdata.Data[kRealY]);
    }

    DEBUG_INFO("CheckTgY", "dy = %f,\t rmsy = %f", dy / fNRawData, rmsy / fNRawData);

    TCanvas *c1, *c2;
    Double_t w = 720, h = 480;
    Int_t splitw = 1, splith = 1;
    if (NPlots == 1) {
        c1 = new TCanvas("YCheck", "YCheck", w, h);
        c2 = new TCanvas("YDiff", "YDiff", w, h);
    } else {
        if (NKines <= 2) {
            w = 720;
            h = 960;
            splitw = 1;
            splith = 2;
        } else if (NKines <= 4) {
            w = 1440;
            h = 960;
            splitw = 2;
            splith = 2;
        } else if (NKines <= 6) {
            w = 1440;
            h = 640;
            splitw = 3;
            splith = 2;
        } else if (NKines <= 9) {
            w = 1440;
            h = 960;
            splitw = 3;
            splith = 3;
        } else if (NKines <= 12) {
            w = 1440;
            h = 720;
            splitw = 4;
            splith = 3;
        }
        c1 = new TCanvas("YCheck", "YCheck", w, h);
        c1->Divide(splitw, splith);
        c2 = new TCanvas("YDiff", "YDiff", w, h);
        c2->Divide(splitw, splith);
    }

    for (Int_t i = 0; i < NPlots; i++) {
        c1->cd(i + 1);

        assert(HTgY[i]);
        HTgYReal[i]->SetLineColor(kRed);
        HTgYReal[i]->Draw();
        HTgY[i]->SetLineColor(kBlack);
        HTgY[i]->Draw("same");

        c2->cd(i + 1);

        assert(HTgYDiff[i]);
        HTgYDiff[i]->Draw();

        Double_t DefResolution = 0.5e-3;
        Double_t FitRangeMultiply = 10;

        TString FitFunc = Form("YDiff%d", i);
        TF1 *f = new TF1(FitFunc, "gaus+[3]", -DefResolution*FitRangeMultiply, DefResolution * FitRangeMultiply);
        f->SetParameter(1, 0.0);
        f->SetParameter(2, DefResolution);
        HTgYDiff[i] -> Fit(FitFunc, "RN0");
        f->SetLineColor(2);
        f->Draw("same");

        TLatex *t = new TLatex(f->GetParameter(1) + DefResolution, f->GetParameter(0) + f->GetParameter(3), Form("\\Delta \\pm \\sigma = (%2.1f \\pm %2.1f) mm", 1000 * f->GetParameter(1), 1000 * (f->GetParameter(2))));
        t->SetTextSize(0.05);
        t->SetTextAlign(12);
        t->SetTextColor(2);
        t->Draw();
    }

    TList* l = new TList();
    l->Add(c1);
    l->Add(c2);

    return l;
}

void OpticsOpt::PrepareDp(void)
{
    // calculate expected dp_kin, dp_kin offsets ....

    // print Central Momentums
    printf("HRSCentralMom[%d] (GeV) = { ", NKines);
    for (Int_t KineID = 0; KineID < NKines; KineID++)
        printf("%f ", HRSCentralMom[KineID]);
    printf("}\n");

    // print radiation loss numbers
    printf("RadiationLossByFoil[%d] (MeV) = { ", NFoils);
    for (Int_t FoilID = 0; FoilID < NFoils; FoilID++)
        printf("%f ", RadiationLossByFoil[FoilID]*1000);
    printf("}\n");

    // print tilt angle
    printf("TiltAngleBySetting[%d] (deg) = { ", NFoils * NKines);
    for (Int_t i = 0; i < NFoils * NKines; i++)
        printf("%f ", TiltAngle[i]*180.0 / TMath::Pi());
    printf("}\n");

    Double_t scatang = 0;
    Double_t dpkinoff = 0, rmsdpkinoff = 0;
    Double_t exttargcorr_dp = 0, rms_exttargcorr_dp = 0;

    for (Int_t i = 0; i < fNRawData; i++) {
        DEBUG_MASSINFO("PrepareDp", "=========== Event %d ===========", i);

        EventData& eventdata = fRawData[i];

        Int_t res = eventdata.Data[kCutID];
        const Int_t KineID = res / (NSieveRow * NSieveCol * NFoils); //starting 0!
        res = res % (NSieveRow * NSieveCol * NFoils);
        const Int_t FoilID = res / (NSieveRow * NSieveCol); //starting 0!
        res = res % (NSieveRow * NSieveCol);
        const Int_t Col = res / (NSieveRow); //starting 0!
        const Int_t Row = res % (NSieveRow); //starting 0!

        assert(KineID < NKines); //check array index size
        assert(FoilID < NFoils); //check array index size

        // write some variables
        eventdata.Data[kKineID] = KineID;
        eventdata.Data[kCentralp] = HRSCentralMom[KineID];

        TVector3 MomDirectionTCS;
        Double_t x_tg;

        if (TargetField) {
            eventdata.Data[kRealTh] = eventdata.Data[kSimOrTh];
            eventdata.Data[kRealPh] = eventdata.Data[kSimOrPh];

            MomDirectionTCS.SetXYZ(eventdata.Data[kRealTh], eventdata.Data[kRealPh], 1);

            x_tg = eventdata.Data[kSimX];
        } else {
            const TVector3 SieveHoleTCS = GetSieveHoleTCS(Col, Row);

            const TVector3 BeamSpotHCS(eventdata.Data[kBeamX], eventdata.Data[kBeamY], targetfoils[FoilID]);
            const TVector3 BeamSpotTCS = fTCSInHCS.Inverse()*(BeamSpotHCS - fPointingOffset);

            MomDirectionTCS = SieveHoleTCS - BeamSpotTCS;

            eventdata.Data[kRealTh] = MomDirectionTCS.X() / MomDirectionTCS.Z();
            eventdata.Data[kRealPh] = MomDirectionTCS.Y() / MomDirectionTCS.Z();

            x_tg = BeamSpotTCS.X() - BeamSpotTCS.Z() * eventdata.Data[kRealTh];
        }

        eventdata.Data[kRealX] = x_tg;

        DEBUG_MASSINFO("PrepareDp", "Reference Angle: th = %f,\t phi = %f", eventdata.Data[kRealTh], eventdata.Data[kRealPh]);

        TVector3 MomDirectionHCS = fTCSInHCS*MomDirectionTCS;
        TVector3 BeamDirection(0, TMath::Tan(TiltAngle[KineID * NFoils + FoilID]), 1);
        const Double_t ScatteringAngle = BeamDirection.Angle(MomDirectionHCS);
        eventdata.Data[kScatterAngle] = ScatteringAngle;
        scatang += ScatteringAngle;

        // calculate difference between dp_kin and dp
        // dp_kin + kDpKinOffsets = dp
        const Double_t DM = ExcitationEnergy[KineID];
        const Double_t Ma = GroundNuclearMass;
        const Double_t P0 = eventdata.Data[kBeamE];
        const Double_t DpKinOffsets = (ScatMom(DM, Ma, P0, ScatteringAngle) - ScatMom(DM, Ma, P0, TMath::Abs(HRSAngle))) / eventdata.Data[kCentralp];
        eventdata.Data[kDpKinOffsets] = DpKinOffsets;

        dpkinoff += DpKinOffsets;
        rmsdpkinoff += DpKinOffsets*DpKinOffsets;

        // calculate kRealDpKin, should be same for same kine settings
        eventdata.Data[kRadiLossDp] = RadiationLossByFoil[FoilID] / eventdata.Data[kCentralp];
        eventdata.Data[kRealDpKin] = ScatMom(DM, Ma, P0, TMath::Abs(HRSAngle)) / eventdata.Data[kCentralp] - 1 - eventdata.Data[kRadiLossDp];

        DEBUG_MASSINFO("PrepareDp", "Reference Dp: pkin = %f,\t p0 = %f,\t radloss = %f,\t ebeam = %f", ScatMom(DM, Ma, P0, TMath::Abs(HRSAngle)), eventdata.Data[kCentralp], eventdata.Data[kRadiLossDp], P0);

        eventdata.Data[kRealDpKinMatrix] = eventdata.Data[kRealDpKin] - x_tg / ExtTarCor_DeltaCorr;

        exttargcorr_dp += x_tg / ExtTarCor_DeltaCorr;
        rms_exttargcorr_dp += (x_tg / ExtTarCor_DeltaCorr)*(x_tg / ExtTarCor_DeltaCorr);

        // calculate expected dp_kin for all other excitation states
        for (Int_t ExcitID = 0; ExcitID < NExcitationStates; ExcitID++) {
            assert(kRealDpKinExcitations + ExcitID < kMaxNEventData); //check array index size
            eventdata.Data[kRealDpKinExcitations + ExcitID] = ScatMom(ExcitationEnergyList[ExcitID], Ma, P0, TMath::Abs(HRSAngle)) / eventdata.Data[kCentralp] - 1;
        }

        DEBUG_MASSINFO("PrepareDp", "Reference Dp: scatang = %f,\t dpkinoff = %f", ScatteringAngle / TMath::Pi()*180, DpKinOffsets);
    }

    DEBUG_INFO("PrepareDp", "Average : scatang = %f", scatang / fNRawData / TMath::Pi()*180);
    DEBUG_INFO("PrepareDp", "Average : dpkinoff = %f, rmsdpkinoff = %f", dpkinoff / fNRawData, TMath::Sqrt(rmsdpkinoff / fNRawData));
    DEBUG_INFO("PrepareDp", "Average Extended Target Corretion: dp = %f,\t rmsdp = %f", exttargcorr_dp / fNRawData, TMath::Sqrt(rms_exttargcorr_dp / fNRawData));

    // make sure kCalcTh, kCalcPh is filled, although not necessary
    SumSquareDTh();
    SumSquareDPhi();
    SumSquareDp();

    DEBUG_INFO("PrepareDp", "Done!");
}

Double_t OpticsOpt::SumSquareDp(Int_t UseFPOff)
{
    // return square sum of diff between calculated dp_kin and expected dp_kin

    Double_t ddp = 0; //Difference
    Double_t rmsdp = 0; //mean square

    static Int_t NCall = 0;
    NCall++;

    Double_t dp, dp_kin;

    for (Int_t i = 0; i < fNRawData; i++) {
        EventData &eventdata = fRawData[i];

        Double_t(*powers)[5] = eventdata.powers;

        if (UseFPOff != 0) {
            Double_t det[4] = {eventdata.Data[kDetX], eventdata.Data[kDetTh], eventdata.Data[kDetY], eventdata.Data[kDetPh]};
            Double_t rot[4] = {0.0, 0.0, 0.0, 0.0};
            DCS2FCS(det, rot);
            eventdata.Data[kRotX] = rot[0];
            eventdata.Data[kRotTh] = rot[1];
            eventdata.Data[kRotY] = rot[2];
            eventdata.Data[kRotPh] = rot[3];

            // calculate the powers we need
            for (Int_t i = 0; i < kPreCalPow; i++) {
                powers[i][0] = pow(rot[0], i);
                powers[i][1] = pow(rot[1], i);
                powers[i][2] = pow(rot[2], i);
                powers[i][3] = pow(rot[3], i);
                powers[i][4] = pow(TMath::Abs(rot[1]), i);
            }
        }

        Double_t x_fp = eventdata.Data[kRotX];

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

        const Int_t KineID = (eventdata.Data[kKineID]);
        assert(KineID < NKines); //check array index size
        const Double_t ArbitaryDpKinShift = fArbitaryDpKinShift[KineID];

        ddp += dp_kin - eventdata.Data[kRealDpKinMatrix] + ArbitaryDpKinShift;
        rmsdp += (dp_kin - eventdata.Data[kRealDpKinMatrix] + ArbitaryDpKinShift)*(dp_kin - eventdata.Data[kRealDpKinMatrix] + ArbitaryDpKinShift);

        DEBUG_MASSINFO("SumSquareDp", "ddp = %f = \t%f - \t%f", dp_kin - eventdata.Data[kRealDpKinMatrix] + ArbitaryDpKinShift, dp_kin, eventdata.Data[kRealDpKinMatrix] + ArbitaryDpKinShift);

        // save the results
        eventdata.Data[kCalcDpKinMatrix] = dp_kin;
        eventdata.Data[kCalcDpKin] = dp_kin + eventdata.Data[kRealX] / ExtTarCor_DeltaCorr;
    }

    DEBUG_INFO("SumSquareDp", "#%d : ddp = %f, rmsdp = %f", NCall, ddp / fNRawData, TMath::Sqrt(rmsdp / fNRawData));

    return rmsdp;
}

TList* OpticsOpt::CheckDp(Int_t PlotType)
{
    // Visualize 1D hitogram of dp_kin
    assert(NKines <= 12);

    const Double_t DpRange = .05;
    const Int_t NDpRange = 5000;

    TH1D * hDpKinAll[NKines];
    memset(hDpKinAll, 0, sizeof (hDpKinAll));
    Double_t RealDpKin[NKines];
    memset(RealDpKin, 0, sizeof (RealDpKin));
    Double_t AverCalcDpKin[NKines];
    memset(AverCalcDpKin, 0, sizeof (AverCalcDpKin));
    Int_t NEvntDpKin[NKines];
    memset(NEvntDpKin, 0, sizeof (NEvntDpKin));
    Double_t RealDpKinAllExcit[NExcitationStates][NKines];
    memset(RealDpKinAllExcit, 0, sizeof (RealDpKinAllExcit));
    Double_t NewArbitaryDpKinShift[NKines];
    memset(NewArbitaryDpKinShift, 0, sizeof (NewArbitaryDpKinShift));

    for (Int_t KineID = 0; KineID < NKines; KineID++) {
        hDpKinAll[KineID] = new TH1D(Form("hDpKinAll%d", KineID), Form("Dp_Kin for Delta Scan Kine. #%d (All Data)", KineID), NDpRange, -DpRange, DpRange);
        assert(hDpKinAll[KineID]); //pointer check
    }

    for (Int_t i = 0; i < fNRawData; i++) {
        const EventData& eventdata = fRawData[i];

        const Int_t KineID = eventdata.Data[kKineID];
        assert(KineID < NKines);

        NEvntDpKin[KineID]++;

        hDpKinAll[KineID]->Fill(eventdata.Data[kCalcDpKin] + eventdata.Data[kRadiLossDp]);
        AverCalcDpKin[KineID] += eventdata.Data[kCalcDpKin] + eventdata.Data[kRadiLossDp];
        RealDpKin[KineID] += eventdata.Data[kRealDpKin] + eventdata.Data[kRadiLossDp];
        for (Int_t ExcitID = 0; ExcitID < NExcitationStates; ExcitID++) {
            assert(kRealDpKinExcitations + ExcitID < kMaxNEventData); //index check
            RealDpKinAllExcit[ExcitID][KineID] += eventdata.Data[kRealDpKinExcitations + ExcitID];
        }
    }

    TCanvas *c1;
    Double_t w = 720, h = 480;
    Int_t splitw = 1, splith = 1;
    if (NKines <= 1) {
        c1 = new TCanvas("CheckDp", "Check Dp Kin Reconstruction", w, h);
    } else {
        if (NKines <= 2) {
            w = 720;
            h = 960;
            splitw = 1;
            splith = 2;
        } else if (NKines <= 4) {
            w = 1440;
            h = 960;
            splitw = 2;
            splith = 2;
        } else if (NKines <= 6) {
            w = 1440;
            h = 640;
            splitw = 3;
            splith = 2;
        } else if (NKines <= 9) {
            w = 1440;
            h = 960;
            splitw = 3;
            splith = 3;
        } else if (NKines <= 12) {
            w = 1440;
            h = 720;
            splitw = 4;
            splith = 3;
        }
        c1 = new TCanvas("CheckDp", "Check Dp Kin Reconstruction", w, h);
        c1->Divide(splitw, splith);
    }
    TCanvas* c2 = new TCanvas("CheckDpGlobal", "Check Dp Kin Reconstruction", 720, 480);

    for (Int_t i = 0; i < NKines; i++) {
        c1->cd(i + 1);
        gPad -> SetLogy();

        AverCalcDpKin[i] /= NEvntDpKin[i];
        RealDpKin[i] /= NEvntDpKin[i];
        DEBUG_MASSINFO("CheckDp", "AverCalcDpKin[%d] = %f", i, AverCalcDpKin[i]);

        // Histograms
        hDpKinAll[i]->SetLineColor(4);
        hDpKinAll[i]->SetFillColor(4);
        hDpKinAll[i]->SetFillStyle(3008);

        const Double_t dpRange = 0.01;
        hDpKinAll[i]->SetAxisRange(RealDpKin[i] - dpRange, RealDpKin[i] + dpRange);
        hDpKinAll[i]->SetXTitle("radiation corrected dp_kin (angular independent dp)");
        hDpKinAll[i]->Draw();

        if (i == 0) {
            TH1D * temp = static_cast<TH1D*> (hDpKinAll[i]->Clone());
            temp->SetAxisRange(-0.04, 0.04);
            c2->cd();
            temp->Draw();
        } else {
            c2->cd();
            hDpKinAll[i]->Draw("same");
        }

        c1->cd(i + 1);
        // expectation lines
        for (Int_t ExcitID = 0; ExcitID < NExcitationStates; ExcitID++) {
            RealDpKinAllExcit[ExcitID][i] /= NEvntDpKin[i];
            const Double_t x = RealDpKinAllExcit[ExcitID][i];
            TLine *l = new TLine(x, 0, x, 20000);
            l->SetLineColor(3);
            l->SetLineWidth(2);
            l->Draw();
        }
        TLine *l = new TLine(RealDpKin[i], 0, RealDpKin[i], 20000);
        l->SetLineColor(6);
        l->SetLineWidth(2);
        l->Draw();
        c2->cd();
        l->Draw();

        c1->cd(i + 1);

        // Fits
        const Double_t DefResolution = 1e-4;
        const Double_t FitRangeMultiply = 5;

        TString FitFunc = Form("DpPeak%d", i);
        TF1 *f = new TF1(FitFunc, "gaus+[3]+[4]*x", AverCalcDpKin[i] - DefResolution*FitRangeMultiply, AverCalcDpKin[i] + DefResolution * FitRangeMultiply);
        f->SetParameter(1, AverCalcDpKin[i]);
        f->SetParameter(2, DefResolution);
        hDpKinAll[i] -> Fit(FitFunc, "RN0");
        f->SetLineColor(2);
        f->Draw("SAME");

        TLatex *t = new TLatex(f->GetParameter(1) + DefResolution, f->GetParameter(0) + f->GetParameter(3) + f->GetParameter(4) * f->GetParameter(1), Form("\\Delta \\pm \\sigma = (%2.1f \\pm %2.1f) \\times 10^{-4}", 10000 * (f->GetParameter(1) - RealDpKin[i]), 10000 * f->GetParameter(2)));
        t->SetTextSize(0.05);
        t->SetTextAlign(12);
        t->SetTextColor(2);
        t->Draw();

        NewArbitaryDpKinShift[i] = f->GetParameter(1) - RealDpKin[i] + fArbitaryDpKinShift[i];
    }

    Info("CheckDp", "New set of arbitary dp shifts:");
    for (Int_t i = 0; i < NKines; i++)
        printf("opt->fArbitaryDpKinShift[%d] = %e;\n", i, NewArbitaryDpKinShift[i]);

    TList* l = new TList();
    l->Add(c1);
    l->Add(c2);

    return l;
}

///////////////////////////////////////////////////////////////////////////////
// declarations for target vertex reconstruction
///////////////////////////////////////////////////////////////////////////////

void OpticsOpt::CalcMatrix(const Double_t x, vector<THaMatrixElement>& matrix)
{
    // calculates the values of the matrix elements for a given location
    // by evaluating a polynomial in x of order it->order with
    // coefficients given by it->poly

    for (vector<THaMatrixElement>::iterator it = matrix.begin();
            it != matrix.end(); it++) {
        it->v = 0.0;

        if (it->order > 0) {
            for (Int_t i = it->order - 1; i >= 1; i--)
                it->v = x * (it->v + it->poly[i]);
            it->v += it->poly[0];
        }
    }
}

Double_t OpticsOpt::CalcTargetVar(const vector<THaMatrixElement>& matrix, const Double_t powers[][5])
{
    // calculates the value of a variable at the target
    // the x-dependence is already in the matrix, so only 1-3 (or np) used

    typedef vector<Int_t>::size_type vsiz_t;

    Double_t retval = 0.0;
    Double_t v = 0;
    for (vector<THaMatrixElement>::const_iterator it = matrix.begin();
            it != matrix.end(); it++)
        if (it->v != 0.0) {
            v = it->v;
            vsiz_t np = it->pw.size(); // generalize for extra matrix elems.

            for (vsiz_t i = 0; i < np; i++)
                v *= powers[it->pw[i]][i + 1];
            retval += v;
        }

    return retval;
}

///////////////////////////////////////////////////////////////////////////////
// Inherited from THaTrackingDetector
///////////////////////////////////////////////////////////////////////////////

Int_t OpticsOpt::Decode(const THaEvData&)
{

    return 0;
}

Int_t OpticsOpt::CoarseTrack(TClonesArray&)
{

    return 0;
}

Int_t OpticsOpt::FineTrack(TClonesArray&)
{

    return 0;
}

THaAnalysisObject::EStatus OpticsOpt::Init(const TDatime&)
{

    return fStatus = kOK;
};

///////////////////////////////////////////////////////////////////////////////
// class for storing matrix element data
///////////////////////////////////////////////////////////////////////////////

THaMatrixElement::THaMatrixElement() : iszero(kTRUE), pw(3), order(0), v(0), poly(OpticsOpt::kPORDER), OptOrder(0)
{
    // Nothing to be done
}

Bool_t THaMatrixElement::match(const THaMatrixElement& rhs) const
{
    // Compare coefficients of this matrix element to another

    if (pw.size() != rhs.pw.size())
        return false;
    for (vector<Int_t>::size_type i = 0; i < pw.size(); i++) {
        if (pw[i] != rhs.pw[i]) return false;
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

ClassImp(OpticsOpt);
