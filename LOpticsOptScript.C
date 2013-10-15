#include <iostream>
#include <map>
#include <cassert>

#include "TROOT.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TVirtualFitter.h"

//#include "LOpticsOpt.h"
#include "SaveCanvas.C"

using namespace std;

class LOpticsOpt;

LOpticsOpt * opt;
UInt_t NPara = 0;
Double_t OldMatrixArray[10000] = {-99}; //NPara
Bool_t free[10000] = {kFALSE}; //NPara
//Double_t OverX = -3.46e-3;
Double_t OverX = -2.70e-3;
//Double_t OverX = -10.0;

UInt_t MaxDataPerGroup = 10000;
//UInt_t MaxDataPerGroup = 100;
TString DataSource = "Sieve.full.f51";
// TString SourceDataBase 	= "OptTh2.1/db_L.vdc.dat.ph2";
// TString DestDataBase 	= "OptTh2.1/db_L.vdc.dat.th2";
// TString DataSource = "Sieve1/Sieve.full.f51";
// TString SourceDataBase 	= "OptTh1.1/db_L.vdc.dat.ph2";
// TString DestDataBase 	= "OptTh1.1/db_L.vdc.dat.th2";
// TString DataSource = "Sieve2/Sieve.x.m4.f51";
// TString SourceDataBase 	= "OptTh2/db_L.vdc.dat.th4";
// TString DestDataBase 	= "OptTh2/db_L.vdc.dat.ph4";
//TString DataSource = "data.p5.0.dat";
//TString SourceDataBase = "OptTest/db_L.vdc.dat.ph2";
//TString DestDataBase = "OptTest/db_L.vdc.dat.th2";

typedef void (*PTRFCN)(Int_t &, Double_t *, Double_t &, Double_t*, Int_t);
PTRFCN myfcn = NULL;

void myfcn1(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
    //compute the sum of squares of dth

    assert(opt);
    assert(opt->fCurrentMatrixElems);

    opt->Array2Matrix(par);
    f = opt->SumSquareDTh();

    return;
}

void myfcn2(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
    //compute the sum of squares of dph

    assert(opt);
    assert(opt->fCurrentMatrixElems);

    opt->Array2Matrix(par);
    f = opt->SumSquareDPhi();

    return;
}

void myfcn12(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
    //compute the sum of squares of dph

    assert(opt);
    assert(opt->fCurrentMatrixElems);

    opt->Array2Matrix(par);
    f = opt->SumSquareDThUBeam();

    return;
}

void myfcn22(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
    //compute the sum of squares of dph

    assert(opt);
    assert(opt->fCurrentMatrixElems);

    opt->Array2Matrix(par);
    f = opt->SumSquareDPhiUBeam();

    return;
}

void DoMinTP(TString SourceDataBase, TString DestDataBase, UInt_t MaxDataPerGroup = 200) {
    // minimize with root

    assert(opt);
    assert(opt->fCurrentMatrixElems);

    opt->LoadDataBase(SourceDataBase);
    NPara = opt->Matrix2Array(OldMatrixArray, free);
    opt->LoadRawData(DataSource, OverX, (UInt_t) - 1, MaxDataPerGroup);
    opt->PrepareSieve();
    //opt->VerifyMatrix_Sieve();

    opt->Print();

    TVirtualFitter::SetDefaultFitter("Minuit2"); //default is Minuit
    TVirtualFitter *fitter = TVirtualFitter::Fitter(NULL, NPara);
    fitter->SetFCN(myfcn);

    for (UInt_t i = 0; i < NPara; i++) {
        Double_t absold = TMath::Abs(OldMatrixArray[i]);
        Double_t abslimit = absold > 0 ? absold * 10000 : 10000;

        fitter->SetParameter(i, Form("TMatrix%03d", i), OldMatrixArray[i], absold > 0 ? absold / 10 : 0.1, -abslimit, abslimit);
        // fitter->SetParameter(1,"asdf",0,0,0,0);

        if (!free[i]) fitter->FixParameter(i);
    }

    fitter->Print();
    cout << fitter->GetNumberFreeParameters() << " Free  / " << fitter->GetNumberTotalParameters() << " Parameters\n";

    assert(opt->fNRawData > 0);
    assert(NPara > 0);
    assert(fitter->GetNumberFreeParameters() > 0);
    assert(fitter->GetNumberTotalParameters() == NPara);

    Double_t arglist[1] = {0};
    fitter->ExecuteCommand("MIGRAD", arglist, 0);

    opt->Print();
    opt->SaveDataBase(DestDataBase);

    opt->SumSquareDTh();
    opt->SumSquareDPhi();
    TCanvas * c1 = opt->CheckSieve();
    c1->Print(DestDataBase + ".Sieve.Opt.png", "png");
    // delete fitter;
}

void PlotDataBase(TString DatabaseFileName, UInt_t MaxDataPerGroup = 1000) {
    opt = new LOpticsOpt();

    assert(opt);

    gStyle->SetOptStat(0);

    opt->LoadDataBase(DatabaseFileName);
    opt->Print();

    opt->LoadRawData(DataSource, OverX, (UInt_t) - 1, MaxDataPerGroup);

    opt->PrepareSieve();
    opt->VerifyMatrix_Sieve();

    opt->SumSquareDTh(kTRUE);
    opt->SumSquareDPhi(kTRUE);

    TCanvas * c1 = opt->CheckSieve();
    SaveCanvas(c1, DatabaseFileName + ".Sieve");

    delete opt;
}

void LOpticsOptScript(TString select, TString SourceDataBase, TString DestDataBase) {
    opt = new LOpticsOpt();

    map<TString, Int_t> selectmap;
    selectmap["theta"] = 1;
    selectmap["phi"] = 2;
    selectmap["delta"] = 3;
    selectmap["y"] = 4;
    selectmap["theta2"] = 12;
    selectmap["phi2"] = 22;

    gStyle->SetOptStat(0);

    switch (selectmap[select]) {
    case 1:
        cout << "Optimizing for Theta\n";
        myfcn = myfcn1;
        opt->fCurrentMatrixElems = &(opt->fTMatrixElems);
        DoMinTP(SourceDataBase, DestDataBase, 500);
        break;
    case 2:
        cout << "Optimizing for Phi\n";
        myfcn = myfcn2;
        opt->fCurrentMatrixElems = &(opt->fPMatrixElems);
        DoMinTP(SourceDataBase, DestDataBase, 500);
        break;
    case 3:
        cout << "Optimizing for Delta\n";
        break;
    case 4:
        cout << "Optimizing for Y\n";
        break;
    case 12:
        cout << "Optimizing for Theta\n";
        myfcn = myfcn12;
        opt->fCurrentMatrixElems = &(opt->fTMatrixElems);
        DoMinTP(SourceDataBase, DestDataBase, 500);
        break;
    case 22:
        cout << "Optimizing for Phi\n";
        myfcn = myfcn22;
        opt->fCurrentMatrixElems = &(opt->fPMatrixElems);
        DoMinTP(SourceDataBase, DestDataBase, 500);
        break;
    default:
        break;
    }

    gSystem->Exec(Form("cp -vf %s %s.source", SourceDataBase.Data(), DestDataBase.Data()));
    //gSystem->Exec(Form("cp -vf LOpticsOptScript.log %s.log", DestDataBase.Data()));

    delete opt;

    return;
}
