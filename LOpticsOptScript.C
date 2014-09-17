#include <iostream>
#include <cassert>

#include "TROOT.h"
#include "TCanvas.h"
#include "TError.h"
#include "TMath.h"
#include "TString.h"
#include "TVirtualFitter.h"

//#include "LOpticsOpt.h"
//#include "SaveCanvas.C"

using namespace std;

class LOpticsOpt;

LOpticsOpt * opt;
UInt_t NPara = 0;
Double_t OldMatrixArray[10000] = {-99}; //NPara
Bool_t free[10000] = {kFALSE}; //NPara
Int_t UseFPOff = 0;

UInt_t MaxDataPerGroup = 100000;
TString DataSource;
TString DataGroup = "1";
TString DataSourceSuf = ".full.f51";

typedef void (*PTRFCN)(Int_t &, Double_t *, Double_t &, Double_t*, Int_t);
PTRFCN myfcn = NULL;

void myfcn1(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t)
{
    // compute the sum of squares of dth

    assert(opt);
    assert(opt->fCurrentMatrixElems);

    opt->Array2Matrix(par, UseFPOff);
    f = opt->SumSquareDTh(UseFPOff);

    return;
}

void myfcn2(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t)
{
    // compute the sum of squares of dph

    assert(opt);
    assert(opt->fCurrentMatrixElems);

    opt->Array2Matrix(par, UseFPOff);
    f = opt->SumSquareDPhi(UseFPOff);

    return;
}

void myfcn3(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t)
{
    // compute the sum of squares of dy

    assert(opt);
    assert(opt->fCurrentMatrixElems);

    opt->Array2Matrix(par, UseFPOff);
    f = opt->SumSquareDY(UseFPOff);

    return;
}

void myfcn4(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t)
{
    // compute the sum of squares of ddp

    assert(opt);
    assert(opt->fCurrentMatrixElems);

    opt->Array2Matrix(par, UseFPOff);
    f = opt->SumSquareDp(UseFPOff);

    return;
}

void myfcn10(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t)
{
    //compute the sum of squares of dph

    assert(opt);
    assert(opt->fCurrentMatrixElems);

    opt->Array2Matrix(par, 1);
    Double_t t = opt->SumSquareDTh(1);
    Double_t p = opt->SumSquareDPhi(1);
    f = sqrt(t * t + p * p);

    return;
}

void DoMinTP(TString SourceDataBase, TString DestDataBase, UInt_t MaxDataPerGroup = 200)
{
    // minimize with root

    assert(opt);
    assert(opt->fCurrentMatrixElems);

    opt->LoadDataBase(SourceDataBase);
    NPara = opt->Matrix2Array(OldMatrixArray, free, UseFPOff);
    opt->LoadRawData(DataSource, (UInt_t) - 1, MaxDataPerGroup);
    opt->PrepareSieve();

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

    opt->SumSquareDTh(UseFPOff);
    opt->SumSquareDPhi(UseFPOff);

    TCanvas * c1 = opt->CheckSieve(1);
    c1->Print(DestDataBase + ".Sieve.Opt.png", "png");

    delete fitter;
}

void DoMinYt(TString SourceDataBase, TString DestDataBase, UInt_t MaxDataPerGroup = 200)
{
    // minimize with root

    assert(opt);
    assert(opt->fCurrentMatrixElems);

    opt->LoadDataBase(SourceDataBase);
    NPara = opt->Matrix2Array(OldMatrixArray, free, UseFPOff);
    opt->LoadRawData(DataSource, (UInt_t) - 1, MaxDataPerGroup);
    opt->PrepareVertex();

    //compensate bias due to dp event selections
    // 	opt->fArbitaryYShift[0] = 2.786177e-05;
    // 	opt->fArbitaryYShift[1] = 8.168538e-05;
    // 	opt->fArbitaryYShift[2] = 5.299596e-05;
    // 	opt->fArbitaryYShift[3] = 3.175602e-05;
    // 	opt->fArbitaryYShift[4] = 9.519830e-05;

    opt->fArbitaryYShift[0] = 0.;
    opt->fArbitaryYShift[1] = 0.;
    opt->fArbitaryYShift[2] = 0.;
    opt->fArbitaryYShift[3] = 0.;
    opt->fArbitaryYShift[4] = 0.;

    opt->Print();

    TVirtualFitter::SetDefaultFitter("Minuit2"); //default is Minuit
    TVirtualFitter *fitter = TVirtualFitter::Fitter(NULL, NPara);
    fitter->SetFCN(myfcn);

    for (UInt_t i = 0; i < NPara; i++) {
        Double_t absold = TMath::Abs(OldMatrixArray[i]);
        Double_t abslimit = absold > 0 ? absold * 10000 : 10000;

        fitter->SetParameter(i, Form("TMatrix%03d", i), OldMatrixArray[i], absold > 0 ? absold / 1000 : 0.1, -abslimit, abslimit);
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

    opt->SumSquareDY(UseFPOff);

    TCanvas * c1 = opt->CheckY();
    c1->Print(DestDataBase + ".Yt.Opt.png", "png");

    delete fitter;
}

void DoMinDp(TString SourceDataBase, TString DestDataBase, UInt_t MaxDataPerGroup = 200)
{
    // minimize with root

    assert(opt);
    assert(opt->fCurrentMatrixElems);

    opt->fCurrentMatrixElems = &(opt->fDMatrixElems);

    opt->LoadDataBase(SourceDataBase);
    NPara = opt->Matrix2Array(OldMatrixArray, free, UseFPOff);
    opt->LoadRawData(DataSource, (UInt_t) - 1, MaxDataPerGroup);
    opt->PrepareDp();

    //compensate bias due to dp event selections
    // 	opt->fArbitaryDpKinShift[0] = 2.786177e-05;
    // 	opt->fArbitaryDpKinShift[1] = 8.168538e-05;
    // 	opt->fArbitaryDpKinShift[2] = 5.299596e-05;
    // 	opt->fArbitaryDpKinShift[3] = 3.175602e-05;
    // 	opt->fArbitaryDpKinShift[4] = 9.519830e-05;

    opt->fArbitaryDpKinShift[0] = 0.;
    opt->fArbitaryDpKinShift[1] = 0.;
    opt->fArbitaryDpKinShift[2] = 0.;
    opt->fArbitaryDpKinShift[3] = 0.;
    opt->fArbitaryDpKinShift[4] = 0.;

    opt->Print();

    TVirtualFitter::SetDefaultFitter("Minuit2"); //default is Minuit
    TVirtualFitter *fitter = TVirtualFitter::Fitter(NULL, NPara);
    fitter->SetFCN(myfcn);

    for (UInt_t i = 0; i < NPara; i++) {
        Double_t absold = TMath::Abs(OldMatrixArray[i]);
        Double_t abslimit = absold > 0 ? absold * 10000 : 10000;

        fitter->SetParameter(i, Form("TMatrix%03d", i), OldMatrixArray[i], absold > 0 ? absold / 1000 : 0.1, -abslimit, abslimit);
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

    opt->SumSquareDp(UseFPOff);

    TCanvas * c1 = opt->CheckDp();
    c1->Print(DestDataBase + ".Dp.Opt.png", "png");

    delete fitter;
}

void PlotDataBase(TString DatabaseFileName, UInt_t MaxDataPerGroup = 1000)
{
    opt = new LOpticsOpt();

    assert(opt);

    gStyle->SetOptStat(0);

    opt->LoadDataBase(DatabaseFileName);
    opt->Print();

    opt->LoadRawData(DataSource, (UInt_t) - 1, MaxDataPerGroup);

    opt->PrepareSieve();
    //opt->PrepareDp();

    // opt->SumSquareDTh(kTRUE);
    // opt->SumSquareDPhi(kTRUE);
    opt->SumSquareDp(kTRUE);

    TCanvas * c1 = opt->CheckSieve();
    c1->Print(DatabaseFileName + ".Sieve.png", "png");
    //TCanvas * c1 = opt->CheckDpGlobal();
    //c1->Print(DatabaseFileName + ".Dp.png","png");

    delete opt;
}

void LOpticsOptScript(TString select, TString SourceDataBase, TString DestDataBase)
{
    opt = new LOpticsOpt();

    Int_t s = 0;
    if (select == "theta") s = 1;
    if (select == "phi") s = 2;
    if (select == "y") s = 3;
    if (select == "delta") s = 4;
    if (select == "off") s = 10;

    gStyle->SetOptStat(0);

    switch (s) {
    case 1:
        Info("LOpticsOptScript", "Optimizing for Theta");
        DataSource = "Sieve" + DataGroup + "/Sieve" + DataSourceSuf;
        myfcn = myfcn1;
        opt->fCurrentMatrixElems = &(opt->fTMatrixElems);
        DoMinTP(SourceDataBase, DestDataBase, 500);
        break;
    case 2:
        Info("LOpticsOptScript", "Optimizing for Phi");
        DataSource = "Sieve" + DataGroup + "/Sieve" + DataSourceSuf;
        myfcn = myfcn2;
        opt->fCurrentMatrixElems = &(opt->fPMatrixElems);
        DoMinTP(SourceDataBase, DestDataBase, 500);
        break;
    case 3:
        Info("LOpticsOptScript", "Optimizing for Y");
        DataSource = "Yt" + DataGroup + "/Yt" + DataSourceSuf;
        myfcn = myfcn3;
        opt->fCurrentMatrixElems = &(opt->fYMatrixElems);
        DoMinYt(SourceDataBase, DestDataBase, 5000);
        break;
    case 4:
        Info("LOpticsOptScript", "Optimizing for Delta");
        DataSource = "Dp" + DataGroup + "/Dp" + DataSourceSuf;
        myfcn = myfcn4;
        opt->fCurrentMatrixElems = &(opt->fDMatrixElems);
        DoMinDp(SourceDataBase, DestDataBase, 5000);
        break;
    case 10:
        Info("LOpticsOptScript", "Optimizing for Focus Plane Offset");
        DataSource = "Sieve" + DataGroup + "/Sieve" + DataSourceSuf;
        myfcn = myfcn10;
        UseFPOff = 0;
        opt->fCurrentMatrixElems = &(opt->fFPMatrixElems);
        DoMinTP(SourceDataBase, DestDataBase, 500);
        break;
    default:
        break;
    }

    gSystem->Exec(Form("cp -vf %s %s.source", SourceDataBase.Data(), DestDataBase.Data()));
    gSystem->Exec(Form("cp -vf log %s.log", DestDataBase.Data()));

    delete opt;

    return;
}
