// Description of G2P optics

#ifndef ROOT_Input
#define ROOT_Input

#include "TROOT.h"
#include "TMath.h"

using namespace std;
const TString InputID = "G2P_LHRS";

/////////////////////////////////////////////////////////////////////////
// Input Sections
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
// HRS Position Inputs
const Double_t D2R = TMath::Pi() / 180.;

// G2P: set as central sieve hole angle
const Double_t HRSAngle = 5.785 * D2R; //WARNING: check

// MissPoint* are in HCS
const Double_t MissPointZ = 0.0; //mm
const Double_t MissPointY = 0.0; //mm

const Double_t BeamShiftX = 0.0;

/////////////////////////////////////////////////////////////////////////
// Sieve Position Inputs
const Double_t SieveYbyCol[] = {3 * 6.1214e-3, 2 * 6.1214e-3, 1 * 6.1214e-3, 0.0, -1 * 4.7752e-3, -2 * 4.7752e-3, -3 * 4.7752e-3, 1e36};
const UInt_t NSieveCol = 7; /*WARNING:7*/

const Double_t SieveXbyRow[] = {-3 * 13.3096e-3, -2 * 13.3096e-3, -1 * 13.3096e-3, 0.0, 1 * 13.3096e-3, 2 * 13.3096e-3, 3 * 13.3096e-3, 1e36};
const UInt_t NSieveRow = 7; /*WARNING:7*/

// SieveOff* are in TCS
const Double_t SieveOffY = 0.; //(2.97*cos(HRSAngle)).
const Double_t SieveOffX = 0.; //-(2.61)
const Double_t ZPos = 799.60 * 1e-3; //800.

/////////////////////////////////////////////////////////////////////////
// Vertex Position Inputs
static const UInt_t NFoils = 1; //WARNING: check
const Double_t targetfoils[] = {-13.6271 * 1e-3, 1e36}; // At end-cap, not target ladder -10.81mm

/////////////////////////////////////////////////////////////////////////
// Excitation State Inputs
const UInt_t NKine = 5; //N Delta Scans

#define DIPOLE_MAG2MOM(Mag) (2.702*(Mag)-1.6e-03*(Mag)*(Mag)*(Mag))
const Double_t HRSCentralMom[] = {
                                  DIPOLE_MAG2MOM((0.8580431 + 0.8580438) / 2.), // -3% run 3181 3182
                                  DIPOLE_MAG2MOM((0.8472928 + 0.8472905) / 2.), // -2% run 3183 3184
                                  DIPOLE_MAG2MOM((0.8328698 + 0.8328679) / 2.), //  0% run 3185 3186
                                  DIPOLE_MAG2MOM((0.8162305 + 0.8162268) / 2.), //  2% run 3187 3188
                                  DIPOLE_MAG2MOM((0.8077355 + 0.8077357) / 2.), //  3% run 3190 3191
                                  0
};

const Double_t GroundNuclearMass = 12.0107 * .931494028; //GeV/c^2  //Carbon Target
const Double_t ExcitationEnergy[] = //selected excitation states for each kinematics
{0., 0., 0., 0., 0.};
//{0.,0.00443891,0.00443891,0.00443891,0.00443891};

const UInt_t NExcitationStates = 1; //C Excitation States
const Double_t ExcitationEnergyList[]
        = {0};

/* const UInt_t NExcitationStates = 8; //C Excitation States */
/* const Double_t ExcitationEnergyList[] */
/* ={0,0.00443891,0.00765420,0.009641,0.010844,0.011160,0.011828,0.012710}; */

/////////////////////////////////////////////////////////////////////////
// Radiation Loss Inputs

// g2p straight through
const Double_t AllLossExceptFoil
        = 1.e-3 * (//in MeV
        // 0.870100332
        1.009711
        );
const Double_t LossEachFoil = 0.501422e-3; // 40mil C12

// Array of FoilID
const Double_t RadiationLossByFoil[] = {
                                        AllLossExceptFoil + LossEachFoil * 1,
                                        0
};
// Warning: these numbers are calculated with small angle approximation

/////////////////////////////////////////////////////////////////////////
// Disable Extended Target Correction

const Double_t ExtTarCor_ThetaCorr = 0;
const Double_t ExtTarCor_DeltaCorr = 1e36;

/////////////////////////////////////////////////////////////////////////
// Database header

const char * Prefix = "L.vdc.";

const char * DatabaseHeader = "\
[ L.global ]   \n\
0.3327 1 0.0 270.2 0.0 -1.6e-03        VDC Angle, Plane Spacing, Gamma Coefficents  \n\
matrix elements   \n\
t 0 0 0  -1.001135e+00 -3.313373e-01 -4.290819e-02  4.470852e-03  0.000000e+00  0.000000e+00  0.000000e+00  0\n\
y 0 0 0  -8.060915e-03  1.071977e-03  9.019102e-04 -3.239615e-04  0.000000e+00  0.000000e+00  0.000000e+00  0\n\
p 0 0 0  -2.861912e-03 -2.469069e-03  8.427172e-03  2.274635e-03  0.000000e+00  0.000000e+00  0.000000e+00  0\n\
";

const char * DatabaseFooter = "\
L 0 0 0 0  25.713\n\
L 1 0 0 0  0.1650\n\
L 2 0 0 0 -0.05\n\
L 0 1 0 0 -11.6554\n\
L 0 2 0 0 -9.4951\n\
L 0 0 1 0  0.0\n\
L 0 0 2 0  0.0\n\
L 0 0 0 1  0.0\n\
L 0 0 0 2  0.0\n\
XF 1 0 0 0 0 -2.181E+00\n\
XF 0 1 0 0 0 -1.980E-01\n\
XF 0 0 0 0 1  1.191E+01\n\
TF 1 0 0 0 0 -1.000E-01\n\
TF 0 1 0 0 0 -4.690E-01\n\
TF 0 0 0 0 1  1.967E+00\n\
PF 0 0 1 0 0  3.630E-01\n\
PF 0 0 0 1 0 -0.902E+00\n\
YF 0 0 1 0 0 -5.950E-01\n\
YF 0 0 0 1 0 -1.274E+00\n\
";
/////////////////////////////////////////////////////////////////////////

#endif
