#include <TApplication.h>
#include <TCanvas.h>
#include <TROOT.h>

#include "Garfield/MediumMagboltz.hh"

using namespace std;
using namespace Garfield;

int main(int argc, char* argv[]) {
    TApplication app("app", &argc, argv);
    /*
    MediumMagboltz gas;
    const size_t nE = 10;
    const double emin = 8000.; // V/cm
    const double emax = 12000.;
    // Flag to request logarithmic spacing.
    constexpr bool useLog = true;
    const size_t nB = 1;
    // Range of magnetic fields [Tesla]
    const double bmin = 0.;
    const double bmax = 0.;
    const size_t nA = 5;
    // Range of angles [rad]
    const double amin = 0.;
    const double amax = HalfPi;
    gas.SetFieldGrid(emin, emax, nE, true);

    gas.SetComposition("cf4", 100.);
    gas.SetTemperature(293.15);
    gas.SetPressure(50.);

    const int ncoll = 10;
    // Run Magboltz to generate the gas table.
    gas.GenerateGasTable(ncoll);

    gas.WriteGasFile("cf4_100.gas");
    */
    
    MediumMagboltz gas;
    gas.LoadGasFile("cf4_100.gas");
    gas.PrintGas();

    TCanvas c1("c1", "", 800, 800);
    gas.PlotTownsend("e", &c1);

    TCanvas c2("c2", "", 800, 800);
    gas.PlotDiffusion("e", &c2);

    app.Run(kTRUE);
}