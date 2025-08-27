#include <TApplication.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/ComponentComsol.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewMedium.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/SolidTube.hh"
#include "Garfield/ViewCell.hh"

using namespace Garfield;
using namespace std;

int main(int argc, char* argv[]) {
    TApplication app("app", &argc, argv);

    // Load field map.
    ComponentComsol fm;
    fm.Initialise("mesh_single_gem_560_210_1cm.mphtxt", "dielectrics.dat", "field_single_gem_560_210_1cm.txt", "mm");
    fm.EnableMirrorPeriodicityX(); // Generating infinite array of cells
    fm.EnableMirrorPeriodicityY(); 
    fm.PrintRange(); 

    // Dimensions of the GEM [cm]
    // We have 1cm drift region on top, then the usual 570um + 2*2um GEM, and 2mm induction region below
    const double rlim = 0.0210 / 2.; // cm, radius of the GEM hole
    const double upperGap = 1.; // cm
    const double lowerGap = 0.2; // cm
    const double copperThick = 0.0002; // cm, thickness of the copper plates
    const double gap = 0.057; // cm, GEM is 570 um thick
    const double pitch = 0.03; // distance between GEM holes
    const double xmin = -10. * pitch; 
    const double xmax = 10. * pitch; // cm
    const double ymin = -10. * pitch;
    const double ymax = 10. * pitch; // cm
    const double zmin = 0.;
    const double zmax = gap + lowerGap + upperGap + (2 * copperThick); // cm

    // Setup the gas.
    MediumMagboltz gas("cf4", 100.); // Material and percentage
    gas.SetTemperature(293.15); // Temperature in Kelvin
    gas.SetPressure(50.); // Pressure in Torr
    gas.EnableCrossSectionOutput();
    gas.Initialise(true);
    // Set gas and see its properties
    fm.SetGas(&gas);
    fm.PrintMaterials();

    // Set up sensor
    Sensor sensor(&fm);
    sensor.SetArea(xmin, ymin, lowerGap + gap, xmax, ymax, zmax);

    AvalancheMicroscopic aval(&sensor); // Microscopic tracking of electrons
    aval.SetCollisionSteps(100);

    AvalancheMC drift(&sensor);
    drift.SetDistanceSteps(2.e-4);




    // NUMBER OF AVALANCHES SIMULATED
    constexpr unsigned int nEvents = 3e4;



    
    // Set initial position [cm], time [ns] and energy [eV] for electrons
    const double x0 = 0.00001, y0 = 0.00001, z0 = zmax - 0.8, t0 = 0.;
    const double e0 = 0.1;
    const double dx0 = 0., dy0 = 0., dz0 = 0.; // vector length == 1 or 0 (direction radomised). 
    // Prepare some constants and plots
    double total_e = 0.;
    double coll_count = 0.; 

    TH2F hColl("hColl", "Entries;x [cm];y [cm]", 200, -3.*pitch, 3.*pitch, 200, -3.*pitch, 3.*pitch);
    
    // Plot the avalanche
    constexpr bool plotAvalanche = true;
    if (plotAvalanche) {
        for (unsigned int i = 0; i < nEvents; i++) {
            std::cout << i << "/" << nEvents << "\n"; // Manual set counter

            aval.AvalancheElectron(x0, y0, z0, t0, e0, dx0, dy0, dz0);

            for (auto & electron : aval.GetElectrons()) {
                auto& p0 = electron.path.front();
                auto& p1 = electron.path.back();

                total_e = total_e + 1;
                
                if (p1.z <= lowerGap + gap + 0.0005 && electron.status != -7) {
                    coll_count = coll_count + 1;
                    hColl.Fill(p1.x, p1.y);
                }
            }
        }

        TCanvas* cColl = new TCanvas("cColl", "Number of electrons collected", 800, 800);
        cColl->SetLogz();
        hColl.SetMinimum(9e-1);
        hColl.Draw("COLZ");
        gPad->SetRightMargin(0.15);
        cColl->Update();

        std::cout << "Total number of electrons: " << total_e << ", Collected electrons: " << coll_count << ". " << std::endl;
    }

    app.Run(kTRUE);
}