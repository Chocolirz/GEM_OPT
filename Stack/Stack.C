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
#include <fstream>
#include <string>
#include <iomanip>
#include <algorithm>

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

    std::ifstream infile("../Ne_r.txt");
    if (!infile) {
        std::cerr << "Error opening file\n";
        return 1;
    }
    
    std::vector<double> startingRadius;
    std::vector<double> probDensity;

    double val1, val2;
    while (infile >> val1 >> val2) {
        startingRadius.push_back(val1 * 1.e-4);
        probDensity.push_back(val2);
    }

    // Load field map.
    ComponentComsol fm;
    fm.Initialise("mesh_triple_LH_510_corrected_inv.mphtxt", "dielectrics.dat", "field_triple_LH_510_corrected_inv.txt", "mm");
    fm.EnableMirrorPeriodicityX(); 
    fm.EnableMirrorPeriodicityY(); 
    fm.PrintRange(); 

    // Dimensions of the GEM [cm]
    // Note that the geometry has changed, in ParallelPlate we have
    // y as the axis perpenticular to the plates, but here we have z. 
    const double rlim = 0.0210 / 2.; // cm, radius of the GEM hole
    const double upperGap = 0.2; // cm
    const double lowerGap = 0.2; // cm
    const double copperThick = 0.0002; // cm, thickness of the copper plates
    const double gap = 0.057; // cm, GEM is 570 um thick
    const double pitch = 0.03; // distance between GEM holes
    const double xmin = -10. * pitch; 
    const double xmax = 10. * pitch; // cm
    const double ymin = -10. * pitch;
    const double ymax = 10. * pitch; // cm
    const double zmin = 0.;
    const double zmax = 3. * gap + 3. * lowerGap + upperGap + (3. * 2. * copperThick); // cm, three GEMs
    const double tsg = 2. * gap + 2. * lowerGap + 4. * copperThick; // Height of the top of the second GEM.
    
    // Generate lattice
    std::vector<std::pair<double, double>> positions;

    // Fill positions with your hex lattice generation logic...
    double a = 0.0280;
    double R_max = 4 * a;
    int max_rows = static_cast<int>(R_max / ((std::sqrt(3)/2) * a)) + 2;
    int max_cols = static_cast<int>(R_max / a) + 2;

    for (int row = -max_rows; row <= max_rows; ++row) {
        double y = row * (std::sqrt(3)/2) * a;
        for (int col = -max_cols; col <= max_cols; ++col) {
            double x = col * a + ((row % 2 != 0) ? a/2 : 0.0);
            double r = std::hypot(x, y);
            if (r <= R_max + 1e-12) {
                positions.emplace_back(x + 0.5*pitch, y);
            }
        }
    }
    
    // Setup the gas.
    MediumMagboltz gas("cf4", 100.); // Material and percentage
    gas.SetTemperature(293.15); // Temperature in Kelvin
    gas.SetPressure(50.); // Pressure in Torr
    gas.LoadIonMobility("IonMobility_CF4+_CF4.txt"); // Ion transport properties
    //gas.SetMaxElectronEnergy(200.); // Electron energy in eV
    gas.EnableCrossSectionOutput();
    gas.Initialise(true);
    // Set gas and see its properties
    fm.SetGas(&gas);
    fm.PrintMaterials();

    // Set up sensor
    Sensor sensor(&fm);
    sensor.SetArea(xmin, ymin, tsg, xmax, ymax, tsg + lowerGap);

    AvalancheMicroscopic aval(&sensor); // Microscopic tracking of electrons
    aval.SetCollisionSteps(100);

    AvalancheMC drift(&sensor);
    drift.SetDistanceSteps(2.e-4);

    // Right now we are not going to consider ions, can add later. 
    
    // Visualising drift lines
    ViewDrift driftView;
    aval.EnablePlotting(&driftView);
    drift.EnablePlotting(&driftView);




    // NUMBER OF AVALANCHES SIMULATED
    constexpr unsigned int nEvents = 3e6;



    
    // Set initial position [cm], time [ns] and energy [eV] for electrons
    double x0, y0, z0 = tsg + lowerGap, t0 = 0., phi0, r0;
    const double e0 = 0.1;
    const double dx0 = 0., dy0 = 0., dz0 = -1; // vector length == 1 or 0 (direction radomised). 
    // Prepare some constants.
    double anode_count = 0, total_e = 0, rz_s, rz_e;

    // Add TH1F plots here
    //TH2F hDisXY("hDisXY", "Entries;x [cm];y [cm]", 100, -5*pitch, 5*pitch, 100, -5*pitch, 5*pitch);
    //TH2F hDisRZ("hDisRZ", "Entries;R [cm];Z [cm]", 100, 0, 5*pitch, 100, tsg, tsg + lowerGap);
    TH1F hR("hR", "R distribution;R [cm];Entries", 220, 0, rlim + 0.002);

    
    // Plot the avalanche
    constexpr bool plotAvalanche = true;
    if (plotAvalanche) {
        for (unsigned int i = 0; i < 105; i++) {
            for (unsigned int j = 0; j <= std::round(nEvents*probDensity[i]); j++) {
                std::cout << "r = " << i << " um: " << j << "/" << std::round(nEvents*probDensity[i]) << "\n"; // Manual set counter

                r0 = startingRadius[i];
                phi0 = 2. * M_PI * Garfield::RndmUniform();
                x0 = r0 * cos(phi0);
                y0 = r0 * sin(phi0);
                std::cout << "x0 = " << x0*1.e4 << " um ; y0 = " << y0*1.e4 << " um. Starting avalanche..." << std::endl;

                aval.AvalancheElectron(x0, y0, z0, t0, e0, dx0, dy0, dz0);

                for (auto & electron : aval.GetElectrons()) {
                    auto& p0 = electron.path.front();
                    auto& p1 = electron.path.back();

                    total_e += 1;

                    
                    double rz_e = rlim * 2;
                    for (const auto& [x, y] : positions) {
                        double dist = sqrt(pow(p1.x - x, 2) + pow(p1.y - y, 2));
                        if (dist < rz_e) {
                            rz_e = dist;
                        }
                    }
                    

                    //double r1 = sqrt(pow((p1.x - 2.5*pitch),2) + pow((p1.y - 1.5*sqrt(3)*pitch),2));
                    //double r2 = sqrt(pow((p1.x + 2.5*pitch),2) + pow((p1.y + 1.5*sqrt(3)*pitch),2));
                    //double r3 = sqrt(pow((p1.x - 2.5*pitch),2) + pow((p1.y + 1.5*sqrt(3)*pitch),2));
                    //double r4 = sqrt(pow((p1.x + 2.5*pitch),2) + pow((p1.y - 1.5*sqrt(3)*pitch),2));
                    //rz_e = std::min({r1, r2, r3, r4});

                    //hDisRZ.Fill(rz_e, p1.z, 1./rz_e);
                    
                    //rz_e = sqrt(pow(p1.x, 2) + pow(p1.y, 2));

                    // Sort the electrons based on their end points
                    if (p1.z < tsg + 0.005 && electron.status != -7 && rz_e < rlim) {
                        anode_count += 1;
                        //hDisXY.Fill(p1.x, p1.y);
                        hR.Fill(rz_e);
                        //hR.Fill(rz_e, 1./rz_e); // We want the actual number, not distribution, don't scale.
                    }
                }
                std::cout << anode_count/total_e << std::endl;
            }
        }

        std::cout << "Generating plots ..." << std::endl;
        /*
        TCanvas* cDisXY = new TCanvas("cDisXY", "Distribution of Electrons at the Exit", 800, 800);
        cDisXY->SetLogz();
        hDisXY.SetMinimum(9e-1);
        hDisXY.Draw("COLZ");
        gPad->SetRightMargin(0.15);
        cDisXY->Update();
        */
        TCanvas* cR = new TCanvas("cR", "Distribution wrt R", 800, 800);
        //cR->SetLogy();
        hR.Draw();
        cR->Update();
        /*
        TCanvas* cDisRZ = new TCanvas("cDisRZ", "Distribution of Electrons at the Exit", 800, 800);
        cDisRZ->SetLogz();
        hDisRZ.SetMinimum(9e-1);
        hDisRZ.Draw("COLZ");
        gPad->SetRightMargin(0.15);
        cDisRZ->Update();
        

        int nbins = hR.GetNbinsX();
        for (int i = 1; i <= nbins; i++) {
            double entries = hR.GetBinContent(i);

            std::cout << "Entries in this bin: " << entries << std::endl;
        }
        */

        std::cout << "-----------------------------------------" << std::endl;
        std::cout << "Number of electrons in total: " << anode_count << std::endl;

        std::cout << "Plots completed. " << std::endl;
    }

    app.Run(kTRUE);
}