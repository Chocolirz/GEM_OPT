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
double totalPhoton = 0; 

std::ofstream outfile;

void userHandle(double x, double y, double z, double t, int type, int level, Medium* m, double e0, double e1, double dx0, double dy0, double dz0, double dx1, double dy1, double dz1) {
    // Check if the collision is an excitation or ionisation
    // Generate a photon using ROOT to a random direction
    // Emit photon in a random direction


    const double theta = M_PI * Garfield::RndmUniform();
    const double phi = 2. * M_PI * Garfield::RndmUniform();
    const double phiDiff = phi - atan(y / x);

    bool scint = (std::rand() % 100) < 33;

    // Check if the photon can reach the anode
    double r = sqrt(x*x + y*y);
    double rDiff = rlim - r;
    double zDiff = z - (tsg + lowerGap);
    double lr = -r*cos(phiDiff) + sqrt((rlim * rlim) - (r*r*sin(phiDiff)*sin(phiDiff)));
    if (scint && level > 0 && level < 8) {
        if (zDiff > 0 && theta > M_PI_2 && tan(theta) > -lr/zDiff) {
            // Write output file if scintillates
            outfile << x << "\t" << y << "\t" << z << "\t" << theta << "\t" << phi << "\n";
        } else if (zDiff < 0 && theta > M_PI_2) {
            // same
            outfile << x << "\t" << y << "\t" << z << "\t" << theta << "\t" << phi << "\n";
        }
        totalPhoton += 1;
    }
}

int main(int argc, char* argv[]) {
    TApplication app("app", &argc, argv);
    
    // Prepare output file
    outfile.open("scint_positions.txt");
    if (!outfile.is_open()) {
        std::cerr << "Failed to open output file!" << std::endl;
        return 1;
    }

    outfile << std::fixed << std::setprecision(12);
    outfile << "# x y z theta phi\n";

    // Load field map.
    ComponentComsol fm;
    fm.Initialise("mesh_triple_LH_510_corrected.mphtxt", "dielectrics.dat", "field_triple_LH_510_corrected.txt", "mm");
    fm.EnableMirrorPeriodicityX(); 
    fm.EnableMirrorPeriodicityY(); 
    fm.PrintRange(); 
    
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
    sensor.SetArea(xmin, ymin, tsg, xmax, ymax, tsg + lowerGap + gap);

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
    constexpr unsigned int nEvents = 1000;



    
    // Set initial position [cm], time [ns] and energy [eV] for electrons
    double x0, y0, z0 = tsg + lowerGap + gap, t0 = 0., phi0, r0;
    const double e0 = 0.1;
    const double dx0 = 0., dy0 = 0., dz0 = -1; // vector length == 1 or 0 (direction radomised). 
    // Prepare some constants.
    double anode_count = 0, total_e = 0, rz_s, rz_e;

    // Add TH1F plots here
    //TH2F hDisXY("hDisXY", "Entries;x [cm];y [cm]", 100, -5*pitch, 5*pitch, 100, -5*pitch, 5*pitch);
    //TH2F hDisRZ("hDisRZ", "Entries;R [cm];Z [cm]", 100, 0, 5*pitch, 100, tsg, tsg + lowerGap);
    //TH1F hR("hR", "R distribution;R [cm];Entries", 220, 0, rlim + 0.002);

    
    // Plot the avalanche
    constexpr bool plotAvalanche = true;
    if (plotAvalanche) {
        for (unsigned int i = 0; i < nEvents; i++) {
            std::cout << i << "/" << nEvents << std::endl;
            r0 = rlim * Garfield::RndmUniform();
            phi0 = 2. * M_PI * Garfield::RndmUniform();
            x0 = r0 * cos(phi0);
            y0 = r0 * sin(phi0);
            std::cout << "x0 = " << x0*1.e4 << " um ; y0 = " << y0*1.e4 << " um. Starting avalanche..." << std::endl;

            aval.SetUserHandleCollision(userHandle);
            aval.AvalancheElectron(x0, y0, z0, t0, e0, dx0, dy0, dz0);
        }

        std::cout << "Generating plots ..." << std::endl;

        std::cout << "-----------------------------------------" << std::endl;
        //std::cout << "Number of electrons in total: " << anode_count << std::endl;

        std::cout << "Plots completed. " << std::endl;
    }

    outfile.close();
    std::cout << "Positions wrote to file. " << std::endl;

    // Things wrote in file are x, y, z, theta, phi for photons that CAN reach the top of second GEM plane. 

    app.Run();
}