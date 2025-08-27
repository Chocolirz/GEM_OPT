#include <TApplication.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>

#include "Garfield/ComponentComsol.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/MediumMagboltz.hh"

using namespace Garfield;
using namespace std;

int main(int argc, char* argv[]) {
    TApplication app("app", &argc, argv);

    // Setup the gas.
    MediumMagboltz gas("cf4", 100.); // Material and percentage
    gas.SetTemperature(293.15); // Temperature in Kelvin
    gas.SetPressure(50.); // Pressure in Torr
    gas.EnableCrossSectionOutput();
    gas.Initialise(true);
    // Set gas and see its properties
    

    // Dimensions of the GEM [cm]
    // Note that the geometry has changed, in ParallelPlate we have
    // y as the axis perpenticular to the plates, but here we have z. 
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
    const double zmax = gap + lowerGap + upperGap + (2 * copperThick); // cm
    const double nSteps = 1.e3; // 1e4 if r, 1e6 if z
    const double dz = zmax / nSteps;

    for (int i = 0, v = 540; v <= 600; v += 10, ++i) {
        std::stringstream ss_mesh, ss_field;
        ss_mesh << "/home/lirui/garfield++/TempTest/GemWithRim/mesh_single_gem_" << v << "_210.mphtxt";
        ss_field << "/home/lirui/garfield++/TempTest/GemWithRim/field_single_gem_" << v << "_210.txt";

        auto fm = new ComponentComsol();
        fm->Initialise(ss_mesh.str(), "/home/lirui/garfield++/TempTest/GemWithRim/dielectrics.dat", ss_field.str(), "mm");
        fm->EnableMirrorPeriodicityX();
        fm->EnableMirrorPeriodicityY();
        fm->SetGas(&gas);

        std::ofstream outfile("efield_z_" + std::to_string(v) + "V_210um.txt");
        //std::ofstream outfile("efield_z_edge_" + std::to_string(v) + "um.txt");
        //std::ofstream outfile("efield_r_" + std::to_string(v) + "um.txt");
        if (!outfile.is_open()) {
            std::cerr << "Failed to open output file!" << std::endl;
            return 1;
        }

        outfile << std::fixed << std::setprecision(12);
        outfile << "# z Ex Ey Ez\n";
        double radius = v * 1.e-4 / (2.*sqrt(2));
        //const double dr = radius / nSteps;

        //for (double r = -radius; r <= radius; r += dr) {
        for (double z = 0.; z <= zmax; z += dz) {
            //auto E = fm->ElectricField(r, r, zmax / 2.);
            auto E = fm->ElectricField(0.000001, 0.000001, z);

            outfile << z << " " << E[0] << " " << E[1] << " " << E[2] << "\n";
            //outfile << r << " " << E[0] << " " << E[1] << " " << E[2] << "\n";
        }

        outfile.close();
        std::cout << "    Electric Field data written to efield_z_580V.txt" << std::endl;
        //std::cout << "    Electric Field data written to efield_z_" << v << "um.txt" << std::endl;
        //std::cout << "    Electric Field data written to efield_r_" << v << "um.txt" << std::endl;
    }

    app.Run();
}