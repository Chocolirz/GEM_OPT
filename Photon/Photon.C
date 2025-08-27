#include <TApplication.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TColor.h>
#include <TGaxis.h>
#include <TF2.h>

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

void SetMyRootStyle() {
    // General font setup (42 = Helvetica, closest to Arial)
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");

    // Font sizes (axis labels and tick labels)
    gStyle->SetLabelSize(0.045, "XYZ");   // Tick labels (≈18 pt)
    gStyle->SetTitleSize(0.05, "XYZ");    // Axis titles (≈21 pt)
    gStyle->SetTitleOffset(1.5, "X");     // Label padding
    gStyle->SetTitleOffset(1.5, "Y");

    // Line and axis settings
    gStyle->SetHistLineWidth(2);          // Line width ≈ 2.5
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameLineWidth(1);

    // Grid and background
    gStyle->SetPadGridX(0);
    gStyle->SetPadGridY(0);
    gStyle->SetPadColor(0);
    gStyle->SetCanvasColor(0);

    // Legend styling
    gStyle->SetLegendBorderSize(0);       // No frame
    gStyle->SetLegendFillColor(0);        // Transparent
    gStyle->SetLegendFillStyle(0);
    gStyle->SetLegendTextSize(0.035);     // ≈18 pt

    // Image color map (approx 'turbo')
    gStyle->SetPalette(kBird);            // Alternatives: kViridis, kInvertedDarkBodyRadiator

    // Disable stat box
    gStyle->SetOptStat(0);

    // Use bold lines and disable default markers
    gStyle->SetLineWidth(2);
    gStyle->SetMarkerStyle(1);            // Optional
}

std::vector<int> tolColors = {
    TColor::GetColor("#4477aa"), // blue
    TColor::GetColor("#ee6677"), // red/pink
    TColor::GetColor("#228833"), // green
    TColor::GetColor("#aa3377"), // purple
    TColor::GetColor("#66ccee"), // cyan
    TColor::GetColor("#ccbb44"), // yellow
    TColor::GetColor("#bbbbbb"), // grey
    kBlack                      // black
};


// Dimensions of the GEM [cm]
// Note that the geometry has changed, in ParallelPlate we have
// y as the axis perpenticular to the plates, but here we have z. 
const double rlim = 0.0190 / 2.; // cm, radius of the GEM hole
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
int anodePhoton;
int totalPhoton;
int level;
double PhotonX;
double PhotonY;

//TH2F hPhXY("hPhXY", "Entries;x [cm];y [cm]", 100, xmin, xmax, 100, ymin, ymax);

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
    double zDiff = z - lowerGap;
    double lr = -r*cos(phiDiff) + sqrt((rlim * rlim) - (r*r*sin(phiDiff)*sin(phiDiff)));
    if (scint && level > 0 && level < 8) {
        if (zDiff > 0 && theta > M_PI_2 && tan(theta) > -lr/zDiff) {
            anodePhoton += 1;
            //PhotonX = x - tan(theta)*cos(phi)*z;
            //PhotonY = y - tan(theta)*sin(phi)*z;
            //hPhXY.Fill(PhotonX, PhotonY, 1.);
        } else if (zDiff < 0 && theta > M_PI_2) {
            anodePhoton += 1;
            //PhotonX = x - tan(theta)*cos(phi)*z;
            //PhotonY = y - tan(theta)*sin(phi)*z;
            //hPhXY.Fill(PhotonX, PhotonY, 1.);
        }
        totalPhoton += 1;
    }
    //std::cout << level;
    //std::cout << e1-e0 << std::endl;
}


int main(int argc, char* argv[]) {
    TApplication app("app", &argc, argv);

    // Load field map.
    ComponentComsol fm;
    // If without holes, run the following line, units should match those in COMSOL files. 
    fm.Initialise("mesh_single_gem_560_190.mphtxt", "dielectrics.dat", "field_single_gem_560_190.txt", "mm");
    // If with copper rims on top and bottom, run the following line:
    //fm.Initialise("mesh_single_gem_560_150_rim.mphtxt", "dielectrics.dat", "field_single_gem_560_150_rim.txt", "mm");
    fm.EnableMirrorPeriodicityX(); // Generating infinite array of cells
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
    gas.PrintGas();
    // Set gas and see its properties
    fm.SetGas(&gas);
    fm.PrintMaterials();

    // Set up sensor
    Sensor sensor(&fm);
    sensor.SetArea(xmin, ymin, zmin, xmax, ymax, zmax);

    AvalancheMicroscopic aval(&sensor); // Microscopic tracking of electrons
    aval.SetCollisionSteps(100);

    AvalancheMC drift(&sensor);
    drift.SetDistanceSteps(2.e-4);




    // NUMBER OF AVALANCHES SIMULATED
    constexpr unsigned int nEvents = 3000;



    
    // Set initial position [cm], time [ns] and energy [eV] for electrons
    const double x0 = 0.00001, y0 = 0.00001, z0 = zmax - upperGap, t0 = 0.;
    const double e0 = 0.1;
    const double dx0 = 0., dy0 = 0., dz0 = 0.; // vector length == 1 or 0 (direction radomised). 
    // Prepare some constants.
    double rz_e;
    double rz_s;
    int total_e;
    double anode_count; // Now we have the anode 2mm below
    int ne, ni;
    TH1F hPhoton("hPhoton", "Number of Photons per Avalanche;Photon Count;Entries", 200, 0., 1300.);
    TH1F hPhTot("hPhTot", "Total Number of Photons per Avalanche;Total Photons;Entries", 200, 0., 5000.);
    //TH2F h2D("h2D", "Photons vs Electrons;Electrons Produced;Photons Produced", 100, 0., 1000, 100, 0., 300);


    /*
    int ngas; 
    int type;
    std::string descr;
    double e;
    */

    // Plot the avalanche
    constexpr bool plotAvalanche = true;
    if (plotAvalanche) {
        for (unsigned int i = 0; i < nEvents; i++) {
            std::cout << i << "/" << nEvents << "\n"; // Manual set counter
            anodePhoton = 0;
            totalPhoton = 0;
            total_e = 0;

            //auto level_str = gas.GetLevel(1, ngas, type, descr, e);
            //std::cout << ngas << "  " << type << "  " << descr << " " << e << std::endl;

            aval.SetUserHandleCollision(userHandle);
            aval.AvalancheElectron(x0, y0, z0, t0, e0, dx0, dy0, dz0);
            
            for (auto & electron : aval.GetElectrons()) {
                total_e += 1;
            }

            std::cout << "Number of Electrons: " << total_e << std::endl;

            if (total_e > 1) {
                hPhoton.Fill(anodePhoton);
                hPhTot.Fill(totalPhoton);
                //h2D.Fill(totalPhoton, anodePhoton);
                std::cout << "Photons at Anode: " << anodePhoton << std::endl;
                std::cout << "Total Photons: " << totalPhoton << std::endl;
            }
        }

        //SetMyRootStyle();

        TCanvas* cPhoton = new TCanvas("cPhoton", "Number of Photons per Avalanche", 800, 800);
        cPhoton->SetLogy();
        hPhoton.Draw();
        cPhoton->Update();
        
        TCanvas* cPhTot = new TCanvas("cPhTot", "Total Number of Photons", 800, 800);
        cPhTot->SetLogy();
        hPhTot.Draw();
        cPhTot->Update();


        /*
        TCanvas* c2D = new TCanvas("c2D", "Photons vs Electrons", 800, 800);
        c2D->SetLogz();
        h2D.Draw("COLZ");
        c2D->Update();
        
        
        TCanvas* cPhXY = new TCanvas("cPhXY", "Photon map at anode", 800, 800);
        cPhXY->SetLogz();
        hPhXY.SetMinimum(9e-1);
        hPhXY.Draw();
        cPhXY->Update();

        
        // Define the 2D Gaussian centered at (0,0)
        TF2* gaus2D = new TF2("gaus2D",
            "[0]*exp(-0.5*((x/[1])**2 + (y/[2])**2))",
            hPhXY.GetXaxis()->GetXmin(), hPhXY.GetXaxis()->GetXmax(),
            hPhXY.GetYaxis()->GetXmin(), hPhXY.GetYaxis()->GetXmax()
        );

        // Set initial guesses: [0]=amplitude, [1]=sigma_x, [2]=sigma_y
        gaus2D->SetParameters(hPhXY.GetMaximum(), 0.1, 0.1);

        // Fit the hPhXYogram
        hPhXY.Fit(gaus2D, "R"); // "R" restricts fit to range of function

        TH2F* residual = (TH2F*)hPhXY.Clone("residual");
        residual->SetTitle("Residuals: Data - Fit");

        // Loop over all bins and compute: data - fit(x, y)
        for (int ix = 1; ix <= hPhXY.GetNbinsX(); ++ix) {
            for (int iy = 1; iy <= hPhXY.GetNbinsY(); ++iy) {
                double x = hPhXY.GetXaxis()->GetBinCenter(ix);
                double y = hPhXY.GetYaxis()->GetBinCenter(iy);
                double data = hPhXY.GetBinContent(ix, iy);
                double model = gaus2D->Eval(x, y);
                double diff = data - model;
                residual->SetBinContent(ix, iy, diff);
            }
        }

        // Draw the residual map
        TCanvas* c2 = new TCanvas("c2", "Residual Map", 800, 800);
        c2->SetLogz();
        residual->SetMinimum(9e-1);
        residual->Draw("COLZ");
        c2->Update();
        */
    }

    app.Run(kTRUE);
}