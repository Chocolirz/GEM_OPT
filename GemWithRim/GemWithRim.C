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
    // If without holes, run the following line, units should match those in COMSOL files. 
    fm.Initialise("mesh_single_gem_600_210.mphtxt", "dielectrics.dat", "field_single_gem_600_210.txt", "mm");
    // If with copper rims on top and bottom, run the following line:
    //fm.Initialise("mesh_single_gem_560_150_rim.mphtxt", "dielectrics.dat", "field_single_gem_560_150_rim.txt", "mm");
    fm.EnableMirrorPeriodicityX(); // Generating infinite array of cells
    fm.EnableMirrorPeriodicityY(); 
    fm.PrintRange(); 

    // Dimensions of the GEM [cm]
    // Note that the geometry has changed, in ParallelPlate we have
    // y as the axis perpenticular to the plates, but here we have z. 
    const double rlim = 0.0075; // cm, radius of the GEM hole
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

    ViewField fieldView(&fm);
    constexpr bool plotPotential = false;
    if (plotPotential) {
        fieldView.SetPlane(0, -1, 0, 0.0001, 0, 0.0001);
        fieldView.SetArea(-pitch, zmin, pitch, zmax);
        
        // 2D plot for potential contours
        TCanvas* cpotential = new TCanvas("cpotential", "Contour lines of potential", 800, 800);
        fieldView.SetCanvas(cpotential);
        fieldView.PlotContour("v");
        cpotential->Update();

        // 2D plot for field contours
        TCanvas* cfield = new TCanvas("cfield", "Contour lines of field", 800, 800);
        fieldView.SetCanvas(cfield);
        fieldView.PlotContour("emag");
        cpotential->Update();

        // Plot of E_z along the axis
        TCanvas* cfield_z = new TCanvas("cfield_z", "Field along z-axis", 800, 800);
        cfield_z->SetLogy();
        fieldView.SetCanvas(cfield_z);
        fieldView.PlotProfile(0.0053033, 0.0053033, zmin, 0.0053033, 0.0053033, zmax, "field"); // slightly off the edge
        cfield_z->Update();

        // Plot of E arrows
        //TCanvas* cfield_arrows = new TCanvas("cfield_arrows", "Field arrows", 800, 800);
        //fieldView.SetCanvas(cfield_arrows);
        //fieldView.Plot("field", "arr");
        //cfield_arrows->Update();
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
    sensor.SetArea(xmin, ymin, zmin, xmax, ymax, zmax);

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
    constexpr unsigned int nEvents = 3000;



    
    // Set initial position [cm], time [ns] and energy [eV] for electrons
    const double x0 = 0.00001, y0 = 0.00001, z0 = zmax - upperGap, t0 = 0.;
    const double e0 = 0.1;
    const double dx0 = 0., dy0 = 0., dz0 = 0.; // vector length == 1 or 0 (direction radomised). 
    // Prepare some constants.
    double rz_e;
    double rz_s;
    double total_e;
    double anode_count; // Now we have the anode 2mm below
    double surv_count; // Count the number of electrons survived to the bottom of GEM
    double wall_count; // Count the electrons hitting the GEM inner cylinder
    double ratio;
    double ratio_s;
    double r1; // Number of electrons ending within the GEM hole.
    double r2; // Number of electrons hitting the GEM glass wall.
    double r3; // Number of electrons hitting the GEM bottom (copper plate).
    double r4; // Number of electrons attachments below the GEM. 

    // Add TH1F plots here
    TH1F hDistIon("hDistIon", "Distance between Ionising Collisions;Distance [cm];Entries", 200, 0., 0.1);
    TH1F hGain("hGain", "Gain Distribution;Anode Count;Entries", 500, 0., 3000.);
    TH1F hEndZ("hEndZ", "Height of Electron End Points;z [cm];Entries", 200, zmin, zmax);
    TH1F hRatio("hRatio", "Gain / Total # of Electrons Created;Ratio;Entries", 200, 0., 0.8);
    TH2F hEndXY("hEndXY", "Entries;x [cm];y [cm]", 500, xmin, xmax, 500, ymin, ymax);
    TH1F hSurv("hSurv", "Ratio of Electrons Survived to Bottom;Ratio;Entries", 200, 0., 1.);
    TH2F hStartRZ("hStartRZ", "Entries;R [cm];Z [cm]", 500, 0., 0.05, 500, zmin, zmax - upperGap);
    TH2F hEndRZ("hEndRZ", "Entries;R [cm];Z [cm]", 500, 0., 0.05, 500, zmin, zmax - upperGap);
    TH1F h1("h1", "Electrons attached in the hole;Ratio;Entries", 200, 0., 0.08);
    TH1F h2("h2", "", 200, 0., 0.8);
    TH1F h3("h3", "", 200, 0., 0.8);
    TH1F h4("h4", "", 200, 0., 0.08);
    TH1F hSEnd("hSend", "Number of Electrons Sent;R [cm];Entries", 200, 0., 0.03); // Number of electrons ends in the region lowergap + copperThick and lowerGap - copperThick
    TH1F hSStart("hSStart", "Number of Electrons Started;R [cm];Entries", 200, 0., 0.03); // Number of electrons starts in the region lowergap + copperThick and lowerGap - copperThick

    aval.EnableDistanceHistogramming(1); // code type == 1, ionization.
    aval.SetDistanceHistogram(&hDistIon);
    
    // Plot the avalanche
    constexpr bool plotAvalanche = true;
    if (plotAvalanche) {
        for (unsigned int i = 0; i < nEvents; i++) {
            std::cout << i << "/" << nEvents << "\n"; // Manual set counter

            aval.AvalancheElectron(x0, y0, z0, t0, e0, dx0, dy0, dz0);

            // Set entries counter
            total_e = 0.;
            anode_count = 0.;
            surv_count = 0.;
            wall_count = 0.;
            r1 = 0.;
            r2 = 0.;
            r3 = 0.;
            r4 = 0.;

            for (auto & electron : aval.GetElectrons()) {
                auto& p0 = electron.path.front();
                auto& p1 = electron.path.back();

                rz_s = sqrt(pow(p0.x, 2.) + pow(p0.y, 2.));
                rz_e = sqrt(pow(p1.x, 2.) + pow(p1.y, 2.));

                total_e = total_e + 1;
                hEndZ.Fill(p1.z);
                /*
                hStartRZ.Fill(rz_s, p0.z, 1./rz_s);
                hEndRZ.Fill(rz_e, p1.z, 1./rz_e);
                // Check the electrons ending at the bottom surface of the GEM
                //if (p1.z <= 0.2 && p1.z > 0.1) {
                //    hEndXY.Fill(p1.x, p1.y);
                //}
                
                // for hSEnd and hSStart
                if (p1.z < lowerGap + copperThick && p1.z > lowerGap - copperThick) {
                    hSEnd.Fill(rz_e);
                }
                if (p0.z < lowerGap + copperThick && p0.z > lowerGap - copperThick) {
                    hSStart.Fill(rz_s);
                }
                */
                // Sort the electrons based on their end points
                if (p1.z < lowerGap) {
                    //surv_count += 1; 
                    // We need to check the ending z-coordinate and the status
                    //if (electron.status == -7) {
                    //    r4 += 1;
                    //} else {
                        if (p1.z < 0.1 && electron.status != -7) {
                            anode_count += 1;
                        } //else {
                            //r3 += 1;
                        //}
                    //}
                } /*else {
                    if (electron.status == -7) {
                        r1 += 1;
                    } else {
                        if (p1.z < lowerGap + copperThick or rz_e > rlim) {
                            r3 += 1;
                        } else {
                            r2 += 1;
                        }
                    }
                }
                */
            }
            /*
            ratio_s = surv_count / (surv_count + wall_count);
            hSurv.Fill(ratio_s);
            */
            // Output the gain
            if (anode_count > 0) {
                hGain.Fill(anode_count);
                //ratio = anode_count / total_e;
                //hRatio.Fill(ratio);
            }
            /*
            //std::cout << "anode_count = " << anode_count << "\n";
            if (r1 > 0) {
                h1.Fill(r1 / total_e);
            }
            if (r2 > 0) {
                h2.Fill(r2 / total_e);
            }
            if (r3 > 0) {
                h3.Fill(r3 / total_e);
            }
            if (r4 > 0) {
                h4.Fill(r4 / total_e);
            }
            */
        }

        TCanvas* hdi = new TCanvas("hdi", "Distribution of Anode Counts per Avalanche", 800, 800);
        hdi->SetLogy();
        hGain.Draw();
        hdi->Update();
        /*
        TCanvas* h_end_z = new TCanvas("h_end_z", "Height of Electron End Points", 800, 800);
        h_end_z->SetLogy();
        hEndZ.Draw();
        h_end_z->Update();

        
        TCanvas* h_ratio = new TCanvas("h_ratio", "Ratio of Electrons Hitting somewhere", 800, 800);
        h_ratio->SetLogy();
        hRatio.SetMinimum(9e-1);
        hRatio.SetLineColor(kRed);
        hRatio.Draw();
        h2.SetMinimum(9e-1);
        h2.SetLineColor(kBlue);
        h2.Draw("same");
        h3.SetMinimum(9e-1);
        h3.SetLineColor(kOrange);
        h3.Draw("same");
        TLegend *legend = new TLegend(0.6, 0.7, 0.88, 0.88);
        legend->AddEntry(&hRatio, "Anode", "l");
        legend->AddEntry(&h2, "Glass", "l");
        legend->AddEntry(&h3, "Cu", "l");
        legend->Draw();
        h_ratio->Update();

        TCanvas* h_ratio_2 = new TCanvas("h_ratio_2", "Ratio of attachments", 800, 800);
        h_ratio_2->SetLogy();
        h1.SetMinimum(9e-1);
        h1.SetLineColor(kRed);
        h1.Draw();
        h4.SetMinimum(9e-1);
        h4.SetLineColor(kGreen + 2);
        h4.Draw("same");
        TLegend *legend1 = new TLegend(0.6, 0.7, 0.88, 0.88);
        legend1->AddEntry(&h1, "Hole_A", "l");
        legend1->AddEntry(&h4, "Drift_A", "l");
        legend1->Draw();
        h_ratio_2->Update();
        
        std::cout << "Anode:                      Mean = " << hRatio.GetMean() << ", Std Dev = " << hRatio.GetRMS() << std::endl;
        std::cout << "Copper at the bottom:       Mean = " << h3.GetMean()     << ", Std Dev = " << h3.GetRMS()     << std::endl;
        std::cout << "Glass wall:                 Mean = " << h2.GetMean()     << ", Std Dev = " << h2.GetRMS()     << std::endl;
        std::cout << "Attachment within GEM hole: Mean = " << h1.GetMean()     << ", Std Dev = " << h1.GetRMS()     << std::endl;
        std::cout << "Attachment in drift region: Mean = " << h4.GetMean()     << ", Std Dev = " << h4.GetRMS()     << std::endl;

        double sum = h1.GetMean() + h2.GetMean() + h3.GetMean() + h4.GetMean() + hRatio.GetMean();
        std::cout << "Sum of means: " << sum << std::endl;
        */

        //TCanvas* h_ratio_s = new TCanvas("h_ratio_s", "Ratio of Electrons Survived to Bottom", 800, 800);
        //h_ratio_s->SetLogy();
        //hSurv.Draw();
        //h_ratio_s->Update();

        //TCanvas* c_end_xy = new TCanvas("c_end_xy", "Electrons Ending at the Bottom Surface of GEM", 800, 800);
        //hEndXY.Draw();
        //c_end_xy->Update();

        /*
        TCanvas* c_end_rz = new TCanvas("c_end_rz", "R-Z relation of endding points", 800, 800);
        c_end_rz->SetLogz();
        hEndRZ.SetMinimum(9e-1);
        hEndRZ.Draw("COLZ");
        gPad->SetRightMargin(0.2);
        c_end_rz->Update();

        TCanvas* c_start_rz = new TCanvas("c_start_rz", "R-Z relation of starting points", 800, 800);
        c_start_rz->SetLogz();
        hStartRZ.SetMinimum(9e-1);
        hStartRZ.Draw("COLZ");
        gPad->SetRightMargin(0.2);
        c_start_rz->Update();
        */

        /*
        // Plot hSEnd and hSStart
        TCanvas* cSEnd = new TCanvas("cSEnd", "Number of Electrons Sent", 800, 800);
        cSEnd->SetLogy();
        hSEnd.SetMinimum(9e-1);
        hSEnd.Draw();
        cSEnd->Update();

        TCanvas* cSStart = new TCanvas("cSStart", "Number of Electrons Started", 800, 800);
        cSStart->SetLogy();
        hSStart.SetMinimum(9e-1);
        hSStart.Draw();
        cSStart->Update();
        */
        // Plot electron paths
        constexpr bool plotDriftLines = false;
        if (plotDriftLines) {
            TCanvas* cd = new TCanvas();
            constexpr bool plotMesh = true; // with mesh, viewed sideways. 
            if (plotMesh) {
                ViewFEMesh* meshView = new ViewFEMesh(&fm);
                meshView->SetArea(-pitch, -pitch, zmin, pitch, pitch, zmax);
                meshView->SetCanvas(cd);
                // x-z projection.
                meshView->SetPlane(0, -1, 0, 0, 0, 0);
                meshView->SetFillMesh(true);
                // Set the color of the kapton and the metal.
                meshView->SetColor(1, kYellow + 3);
                meshView->SetColor(2, kGray);
                meshView->EnableAxes();
                meshView->SetViewDrift(&driftView);
                meshView->Plot();
            } else { // without mesh, a 3D trajectory
                driftView.SetArea(-pitch, -pitch, zmin, pitch, pitch, zmax);
                driftView.SetCanvas(cd);
                constexpr bool twod = false;
                driftView.Plot(twod);
            }
        }
    }

    app.Run(kTRUE);
}