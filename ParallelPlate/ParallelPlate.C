#include <TApplication.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>

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
#include "Garfield/Sensor.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewMedium.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/SolidTube.hh"
#include "Garfield/ViewCell.hh"
#include "Garfield/SolidHole.hh"

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
    gStyle->SetHistLineWidth(2);
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

// Set up x-range.
const double xmin = -0.5; // define x range
const double xmax = 0.5;
// Define the geometry
const double gap = 0.057; // in cm
const double vdrift = -560.; // in volts
const double rlim = 0.0170 / 2.; // in cm
const double rbot = 0.0185 / 2.; // in cm
int anodePhoton;
int totalPhoton;
bool scint;

void userHandle(double x, double y, double z, double t, int type, int level, Medium* m, double e0, double e1, double dx0, double dy0, double dz0, double dx1, double dy1, double dz1) {
    // Check if the collision is an excitation or ionisation
    // Generate a photon using ROOT to a random direction
    // Emit photon in a random direction
    const double theta = M_PI * Garfield::RndmUniform();
    const double phi = 2. * M_PI * Garfield::RndmUniform();
    const double phiDiff = phi - atan(y / x);

    scint = (std::rand() % 100) < 33;

    // Check if the photon can reach the anode
    double r = sqrt(x*x + y*y);
    double rDiff = rlim - r;
    double zDiff = z;
    double lr = -r*cos(phiDiff) + sqrt((rlim * rlim) - (r*r*sin(phiDiff)*sin(phiDiff)));
    if (scint && level > 0 && level < 8) {
      totalPhoton += 1;
    }
    //std::cout << level;
    //std::cout << e1-e0 << std::endl;
}


int main(int argc, char* argv[]) {
 
  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();
 
  // Setup the gas.
  MediumMagboltz gas("cf4", 100.); // Material and percentage
  //MediumMagboltz gas("cf4", 70., "ar", 30.);
  gas.SetTemperature(293.15); // Temperature in Kelvin
  gas.SetPressure(50.); // Pressure in Torr
  //gas.LoadIonMobility("IonMobility_CF4+_CF4.txt"); // Ion transport properties
  //gas.SetMaxElectronEnergy(200.); // Electron energy in eV
  gas.EnableCrossSectionOutput();
  gas.Initialise(true);


  
  //SolidTube tube(0., 0.5 * gap, 0., rlim, 0.5 * gap, 0., 1., 0.);
  //GeometrySimple geo;
  //geo.AddSolid(&tube, &gas);

  //GeometrySimple geo;
  //for (unsigned int i = 0; i < 100; i++) {
  //  SolidTube* tube = new SolidTube(0., (0.995 * gap) - (0.01 * i * gap), 0., rlim + (0.01*i*(rbot - rlim)), 0.01 * gap, 0., 1., 0.);
  //  geo.AddSolid(tube, &gas);
  //}
  
  // Setup the cell.
  ComponentAnalyticField cmp;
  cmp.AddPlaneY(0., 0.); // Position and voltage
  cmp.AddPlaneY(gap, vdrift);
  cmp.SetMedium(&gas);
  //cmp.SetGeometry(&geo);
 


  // Plot the potential.
  constexpr bool plotPotential = false;
  if (plotPotential) {
    ViewField fieldView(&cmp);
    fieldView.SetArea(xmin, 0., xmax, gap);

    TCanvas* c1 = new TCanvas("c1", "", 600, 600);
    fieldView.SetCanvas(c1);
    fieldView.PlotContour();

    ViewCell cellView(&cmp);
    cellView.SetCanvas(fieldView.GetCanvas());
    cellView.SetArea(xmin, 0., xmax, gap);
    cellView.Plot2d();
  }

  // Set up sensor
  Sensor sensor(&cmp);
  sensor.SetArea(xmin, 0., xmin, xmax, gap, xmax);
 
  AvalancheMicroscopic aval(&sensor); // Microscopic tracking of electrons
  aval.SetCollisionSteps(100);
 
  AvalancheMC drift;
  drift.SetSensor(&sensor); // Track ions using MC technique
  drift.SetDistanceSteps(2.e-4);

 

  // Visualising drift lines
  ViewDrift driftView;
  driftView.SetPlane(0, 0, 1, 0, 0, 0);
  driftView.SetArea(-0.05, 0., 0.05, gap);
  //driftView.SetArea(xmin, xmin, 0., xmax, xmax, gap); // for 3D
  aval.EnablePlotting(&driftView);
  drift.EnablePlotting(&driftView);




  //NUMBER OF AVALANCHES SIMULATED
  constexpr unsigned int nEvents = 1;

 


  // Count the total number of ions produced the back-flowing ions.
  unsigned int nTotal = 0;
  unsigned int nBF = 0;

  // Set initial position [cm] and time [ns] for electron.
  const double x0 = 0., y0 = gap, z0 = 0., t0 = 0.;
  const double e0 = 0.1; // initial energy [eV]
  const double dx0 = 0., dy0 = -1., dz0 = 0.; // randomised if null, velocity direction
 
  double rz_e;
  double rz_s;
  int total_e;
  double anode_count;
  double ratio;

  // Plot the radius at the anode for electrons
  TH1F rend("rend","End;R [cm];Entries", 500, 0., 0.07); // rlim + 0.0015 if cylinder is added
  // and the r-z dots for each electron
  TH2F xyend("xyend","End ; R [cm]; Z [cm]", 300, 0., 0.07, 200, 0., gap);
  // Mean distance between collisions causing ionisation
  TH1F hDistIon("hDistIon", "Distance between Ionising Collisions;Distance [cm];Entries", 200, 0., 0.04);
  // Energy Histogram
  TH1F hEnergy("hEnergy", "Electron Energy Distribution", 200, 0., 130.);
  // Gain Histogram
  TH1F hGain("hGain", "Gain Distribution (polya)", 400, 0., 2000.);
  // Plot number of electrons reaching the anode / number of electrons created
  TH1F hRatio("hRatio", "Gain / Total No. e. Created;Ratio;Entries", 200, 0., 1.);
  // Plot y coordinates of electron endpoints.
  TH1F hEndY("hEndY", "Height of Electron Endpoints;Height [cm];Entries", 200, 0., gap + 0.003);

  TH1F hAtt("hAtt", "Attachments;Number;Entries", 200, 0., 20000.);
  //TH1F hPhoton("hPhoton", "Number of Photons per Avalanche;Photon Count;Entries", 200, 0., 1300.);
  //TH1F hPhTot("hPhTot", "Total Number of Photons per Avalanche;Total Photons;Entries", 200, 0., 5000.);
  //TH2F h2D("h2D", "Photons vs Electrons;Electrons Produced;Photons Produced", 200, 0., 150000, 200, 0., 60000);
  /*
  aval.EnableElectronEnergyHistogramming(&hEnergy);
  aval.EnableDistanceHistogramming(1); //ionization
  aval.SetDistanceHistogram(&hDistIon);
  */
  // Plot the avalanche
  constexpr bool plotAvalanche = true;
  if (plotAvalanche) {  
    // Run calculation of avalanche
    for (unsigned int i = 0; i < nEvents; i++){
     
      std::cout << i << "/" << nEvents << "\n";
      total_e = 0;
      totalPhoton = 0;

      //aval.SetUserHandleCollision(userHandle);
      aval.AvalancheElectron(x0, y0, z0, t0, e0, dx0, dy0, dz0);

      for (auto & electron : aval.GetElectrons()) {
        total_e += 1;
      }
      /*
      std::cout << "Number of Electrons: " << total_e << std::endl;

      if (totalPhoton != 0) {
        hPhoton.Fill(anodePhoton);
        hPhTot.Fill(totalPhoton);
        h2D.Fill(total_e, totalPhoton);
        std::cout << "Photons at Anode: " << anodePhoton << std::endl;
        std::cout << "Total Photons: " << totalPhoton << std::endl;
      }
      */
      
      int ne = 0, ni = 0;
      aval.GetAvalancheSize(ne, ni);
      //hGain.Fill(ne); 

      double e_tracks = aval.GetNumberOfElectronEndpoints();
      double attachment = e_tracks - ne; // account for electrons attached by df4
      hAtt.Fill(attachment);
      std::cout << ne << "  " << attachment << std::endl;
      /*
      // set electron entries counter
      total_e = 0.;
      anode_count = 0.;
      
      for (auto & electron : aval.GetElectrons()){
        auto& p0 = electron.path.front();
        auto& p1 = electron.path.back();

        rz_s = pow(pow(p0.x,2)+pow(p0.z,2), 0.5);
        rz_e = pow(pow(p1.x,2)+pow(p1.z,2), 0.5);
        hEndY.Fill(p1.y);

        total_e = total_e + 1;
        xyend.Fill(rz_e,p1.y,1/rz_e);

        if (p1.y < 0.005 && electron.status != -7 && rz_e < rbot) {
          rend.Fill(rz_e,1/rz_e);
          anode_count += 1;
        }
      }
      if (anode_count != 0) {
        hGain.Fill(anode_count);
        ratio = anode_count / total_e;
        hRatio.Fill(ratio);
        std::cout << "anode count "<< anode_count << ", total_e" << total_e << ", ratio" << ratio << ")\n";
      }
      */
      // Run if ion path is wanted.
      //for (const auto& electron : aval.GetElectrons()) {
      //  const auto& p0 = electron.path[0];
      //  drift.DriftIon(p0.x, p0.y, p0.z, p0.t);
      //  ++nTotal;
      //  const auto& endpoint = drift.GetIons().front().path.back();
      //  if (endpoint.z <= 0.) ++nBF;
      //}
    }

    // If 3D/2D trajectory is wanted, say true.
    constexpr bool plotDriftLines = true;
    if (plotDriftLines) {
      TCanvas* cD = new TCanvas("cD", "", 600, 600);
      driftView.SetCanvas(cD);
      constexpr bool twod = true;
      driftView.Plot(twod);
      //driftView.Plot(); // for 3D
      cD->Update();
    }

    SetMyRootStyle();

    //TCanvas* y_end = new TCanvas("y_end", "Height at Endpoint", 800, 800);
    //y_end->SetLogy();
    //hEndY.Draw();
    //y_end->Update();

    //TCanvas* r_end = new TCanvas("r_end", "Radius at Anode", 5000, 0, 500, 500);
    //r_end->SetLogy();
    //rend.Draw();
    //r_end->Update();

    //TCanvas* cRatio = new TCanvas("cRatio", "rlim Dependence", 800, 800);
    //cRatio->SetLogy();
    //hRatio.Draw();
    //cRatio->Update();

    TCanvas* hdi = new TCanvas("hdi", "hdi",800,800);
    hdi->SetLogy();
    //hGain.SetLineColor(kRed);
    //hGain.Draw();
    hAtt.SetLineColor(kBlue);
    hAtt.Draw();
    hdi->Update();
    //std::cout << rlim*2.*10000. << "    " << rbot*2.*10000. << "    " << hGain.GetMean() << "\n";

    //TCanvas* hxyend = new TCanvas("hxyend", "hxyend", 800, 800);
    //hxyend->SetLogz();
    //xyend.Draw("colz");
    //xyend.SetLeftMargin(0.2);
    //hxyend->Update();

    //TCanvas* cDist = new TCanvas("cDist", "Ionising Collision Distances", 800, 600);
    //cDist->SetLogy();
    //hDistIon.Draw();
    //cDist->Update();

    //TCanvas* cEnergy = new TCanvas("cEnergy", "Electron Energy Distribution", 600, 600);
    //cEnergy->SetLogy();
    //hEnergy.Draw();
    //cEnergy->Update();

    //TCanvas* c2D = new TCanvas("c2D", "Photons vs Electrons", 800, 800);
    //c2D->SetLogz();
    //h2D.Draw();
    //c2D->Update();
  }
 
  app.Run(kTRUE);
}