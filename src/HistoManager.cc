//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
// $Id: HistoManager.cc 95740 2016-02-23 09:34:37Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "DetectorConstruction.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("Hadr07")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);

  // Define histograms start values
  const G4int kMaxHisto = 16;
  const G4String id[] = { "f0", "f1", "f2", "f3", "f4", "f5", "f6", "f7", "f8", "f9",
                         "f10","f11","f12","f13","f14","f15","f16","f17","f18","f19",
                         "f20","f21","f22"};

  const G4String title[] = 
                { "dummy",                                        //0
                  "PhotonsProducedSpectrum",         //1
                  "PhotonReachingPMTSpectrum",         //2
                  "PhotonsDetectedSpectrum",         //3
                  "Energy of Gamma produced in Pb glass counter 4 (MeV)",         //4
                  "Energy of Gamma produced in Pb glass counter 5 (MeV)",         //5
                  "Energy of Gamma produced in Pb glass counter 6 (MeV)",         //6
                  "Energy of Gamma produced in Pb glass counter 7 (MeV)",         //7
                  "Energy of Gamma produced in Pb glass counter 8 (MeV)",         //8
                  "Energy of Gamma produced in Pb glass counter 9 (MeV)",         //9
                  "Edep (MeV/mm) along calorimeter",                 //10
                  "Energy of photons (MeV) emerging the Pb glass", //11
                  "Number of photons produced per event", //12
                  //"Number of photons emerging per event", //13
                  "Photoelectrons", //13
                  "Energy Spectrum", //14
                  "Energy detected", //15
                 };

  // Default values (to be reset via /analysis/h1/set command)
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

 /* // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
   analysisManager->SetH1Activation(ih, false);
  }*/

G4int ih = analysisManager->CreateH1(id[0], title[0], nbins, vmin, vmax);
   analysisManager->SetH1Activation(ih, false);
for (G4int k=1; k<kMaxHisto; k++) {
    ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
   analysisManager->SetH1Activation(ih, true);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
