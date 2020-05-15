#include "EventAction.hh"

#include "Run.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* det)
:G4UserEventAction(), fDetector(det)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  //energy deposited per event
  for (G4int k=0; k<kMaxAbsor; k++) { fEdepAbsor[k] = 0.0; }
  fPhotonCount = 0.0;
  fPhotonCountEnd = 0.0;
  TotalEdepEvent = 0.0;
  qenergy = 0.0;
  nhits = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  
  //Energy deposited per event
  //
  G4double TotalEdep(0.);
  for (G4int k=1; k<=fDetector->GetNbOfAbsor(); k++) {
    if (fEdepAbsor[k] > 0.) {
      run->AddEdep(k,fEdepAbsor[k]);
      TotalEdep += fEdepAbsor[k];
    }
  }
  
  if (TotalEdep > 0.) {
    run->AddTotEdep(TotalEdep);
  }


   // G4double qe = 0.3;                      // quantum efficiency of photocathode
    G4double p2c = 0.38;                 // photon-to-charge ratio (gain)
    G4double q_dep = fPhotonCount * p2c; //* qe;      // deposited charge


  G4AnalysisManager::Instance()->FillH1(12, q_dep);
  //G4AnalysisManager::Instance()->FillH1(13, fPhotonCountEnd);
  G4AnalysisManager::Instance()->FillH1(14, TotalEdepEvent);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

