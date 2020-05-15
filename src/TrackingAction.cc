#include "TrackingAction.hh"

#include "Run.hh"
#include "HistoManager.hh"
#include "EventAction.hh"

#include "G4Track.hh"
#include "G4StepStatus.hh"
#include "G4RunManager.hh"
#include "G4Material.hh"
#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* det, EventAction* event)
:G4UserTrackingAction(),fDetector(det),fEventAction(event)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{  
  //count secondary particles
  if (track->GetTrackID() == 1) return;
  G4int iabs = track->GetTouchableHandle()->GetCopyNumber();
  G4String name   = track->GetDefinition()->GetParticleName();
  G4double energy = track->GetKineticEnergy();
  Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());    
  if(iabs > 0){
     run->ParticleCount(iabs,name,energy);
         if(name == "opticalphoton"){
          fEventAction->AddPhoton();
          G4AnalysisManager::Instance()->FillH1(1, energy);
          }
     }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
  // Get Run
  Run* run 
    = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  
  // Find location
  G4StepStatus status = track->GetStep()->GetPostStepPoint()->GetStepStatus();
  
  //Primary particle : absorbed or transmited?
  if (track->GetTrackID() == 1) {
    G4int flag = 0;
    if (status == fWorldBoundary) flag = 1;
    run->AddTrackStatus(flag);
  }

  // keep only emerging particles
  if (status != fWorldBoundary) return;

  // count particles
  const G4ParticleDefinition* particle = track->GetParticleDefinition();
  G4String name   = particle->GetParticleName();
  G4double energy = track->GetKineticEnergy();
  run->ParticleCount(0,name,energy);

G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

