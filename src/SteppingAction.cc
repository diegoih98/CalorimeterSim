#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "Run.hh"
#include "HistoManager.hh"

#include "G4SteppingManager.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* event)
:G4UserSteppingAction(),fDetector(det), fEventAction(event)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // Count process
  
  const G4StepPoint* endPoint = step->GetPostStepPoint();
  const G4VProcess* process   = endPoint->GetProcessDefinedStep();
  Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->CountProcesses(process);
  
 G4double edep = step->GetTotalEnergyDeposit();
 if (edep <= 0.) return;
 
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();     
 
 //Longitudinal profile of deposited energy

 G4StepPoint* prePoint  = step->GetPreStepPoint();
 G4StepPoint* postPoint = step->GetPostStepPoint(); 
 G4ThreeVector P1 = prePoint ->GetPosition();
 G4ThreeVector P2 = postPoint->GetPosition();
 G4ThreeVector point = P1 + G4UniformRand()*(P2 - P1);  //Randomize point of energy deposition
 G4double z = point.z();
 G4double zshifted = z + 0.5*fDetector->GetAbsorSizeZ();  
 analysisManager->FillH1(10, zshifted, edep);

 fEventAction->AddTEdep(edep);

 //Total energy deposit in calorimeter
  
  const G4VPhysicalVolume* pv = prePoint->GetPhysicalVolume();
    G4int copyNo = pv->GetCopyNo();
   fEventAction->AddEdep(copyNo, edep/MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

