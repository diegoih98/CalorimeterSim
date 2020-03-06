//#include "LGSimRunAction.hh"
#include "SensitiveDetector.hh"
#include "EventAction.hh"
#include "G4RunManager.hh"
#include "HistoManager.hh"
#include "G4SystemOfUnits.hh"

SensitiveDetector::SensitiveDetector(G4String name) 
//:G4VSensitiveDetector(name),  nhits(0) {}
:G4VSensitiveDetector(name), qe_nhits(0.), p2c(3.8e-13) {}

SensitiveDetector::~SensitiveDetector() {}


G4double SensitiveDetector::GetQE(G4double energy)
{
    G4double qetable[] = {0., 0.04, 0.2, 0.22, 0.3, 0.3, 0.15, 0.03, 0.};
    G4double eVtable[] = {1.77*eV, 2.07*eV, 2.48*eV, 2.76*eV, 3.10*eV, 3.54*eV, 4.13*eV, 4.43*eV, 4.77*eV};
    
    G4int tableSize = 9;
    
    G4int i = -1;
    while(energy >= eVtable[i+1]){
        i++;
    }
    
    if(i < 0 || i >= tableSize - 1)
        return 0.;
    else
        return (qetable[i+1] - qetable[i]) * (energy - eVtable[i]) 
                / (eVtable[i+1] - eVtable[i]) + qetable[i];
    
}






void SensitiveDetector::Initialize(G4HCofThisEvent*)
{
  //  nhits = 0;
  qe_nhits = 0;
    p2c = 8.4e-13;                 // photon-to-charge ratio (gain)

}

G4bool SensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{   
    auto aTrack = aStep->GetTrack();
    auto aParticle = aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName();
    
    if(aParticle == "opticalphoton" && aStep->IsFirstStepInVolume()){
        qe_nhits += GetQE(aTrack->GetKineticEnergy());
    }
    
    return true;
}

void SensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{


  G4double q_dep = qe_nhits * p2c;      // deposited charge
    
    if(q_dep > 0){
       // G4cout << "Hit: " << qe_nhits << " photons" << G4endl;
        //G4cout << "Deposited charge: " << q_dep * 1e12 << " pC" << G4endl;

          G4AnalysisManager::Instance()->FillH1(15, q_dep * 1e12);

    }
    else
        G4cout << "No hit detected." << G4endl;
        
}
