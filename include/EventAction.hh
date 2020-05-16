#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction(DetectorConstruction*);
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    
    void AddEdep(G4int k, G4double edep) { fEdepAbsor[k] += edep; }
    void AddPhoton() { fPhotonCount += 1; }
    void AddPhotonEnd() { fPhotonCountEnd += 1; }
    void AddHits() { nhits += 1; }
    G4int GetHits() {return nhits;}
    void AddCharge(G4double q) { qenergy += q; }

    void AddTEdep(G4double edep) { TotalEdepEvent += edep; }
    
  private:
    DetectorConstruction* fDetector;

    G4double  fEdepAbsor[kMaxAbsor];
    G4double  fPhotonCount;
    G4double  fPhotonCountEnd;

    G4int nhits;
    G4double TotalEdepEvent;

    G4double qenergy;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
