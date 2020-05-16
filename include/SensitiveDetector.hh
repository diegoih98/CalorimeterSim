
#ifndef SensitiveDetector_HH
#define SensitiveDetector_HH 1

#include "G4VSensitiveDetector.hh"

//class EventAction;

class SensitiveDetector : public G4VSensitiveDetector
{
    public:
        SensitiveDetector(G4String name);
        ~SensitiveDetector();
        
        G4double GetQE(G4double energy);    
    
        void Initialize(G4HCofThisEvent*);
        G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);
        void EndOfEvent(G4HCofThisEvent*);
    
    private:
        G4double qe_nhits;
        G4double p2c;
        G4double qe; 
};

#endif
