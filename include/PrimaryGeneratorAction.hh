#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include "g4root.hh"
#include "G4GeneralParticleSource.hh"

class G4Event;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();    
   ~PrimaryGeneratorAction();

  public:
    void SetDefaultKinematic();  
    virtual void GeneratePrimaries(G4Event*);
         
    G4GeneralParticleSource* GetParticleGun() {return gun;}
    
  private:
    G4GeneralParticleSource* gun;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

