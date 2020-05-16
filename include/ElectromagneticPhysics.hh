#ifndef ElectromagneticPhysics_h
#define ElectromagneticPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class G4Cerenkov;
class G4Scintillation;
class G4OpAbsorption;
class G4OpRayleigh;
class G4OpMieHG;
class G4OpBoundaryProcess;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ElectromagneticPhysics : public G4VPhysicsConstructor
{
  public: 
    ElectromagneticPhysics(const G4String& name = "standard");
   ~ElectromagneticPhysics();

  public: 

    virtual void ConstructParticle() {};

    virtual void ConstructProcess();

  private:
      G4Cerenkov*          theCerenkovProcess;
      G4Scintillation*     theScintillationProcess;
      G4OpAbsorption*      theAbsorptionProcess;
      G4OpRayleigh*        theRayleighScatteringProcess;
      G4OpMieHG*           theMieHGScatteringProcess;
      G4OpBoundaryProcess* theBoundaryProcess;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
