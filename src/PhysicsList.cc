#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"


#include "G4MTRunManager.hh"
#include "G4RunManager.hh"




#include "HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
//#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4IonPhysics.hh"
#include "G4IonINCLXXPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "GammaNuclearPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"

#include "ElectromagneticPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "StepMaxBuilder.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
:G4VModularPhysicsList()
{


  G4int verb = 1;
  SetVerboseLevel(verb);
  
  //add new units
  //
  new G4UnitDefinition( "millielectronVolt", "meV", "Energy", 1.e-3*eV);   
  new G4UnitDefinition( "mm2/g",  "mm2/g", "Surface/Mass", mm2/g);
  new G4UnitDefinition( "um2/mg", "um2/mg","Surface/Mass", um*um/mg);
    


 // Hadron Elastic scattering
  RegisterPhysics( new HadronElasticPhysicsHP(verb) );
  
  // Hadron Inelastic Physics
  //RegisterPhysics( new G4HadronPhysicsQGSP_BERT(verb));
  RegisterPhysics( new G4HadronPhysicsFTFP_BERT_HP(verb));
  ////RegisterPhysics( new G4HadronPhysicsQGSP_BIC_HP(verb));
  ////RegisterPhysics( new G4HadronInelasticQBBC(verb));        
  ////RegisterPhysics( new G4HadronPhysicsINCLXX(verb));
  
  // Ion Physics
  RegisterPhysics( new G4IonPhysics(verb));
  ////RegisterPhysics( new G4IonINCLXXPhysics(verb));
  
  // stopping Particles
  //RegisterPhysics( new G4StoppingPhysics(verb)); //should not been here
      
  // Gamma-Nuclear Physics
  RegisterPhysics( new GammaNuclearPhysics("gamma"));
  
  // EM physics
  RegisterPhysics(new ElectromagneticPhysics());
  //RegisterPhysics(new G4EmStandardPhysics());
  ////RegisterPhysics(new G4EmStandardPhysics_option3());
  

 //Optical Photons

  RegisterPhysics(new G4OpticalPhysics());

  // Decay
  RegisterPhysics(new G4DecayPhysics());

  // Radioactive decay
  RegisterPhysics(new G4RadioactiveDecayPhysics());

  // Step Max
  RegisterPhysics(new StepMaxBuilder());



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
  SetCutValue(0*mm, "proton");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
