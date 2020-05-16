#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Cache.hh"
#include "G4OpticalSurface.hh"
#include "G4Tubs.hh"

//class G4Tubs;
//class G4Sphere;
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;
//class EventAction;

const G4int kMaxAbsor = 10;                        // 0 + 9

class G4GlobalMagFieldMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  DetectorConstruction();
 ~DetectorConstruction();

public:

  G4Material* 
  MaterialWithSingleIsotope(G4String, G4String, G4double, G4int, G4int);
  
  void SetNbOfAbsor     (G4int);      
  void SetAbsorMaterial (const G4String&);     
                
     
  virtual G4VPhysicalVolume* Construct();
     
public:

  G4int       GetNbOfAbsor()             {return fNbOfAbsor;}
  G4Material* GetAbsorMaterial (G4int i) {return fAbsorMaterial[i];};
  G4double    GetAbsorThickness(G4int i) {return fAbsorThickness[i];};      
  G4Material* GetWorldMaterial()         {return fDefaultMaterial;}
            
  G4double GetAbsorSizeZ()               {return fAbsorSizeZ;}; 
  G4double GetAbsorSizeXY()              {return fAbsorSizeXY;};
  
  G4double GetWorldSizeZ()               {return fWorldSizeZ;}; 
  G4double GetWorldSizeXY()              {return fWorldSizeXY;}; 

 // G4LogicalVolume* GetLogPhotoCath() {return fPhotocath_log;}

  
  void PrintParameters();
   
private:

  G4int              fNbOfAbsor;
  G4Material*        fAbsorMaterial [kMaxAbsor];
  G4double           fAbsorThickness[kMaxAbsor];
  G4double           fXfront[kMaxAbsor];  

  G4double           fAbsorSizeZ;
  G4double           fAbsorSizeXY;
  
  G4double           fWorldSizeZ;
  G4double           fWorldSizeXY;  
  G4Material*        fDefaultMaterial;  
  
  G4VPhysicalVolume* fPhysiWorld;


//0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0

  G4double fEcalLength;
  G4double fEcalWidth;
  G4double fVertexLength;
  G4double fPadLength;
  G4double fPadWidth;
  G4double fAbsLength;
  G4double fWorldZ;
  G4double fractionmass;

  G4Material* fCalMaterial;
  G4Material* fVertMaterial;
  G4Material* fAbsMaterial;
  G4Material* fWorldMaterial;
  G4Material* fYorkMaterial;

  G4Material* fSiliconDioxide; 
  G4Material* fLeadOxide; 
  G4Material* fPotassiumOxide;
  G4Material* fSodiumMonoxide;

  G4Material* fGlass;

  G4Material* fAntimonyTrioxideMaterial;

  G4LogicalVolume* fLogicWorld;
  G4LogicalVolume* fLogicCal;
  G4LogicalVolume* fLogicAl;

    G4LogicalVolume* cathodeLV;
    G4LogicalVolume* alPlateLV;
    G4LogicalVolume* lgTubeLV;
    G4LogicalVolume* alShieldingLV;
    G4LogicalVolume* lgLV;

  G4VPhysicalVolume* lgPV;
  G4VPhysicalVolume* alPlatePV;
  G4VPhysicalVolume* alCoverPV;
  G4VPhysicalVolume* lgTubePV;
  G4VPhysicalVolume* modulePV[9]; // 0 + 9
//0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0

  DetectorMessenger* fDetectorMessenger;
  G4Cache<G4GlobalMagFieldMessenger*> fFieldMessenger;


private:

  void DefineMaterials();
  void ComputeParameters();
  void SurfaceProperties();
  void ConstructSDandField();
  G4VPhysicalVolume* ConstructVolumes();
  DetectorConstruction & operator=(const DetectorConstruction &right);
  DetectorConstruction(const DetectorConstruction&);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
