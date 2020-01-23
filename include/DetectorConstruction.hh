//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// $Id: DetectorConstruction.hh 78560 2014-01-07 10:06:52Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

const G4int kMaxAbsor = 10;                        // 0 + 8

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
  virtual void ConstructSDandField();
     
public:

  G4int       GetNbOfAbsor()             {return fNbOfAbsor;}
  G4Material* GetAbsorMaterial (G4int i) {return fAbsorMaterial[i];};
  G4double    GetAbsorThickness(G4int i) {return fAbsorThickness[i];};      
  G4double    GetXfront        (G4int i) {return fXfront[i];};
  G4Material* GetWorldMaterial()         {return fDefaultMaterial;}
            
  G4double GetAbsorSizeZ()               {return fAbsorSizeZ;}; 
  G4double GetAbsorSizeXY()              {return fAbsorSizeXY;};
  
  G4double GetWorldSizeZ()               {return fWorldSizeZ;}; 
  G4double GetWorldSizeXY()              {return fWorldSizeXY;}; 

  G4LogicalVolume* GetLogPhotoCath() {return fPhotocath_log;}

  
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
  G4LogicalVolume* fLogicA1;
  G4LogicalVolume* fLogicA2;
  G4LogicalVolume* fLogicA3;
  G4LogicalVolume* fLogicA4;

   G4VPhysicalVolume* physAl; 
   G4VPhysicalVolume* physPTM;



    G4Box* fScint_box;
    G4Box* fHousing_box;
  //  G4Tubs* fPmt;
  //  G4Tubs* fPhotocath;
    //G4Sphere* fSphere;
    G4LogicalVolume* fScint_log;
    G4LogicalVolume* fHousing_log;
    G4LogicalVolume* fPmt_log;
    G4LogicalVolume* fPhotocath_log;
    G4LogicalVolume* fSphere_log;
//0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0

  DetectorMessenger* fDetectorMessenger;
  G4Cache<G4GlobalMagFieldMessenger*> fFieldMessenger;


private:

  void DefineMaterials();
  void ComputeParameters();
  void SurfaceProperties();
  G4VPhysicalVolume* ConstructVolumes();
  DetectorConstruction & operator=(const DetectorConstruction &right);
  DetectorConstruction(const DetectorConstruction&);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

