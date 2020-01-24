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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id: DetectorConstruction.cc 77288 2013-11-22 10:52:58Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4RunManager.hh"
#include <iomanip>

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4Trap.hh"

#include "G4ios.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4MaterialPropertiesTable.hh"

#include "G4LogicalSkinSurface.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),fDefaultMaterial(0),fPhysiWorld(0),
 fDetectorMessenger(0)
{
  // default parameter values of the absorbers
  fNbOfAbsor = 9;
  fAbsorThickness[0] = 340*mm;        //dummy, for initialization   
  fAbsorThickness[1] = 340*mm;
  fAbsorThickness[2] = 340*mm;  
  fAbsorThickness[3] = 340*mm;  
  fAbsorThickness[4] = 340*mm;  
  fAbsorThickness[5] = 340*mm;  
  fAbsorThickness[6] = 340*mm;  
  fAbsorThickness[7] = 340*mm;  
  fAbsorThickness[8] = 340*mm; 
  fAbsorThickness[9] = 340*mm;       
  fAbsorSizeXY       = 61*mm;


  ComputeParameters();

  // materials
  DefineMaterials();
  SetAbsorMaterial("CalMaterial");

  // create commands for interactive definition of the calorimeter
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Default materials

  G4NistManager* man = G4NistManager::Instance();
  //  man->SetVerbose(1);
  fWorldMaterial = man->FindOrBuildMaterial("G4_AIR");
  fAbsMaterial   = man->FindOrBuildMaterial("G4_Al");

  fSiliconDioxide = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  fLeadOxide = man->FindOrBuildMaterial("G4_LEAD_OXIDE");
  fPotassiumOxide = man->FindOrBuildMaterial("G4_POTASSIUM_OXIDE");
  fSodiumMonoxide = man->FindOrBuildMaterial("G4_SODIUM_MONOXIDE");
 
 //Oxygen
   G4double A = 16.0 * g/mole;
   G4double Z = 8;
   G4Element* fOxigen = new G4Element ("Oxygen", "O", Z, A);

  //Antimony  
   A = 121.76 * g/mole;
   Z = 51;
   G4Element* fAntimony = new G4Element ("Antimony", "Sb", Z, A);

  //Hydrogen 
   A = 11.01*g/mole;
   Z = 1.;
  G4Element* fH = new G4Element("H", "H", Z, A);

  //Carbon  
   A = 12.01*g/mole;
   Z = 6.;
  G4Element* fC = new G4Element("C", "C", Z, A);  

//Glass

  fGlass = new G4Material("Glass", 1.032*g/cm3,2);
  fGlass->AddElement(fC,91.533*perCent);
  fGlass->AddElement(fH,8.467*perCent);


  fAntimonyTrioxideMaterial= new G4Material("Sb2O3", 5.2*g/cm3, 2); 
  fAntimonyTrioxideMaterial->AddElement(fAntimony, 2);
  fAntimonyTrioxideMaterial->AddElement(fOxigen, 3);

  fCalMaterial = new G4Material("fCalMaterial",5.2*g/cm3,5);
  fCalMaterial->AddMaterial(fSiliconDioxide, fractionmass=27.3*perCent);
  fCalMaterial->AddMaterial(fLeadOxide, fractionmass=70.9*perCent);
  fCalMaterial->AddMaterial(fPotassiumOxide, fractionmass=0.9*perCent);
  fCalMaterial->AddMaterial(fSodiumMonoxide, fractionmass=0.6*perCent);
  fCalMaterial->AddMaterial(fAntimonyTrioxideMaterial, fractionmass=0.3*perCent);

  fDefaultMaterial = fWorldMaterial;

  G4double cal_Energy[]    = { 7.0*eV , 7.07*eV, 7.14*eV };
  const G4int calnum = sizeof(cal_Energy)/sizeof(G4double);

  G4double cal_RIND[]  = { 1.8 , 1.8, 1.8 };
  assert(sizeof(cal_RIND) == sizeof(cal_Energy));
  G4MaterialPropertiesTable* fCalMaterial_mt = new G4MaterialPropertiesTable();
  fCalMaterial_mt->AddProperty("RINDEX", cal_Energy, cal_RIND,  calnum);
  fCalMaterial->SetMaterialPropertiesTable(fCalMaterial_mt);

 
  G4double glass_RIND[]={1.49,1.49,1.49};
  assert(sizeof(glass_RIND) == sizeof(cal_Energy));
  G4double glass_AbsLength[]={420.*cm,420.*cm,420.*cm};
  assert(sizeof(glass_AbsLength) == sizeof(cal_Energy));
  G4MaterialPropertiesTable *glass_mt = new G4MaterialPropertiesTable();
  glass_mt->AddProperty("ABSLENGTH",cal_Energy,glass_AbsLength,calnum);
  glass_mt->AddProperty("RINDEX",cal_Energy,glass_RIND,calnum);
  fGlass->SetMaterialPropertiesTable(glass_mt);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::MaterialWithSingleIsotope( G4String name,
                           G4String symbol, G4double density, G4int Z, G4int A)
{
 // define a material from an isotope
 //
 G4int ncomponents;
 G4double abundance, massfraction;

 G4Isotope* isotope = new G4Isotope(symbol, Z, A);
 
 G4Element* element  = new G4Element(name, symbol, ncomponents=1);
 element->AddIsotope(isotope, abundance= 100.*perCent);
 
 G4Material* material = new G4Material(name, density, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);

 return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ComputeParameters()
{
  fAbsorSizeZ = 360*mm;
  fWorldSizeZ  = 100.*cm;
  fWorldSizeXY = 50.*cm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // complete the Calor parameters definition
  ComputeParameters();

  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  G4LogicalSkinSurface::CleanSurfaceTable();

  //
  // World
  //
  G4Box* solidWorld =
    new G4Box("World",                                             //name
               fWorldSizeXY/2,fWorldSizeXY/2,fWorldSizeZ/2);       //size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,              //solid
                        fDefaultMaterial,        //material
                        "World");                //name

  fPhysiWorld = 
    new G4PVPlacement(0,                        //no rotation
                        G4ThreeVector(),        //at (0,0,0)
                      logicWorld,               //logical volume
                      "World",                  //name
                       0,                       //mother volume
                       false,                   //no boolean operation
                       0);                      //copy number
                                 
  //
  // Absorbers
  //
  for (G4int k=1; k<=fNbOfAbsor; k++) {
    fAbsorMaterial[k] = fCalMaterial;
    G4String matname = fAbsorMaterial[k]->GetName();
         }


  fEcalLength   = 40.*cm;

  G4Box* solidE = new G4Box("VolE",fWorldSizeXY*.4,fWorldSizeXY*.4,fEcalLength/2);
  G4LogicalVolume* logicECal = 
    new G4LogicalVolume( solidE,fDefaultMaterial,"VolE");
  G4VPhysicalVolume* physE = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),
                                               "VolE",logicECal,fPhysiWorld,false,0);



 //PMTs
  G4double innerRadius_pmt = 0*mm;
  G4double fOuterRadius_pmt = 34*mm;
  G4double height_pmt = 70*mm;
  G4double startAngle_pmt = 0*deg;
  G4double spanningAngle_pmt = 360*deg;
 
  G4Tubs* fPmt = new G4Tubs("pmt_tube", innerRadius_pmt, fOuterRadius_pmt, height_pmt, startAngle_pmt, spanningAngle_pmt);
 
  G4Tubs* fPhotocath = new G4Tubs("photocath_tube",innerRadius_pmt,fOuterRadius_pmt,height_pmt/2,startAngle_pmt,spanningAngle_pmt);
 
 fPmt_log = new G4LogicalVolume(fPmt,G4Material::GetMaterial("Glass"),
                                 "pmt_log");
  fPhotocath_log = new G4LogicalVolume(fPhotocath,
                                       G4Material::GetMaterial("Al"),
                                       "photocath_log");

fPmt_log = new G4LogicalVolume(fPmt,fGlass,
                                 "pmt_log");
  fPhotocath_log = new G4LogicalVolume(fPhotocath,
                                       fAbsMaterial,
                                       "photocath_log");
 
  new G4PVPlacement(0,G4ThreeVector(0,0,-height_pmt/2),fPhotocath_log,"photocath",fPmt_log,false,0);


//Aluminium covered

G4double pDz = 170*mm;
G4double pPhi = 0*deg;
G4double pDx1 = 56.5*mm;
G4double pAlp1 = 0 *deg;
G4double pDx3 = 67.5*mm;
G4double pAlp2 =  0*deg;
G4double pTheta = 1.8*deg;
G4double pDy1 = 61*mm;
G4double pDx2 = 56.5*mm;
G4double pDy2 = 61*mm;
G4double pDx4 = 67.5*mm;
G4double width = 113*mm;
G4double height = 122*mm;
G4double gap = 0.2*mm;



G4Trap* solidAl = new G4Trap("ShapeAl",pDz + gap,pTheta,pPhi,pDy1 + gap,pDx1 + gap,pDx2 + gap,pAlp1,pDy2 + gap,pDx3 + gap,pDx4 + gap,pAlp2);


  fLogicAl = new G4LogicalVolume(solidAl,fAbsMaterial,"EcalAl");


// Crystals

G4Trap* solidC = new G4Trap("Shape2",pDz,pTheta,pPhi,pDy1,pDx1,pDx2,pAlp1,pDy2,pDx3,pDx4,pAlp2);


  fLogicCal = new G4LogicalVolume( solidC,fCalMaterial,"Ecal");


// Pb Counters (Crystal + Al cover + PMTs)

  G4double T = 3.6;
  G4double angle = std::sin(T*deg);
  G4double angle2 = std::cos(T*deg);


  G4double x0 = (-(width + gap)*1.0) + (-(pDz+gap)*angle);
  G4double y  = -(height + gap);
  G4double zheta0 = -(pDz+gap)*(1-angle2);
  G4double zheta;
  G4double x;
  G4int k = 1;
  G4int i,j;

    G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateY(-3.6*deg);

 for (i=0; i<3; i++) {
    x  = x0;
    zheta = zheta0;  
    for (j=0; j<3; j++) {
    G4ThreeVector pos2 = G4ThreeVector(x,y,zheta);
    G4Transform3D transform = G4Transform3D(rotm,pos2);
    G4ThreeVector pos3 = G4ThreeVector(x,y,zheta+pDz+height_pmt-2*gap);
    G4Transform3D transform3 = G4Transform3D(rotm,pos3);
    
      physAl = new G4PVPlacement(transform, "EcalAl", fLogicAl, physE, false, k);

      physPTM = new G4PVPlacement(transform3, "PTM", fPmt_log, physE, false, k);

    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),"Ecal",fLogicCal,physAl,false,k);
      k++;
     // T += 3.6;
      angle = std::sin(T*deg);
      angle2 = std::cos(T*deg);
      x += (width + 2*gap) +((pDz+gap)*angle);
      zheta += (pDz+gap)*(1-angle2);

       rotm.rotateY(3.6*deg);
    }
    y += height + 2*gap;
    rotm.rotateY(-(3.6*j)*deg);

  }


  // surface propoerties

   SurfaceProperties();

  // color regions

 logicECal-> SetVisAttributes(G4VisAttributes::GetInvisible());

  G4VisAttributes* regCcolor = new G4VisAttributes(G4Colour(0., 0.7, 0.3));
  fLogicCal->SetVisAttributes(regCcolor);
 // logicECal->SetVisAttributes(regCcolor);

  G4VisAttributes* Ccolor = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
  Ccolor -> SetForceAuxEdgeVisible(true);
  fPmt_log->SetVisAttributes(Ccolor);
  fLogicAl->SetVisAttributes(Ccolor);
  PrintParameters();

  //always return the physical World
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n-------------------------------------------------------------"
         << "\n ---> The Calorimeter is " << fNbOfAbsor << " crystals of:";
  for (G4int i=1; i<=fNbOfAbsor; i++)
     {
      G4cout << "\n \t" << std::setw(12) << fAbsorMaterial[i]->GetName() <<": "
              << std::setw(6) << G4BestUnit(fAbsorThickness[i],"Length");
     }
  G4cout << "\n-------------------------------------------------------------\n"
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void DetectorConstruction::SurfaceProperties(){
  G4double ephoton[] = {7.0*eV, 7.14*eV};
  const G4int num = sizeof(ephoton)/sizeof(G4double);

  //**Photocathode surface properties
  G4double photocath_EFF[]={1.,1.}; //Enables 'detection' of photons
  assert(sizeof(photocath_EFF) == sizeof(ephoton));
  G4double photocath_ReR[]={1.92,1.92};
  assert(sizeof(photocath_ReR) == sizeof(ephoton));
  G4double photocath_ImR[]={1.69,1.69};
  assert(sizeof(photocath_ImR) == sizeof(ephoton));
  G4MaterialPropertiesTable* photocath_mt = new G4MaterialPropertiesTable();
  photocath_mt->AddProperty("EFFICIENCY",ephoton,photocath_EFF,num);
  photocath_mt->AddProperty("REALRINDEX",ephoton,photocath_ReR,num);
  photocath_mt->AddProperty("IMAGINARYRINDEX",ephoton,photocath_ImR,num);
  G4OpticalSurface* photocath_opsurf=
    new G4OpticalSurface("photocath_opsurf",glisur,polished,
                         dielectric_metal);
  photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);

   //**Al surface properties
   G4double refl_Al[] = {1.,1.};
   assert(sizeof(refl_Al) == sizeof(ephoton));
   G4double effi_Al[] = {0, 0};
   assert(sizeof(effi_Al) == sizeof(ephoton));
   G4MaterialPropertiesTable* AlSurfaceProperty = new G4MaterialPropertiesTable();
   AlSurfaceProperty->AddProperty("REFLECTIVITY",ephoton,refl_Al,num);
   AlSurfaceProperty->AddProperty("EFFICIENCY",ephoton,effi_Al,num);
   G4OpticalSurface* AlSurface = 
    // new G4OpticalSurface("AlSurface",glisur,ground,dielectric_metal,1.);
    new G4OpticalSurface("AlSurface",glisur,polished,dielectric_metal);
   AlSurface -> SetMaterialPropertiesTable(AlSurfaceProperty);

  //**Create logical skin surfaces

  new G4LogicalSkinSurface("photocath_surf",fPhotocath_log,photocath_opsurf);
  new G4LogicalSkinSurface("Al_Surface", fLogicAl, AlSurface);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfAbsor(G4int ival)
{
  // set the number of Crystals 
  //
  if (ival < 1 || ival > (kMaxAbsor-1))
    { G4cout << "\n ---> warning from SetfNbOfAbsor: "
             << ival << " must be at least 1 and and most " << kMaxAbsor-1
             << ". Command refused" << G4endl;
      return;
    }
  fNbOfAbsor = ival;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorMaterial(const G4String& material)
{
      for(G4int i=1; i<=fNbOfAbsor; i++)
      {fAbsorMaterial[i] =fCalMaterial;}
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

