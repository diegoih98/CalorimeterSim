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
#include "SensitiveDetector.hh"
#include "EventAction.hh"

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
#include "G4SubtractionSolid.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4SDManager.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction()/*, fEventAction(0)*/,fDefaultMaterial(0),fPhysiWorld(0),
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
 // fEventAction = new EventAction(this);
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



  //Lead Glass

    G4double refEnergy[] = {1.24*eV, 1.27*eV, 1.29*eV, 1.32*eV, 1.35*eV, 
                            1.38*eV, 1.41*eV, 1.44*eV, 1.48*eV, 1.51*eV, 
                            1.55*eV, 1.59*eV, 1.63*eV, 1.68*eV, 1.72*eV, 
                            1.77*eV, 1.82*eV, 1.88*eV, 1.94*eV, 2.0*eV, 
                            2.07*eV, 2.14*eV, 2.21*eV, 2.3*eV, 2.38*eV, 
                            2.48*eV, 2.58*eV, 2.7*eV, 2.82*eV, 2.95*eV, 
                            3.1*eV, 3.26*eV, 3.44*eV, 3.65*eV, 3.87*eV,
                            4.13*eV, 4.43*eV, 4.77*eV, 5.17*eV, 5.64*eV, 6.2*eV};
    G4double refIndex[] = {1.77568, 1.77635, 1.77705, 1.77779, 1.77858, 
                            1.77941, 1.78029, 1.78123, 1.78223, 1.7833, 
                            1.78445, 1.78569, 1.78703, 1.78847, 1.79004, 
                            1.79176, 1.79363, 1.79569, 1.79796, 1.80048, 
                            1.80328, 1.80641, 1.80993, 1.81391, 1.81844, 
                            1.82365, 1.82968, 1.83673, 1.84508, 1.85511, 
                            1.86739, 1.88274, 1.9025, 1.92879, 1.96527,
                            2.01844, 2.10029, 2.23367, 2.4627, 2.87238, 3.62362};           
    const G4int refNbins = sizeof(refEnergy) / sizeof(G4double);
    assert(sizeof(refEnergy) == sizeof(refIndex));
    

    G4double absEnergy[] = {0.62*eV, 0.83*eV, 1.17*eV, 1.55*eV, 1.77*eV, 
                            1.91*eV, 2.07*eV, 2.25*eV, 2.48*eV, 2.58*eV, 
                            2.7*eV, 2.82*eV, 2.95*eV, 3.1*eV, 3.18*eV, 
                            3.26*eV, 3.35*eV, 3.44*eV};
    G4double absLength[] = {79.32*cm, 141.01*cm, 188.99*cm, 572.76*cm, 572.76*cm, 
                            284.94*cm, 284.94*cm, 284.94*cm, 188.99*cm, 188.99*cm, 
                            112.23*cm, 93.03*cm, 54.64*cm, 25.80*cm, 14.93*cm, 
                            7.41*cm, 3.14*cm, 1.04*cm};
    const G4int absNbins = sizeof(absEnergy) / sizeof(G4double);
    assert(sizeof(absEnergy) == sizeof(absLength));

    G4MaterialPropertiesTable* lgMPT = new G4MaterialPropertiesTable();
    lgMPT->AddProperty("RINDEX", refEnergy, refIndex, refNbins)->SetSpline(true);
    lgMPT->AddProperty("ABSLENGTH", absEnergy, absLength, absNbins)->SetSpline(true);
    fCalMaterial->SetMaterialPropertiesTable(lgMPT);

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
  fWorldSizeZ  = (3*100.)*cm;
  fWorldSizeXY = (5*50.)*cm;

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
//G4double gap = 0.2*mm;
//G4double gap = 3.0*mm;
G4double gap = 2.0*mm;
G4double gap2 = 3.0*mm;

G4double TubeG =18.85*cm;
G4double CathG =20.8*cm;
G4double Lateral = 135*mm;
G4double PlateG = 17.1*cm;


// Crystals

    G4Trap* solidC = new G4Trap("Shape2",pDz,pTheta,pPhi,pDy1,pDx1,pDx2,pAlp1,pDy2,pDx3,pDx4,pAlp2);
    fLogicCal = new G4LogicalVolume( solidC,fCalMaterial,"Ecal");


    G4Tubs *lgTubeSolid = new G4Tubs("LGTubeSolid", 0, 4.3*cm, 1.85*cm, 0, 360.*deg);
   // G4LogicalVolume* lgTubeLV = new G4LogicalVolume(lgTubeSolid, fCalMaterial, "LGTubeLV");
   lgTubeLV = new G4LogicalVolume(lgTubeSolid, fCalMaterial, "LGTubeLV");


    // Al plate
    G4Box* alBoxSolid = new G4Box("AlBoxSolid", 6.75*cm, 6.1*cm, 0.1*cm);
    G4Tubs* alTubeSolid = new G4Tubs("AlTubeSolid", 0, 4.3*cm, 0.1*cm, 0, 360.*deg);
    G4VSolid* alPlateSolid = new G4SubtractionSolid("AlPlateSolid", alBoxSolid, alTubeSolid, G4Transform3D());
    alPlateLV = new G4LogicalVolume(alPlateSolid, fAbsMaterial, "AlPlateLV");  



  /*  // Lead-glass mirror surface
    G4OpticalSurface* lgMirrorSurface = new G4OpticalSurface("LGMirrorSurface");
    lgMirrorSurface->SetType(dielectric_LUT);
    lgMirrorSurface->SetModel(LUT);
    lgMirrorSurface->SetFinish(polishedlumirrorair);*/


    // Lead-glass taped surface
    G4OpticalSurface* lgTapedSurface = new G4OpticalSurface("LGTapedSurface");
    lgTapedSurface->SetType(dielectric_LUT);
    lgTapedSurface->SetModel(LUT);
    lgTapedSurface->SetFinish(groundfrontpainted);



    // PMT Photocathode as SD
    G4Tubs *cathodeSolid= new G4Tubs("CathodeSolid", 0, 4.3*cm, 0.1*cm, 0, 360.*deg);
    cathodeLV = new G4LogicalVolume(cathodeSolid, fCalMaterial, "CathodeLV");
   // cathodeLV = new G4LogicalVolume(cathodeSolid, fAbsMaterial, "CathodeLV");





// Pb Counters (Crystal + Al cover + PMTs)

  G4double T = 3.6;
  G4double angle = std::sin(T*deg);
  G4double angle2 = std::cos(T*deg);
  G4double dis = std::tan(T*deg)*pDz;

  G4double x0 = (-(width + gap)*1.0) + (-(pDz+gap)*angle);
  G4double x01 = (-(Lateral)*1.0 + gap2/2 +4*gap2) + (-(pDz+gap2)*angle);

  G4double y  = -(height + gap);


  G4double zheta0 = -(pDz+gap)*(1-angle2);
  G4double zheta01 = -(pDz+gap+PlateG)*(1-angle2);

  G4double zheta;
  G4double zheta1;

  G4double x;
  G4double x1;
  G4int k = 1;
  G4int i,j;

    G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateY(-3.6*deg);

 for (i=0; i<3; i++) {
    x  = x0;
    x1  = x01;
    zheta = zheta0;  
    zheta1 = zheta01;  
    for (j=0; j<3; j++) {

    G4ThreeVector pos2 = G4ThreeVector(x,y,zheta);
    G4Transform3D transform = G4Transform3D(rotm,pos2);

   // G4ThreeVector pos3 = G4ThreeVector(x1,y,zheta+TubeG);
    G4ThreeVector pos3 = G4ThreeVector(x1,y,zheta+TubeG-3.0*mm);
    G4Transform3D transform3 = G4Transform3D(rotm,pos3);

    G4ThreeVector pos4 = G4ThreeVector(x,y,zheta-(gap2/2));
    G4Transform3D transform4 = G4Transform3D(rotm,pos4);

    G4ThreeVector pos6 = G4ThreeVector(x1,y,zheta+17.1*cm -2.0*mm);
    G4Transform3D transform6 = G4Transform3D(rotm,pos6);

    //G4ThreeVector pos5 = G4ThreeVector(x1,y,zheta+CathG);
    G4ThreeVector pos5 = G4ThreeVector(x1,y,zheta+CathG-3.0*mm);
    G4Transform3D transform5 = G4Transform3D(rotm,pos5);


    new G4PVPlacement(transform3, "LGTubePV", lgTubeLV, fPhysiWorld, false, k);
    alPlatePV[k] = new G4PVPlacement(transform6, "AlPlatePV", alPlateLV, fPhysiWorld, false, k);
    new G4PVPlacement(transform5, "CathodePV", cathodeLV, fPhysiWorld, false, k);
    physCal[k] = new G4PVPlacement(transform,"Ecal",fLogicCal,fPhysiWorld,false,k);

      k++;
      angle = std::sin(T*deg);
      angle2 = std::cos(T*deg);
      x += (width + 2*gap) +((pDz+gap)*angle);
      x1 += (Lateral + -4*gap2 + 2*gap2) +((pDz+gap2)*angle);
      zheta += (pDz+gap)*(1-angle2);
      zheta1 += (pDz+gap+PlateG)*(1-angle2);

       rotm.rotateY(3.6*deg);
    }
    y += height + 2*gap;
    rotm.rotateY(-(3.6*j)*deg);

  }


  // surface propoerties

 for (G4int l=1; l<=9; l++) {
   // new G4LogicalBorderSurface("LGMirrorBorderSurface", physCal[k], fPhysiWorld, lgMirrorSurface);
    new G4LogicalBorderSurface("LGTapedBorderSurface", physCal[k], alPlatePV[k], lgTapedSurface);
}
  SurfaceProperties();

  // color regions

   G4Colour myColour(G4Colour::White());
   G4VisAttributes* magnetVisAtt= new G4VisAttributes(myColour);
   magnetVisAtt ->SetForceSolid(true);
    //magnetVisAtt ->SetForceAuxEdgeVisible(true);
    fLogicCal ->SetVisAttributes(magnetVisAtt);

 /* G4VisAttributes* regCcolor = new G4VisAttributes(G4Colour(0., 0.7, 0.3));
  alShieldingLV->SetVisAttributes(regCcolor);*/

   G4Colour myColourB(G4Colour::Green());
  G4VisAttributes* Ccolor = new G4VisAttributes(myColourB);
 // Ccolor -> SetForceAuxEdgeVisible(true);
  Ccolor -> SetForceSolid(true);
  alPlateLV->SetVisAttributes(Ccolor);

  G4VisAttributes* color = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
 // color -> SetForceAuxEdgeVisible(true);
  color -> SetForceSolid(true);
  lgTubeLV->SetVisAttributes(color);

  G4VisAttributes* colorr = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
  //colorr -> SetForceAuxEdgeVisible(true);
  colorr -> SetForceSolid(true);
   cathodeLV->SetVisAttributes(colorr);

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

    G4double ephoton[] = {  0.62*eV, 0.83*eV, 1.17*eV, 1.24*eV, 1.27*eV, 
                            1.29*eV, 1.32*eV, 1.35*eV, 1.38*eV, 1.41*eV,
                            1.44*eV, 1.48*eV, 1.51*eV, 1.55*eV, 1.59*eV, 
                            1.63*eV, 1.68*eV, 1.72*eV, 1.77*eV, 1.82*eV, 
                            1.88*eV, 1.91*eV, 1.94*eV, 2.0*eV,  2.07*eV, 
                            2.14*eV, 2.21*eV, 2.25*eV, 2.3*eV,  2.38*eV, 
                            2.48*eV, 2.58*eV, 2.7*eV,  2.82*eV, 2.95*eV, 
                            3.1*eV,  3.18*eV, 3.26*eV, 3.35*eV, 3.44*eV,
                            3.65*eV, 3.87*eV, 4.13*eV, 4.43*eV, 4.77*eV, 
                            5.17*eV, 5.64*eV, 6.2*eV};
  const G4int num = sizeof(ephoton)/sizeof(G4double);
    
   //Al surface properties 
   G4double refl_Al[] = {1.,1.,1.,1.,1.,
                         1.,1.,1.,1.,1.,
                         1.,1.,1.,1.,1.,
                         1.,1.,1.,1.,1.,
                         1.,1.,1.,1.,1.,
                         1.,1.,1.,1.,1.,
                         1.,1.,1.,1.,1.,
                         1.,1.,1.,1.,1.,
                         1.,1.,1.,1.,1.,
                         1.,1.,1.};
   assert(sizeof(refl_Al) == sizeof(ephoton));
   G4double effi_Al[] = {0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0};
   assert(sizeof(effi_Al) == sizeof(ephoton));
   G4MaterialPropertiesTable* AlSurfaceProperty = new G4MaterialPropertiesTable();
   AlSurfaceProperty->AddProperty("REFLECTIVITY",ephoton,refl_Al,num);
   AlSurfaceProperty->AddProperty("EFFICIENCY",ephoton,effi_Al,num);
   G4OpticalSurface* AlSurface = 
    new G4OpticalSurface("AlSurface",glisur,polished,dielectric_metal);
   AlSurface -> SetMaterialPropertiesTable(AlSurfaceProperty);



   //Taped Surface
   G4double refl_Taped[] = {0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0,
                            0, 0, 0};
   assert(sizeof(refl_Taped) == sizeof(ephoton));
   //G4double rindex[] = {0.,0.};
   //assert(sizeof(effi_Al) == sizeof(ephoton));
   G4MaterialPropertiesTable* TapedSurfaceProperty = new G4MaterialPropertiesTable();
   TapedSurfaceProperty->AddProperty("REFLECTIVITY",ephoton,refl_Taped,num);
   //TapedSurfaceProperty->AddProperty("EFFICIENCY",ephoton,effi_Al,num);
   G4OpticalSurface* TapedSurface = 
    new G4OpticalSurface("TapedSurface",unified,polished,dielectric_dielectric);
   TapedSurface -> SetMaterialPropertiesTable(TapedSurfaceProperty);

  //Photocathode surface properties
 /* G4double photocath_EFF[]={1.,1.}; //Enables 'detection' of photons
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
  photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);*/

  //**Create logical skin surfaces
  //new G4LogicalSkinSurface("Al_Surface", alShieldingLV, AlSurface);
 

  new G4LogicalSkinSurface("Al_Surface", fLogicCal, AlSurface);
  //new G4LogicalSkinSurface("Taped_Surface", alPlateLV, TapedSurface); //(This is good)
  //new G4LogicalSkinSurface("photocath_surf",lgTubeLV,TapedSurface);
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
{ 
   /* SensitiveDetector* pmtSD = new SensitiveDetector*("PMTSD");
    G4SDManager::GetSDMpointer()->AddNewDetector(pmtSD);
    SetSensitiveDetector("LGTubeLV", pmtSD);*/



 SensitiveDetector* aSD = new SensitiveDetector("SimSD");
 // cathodeLV->SetSensitiveDetector(aSD);
  lgTubeLV->SetSensitiveDetector(aSD);
  G4SDManager* sdManager = G4SDManager::GetSDMpointer();
  sdManager->AddNewDetector(aSD);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

