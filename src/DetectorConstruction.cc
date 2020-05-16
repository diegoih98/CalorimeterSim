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

#include "G4UnionSolid.hh"
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

  G4NistManager* man = G4NistManager::Instance();

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



//Lead Glass Optical Propeties

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
{}

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

  ComputeParameters();

// Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  G4LogicalSkinSurface::CleanSurfaceTable();


  for (G4int k=1; k<=fNbOfAbsor; k++) {
    fAbsorMaterial[k] = fCalMaterial;
    G4String matname = fAbsorMaterial[k]->GetName();
         }

// World

  G4Box* worldSolid = new G4Box("World", fWorldSizeXY/2,fWorldSizeXY/2,fWorldSizeZ/2);   //size
  G4LogicalVolume* worldLV =  new G4LogicalVolume(worldSolid, fDefaultMaterial, "WorldLV");
  fPhysiWorld = 
    new G4PVPlacement(0,                        //no rotation
                        G4ThreeVector(),        //at (0,0,0)
                      worldLV,               //logical volume
                      "World",                  //name
                       0,                       //mother volume
                       false,                   //no boolean operation
                       0);                      //copy number



// Module  (moduleLV)
  //G4Trap* moduleTrapSolid = new G4Trap("ModuleTrapSolid", 122.4*mm, 340.2*mm, 135.4*mm, 113.4*mm);
  G4Trap* moduleTrapSolid = new G4Trap("ModuleTrapSolid", 126.*mm, 344.*mm, 139.*mm, 117.*mm);
  G4Box* moduleBoxSolid = new G4Box("ModuleBoxSolid", 69.5*mm, 100.*mm, 63*mm/*67.7*mm, 20.*mm, 61*mm*/);
  G4UnionSolid* moduleSolid = new G4UnionSolid("ModuleSolid", moduleTrapSolid, moduleBoxSolid,
                              G4Transform3D(G4RotationMatrix(), G4ThreeVector(5.5*mm, -272.*mm /*-190.1*mm*/, 0)));
  G4LogicalVolume* moduleLV = new G4LogicalVolume(moduleSolid, fDefaultMaterial, "ModuleLV");

  // Lead-glass (lgLV)
  G4Trap* lgTrapSolid = new G4Trap("LGTrapSolid", 122.*mm, 340.*mm, 135.*mm, 113.*mm);
  G4Tubs* lgTubeSolid = new G4Tubs("LGTubeSolid", 0, 38*mm, /*30.*mm*/15.*mm, 0, 360.*deg);
  G4UnionSolid* lgSolid = new G4UnionSolid("LGSolid", lgTrapSolid, lgTubeSolid,
                         G4Transform3D(G4RotationMatrix().rotateX(90.*deg), G4ThreeVector(5.5*mm, /*-200.*mm*/ -185.*mm, 0)));
  lgLV = new G4LogicalVolume(lgSolid,fCalMaterial, "LGLV");
  //lgLV = new G4LogicalVolume(lgTrapSolid,fCalMaterial, "LGLV");
  
  //Al plate "flange" (alPlateLV)
  G4Box* PlateSolid = new G4Box("PlateSolid", 67.5*mm, 0.2*mm, 61*mm);
  G4Tubs* TubeSolid = new G4Tubs("TubeSolid", 0, 38.*mm, 35.*mm, 0, 360.*deg);
  G4SubtractionSolid* alPlateSolid = new G4SubtractionSolid("AlPlateSolid", PlateSolid, TubeSolid, G4Transform3D(G4RotationMatrix().rotateX(90.*deg), G4ThreeVector(0,0,0)));
              alPlateLV = new G4LogicalVolume(alPlateSolid, fAbsMaterial, "AlPlateLV");

  //PMT Photocathode (cathodeLV)
  G4Tubs *cathodeSolid= new G4Tubs("CathodeSolid", 0, 38.*mm, /*1.1*mm*/ /*0.2*mm*/ 0.1*mm  , 0, 360.*deg);
  cathodeLV = new G4LogicalVolume(cathodeSolid, fCalMaterial, "CathodeLV");

  
  
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Placement

  G4Transform3D alPlatePlacement = G4Transform3D(G4RotationMatrix(), G4ThreeVector(5.5*mm, -170.3*mm /*-180.*mm*/, 0));
  G4Transform3D lgPlacement = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,-0.1*mm, 0));
  G4Transform3D cathodePlacement = G4Transform3D(G4RotationMatrix().rotateX(90.*deg), G4ThreeVector(5.5*mm, /*-170.3*mm*/ -200.2*mm, 0));

  G4bool bCheckOverlaps = false;
            alPlatePV = new G4PVPlacement(alPlatePlacement, alPlateLV, "ECalAlPlate", moduleLV, false, 0, bCheckOverlaps);
            lgPV = new G4PVPlacement(lgPlacement, lgLV, "ECalPbGlass", moduleLV, false, 1/*0*/, bCheckOverlaps);
  new G4PVPlacement(cathodePlacement, cathodeLV, "ECalPMTCathode", moduleLV, false, 0, bCheckOverlaps);

  G4RotationMatrix moduleRot = G4RotationMatrix().rotateX(-90.*deg);
  G4RotationMatrix columnRot;

  G4double tiltAngle = -3.70221285*deg;
  G4ThreeVector pos;

  for(G4int i = -1; i <= 1; i++){
    for(G4int j = -1; j <= 1; j++){
      G4int copyNo = 3*i + j + 5;
      
      columnRot = G4RotationMatrix().rotateY(j*tiltAngle);
      
      if(j < 0) pos = G4ThreeVector(124.253698*mm, i*122.4*mm, -4.37128422*mm);
      else if(j > 0) pos = G4ThreeVector(-124.277653*mm, i*122.4*mm, -3.66100488*mm);
      else pos = G4ThreeVector(0, i*122.4*mm, 0);
      
      modulePV[copyNo - 1] = new G4PVPlacement(G4Transform3D(columnRot*moduleRot, pos), 
                            moduleLV, "ECalWorld", worldLV, true, copyNo, bCheckOverlaps);
    }
  }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Surface Optical Porperties

     SurfaceProperties();
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Visualisation 

  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  moduleLV->SetVisAttributes(G4VisAttributes::GetInvisible());


   G4Colour myColour(G4Colour::White());
   G4VisAttributes* magnetVisAtt= new G4VisAttributes(myColour);
   magnetVisAtt ->SetForceSolid(true);
    //magnetVisAtt ->SetForceAuxEdgeVisible(true);
    lgLV ->SetVisAttributes(magnetVisAtt);

   G4Colour myColourB(G4Colour::Green());
  G4VisAttributes* Ccolor = new G4VisAttributes(myColourB);
  //Ccolor -> SetForceAuxEdgeVisible(true);
  Ccolor -> SetForceSolid(true);
  alPlateLV->SetVisAttributes(Ccolor);

  G4VisAttributes* colorr = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
  //colorr -> SetForceAuxEdgeVisible(true);
  colorr -> SetForceSolid(true);
   cathodeLV->SetVisAttributes(colorr);


PrintParameters();



  //always return the physical World
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

    G4double ephoton[] = {  1.379*eV, 1.551*eV, 1.909*eV, 1.960*eV, 2.332*eV, 
                            2.433*eV, 2.612*eV, 3.102*eV, 5.003*eV, 6.706*eV
                           };

  const G4int num = sizeof(ephoton)/sizeof(G4double);
    

   //Taped surface properties 

   G4double refl_Al[num];
  for(int i=0; i<num; i++) refl_Al[i] = 1.;
   assert(sizeof(refl_Al) == sizeof(ephoton));

  G4MaterialPropertiesTable* AlSurfaceProperty = new G4MaterialPropertiesTable();
   AlSurfaceProperty->AddProperty("REFLECTIVITY",ephoton,refl_Al,num);
   G4OpticalSurface* AlSurface = 
    new G4OpticalSurface("AlSurface",unified,polishedfrontpainted,dielectric_dielectric);
   AlSurface -> SetMaterialPropertiesTable(AlSurfaceProperty);


   //Al Surface
   G4double refl_Taped[num];
  for(int i=0; i<num; i++) refl_Taped[i] = 1.;
   assert(sizeof(refl_Taped) == sizeof(ephoton));

   G4MaterialPropertiesTable* TapedSurfaceProperty = new G4MaterialPropertiesTable();
   TapedSurfaceProperty->AddProperty("REFLECTIVITY",ephoton,refl_Taped,num);
   G4OpticalSurface* TapedSurface = 
    new G4OpticalSurface("TapedSurface",unified,groundfrontpainted,dielectric_dielectric);
   TapedSurface -> SetMaterialPropertiesTable(TapedSurfaceProperty);


  //Logical skin surfaces

    new G4LogicalBorderSurface("ECalMirrorBorderSurface", lgPV, alPlatePV, TapedSurface);
  for(G4int copyNo = 1; copyNo <= 9; copyNo++)
    new G4LogicalBorderSurface("LGTapedBorderSurface", lgPV, modulePV[copyNo - 1], AlSurface);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfAbsor(G4int ival)
{
  // set the number of Counters 
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
SensitiveDetector* aSD = new SensitiveDetector("SimSD");
  cathodeLV->SetSensitiveDetector(aSD);
  G4SDManager* sdManager = G4SDManager::GetSDMpointer();
  sdManager->AddNewDetector(aSD);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


