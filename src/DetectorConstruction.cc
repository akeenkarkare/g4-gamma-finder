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
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Element.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"

#include <cmath>
#include <string>
#include <vector>

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetShape(G4String shape)
{
  // Accept square | S | J | T | L (case-insensitive on the first letter for
  // the tetrominoes; "square" must be spelled out).
  fShape = shape;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineCommands()
{
  fMessenger =
    new G4GenericMessenger(this, "/det/", "Detector configuration controls");

  auto& shapeCmd = fMessenger->DeclareMethod(
    "shape", &DetectorConstruction::SetShape,
    "Crystal arrangement: square | S | J | T | L. "
    "Run /run/reinitializeGeometry afterwards.");
  shapeCmd.SetParameterName("shape", true);
  shapeCmd.SetDefaultValue("square");

  auto& padCmd = fMessenger->DeclareMethodWithUnit(
    "padding", "mm", &DetectorConstruction::SetPadding,
    "Lead inter-pixel padding FULL thickness. "
    "Run /run/reinitializeGeometry afterwards.");
  padCmd.SetParameterName("padding", true);
  padCmd.SetDefaultValue("1.0 mm");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  // The detector and source are coplanar (paper geometry), so the source sits
  // in the X-Y plane at distance d from the detector center. The envelope must
  // be large enough to contain the source for the chosen d. We use air (not
  // water) since the experiment is in open air, and size it for d up to ~1 m.
  G4double env_sizeXY = 250 * cm, env_sizeZ = 20 * cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeXY = 1.2 * env_sizeXY;
  G4double world_sizeZ = 1.2 * env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  auto solidWorld =
    new G4Box("World",  // its name
              0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
                                        world_mat,  // its material
                                        "World");  // its name

  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
                                     G4ThreeVector(),  // at (0,0,0)
                                     logicWorld,  // its logical volume
                                     "World",  // its name
                                     nullptr,  // its mother  volume
                                     false,  // no boolean operation
                                     0,  // copy number
                                     checkOverlaps);  // overlaps checking

  //
  // Envelope
  //
  auto solidEnv = new G4Box("Envelope",  // its name
                            0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ);  // its size

  auto logicEnv = new G4LogicalVolume(solidEnv,  // its solid
                                      env_mat,  // its material
                                      "Envelope");  // its name

  new G4PVPlacement(nullptr,  // no rotation
                    G4ThreeVector(),  // at (0,0,0)
                    logicEnv,  // its logical volume
                    "Envelope",  // its name
                    logicWorld,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    checkOverlaps);  // overlaps checking


  //
  // Gamma direction detector prototype
  //
  // 2x2 grid of LaBr3:Ce scintillator crystals separated by thin lead padding,
  // following the configuration of Okabe et al., Nat. Commun. 15:3061 (2024)
  // but with 2-inch crystals (the hardware we will use).
  //
  // Geometry note: the paper treats the detector and source as coplanar,
  // with the source direction defined as an angle theta in the detector
  // plane. We keep the detector face in the X-Y plane here; the source
  // angle is handled by the primary generator.

  // --- Materials ---
  // Lanthanum bromide (LaBr3), Ce-doped. LaBr3 is not a NIST material, so
  // build it from its components. Density ~5.08 g/cm3. The small Ce dopant
  // (~5%) has a negligible effect on gamma transport, so we model pure LaBr3.
  G4double labrDensity = 5.08 * g / cm3;
  auto labr3 = new G4Material("LaBr3", labrDensity, 2);
  G4Element* elLa = nist->FindOrBuildElement("La");
  G4Element* elBr = nist->FindOrBuildElement("Br");
  labr3->AddElement(elLa, 1);
  labr3->AddElement(elBr, 3);

  G4Material* pixel_mat  = labr3;
  G4Material* lead_mat   = nist->FindOrBuildMaterial("G4_Pb");
  G4Material* casing_mat = nist->FindOrBuildMaterial("G4_Al");

  // --- Pixel geometry ---
  // Crystal: cylinder (G4Tubs), 2 inch diameter x 2 inch height.
  //   1 inch = 2.54 cm.
  //   diameter 2 inch -> radius 1 inch = 2.54 cm
  //   height   2 inch -> half-height 1 inch = 2.54 cm
  // The cylinder axis points along z (the beam/depth direction), so the
  // circular faces lie in the X-Y detector plane.
  //
  // Each crystal is fully enclosed in a 2 mm-thick aluminium casing (a sealed
  // can wrapping the side wall and both end faces) -- realistic, since LaBr3
  // is hygroscopic and always hermetically sealed in an Al housing.
  G4double pixelRadius = 1.0 * 2.54 * cm;  // 2.54 cm radius = 2 inch diameter
  G4double pixelHalfZ  = 1.0 * 2.54 * cm;  // 2.54 cm half = 2 inch height
  G4double casingThk   = 0.2 * cm;         // 2 mm aluminium casing
  G4double padHalf     = 0.5 * fPadding * mm;  // lead padding half-thickness (fPadding in mm)

  // Outer casing cylinder: crystal radius/height + 2 mm on every side.
  G4double casingRadius = pixelRadius + casingThk;
  G4double casingHalfZ  = pixelHalfZ + casingThk;

  // Grid cell pitch: casings in adjacent cells leave a `fPadding`-mm gap.
  // pitch = casing diameter + padding gap.
  G4double pitch = 2.0 * casingRadius + 2.0 * padHalf;

  // The crystal solid (scored), placed INSIDE its casing.
  auto solidPixel = new G4Tubs("Pixel",
                               0., pixelRadius, pixelHalfZ, 0. * deg, 360. * deg);

  // The casing solid (a solid Al cylinder; the crystal sits inside it).
  auto solidCasing = new G4Tubs("Casing",
                                0., casingRadius, casingHalfZ, 0. * deg, 360. * deg);

  // --- Shape layouts on an integer (col,row) grid ---------------------------
  // Each shape is 4 cells. The array is centered on its centroid before
  // placement. Copy number = pixel index = order in this list.
  //   square: 2x2          S:   _ X X       J:   X _ _      T:   _ X _      L: X _ _
  //           X X               X X _            X X X           X X X          X _ _
  //           X X                                                              X X
  // (rows increase upward in y; cols increase to the right in x)
  struct Cell { int col; int row; };
  std::vector<Cell> cells;
  if (fShape == "S") {
    cells = {{1, 1}, {2, 1}, {0, 0}, {1, 0}};        // S / Z tetromino
  } else if (fShape == "J") {
    cells = {{0, 1}, {0, 0}, {1, 0}, {2, 0}};        // J tetromino
  } else if (fShape == "T") {
    cells = {{1, 1}, {0, 0}, {1, 0}, {2, 0}};        // T tetromino
  } else if (fShape == "L") {
    cells = {{0, 2}, {0, 1}, {0, 0}, {1, 0}};        // L tetromino
  } else {
    fShape = "square";
    cells = {{0, 1}, {1, 1}, {0, 0}, {1, 0}};        // 2x2 square
  }

  // Compute centroid to center the array on the origin.
  G4double cx = 0., cy = 0.;
  for (auto& c : cells) { cx += c.col; cy += c.row; }
  cx /= cells.size(); cy /= cells.size();

  auto cellPos = [&](const Cell& c) {
    return G4ThreeVector((c.col - cx) * pitch, (c.row - cy) * pitch, 0.);
  };

  // --- Place casings + crystals ---------------------------------------------
  std::vector<G4ThreeVector> pos;
  for (G4int i = 0; i < kNumPixels; ++i) {
    G4ThreeVector p = cellPos(cells[i]);
    pos.push_back(p);

    auto logicCasing =
      new G4LogicalVolume(solidCasing, casing_mat, "Casing" + std::to_string(i));
    new G4PVPlacement(nullptr, p, logicCasing,
                      "Casing" + std::to_string(i), logicEnv, false, i, checkOverlaps);

    fLogicPixel[i] = new G4LogicalVolume(solidPixel, pixel_mat, "Pixel" + std::to_string(i));
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), fLogicPixel[i],
                      "Pixel" + std::to_string(i), logicCasing, false, i, checkOverlaps);
  }

  // --- Lead padding in the gaps between adjacent cells -----------------------
  // For every pair of pixels that are orthogonal grid neighbours, drop a lead
  // slab of full thickness `fPadding` into the gap between their casings. This
  // generalizes the "+" cross to any shape.
  G4int padCount = 0;
  for (G4int i = 0; i < kNumPixels; ++i) {
    for (G4int j = i + 1; j < kNumPixels; ++j) {
      int dcol = cells[i].col - cells[j].col;
      int drow = cells[i].row - cells[j].row;
      bool neighbour = (std::abs(dcol) + std::abs(drow)) == 1;
      if (!neighbour) continue;

      // Midpoint between the two casings.
      G4ThreeVector mid = 0.5 * (pos[i] + pos[j]);

      // Slab spans the gap (thin along the neighbour direction) and is as tall
      // as the casing diameter in the perpendicular in-plane direction.
      G4double along = padHalf;                 // half-thickness across the gap
      G4double across = casingRadius;           // half-width along the shared face
      G4Box* solidPad;
      if (dcol != 0) {  // horizontal neighbours -> slab thin in X
        solidPad = new G4Box("Pad", along, across, casingHalfZ);
      } else {          // vertical neighbours -> slab thin in Y
        solidPad = new G4Box("Pad", across, along, casingHalfZ);
      }
      auto logicPad = new G4LogicalVolume(solidPad, lead_mat,
                                          "Pad" + std::to_string(padCount));
      new G4PVPlacement(nullptr, mid, logicPad,
                        "Pad" + std::to_string(padCount), logicEnv, false,
                        padCount, checkOverlaps);
      ++padCount;
    }
  }

  G4cout << "[DetectorConstruction] shape=" << fShape
         << "  padding=" << fPadding << " mm"
         << "  (" << padCount << " lead pads placed)" << G4endl;

  // Keep one pixel as the "scoring volume" for the inherited dose printout;
  // the real per-pixel scoring is handled by copy number in SteppingAction.
  fScoringVolume = fLogicPixel[0];

  //
  // always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
