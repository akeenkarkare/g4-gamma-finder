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
#include "G4RotationMatrix.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"

#include <cmath>
#include <string>
#include <vector>

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Active crystal count, set when the shape is built. Default 4 until Construct().
G4int DetectorConstruction::fNumPixels = 4;

G4int DetectorConstruction::GetNumPixels() { return fNumPixels; }

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
    "Crystal arrangement: square|S|J|T|L (4) | P|F|W (5) | grid6|A6 (6). "
    "Run /run/reinitializeGeometry afterwards.");
  shapeCmd.SetParameterName("shape", true);
  shapeCmd.SetDefaultValue("square");

  auto& padCmd = fMessenger->DeclareMethodWithUnit(
    "padding", "mm", &DetectorConstruction::SetPadding,
    "Lead inter-pixel padding FULL thickness. "
    "Run /run/reinitializeGeometry afterwards.");
  padCmd.SetParameterName("padding", true);
  padCmd.SetDefaultValue("1.0 mm");

  auto& asymCmd = fMessenger->DeclareMethodWithUnit(
    "asymCasing", "mm", &DetectorConstruction::SetAsymCasing,
    "Extra aluminium thickness on the -y half of each crystal casing (mm). "
    "0 = uniform casing. Breaks rotational symmetry to make a single crystal "
    "directionally sensitive.");
  asymCmd.SetParameterName("thk", true);
  asymCmd.SetDefaultValue("0 mm");

  auto& asymMatCmd = fMessenger->DeclareMethod(
    "asymMaterial", &DetectorConstruction::SetAsymMaterial,
    "Material of the asymmetric directional shield: Al | Pb | W.");
  asymMatCmd.SetParameterName("mat", true);
  asymMatCmd.SetDefaultValue("Al");
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

  // Optional ASYMMETRIC casing shaping: a thick Al half-arc added on the -y
  // side of each crystal. This breaks the casing's rotational symmetry so a
  // single crystal's count rate depends on the source direction (gammas from
  // +y pass thin Al; from -y they cross the thick arc). Built as a 180-deg
  // phi-segment ring of extra Al hugging the outside of the casing.
  G4LogicalVolume* logicAsym = nullptr;
  if (fAsymCasingThk > 0.) {
    G4double extra = fAsymCasingThk;  // already in Geant4 length units (mm) from the messenger
    // phi=0 is +x, 90 +y, 180 -x, 270 -y. Span [180,360] deg = the -y half
    // (the lower half-plane, centered on -y at 270 deg).
    auto solidAsym = new G4Tubs("AsymShield",
                                casingRadius,                 // inner = casing outer
                                casingRadius + extra,         // outer = + extra Al
                                casingHalfZ,
                                180. * deg, 180. * deg);      // start 180, span 180 -> -y half
    G4Material* asym_mat = casing_mat;  // default Al
    if (fAsymMaterial == "Pb") asym_mat = lead_mat;
    else if (fAsymMaterial == "W") asym_mat = nist->FindOrBuildMaterial("G4_W");
    logicAsym = new G4LogicalVolume(solidAsym, asym_mat, "AsymShield");
  }

  // --- Shape layouts on an integer (col,row) grid ---------------------------
  // Each shape is 4 cells. The array is centered on its centroid before
  // placement. Copy number = pixel index = order in this list.
  //   square: 2x2          S:   _ X X       J:   X _ _      T:   _ X _      L: X _ _
  //           X X               X X _            X X X           X X X          X _ _
  //           X X                                                              X X
  // (rows increase upward in y; cols increase to the right in x)
  struct Cell { int col; int row; };
  std::vector<Cell> cells;
  // --- single crystal (for casing-asymmetry tests) ---
  if (fShape == "single") {
    cells = {{0, 0}};
  // --- 3-crystal arrays (for the 'fewer shaped detectors' study) ---
  } else if (fShape == "tri") {
    cells = {{0, 0}, {2, 0}, {1, 1}};   // triangle (compact, asymmetric)
  } else if (fShape == "L3") {
    cells = {{0, 1}, {0, 0}, {1, 0}};   // L-tromino
  } else if (fShape == "tri3") {
    // Equilateral-ish triangle: 3 crystals ~120 deg apart around the centre, so
    // their outward-facing (shielded-inward) sectors partition the circle.
    // Positions on a ring of radius ~1.15 cells: (0,2),(2,-1),(-2,-1)/scaled.
    cells = {{0, 2}, {2, -1}, {-2, -1}};
  // --- 4-crystal (tetromino) shapes ---
  } else if (fShape == "S") {
    cells = {{1, 1}, {2, 1}, {0, 0}, {1, 0}};        // S / Z tetromino
  } else if (fShape == "J") {
    cells = {{0, 1}, {0, 0}, {1, 0}, {2, 0}};        // J tetromino
  } else if (fShape == "T") {
    cells = {{1, 1}, {0, 0}, {1, 0}, {2, 0}};        // T tetromino
  } else if (fShape == "L") {
    cells = {{0, 2}, {0, 1}, {0, 0}, {1, 0}};        // L tetromino

  // --- 5-crystal: all 12 free pentominoes ---
  } else if (fShape == "pF") {
    cells = {{1, 2}, {2, 2}, {0, 1}, {1, 1}, {1, 0}};
  } else if (fShape == "pI") {
    cells = {{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}};
  } else if (fShape == "pL") {
    cells = {{0, 0}, {0, 1}, {0, 2}, {0, 3}, {1, 0}};
  } else if (fShape == "pN") {
    cells = {{0, 0}, {0, 1}, {1, 1}, {1, 2}, {1, 3}};
  } else if (fShape == "pP") {
    cells = {{0, 0}, {0, 1}, {1, 0}, {1, 1}, {1, 2}};
  } else if (fShape == "pT") {
    cells = {{0, 2}, {1, 2}, {2, 2}, {1, 1}, {1, 0}};
  } else if (fShape == "pU") {
    cells = {{0, 0}, {0, 1}, {1, 0}, {2, 0}, {2, 1}};
  } else if (fShape == "pV") {
    cells = {{0, 0}, {0, 1}, {0, 2}, {1, 0}, {2, 0}};
  } else if (fShape == "pW") {
    cells = {{0, 2}, {0, 1}, {1, 1}, {1, 0}, {2, 0}};
  } else if (fShape == "pX") {
    cells = {{1, 2}, {0, 1}, {1, 1}, {2, 1}, {1, 0}};
  } else if (fShape == "pY") {
    cells = {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {0, 1}};
  } else if (fShape == "pZ") {
    cells = {{0, 2}, {1, 2}, {1, 1}, {1, 0}, {2, 0}};

  // --- 6-crystal hexominoes (asymmetric-weighted + 1 symmetric baseline) ---
  } else if (fShape == "hGrid") {       // 2x3 rectangle (symmetric baseline)
    cells = {{0, 1}, {1, 1}, {2, 1}, {0, 0}, {1, 0}, {2, 0}};
  } else if (fShape == "hA") {          // staggered, no symmetry axis
    cells = {{1, 2}, {0, 1}, {1, 1}, {2, 1}, {2, 0}, {0, 0}};
  } else if (fShape == "hL") {          // long L (asymmetric)
    cells = {{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 0}};
  } else if (fShape == "hY") {          // Y-hexomino (asymmetric)
    cells = {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {0, 2}};
  } else if (fShape == "hN") {          // extended N / S (asymmetric)
    cells = {{0, 0}, {0, 1}, {0, 2}, {1, 2}, {1, 3}, {1, 4}};
  } else if (fShape == "hZ") {          // extended Z (point asymmetric)
    cells = {{0, 0}, {1, 0}, {1, 1}, {1, 2}, {1, 3}, {2, 3}};
  } else if (fShape == "hW") {          // staircase (asymmetric)
    cells = {{0, 0}, {1, 0}, {1, 1}, {2, 1}, {2, 2}, {3, 2}};
  } else if (fShape == "hF") {          // F-like extension (asymmetric)
    cells = {{1, 3}, {2, 3}, {1, 2}, {0, 1}, {1, 1}, {1, 0}};

  } else {
    fShape = "square";
    cells = {{0, 1}, {1, 1}, {0, 0}, {1, 0}};        // 2x2 square
  }

  // Record the active crystal count for the rest of the framework.
  fNumPixels = static_cast<G4int>(cells.size());

  // Compute centroid to center the array on the origin.
  G4double cx = 0., cy = 0.;
  for (auto& c : cells) { cx += c.col; cy += c.row; }
  cx /= cells.size(); cy /= cells.size();

  auto cellPos = [&](const Cell& c) {
    return G4ThreeVector((c.col - cx) * pitch, (c.row - cy) * pitch, 0.);
  };

  // --- Place casings + crystals ---------------------------------------------
  std::vector<G4ThreeVector> pos;
  for (G4int i = 0; i < fNumPixels; ++i) {
    G4ThreeVector p = cellPos(cells[i]);
    pos.push_back(p);

    auto logicCasing =
      new G4LogicalVolume(solidCasing, casing_mat, "Casing" + std::to_string(i));
    new G4PVPlacement(nullptr, p, logicCasing,
                      "Casing" + std::to_string(i), logicEnv, false, i, checkOverlaps);

    fLogicPixel[i] = new G4LogicalVolume(solidPixel, pixel_mat, "Pixel" + std::to_string(i));
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), fLogicPixel[i],
                      "Pixel" + std::to_string(i), logicCasing, false, i, checkOverlaps);

    // Optional asymmetric directional shield arc. The shield is placed on the
    // crystal's INWARD side (facing the array centre) so each crystal is OPEN to
    // its own outward sector. This partitions the circle into distinct sectors
    // and breaks the front/back degeneracy a plain array suffers.
    // The shield solid is a half-arc centred on -y (270 deg); we rotate it so
    // its centre points from the crystal toward the array centroid (origin).
    if (logicAsym) {
      G4double inwardAng = std::atan2(-p.y(), -p.x());  // direction crystal->centre
      // Default arc centre is at 270 deg (-y). Rotate so the arc centre aligns
      // with the inward direction.
      G4double rotZ = inwardAng - (-90. * deg);  // -90 deg = -y in atan2 convention
      auto rot = new G4RotationMatrix();
      rot->rotateZ(-rotZ);
      new G4PVPlacement(rot, p, logicAsym,
                        "AsymShield" + std::to_string(i), logicEnv, false, i, checkOverlaps);
    }
  }

  // --- Lead padding in the gaps between adjacent cells -----------------------
  // For every pair of pixels that are orthogonal grid neighbours, drop a lead
  // slab of full thickness `fPadding` into the gap between their casings. This
  // generalizes the "+" cross to any shape.
  G4int padCount = 0;
  for (G4int i = 0; i < fNumPixels; ++i) {
    for (G4int j = i + 1; j < fNumPixels; ++j) {
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
