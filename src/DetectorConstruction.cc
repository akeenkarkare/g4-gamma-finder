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

#include <cctype>
#include <cmath>
#include <sstream>
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

  auto& asymModeCmd = fMessenger->DeclareMethod(
    "asymMode", &DetectorConstruction::SetAsymMode,
    "Shield geometry: 'arc' (open half + shielded half) or "
    "'window' (near-full lead ring with a hole/aperture facing outward).");
  asymModeCmd.SetParameterName("mode", true);
  asymModeCmd.SetDefaultValue("arc");

  auto& zStagCmd = fMessenger->DeclareMethodWithUnit(
    "zStagger", "cm", &DetectorConstruction::SetZStagger,
    "Per-crystal vertical (z) step. Lifts crystal i by i*step to make the "
    "array 3D/chiral and break the planar front/back degeneracy. 0 = coplanar.");
  zStagCmd.SetParameterName("z", true);
  zStagCmd.SetDefaultValue("0 cm");

  auto& cellsCmd = fMessenger->DeclareMethod(
    "cells", &DetectorConstruction::SetCells,
    "Custom polyomino as 'col,row;col,row;...' (integer grid cells). Sets shape "
    "to 'custom'. Lets any arrangement (e.g. heptominoes) be built without "
    "recompiling. Up to kMaxPixels crystals.");
  cellsCmd.SetParameterName("spec", true);
  cellsCmd.SetDefaultValue("");
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
  G4double env_sizeXY = 250 * cm, env_sizeZ = 60 * cm;  // tall enough for z-staggered arrays
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

  // Tungsten heavy alloy (W-Ni-Fe), a dense, machinable, non-toxic shielding
  // material -- denser than lead (~18 vs 11.3 g/cm3) so it shields more per mm.
  // Typical "90W" composition: 90% W, 6% Ni, 4% Fe by mass.
  G4double wnifeDensity = 18.0 * g / cm3;
  auto wnife = new G4Material("WNiFe", wnifeDensity, 3);
  wnife->AddElement(nist->FindOrBuildElement("W"),  0.90);
  wnife->AddElement(nist->FindOrBuildElement("Ni"), 0.06);
  wnife->AddElement(nist->FindOrBuildElement("Fe"), 0.04);

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
  // Build the directional shield. Convention: the shield's SENSITIVE direction
  // (where gammas can reach the crystal) is at +y (90 deg) by construction; the
  // placement code then rotates each crystal's shield to point its sensitive
  // direction OUTWARD (away from the array centre).
  //   - "arc":    shielded half on -y, open half on +y  (broad acceptance)
  //   - "window": near-full lead ring with a small HOLE on +y (narrow aperture,
  //               solid lead on the opposite side -> breaks front/back)
  G4LogicalVolume* logicAsym = nullptr;
  if (fAsymCasingThk > 0.) {
    G4double extra = fAsymCasingThk;  // Geant4 length units (mm) from the messenger
    G4Material* asym_mat = casing_mat;  // default Al
    if (fAsymMaterial == "Pb") asym_mat = lead_mat;
    else if (fAsymMaterial == "W") asym_mat = nist->FindOrBuildMaterial("G4_W");
    else if (fAsymMaterial == "WNiFe") asym_mat = wnife;

    G4Tubs* solidAsym;
    if (fAsymMode == "window") {
      // Lead ring covering 300 deg, leaving a 60-deg HOLE centred on +y (90 deg).
      // Ring spans phi [120, 420] = [120, 360]+[0,60] -> hole is [60,120] around +y.
      // G4Tubs takes (start, span): start 120 deg, span 300 deg.
      solidAsym = new G4Tubs("AsymShield",
                             casingRadius, casingRadius + extra, casingHalfZ,
                             120. * deg, 300. * deg);
    } else {
      // "arc": shielded half on -y (phi [180,360]); open half on +y.
      solidAsym = new G4Tubs("AsymShield",
                             casingRadius, casingRadius + extra, casingHalfZ,
                             180. * deg, 180. * deg);
    }
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
  // --- custom polyomino from /det/cells "col,row;col,row;..." ---
  // Lets any arrangement (e.g. all 108 heptominoes) be specified at runtime
  // without recompiling. Parsed into integer (col,row) cells.
  if (fShape == "custom" && !fCellSpec.empty()) {
    std::stringstream ss(fCellSpec);
    std::string tok;
    while (std::getline(ss, tok, ';')) {
      if (tok.empty()) continue;
      auto comma = tok.find(',');
      if (comma == std::string::npos) continue;
      int col = std::stoi(tok.substr(0, comma));
      int row = std::stoi(tok.substr(comma + 1));
      cells.push_back({col, row});
    }
  // --- single crystal (for casing-asymmetry tests) ---
  } else if (fShape == "single") {
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
  // --- 3-crystal arrays with UNEQUAL spacing (tests the spacing lever) ---
  } else if (fShape == "tri_uneq1") {
    // Collinear but unequal gaps: close pair + far one. Strong 1/r^2 asymmetry.
    cells = {{0, 0}, {1, 0}, {4, 0}};
  } else if (fShape == "tri_uneq2") {
    // Right-angle, unequal arms: one neighbour near, one far. Breaks more symmetry.
    cells = {{0, 0}, {1, 0}, {0, 3}};
  } else if (fShape == "tri_uneq3") {
    // Scalene triangle: all three pairwise gaps different.
    cells = {{0, 0}, {3, 1}, {1, 3}};
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

  // --- 6-crystal: all 35 free hexominoes (named h01..h35) ---
  // Each entry is a canonical 6-cell polyomino. Used for the crystal-count knee
  // study so the BEST 6-shape is found by exhaustive search, not guessed.
  } else if (fShape.rfind("h", 0) == 0 && fShape.size() == 3 &&
             std::isdigit(fShape[1]) && std::isdigit(fShape[2])) {
    static const std::vector<std::vector<Cell>> HEX = {
      /*h01 I*/ {{0,0},{1,0},{2,0},{3,0},{4,0},{5,0}},
      /*h02*/   {{0,0},{1,0},{2,0},{3,0},{4,0},{0,1}},
      /*h03*/   {{0,0},{1,0},{2,0},{3,0},{4,0},{1,1}},
      /*h04*/   {{0,0},{1,0},{2,0},{3,0},{4,0},{2,1}},
      /*h05*/   {{0,0},{1,0},{2,0},{3,0},{0,1},{1,1}},
      /*h06*/   {{0,0},{1,0},{2,0},{3,0},{0,1},{2,1}},
      /*h07*/   {{0,0},{1,0},{2,0},{3,0},{0,1},{3,1}},
      /*h08*/   {{0,0},{1,0},{2,0},{3,0},{1,1},{2,1}},
      /*h09*/   {{0,0},{1,0},{2,0},{3,0},{1,1},{3,1}},
      /*h10*/   {{0,0},{1,0},{2,0},{3,0},{2,1},{3,1}},
      /*h11*/   {{1,0},{2,0},{3,0},{0,1},{1,1},{2,1}},
      /*h12*/   {{0,0},{1,0},{2,0},{1,1},{2,1},{3,1}},
      /*h13*/   {{0,0},{1,0},{2,0},{0,1},{1,1},{2,1}}, // 3x2 block (alt)
      /*h14*/   {{0,1},{1,1},{2,1},{2,0},{3,0},{4,0}},
      /*h15*/   {{0,0},{1,0},{2,0},{2,1},{3,1},{3,2}},
      /*h16*/   {{0,0},{1,0},{2,0},{0,1},{0,2},{0,3}}, // L
      /*h17*/   {{0,0},{1,0},{2,0},{2,1},{2,2},{2,3}},
      /*h18*/   {{0,3},{0,2},{0,1},{0,0},{1,0},{2,0}},
      /*h19*/   {{0,0},{0,1},{0,2},{0,3},{1,3},{2,3}},
      /*h20*/   {{0,0},{0,1},{0,2},{0,3},{1,0},{1,3}}, // U-ish tall
      /*h21*/   {{0,0},{1,0},{0,1},{0,2},{0,3},{0,4}}, // Y-ish
      /*h22*/   {{0,0},{0,1},{0,2},{0,3},{0,4},{1,2}},
      /*h23*/   {{0,0},{1,0},{1,1},{1,2},{2,2},{2,3}}, // W-ish
      /*h24*/   {{0,0},{1,0},{1,1},{2,1},{2,2},{3,2}}, // staircase
      /*h25*/   {{0,0},{1,0},{1,1},{1,2},{1,3},{2,3}}, // N/Z extended
      /*h26*/   {{0,0},{0,1},{1,1},{1,2},{2,2},{2,3}},
      /*h27*/   {{1,0},{1,1},{0,1},{0,2},{1,2},{1,3}}, // F-ish
      /*h28*/   {{1,0},{0,1},{1,1},{2,1},{1,2},{1,3}}, // plus-tail
      /*h29*/   {{1,0},{0,1},{1,1},{2,1},{0,2},{2,2}},
      /*h30*/   {{0,0},{2,0},{0,1},{1,1},{2,1},{1,2}},
      /*h31*/   {{1,0},{0,1},{1,1},{2,1},{1,2},{0,2}},
      /*h32*/   {{0,0},{1,0},{2,0},{1,1},{1,2},{2,2}},
      /*h33*/   {{0,0},{1,0},{2,0},{0,1},{1,1},{1,2}},
      /*h34*/   {{0,0},{1,0},{1,1},{2,1},{0,1},{1,2}},
      /*h35*/   {{0,1},{1,1},{2,1},{1,0},{1,2},{0,2}},
    };
    int idx = (fShape[1]-'0')*10 + (fShape[2]-'0') - 1;  // h01->0 .. h35->34
    if (idx >= 0 && idx < (int)HEX.size()) cells = HEX[idx];
    else { fShape = "square"; cells = {{0,1},{1,1},{0,0},{1,0}}; }

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
  // Optional z-stagger: crystal i is lifted by i*fZStagger. Combined with the
  // outward arc shields, this makes inter-crystal shadowing front/back
  // ASYMMETRIC (a back-side gamma must cross a neighbour at a different height),
  // which is the only thing that can break the planar mirror degeneracy.
  std::vector<G4ThreeVector> pos;
  G4double zMid = 0.5 * (fNumPixels - 1) * fZStagger;  // centre the stack in z
  for (G4int i = 0; i < fNumPixels; ++i) {
    G4ThreeVector p = cellPos(cells[i]);
    p.setZ(i * fZStagger - zMid);
    pos.push_back(p);

    auto logicCasing =
      new G4LogicalVolume(solidCasing, casing_mat, "Casing" + std::to_string(i));
    new G4PVPlacement(nullptr, p, logicCasing,
                      "Casing" + std::to_string(i), logicEnv, false, i, checkOverlaps);

    fLogicPixel[i] = new G4LogicalVolume(solidPixel, pixel_mat, "Pixel" + std::to_string(i));
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), fLogicPixel[i],
                      "Pixel" + std::to_string(i), logicCasing, false, i, checkOverlaps);

    // Place the directional shield so its SENSITIVE direction (open half / hole,
    // built at +y) points OUTWARD -- away from the array centre. Each crystal is
    // thus most sensitive to its own outward sector; the solid-lead opposite
    // side blocks the back, breaking the front/back degeneracy.
    if (logicAsym) {
      // Outward direction for this crystal (centre->crystal). For a single
      // crystal at the origin, default outward = +y.
      G4double outwardAng = (p.mag() > 1e-6) ? std::atan2(p.y(), p.x()) : 90. * deg;
      // Shield's sensitive dir is built at +y = 90 deg; rotate to outwardAng.
      G4double rotZ = outwardAng - 90. * deg;
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
