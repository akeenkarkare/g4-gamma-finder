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
/// \brief Definition of the B1::DetectorConstruction class

#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4GenericMessenger;

namespace B1
{

/// Detector construction class to define materials and geometry.
///
/// The detector is a 4-crystal array whose ARRANGEMENT (shape) and lead
/// padding thickness are configurable at run time via /det/ UI commands:
///   /det/shape   square | S | J | T | L
///   /det/padding <thickness> mm
/// Re-issue /run/reinitializeGeometry after changing either.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

    // Maximum supported crystals (sizes fixed arrays). The ACTIVE count for the
    // current shape is fNumPixels (4, 5, or 6), available via GetNumPixels().
    static constexpr G4int kMaxPixels = 6;
    static G4int GetNumPixels();              // active count for the current shape

    // Configuration setters (used by the messenger).
    void SetShape(G4String shape);
    void SetPadding(G4double pad) { fPadding = pad; }
    void SetAsymCasing(G4double t) { fAsymCasingThk = t; }  // extra shield thickness on one side; 0 = uniform
    void SetAsymMaterial(G4String m) { fAsymMaterial = m; } // "Al" | "Pb" | "W"
    void SetAsymMode(G4String m) { fAsymMode = m; }         // "arc" (half-shield) | "window" (near-full ring w/ hole)
    void SetZStagger(G4double z) { fZStagger = z; }         // per-crystal z step (breaks planar front/back symmetry)
    G4String GetShape() const { return fShape; }
    G4double GetPadding() const { return fPadding; }

  protected:
    void DefineCommands();

    G4LogicalVolume* fScoringVolume = nullptr;
    G4LogicalVolume* fLogicPixel[kMaxPixels] = {nullptr};
    G4GenericMessenger* fMessenger = nullptr;

    G4String fShape = "square";       // square|S|J|T|L (4) | P|F (5) | H|U (6)
    G4double fPadding = 1.0;          // lead padding FULL thickness in mm (default 1 mm)
    G4double fAsymCasingThk = 0.;     // extra shield on the -y half of each casing; 0 = uniform
    G4String fAsymMaterial = "Al";    // material of the asymmetric shield arc
    G4String fAsymMode = "arc";       // "arc" = open half + shielded half; "window" = near-full ring with a hole
    G4double fZStagger = 0.;          // per-crystal z step (i * fZStagger); 0 = coplanar

    static G4int fNumPixels;          // active crystal count for the built shape
};

}  // namespace B1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
