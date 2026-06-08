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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the B1::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "Randomize.hh"

#include <cmath>
#include <vector>

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName = "gamma");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(0.5 * MeV);  // paper's simulation energy

  // Default source placement: 30 deg, 50 cm away, coplanar with detector.
  fAngle = 30. * deg;
  fDistance = 50. * cm;

  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::DefineCommands()
{
  // UI commands under /det/ to move the source without recompiling.
  //   /det/sourceAngle    30 deg
  //   /det/sourceDistance 50 cm
  fMessenger =
    new G4GenericMessenger(this, "/det/", "Gamma source placement controls");

  auto& angleCmd = fMessenger->DeclareMethodWithUnit(
    "sourceAngle", "deg", &PrimaryGeneratorAction::SetSourceAngle,
    "Source angle, clockwise from detector front (+Y).");
  angleCmd.SetParameterName("angle", false);
  angleCmd.SetDefaultValue("30. deg");

  auto& distCmd = fMessenger->DeclareMethodWithUnit(
    "sourceDistance", "cm", &PrimaryGeneratorAction::SetSourceDistance,
    "Source distance from the detector center.");
  distCmd.SetParameterName("distance", false);
  distCmd.SetDefaultValue("50. cm");

  fMessenger->DeclareMethod(
    "randomizeSource", &PrimaryGeneratorAction::RandomizeSource,
    "Pick a new random source: distance in [20,500] cm, angle in [0,360) deg. "
    "Call once before each beamOn to generate one training configuration.");

  auto& typeCmd = fMessenger->DeclareMethod(
    "sourceType", &PrimaryGeneratorAction::SetSourceType,
    "Source energy spectrum: 'mono' (single /gun/energy line) or "
    "'euba' (Eu-152 + Ba-133 mixed calibration source, discrete lines).");
  typeCmd.SetParameterName("type", true);
  typeCmd.SetDefaultValue("mono");

  auto& bgCmd = fMessenger->DeclareMethod(
    "bgFraction", &PrimaryGeneratorAction::SetBgFraction,
    "Fraction (0..1) of events drawn as 138La internal background, "
    "interleaved with the external source to make a realistic mixed config.");
  bgCmd.SetParameterName("frac", true);
  bgCmd.SetDefaultValue("0.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::RandomizeSource()
{
  // Matches the paper's training distribution (gen_data_square.py):
  //   distance in [20, 500] cm, angle in [0, 360) deg.
  G4double dmin = 20. * cm, dmax = 500. * cm;
  fDistance = dmin + (dmax - dmin) * G4UniformRand();
  fAngle = 360. * deg * G4UniformRand();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  // --- Internal LaBr3 radioactivity (138La) ---
  // Instead of an external source, the primary is generated at a random point
  // INSIDE a randomly chosen crystal and emitted isotropically. 138La decays:
  //   ~66% electron capture -> 1436 keV gamma (+ ~32 keV Ba X-ray, omitted)
  //   ~34% beta-minus       ->  789 keV gamma
  // This is the dominant intrinsic gamma background of LaBr3:Ce. It is
  // isotropic and originates within the crystals, so it adds a near-uniform
  // baseline across pixels -- distinguishable from a directional source both
  // by energy (789/1436 keV vs the low-energy Eu/Ba lines) and by symmetry.
  // Pure internal mode, OR a mixed config where a fraction fBgFraction of the
  // events are interleaved 138La internal decays alongside the external source.
  if (fSourceType == "internal" ||
      (fBgFraction > 0. && G4UniformRand() < fBgFraction)) {
    GenerateInternalDecay(event);
    return;
  }

  // Called at the beginning of each event (external source).
  //
  // Paper geometry: the source and detector are coplanar (X-Y plane). The
  // source direction is an angle theta measured clockwise from the detector
  // "front" (+Y axis). We place a point source at distance fDistance in that
  // direction and aim each emitted gamma straight back at the detector center
  // (origin). This is the directional configuration the detector is trained on.
  //
  //   x = d * sin(theta)
  //   y = d * cos(theta)   (theta = 0 => directly in front, +Y)
  G4double x0 = fDistance * std::sin(fAngle);
  G4double y0 = fDistance * std::cos(fAngle);
  G4double z0 = 0.;
  G4ThreeVector sourcePos(x0, y0, z0);

  fParticleGun->SetParticlePosition(sourcePos);

  // Aim at the detector center (origin).
  G4ThreeVector dir = (G4ThreeVector(0., 0., 0.) - sourcePos).unit();
  fParticleGun->SetParticleMomentumDirection(dir);

  // Pick this primary's energy according to the source type.
  fParticleGun->SetParticleEnergy(SampleEnergy());

  fParticleGun->GeneratePrimaryVertex(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GenerateInternalDecay(G4Event* event)
{
  // Locate the four crystal placements and their cylinder dimensions from the
  // geometry store (works for any shape/spacing configuration).
  auto pvStore = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4ThreeVector> centers;
  G4double rad = 0., halfz = 0.;
  for (G4int i = 0; i < DetectorConstruction::GetNumPixels(); ++i) {
    G4String name = "Pixel" + std::to_string(i);
    G4VPhysicalVolume* pv = pvStore->GetVolume(name, false);
    if (!pv) continue;
    // Global position = casing translation (crystal sits at casing-local origin).
    G4VPhysicalVolume* casing = pvStore->GetVolume("Casing" + std::to_string(i), false);
    G4ThreeVector c = casing ? casing->GetTranslation() : pv->GetTranslation();
    centers.push_back(c);
    auto tubs = dynamic_cast<G4Tubs*>(pv->GetLogicalVolume()->GetSolid());
    if (tubs) { rad = tubs->GetOuterRadius(); halfz = tubs->GetZHalfLength(); }
  }
  if (centers.empty()) return;

  // Pick a random crystal and a uniform random point inside its cylinder.
  G4int k = std::min((G4int)(G4UniformRand() * centers.size()),
                     (G4int)centers.size() - 1);
  G4double rr = rad * std::sqrt(G4UniformRand());      // uniform in disk
  G4double phi = 2. * CLHEP::pi * G4UniformRand();
  G4double dz = halfz * (2. * G4UniformRand() - 1.);
  G4ThreeVector pos = centers[k] +
    G4ThreeVector(rr * std::cos(phi), rr * std::sin(phi), dz);
  fParticleGun->SetParticlePosition(pos);

  // Isotropic emission direction.
  G4double ct = 2. * G4UniformRand() - 1.;
  G4double st = std::sqrt(1. - ct * ct);
  G4double ph = 2. * CLHEP::pi * G4UniformRand();
  fParticleGun->SetParticleMomentumDirection(
    G4ThreeVector(st * std::cos(ph), st * std::sin(ph), ct));

  // 138La gamma line: 66% -> 1436 keV (EC branch), 34% -> 789 keV (beta branch).
  G4double e = (G4UniformRand() < 0.66) ? 1436.0 * keV : 789.0 * keV;
  fParticleGun->SetParticleEnergy(e);

  fParticleGun->GeneratePrimaryVertex(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PrimaryGeneratorAction::SampleEnergy()
{
  if (fSourceType != "euba") {
    // "mono": keep whatever energy the gun is set to (default 0.5 MeV or
    // whatever /gun/energy specified).
    return fParticleGun->GetParticleEnergy();
  }

  // Eu-152 + Ba-133 mixed source: discrete gamma lines, sampled in proportion
  // to their emission intensities (gammas per 100 decays). Values are the main
  // lines; minor lines are omitted for clarity.
  //   {energy in keV, relative intensity}
  static const G4int N = 12;
  static const G4double eKeV[N] = {
    // Ba-133
    81.0, 276.4, 302.9, 356.0, 383.8,
    // Eu-152
    121.8, 244.7, 344.3, 778.9, 964.1, 1112.1, 1408.0};
  static const G4double inten[N] = {
    34.1, 7.16, 18.3, 62.0, 8.94,        // Ba-133
    28.4, 7.55, 26.6, 12.9, 14.6, 13.6, 21.0};  // Eu-152

  // Cumulative-sum sampling.
  static G4double cum[N];
  static G4bool init = false;
  static G4double tot = 0.;
  if (!init) {
    for (G4int i = 0; i < N; ++i) { tot += inten[i]; cum[i] = tot; }
    init = true;
  }
  G4double r = G4UniformRand() * tot;
  for (G4int i = 0; i < N; ++i) {
    if (r <= cum[i]) return eKeV[i] * keV;
  }
  return eKeV[N - 1] * keV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
