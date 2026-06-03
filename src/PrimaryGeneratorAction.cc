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

#include "G4Box.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <cmath>

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
  // Called at the beginning of each event.
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

  fParticleGun->GeneratePrimaryVertex(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
