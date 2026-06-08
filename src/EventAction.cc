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
/// \file EventAction.cc
/// \brief Implementation of the B1::EventAction class

#include "EventAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction) : fRunAction(runAction)
{
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SetRoiLow(G4double lo)  { fRoiLow = lo / keV; fRoiEnabled = true; }
void EventAction::SetRoiHigh(G4double hi) { fRoiHigh = hi / keV; fRoiEnabled = true; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::DefineCommands()
{
  fMessenger = new G4GenericMessenger(this, "/det/", "Analysis controls");
  auto& loCmd = fMessenger->DeclareMethodWithUnit(
    "roiLow", "keV", &EventAction::SetRoiLow,
    "Lower edge of the energy ROI applied to each pixel's per-event deposit.");
  loCmd.SetParameterName("lo", true);
  loCmd.SetDefaultValue("40 keV");
  auto& hiCmd = fMessenger->DeclareMethodWithUnit(
    "roiHigh", "keV", &EventAction::SetRoiHigh,
    "Upper edge of the energy ROI applied to each pixel's per-event deposit.");
  hiCmd.SetParameterName("hi", true);
  hiCmd.SetDefaultValue("450 keV");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  for (G4int i = 0; i < kMaxPixels; ++i) fEdep[i] = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{
  auto analysisManager = G4AnalysisManager::Instance();
  const G4int n = DetectorConstruction::GetNumPixels();

  // Fill the FULL per-crystal energy spectrum (no ROI cut) so the histograms
  // always show the complete spectrum for inspection.
  for (G4int i = 0; i < n; ++i) {
    if (fEdep[i] > 0.) analysisManager->FillH1(i, fEdep[i] / keV);
  }

  // Build the aggregated readout vector. If an ROI is set, a pixel only
  // contributes when its deposit falls in the window -- this is the photopeak
  // energy cut, applied here exactly as it would be to real measured data.
  G4double gated[kMaxPixels] = {0.};
  G4double total = 0.;
  for (G4int i = 0; i < n; ++i) {
    G4double e = fEdep[i] / keV;
    G4bool inRoi = !fRoiEnabled || (e >= fRoiLow && e <= fRoiHigh);
    gated[i] = (e > 0. && inRoi) ? fEdep[i] : 0.;
    total += gated[i];
  }
  fRunAction->AddPixelVector(gated);
  fRunAction->AddEdep(total);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
