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
/// \file RunAction.cc
/// \brief Implementation of the B1::RunAction class

#include "RunAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4AccumulableManager.hh"
#include "G4AnalysisManager.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include <string>

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
  // add new units for dose
  //
  const G4double milligray = 1.e-3 * gray;
  const G4double microgray = 1.e-6 * gray;
  const G4double nanogray = 1.e-9 * gray;
  const G4double picogray = 1.e-12 * gray;

  new G4UnitDefinition("milligray", "milliGy", "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy", "Dose", microgray);
  new G4UnitDefinition("nanogray", "nanoGy", "Dose", nanogray);
  new G4UnitDefinition("picogray", "picoGy", "Dose", picogray);

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Register(fEdep);
  accumulableManager->Register(fEdep2);

  // --- Set up the analysis manager / ntuple for the aggregated readout ---
  // One row per CONFIG (per beamOn): the summed per-pixel energy plus the
  // source angle/distance label and the event count. This is the ML training
  // sample format -- N configs -> N rows.
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("csv");
  analysisManager->SetNtupleMerging(true);  // merge per-thread output into one file
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetH1Activation(true);   // enable histogram filling

  // Column layout (fixed so old/new files are parseable by name):
  //   e0..e5            : per-crystal total ROI energy (cols 0..5)
  //   c0_b0,c0_b1,...   : per-crystal per-band COUNTS (cols 6 .. 6+2*kMaxPixels-1)
  //   angle_deg, dist_cm, nEvents
  // Shapes with fewer crystals leave the unused columns at 0. The band columns
  // give the ML the spectral shape (peak vs Compton) per crystal, not just total.
  analysisManager->CreateNtuple("pixels", "Aggregated per-config readout");
  for (G4int i = 0; i < kMaxPixels; ++i)
    analysisManager->CreateNtupleDColumn("e" + std::to_string(i) + "_keV");
  for (G4int i = 0; i < kMaxPixels; ++i)
    for (G4int b = 0; b < kBands; ++b)
      analysisManager->CreateNtupleDColumn("c" + std::to_string(i) + "_b" + std::to_string(b));
  analysisManager->CreateNtupleDColumn("angle_deg");
  analysisManager->CreateNtupleDColumn("dist_cm");
  analysisManager->CreateNtupleIColumn("nEvents");
  analysisManager->FinishNtuple();

  // --- Per-crystal energy spectra (histograms) ---
  // One H1 per crystal (up to kMaxPixels): counts vs deposited energy. Range
  // 0-1600 keV captures the low-energy source photopeaks AND the high-energy
  // 138La intrinsic lines (789, 1436 keV). Unused crystals stay empty.
  for (G4int i = 0; i < kMaxPixels; ++i) {
    analysisManager->CreateH1("spectrum" + std::to_string(i),
                              "Energy spectrum pixel " + std::to_string(i) + " (keV)",
                              1600, 0., 1600.);
  }

  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::DefineCommands()
{
  fMessenger =
    new G4GenericMessenger(this, "/det/", "Dataset output controls");

  auto& openCmd = fMessenger->DeclareMethod(
    "openFile", &RunAction::OpenDataFile,
    "Open a dataset file (base name, no extension). Keeps it open across "
    "multiple beamOn runs so a whole angle scan goes into one CSV.");
  openCmd.SetParameterName("name", true);
  openCmd.SetDefaultValue("pixels");

  fMessenger->DeclareMethod(
    "writeFile", &RunAction::WriteDataFile,
    "Write and close the dataset file. Call once at the end of a scan.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::OpenDataFile(G4String name)
{
  fFileName = name;
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->OpenFile(fFileName);
  fFileOpen = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::WriteDataFile()
{
  if (!fFileOpen) return;
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
  fFileOpen = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddPixelVector(const G4double e[kMaxPixels])
{
  for (G4int i = 0; i < kMaxPixels; ++i) fPixelSum[i] += e[i];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddBandCount(G4int pixel, G4int band)
{
  if (pixel >= 0 && pixel < kMaxPixels && band >= 0 && band < kBands)
    fBandSum[pixel][band] += 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

  // Reset this config's per-pixel totals and per-band counts.
  for (G4int i = 0; i < kMaxPixels; ++i) {
    fPixelSum[i] = 0.;
    for (G4int b = 0; b < kBands; ++b) fBandSum[i][b] = 0.;
  }

  // If the user did not explicitly open a dataset file (e.g. a simple
  // single-run macro), open one automatically so output is never lost.
  if (!fFileOpen) {
    OpenDataFile(fFileName);
    fAutoFile = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();

  G4double rms = edep2 - edep * edep / nofEvents;
  if (rms > 0.)
    rms = std::sqrt(rms);
  else
    rms = 0.;

  const auto detConstruction = static_cast<const DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detConstruction->GetScoringVolume()->GetMass();
  G4double dose = edep / mass;
  G4double rmsDose = rms / mass;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const auto generatorAction = static_cast<const PrimaryGeneratorAction*>(
    G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction) {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy, "Energy");
  }

  // Print
  //
  if (IsMaster()) {
    G4cout << G4endl << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout << G4endl << "--------------------End of Local Run------------------------";
  }

  G4cout << G4endl << " The run is " << nofEvents << " " << runCondition << G4endl << G4endl;
  G4cout << "  --> cumulated edep per run in scoring volume = " << G4BestUnit(edep, "Energy") 
         << " = " << edep/joule << " joule" << G4endl;  
  G4cout << "  --> mass of scoring volume = " << G4BestUnit(mass, "Mass") << G4endl << G4endl; 
  G4cout << " Absorbed dose per run in scoring volume = edep/mass = " << G4BestUnit(dose, "Dose")
         << "; rms = " << G4BestUnit(rmsDose, "Dose") << G4endl
         << "------------------------------------------------------------" << G4endl << G4endl;

  // --- Write ONE aggregated ntuple row for this config (this run) ---
  if (fFileOpen) {
    G4double angleDeg = 0., distCm = 0.;
    if (generatorAction) {
      angleDeg = generatorAction->GetSourceAngle() / deg;
      distCm = generatorAction->GetSourceDistance() / cm;
    }
    auto analysisManager = G4AnalysisManager::Instance();
    G4int col = 0;
    for (G4int i = 0; i < kMaxPixels; ++i)
      analysisManager->FillNtupleDColumn(col++, fPixelSum[i] / keV);   // energy cols
    for (G4int i = 0; i < kMaxPixels; ++i)
      for (G4int b = 0; b < kBands; ++b)
        analysisManager->FillNtupleDColumn(col++, fBandSum[i][b]);     // band-count cols
    analysisManager->FillNtupleDColumn(col++, angleDeg);
    analysisManager->FillNtupleDColumn(col++, distCm);
    analysisManager->FillNtupleIColumn(col++, nofEvents);
    analysisManager->AddNtupleRow();
  }

  // If the file was opened automatically for a single-run macro, write and
  // close it here. For an explicit scan (/det/openFile ... /det/writeFile),
  // leave it open so subsequent beamOn runs keep appending.
  if (fAutoFile) {
    WriteDataFile();
    fAutoFile = false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddEdep(G4double edep)
{
  fEdep += edep;
  fEdep2 += edep * edep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
