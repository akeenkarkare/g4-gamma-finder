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
/// \file RunAction.hh
/// \brief Definition of the B1::RunAction class

#ifndef B1RunAction_h
#define B1RunAction_h 1

#include "G4UserRunAction.hh"

#include "G4Accumulable.hh"
#include "globals.hh"

class G4Run;
class G4GenericMessenger;

namespace B1
{

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    ~RunAction() override;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

    void AddEdep(G4double edep);

    // Per-config aggregation: EventAction pushes each event's 4-pixel vector
    // here; EndOfRunAction writes ONE aggregated row (the config readout) so a
    // dataset of N configs is N rows, not N x events rows.
    static constexpr G4int kMaxPixels = 6;
    void AddPixelVector(const G4double e[kMaxPixels]);

    // Dataset-file control (driven by /det/ UI commands). Open the ntuple
    // file once, run many beamOn's (configs), then write once.
    void OpenDataFile(G4String name);
    void WriteDataFile();

  private:
    void DefineCommands();

    G4Accumulable<G4double> fEdep = 0.;
    G4Accumulable<G4double> fEdep2 = 0.;
    G4GenericMessenger* fMessenger = nullptr;
    G4String fFileName = "pixels";
    G4bool fFileOpen = false;
    G4bool fAutoFile = false;  // true if opened automatically (single-run macro)

    G4double fPixelSum[kMaxPixels] = {0.};  // per-run pixel totals
};

}  // namespace B1

#endif
