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
/// \file EventAction.hh
/// \brief Definition of the B1::EventAction class

#ifndef B1EventAction_h
#define B1EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;
class G4GenericMessenger;

namespace B1
{

class RunAction;

/// Event action class

class EventAction : public G4UserEventAction
{
  public:
    static constexpr G4int kNumPixels = 4;

    EventAction(RunAction* runAction);
    ~EventAction() override;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    // Add energy deposit to a specific pixel (0..3).
    void AddEdep(G4int pixel, G4double edep)
    {
      if (pixel >= 0 && pixel < kNumPixels) fEdep[pixel] += edep;
    }

    // Energy region-of-interest. When fRoiEnabled, a pixel's deposit
    // contributes to the aggregated readout only if it falls in [lo, hi] --
    // mimicking the photopeak energy cut applied to real measured data.
    // The messenger passes values in Geant4 energy units; we store keV.
    void SetRoiLow(G4double lo);
    void SetRoiHigh(G4double hi);

  private:
    void DefineCommands();

    RunAction* fRunAction = nullptr;
    G4double fEdep[kNumPixels] = {0., 0., 0., 0.};

    G4GenericMessenger* fMessenger = nullptr;
    G4bool fRoiEnabled = false;
    G4double fRoiLow = 0.;     // keV
    G4double fRoiHigh = 1e9;   // keV
};

}  // namespace B1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
