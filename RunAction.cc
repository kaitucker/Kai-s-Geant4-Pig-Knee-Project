#include "RunAction.hh"
#include "DetectorConstruction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

RunAction::RunAction()
 : G4UserRunAction(),
   fEdepHCID(-1)
{}

RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction(const G4Run* run)
{
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;

  // Get the Hits Collection ID for our energy deposit primitive scorer.
  // The names are defined in DetectorConstruction.cc ("TendonMFD"/"Edep")
  if (fEdepHCID < 0) {
    auto sdManager = G4SDManager::GetSDMpointer();
    fEdepHCID = sdManager->GetCollectionID("TendonMFD/Edep");
  }
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Get the hits collection for this run
  auto hce = run->GetHCofThisEvent();
  if (!hce) return;

  // Get the hits map for our energy deposit scorer
  auto hitsMap = static_cast<G4THitsMap<G4double>*>(hce->GetHC(fEdepHCID));

  G4double totalEdep = 0.;
  if (hitsMap && !hitsMap->GetMap()->empty()) {
      // The map contains the total energy deposited in the volume.
      // For a single volume, it's the first (and only) entry.
      totalEdep = *(hitsMap->GetMap()->begin()->second);
  }

  // Get the mass of the scoring volume from DetectorConstruction
  const auto detConstruction = static_cast<const DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detConstruction->GetScoringVolume()->GetMass();

  // Calculate the dose
  G4double dose = 0.;
  if (mass > 0.) {
    dose = totalEdep / mass;
  }

  // Print the final, merged results
  if (IsMaster()) {
    G4cout << "--------------------End of Global Run-----------------------" << G4endl;
    G4cout << " The run consisted of " << nofEvents << " events." << G4endl;
    G4cout << " The mass of the scoring volume (PatellarTendonLog) is: "
           << G4BestUnit(mass, "Mass") << G4endl;
    G4cout << " Total energy deposited in scoring volume: "
           << G4BestUnit(totalEdep, "Energy") << G4endl;
    G4cout << " Mean dose to scoring volume: "
           << G4BestUnit(dose, "Dose") << G4endl;
    G4cout << "------------------------------------------------------------" << G4endl << G4endl;
  }
}
