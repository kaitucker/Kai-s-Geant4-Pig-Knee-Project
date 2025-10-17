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
  // The names are defined in DetectorConstruction.cc ("TendonMFD" and "Edep").
  if (fEdepHCID < 0) {
    auto sdManager = G4SDManager::GetSDMpointer();
    fEdepHCID = sdManager->GetCollectionID("TendonMFD/Edep");
  }
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // The 'GetHCofThisEvent()' method was incorrect.
  // We get the merged hits collection directly from the G4Run object using its ID.
  auto hitsMap = static_cast<G4THitsMap<G4double>*>(run->GetHC(fEdepHCID));

  G4double totalEdep = 0.;
  if (hitsMap && !hitsMap->GetMap()->empty()) {
      // The map contains the total energy deposited.
      // For a single scoring volume, it's the first (and only) entry.
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

  // Print the final, merged results for the entire run
  if (IsMaster()) {
    G4cout << G4endl;
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
