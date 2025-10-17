#include "RunAction.hh"
#include "DetectorConstruction.hh" // <-- ADD

#include "G4Run.hh"
#include "G4RunManager.hh" // <-- ADD
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4AccumulableManager.hh" // <-- ADD

#include <iomanip> // <-- ADD for formatting output

RunAction::RunAction()
 : G4UserRunAction(),
   fTotalEdep(0.),
   fTotalEdep2(0.)
{
  // Register the accumulables with the manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Register(fTotalEdep);
  accumulableManager->Register(fTotalEdep2);
}

RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction(const G4Run*)
{
  // Reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

  G4cout << "### Run start" << G4endl;
}

void RunAction::AddEdep(G4double edep)
{
  fTotalEdep += edep;
  fTotalEdep2 += edep*edep;
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables from all threads
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Get the final values
  G4double edep  = fTotalEdep.GetValue();
  G4double edep2 = fTotalEdep2.GetValue();

  // Calculate the standard deviation
  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;

  // Get the mass of the scoring volume from DetectorConstruction
  const auto detConstruction = static_cast<const DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detConstruction->GetScoringVolume()->GetMass();

  // Calculate the dose and its uncertainty
  G4double dose = 0.;
  G4double rmsDose = 0.;
  if (mass > 0.) {
    dose = edep / mass;
    rmsDose = rms / mass;
  }

  // Print the results
  G4cout << "### Run end" << G4endl;
  G4cout << "--------------------End of Run-----------------------" << G4endl;
  G4cout << " The run consisted of " << nofEvents << " events." << G4endl;
  G4cout << " The mass of the scoring volume (PatellarTendonLog) is: "
         << G4BestUnit(mass, "Mass") << G4endl;
  G4cout << " Total energy deposited in scoring volume: "
         << G4BestUnit(edep, "Energy") << " +/- " << G4BestUnit(rms, "Energy") << G4endl;
  G4cout << " Mean dose to scoring volume: "
         << G4BestUnit(dose, "Dose") << " +/- " << G4BestUnit(rmsDose, "Dose") << G4endl;
  G4cout << "------------------------------------------------------------" << G4endl << G4endl;
}
