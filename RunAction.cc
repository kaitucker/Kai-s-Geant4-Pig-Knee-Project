#include "RunAction.hh"
#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include <cmath>
#include "G4ios.hh"

RunAction::RunAction() : G4UserRunAction(), fTotalEdep(0.) {}
RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run*) {
  fTotalEdep = 0.;
  G4cout << "### Run start" << G4endl;
}

void RunAction::EndOfRunAction(const G4Run*) {
  G4cout << "### Run end" << G4endl;
  G4cout << "Total energy deposited (internal units): " << fTotalEdep << G4endl;
  G4cout << "Total energy deposited: " << G4BestUnit(fTotalEdep, "Energy") << G4endl;

  // Approximate pig knee as ellipsoid with same semi-axes as in DetectorConstruction:
  G4double a = 45*mm, b = 30*mm, c = 60*mm;
  G4double volume = (4.0/3.0) * M_PI * a * b * c; // volume in CLHEP length^3 units
  // density of water (material used for pig knee)
  auto nist = G4NistManager::Instance();
  G4Material* water = nist->FindOrBuildMaterial("G4_WATER");
  G4double density = water->GetDensity(); // mass/volume in internal units
  G4double mass = density * volume; // mass in internal mass units
  G4double mass_kg = mass / kg;

  G4cout << "Pig knee approx volume: " << G4BestUnit(volume, "Volume") << G4endl;
  G4cout << "Pig knee material density: " << G4BestUnit(density, "Volumic Mass") << G4endl;
  G4cout << "Pig knee mass (kg): " << mass_kg << G4endl;

  if (mass_kg > 0) {
    G4double dose_Gy = (fTotalEdep / joule) / mass_kg; // J/kg = Gy
    G4cout << "Approx dose (Gy) in pig knee: " << dose_Gy << G4endl;
  } else {
    G4cout << "Cannot compute mass (mass_kg==0)" << G4endl;
  }
}
