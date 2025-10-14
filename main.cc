#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
int main(int argc, char** argv) {
// Make run manager
auto* runManager = G4RunManagerFactory::CreateRunManager();
// Set user initialization classes
runManager->SetUserInitialization(new DetectorConstruction());
runManager->SetUserInitialization(new PhysicsList());
runManager->SetUserInitialization(new ActionInitialization());
// Initialize Geant4 kernel
runManager->Initialize();
// Visualization
G4VisManager* visManager = new G4VisExecutive;
visManager->Initialize();
G4UImanager* UImanager = G4UImanager::GetUIpointer();
if (argc == 1) {
// interactive
G4UIExecutive* ui = new G4UIExecutive(argc, argv);
UImanager->ApplyCommand("/control/execute init_vis.mac");
ui->SessionStart();
delete ui;
} else {
// batch macro provided as argument
G4String command = "/control/execute ";
G4String fileName = argv[1];
UImanager->ApplyCommand(command + fileName);
}
delete visManager;
delete runManager;
return 0;
}
