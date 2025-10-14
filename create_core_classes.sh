cat > include/DetectorConstruction.hh <<EOL
#ifndef DETECTORCONSTRUCTION_H
#define DETECTORCONSTRUCTION_H

#include "G4VUserDetectorConstruction.hh"

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();
    virtual G4VPhysicalVolume* Construct() override;
};

#endif
EOL

cat > src/DetectorConstruction.cc <<EOL
#include "DetectorConstruction.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"

DetectorConstruction::DetectorConstruction() {}
DetectorConstruction::~DetectorConstruction() {}

G4VPhysicalVolume* DetectorConstruction::Construct() {
    G4NistManager* nist = G4NistManager::Instance();

    // World
    G4double world_size = 1.0*m;
    G4Box* solidWorld = new G4Box("World", world_size/2, world_size/2, world_size/2);
    G4Material* worldMat = nist->FindOrBuildMaterial("G4_AIR");
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
    G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0);

    // Placeholder voxel box
    G4Box* solidVoxel = new G4Box("Voxel", 1.0*cm, 1.0*cm, 1.0*cm);
    G4Material* voxelMat = nist->FindOrBuildMaterial("G4_WATER");
    new G4PVPlacement(0, G4ThreeVector(0,0,0), new G4LogicalVolume(solidVoxel, voxelMat, "Voxel"), "Voxel", logicWorld, false, 0);

    return physWorld;
}
EOL

# PrimaryGeneratorAction
cat > include/PrimaryGeneratorAction.hh <<EOL
#ifndef PRIMARYGENERATORACTION_H
#define PRIMARYGENERATORACTION_H

#include "G4VUserPrimaryGeneratorAction.hh"

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();
    virtual void GeneratePrimaries(G4Event* event) override;
};

#endif
EOL

cat > src/PrimaryGeneratorAction.cc <<EOL
#include "PrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Event.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction() {}
PrimaryGeneratorAction::~PrimaryGeneratorAction() {}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    G4ParticleGun* particleGun = new G4ParticleGun(1);
    particleGun->SetParticleDefinition(G4Gamma::Definition());
    particleGun->SetParticleEnergy(100*keV);
    particleGun->SetParticlePosition(G4ThreeVector(0,0,0));
    particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    particleGun->GeneratePrimaryVertex(event);
}
EOL

# PhysicsList
cat > include/PhysicsList.hh <<EOL
#ifndef PHYSICSLIST_H
#define PHYSICSLIST_H

#include "G4VModularPhysicsList.hh"

class PhysicsList : public G4VModularPhysicsList {
public:
    PhysicsList();
    virtual ~PhysicsList();
};

#endif
EOL

cat > src/PhysicsList.cc <<EOL
#include "PhysicsList.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"

PhysicsList::PhysicsList() {
    RegisterPhysics(new G4EmLivermorePhysics());
    G4LossTableManager::Instance()->SetAtomDeexcitation(new G4UAtomicDeexcitation());
}

PhysicsList::~PhysicsList() {}
EOL

# RunAction
cat > include/RunAction.hh <<EOL
#ifndef RUNACTION_H
#define RUNACTION_H

#include "G4UserRunAction.hh"

class RunAction : public G4UserRunAction {
public:
    RunAction();
    virtual ~RunAction();
    virtual void EndOfRunAction(const G4Run*) override;
};

#endif
EOL

cat > src/RunAction.cc <<EOL
#include "RunAction.hh"
#include <fstream>
#include <iostream>

RunAction::RunAction() {}
RunAction::~RunAction() {}

void RunAction::EndOfRunAction(const G4Run*) {
    // Placeholder: write dummy CSV
    std::ofstream file("dose_output.csv");
    file << "voxel_id,energy_deposit_J\n";
    file << "0,0.0\n";
    file.close();
    std::cout << "Run completed. Dummy dose CSV written.\n";
}
EOL

echo "Minimal Geant4 core classes created in include/ and src/"
