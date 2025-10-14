#include "PhysicsList.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
    SetVerboseLevel(1);
    
    // Use Livermore physics for better low-energy photon interactions
    RegisterPhysics(new G4EmLivermorePhysics());
    RegisterPhysics(new G4DecayPhysics());
    RegisterPhysics(new G4RadioactiveDecayPhysics());
    
    // Enable atomic deexcitation for accurate fluorescence
    G4LossTableManager::Instance()->SetAtomDeexcitation(new G4UAtomicDeexcitation());
}

PhysicsList::~PhysicsList() {}

void PhysicsList::SetCuts()
{
    G4VModularPhysicsList::SetCuts();
    
    G4cout << "PhysicsList::SetCuts: Applying special cuts to tendon region." << G4endl;
    
    // Find tendon logical volume
    auto tendonLogical = G4LogicalVolumeStore::GetInstance()->GetVolume("PatellarTendonLog", false);
    
    if (tendonLogical) {
        G4cout << "Found PatellarTendonLog. Setting up special region with tight cuts." << G4endl;
        
        // Create region for tendon
        auto tendonRegion = new G4Region("TendonRegion");
        tendonRegion->AddRootLogicalVolume(tendonLogical);
        
        // Set very fine production cuts for accurate energy deposition tracking
        auto tendonCuts = new G4ProductionCuts();
        tendonCuts->SetProductionCut(10*um);  // 10 micrometer cut
        
        tendonRegion->SetProductionCuts(tendonCuts);
        tendonRegion->DumpProductionCuts();
        
    } else {
        G4cout << "WARNING: PatellarTendonLog not found! Detector scoring may not work." << G4endl;
    }
}
