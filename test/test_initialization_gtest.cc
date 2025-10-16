#include <gtest/gtest.h>
#include "G4RunManagerFactory.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

TEST(Geant4Setup, InitializesSuccessfully) {
    auto* runManager = G4RunManagerFactory::CreateRunManager();
    ASSERT_NE(runManager, nullptr);

    EXPECT_NO_THROW({
        runManager->SetUserInitialization(new DetectorConstruction());
        runManager->SetUserInitialization(new PhysicsList());
        runManager->SetUserInitialization(new ActionInitialization());
        runManager->Initialize();
    });

    delete runManager;
}

TEST(Geant4Setup, DetectsTendonLogicalVolume) {
    DetectorConstruction det;
    auto world = det.Construct();
    ASSERT_NE(world, nullptr);
    // Check that the tendon logical volume exists in the logical volume store
    auto store = G4LogicalVolumeStore::GetInstance();
    auto tendon = store->GetVolume("PatellarTendonLog", false);
    EXPECT_NE(tendon, nullptr);
}
