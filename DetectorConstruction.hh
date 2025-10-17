#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1
#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();
    virtual G4VPhysicalVolume* Construct() override;
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

    private:
      // ADD THIS MEMBER VARIABLE
      G4LogicalVolume* fScoringVolume = nullptr;
};

#endif
