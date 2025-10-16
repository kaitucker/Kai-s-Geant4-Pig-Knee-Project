#ifndef ACTIONINITIALIZATION_HH
#define ACTIONINITIALIZATION_HH 1

#include "G4VUserActionInitialization.hh"

class ActionInitialization : public G4VUserActionInitialization {
public:
    ActionInitialization();
    virtual ~ActionInitialization();
    virtual void BuildForMaster() const override;
    virtual void Build() const override;
};

#endif
