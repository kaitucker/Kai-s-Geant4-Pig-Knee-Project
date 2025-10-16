#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class ActionInitialization : public G4VUserActionInitialization
{
public:
  ActionInitialization();
  virtual ~ActionInitialization();

  // master thread actions
  virtual void BuildForMaster() const override;

  // worker threads actions
  virtual void Build() const override;
};

#endif
