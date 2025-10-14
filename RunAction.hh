#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"

class G4Run;

class RunAction : public G4UserRunAction {
public:
  RunAction();
  virtual ~RunAction();
  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);
  void AddEdep(G4double edep) { fTotalEdep += edep; }
private:
  G4double fTotalEdep;
};
#endif
