#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "EventAction.hh"

ActionInitialization::ActionInitialization()
: G4VUserActionInitialization()
{}
ActionInitialization::~ActionInitialization()
{}
void ActionInitialization::BuildForMaster() const
{
// Only the RunAction is needed on the master for summaries
SetUserAction(new RunAction());
}
//ONLY CARE ABOUT THIS:
void ActionInitialization::Build() const
{
  RunAction* runAction = new RunAction();
  SetUserAction(runAction);

  SetUserAction(new PrimaryGeneratorAction());

  auto eventAction = new EventAction(runAction);
  SetUserAction(eventAction);
}
