#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
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
// create the RunAction first so pointers can be passed to other actions
RunAction* runAction = new RunAction();
SetUserAction(new PrimaryGeneratorAction()); // must be registered
SetUserAction(runAction); // RunAction for this worker
// Register SteppingAction and pass runAction pointer so it can call AddEdep(...)
SetUserAction(new SteppingAction(runAction, "PatellarTendonLog"));

}
