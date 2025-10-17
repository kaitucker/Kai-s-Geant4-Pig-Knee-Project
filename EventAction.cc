#include "EventAction.hh"
#include "RunAction.hh"
#include "G4Event.hh"

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.)
{}

EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction(const G4Event*)
{
  // Reset the energy deposit accumulator at the start of each event
  fEdep = 0.;
}

void EventAction::EndOfEventAction(const G4Event*)
{
  // At the end of the event, add the total accumulated energy
  // to the RunAction for this event.
  fRunAction->AddEdep(fEdep);
}

void EventAction::AddEdep(G4double edep)
{
  // Accumulate energy deposited in a step
  fEdep += edep;
}
