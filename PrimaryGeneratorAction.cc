#include "PrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>
PrimaryGeneratorAction::PrimaryGeneratorAction() {
fParticleGun = new G4ParticleGun(1);
auto particleTable = G4ParticleTable::GetParticleTable();
auto gamma = particleTable->FindParticle("gamma");
fParticleGun->SetParticleDefinition(gamma);
// Default energy (monoenergetic)
fParticleGun->SetParticleEnergy(150*keV);
// Default source position, simulating an X-ray tube anode.
// This is positioned far from the phantom, centered on the collimator.
// X = 15.5 cm (center of plexi box)
// Z = 100 cm (typical source position)
fParticleGun->SetParticlePosition(G4ThreeVector(15.5*cm, 0., 100.0*cm));
}
PrimaryGeneratorAction::~PrimaryGeneratorAction() { delete fParticleGun; }
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
// To simulate an X-ray tube, the source should emit photons in a wide cone
// that fully illuminates the collimator jaws. The jaws then shape the beam.
G4double coneHalfAngleDeg = 45.0; // Wide cone to cover the collimator
G4double cosMax = std::cos(coneHalfAngleDeg * CLHEP::deg);
G4double u = G4UniformRand();
G4double cosTheta = cosMax + (1.0 - cosMax) * u; // uniform in cos(theta)
G4double sinTheta = std::sqrt(1.0 - cosTheta*cosTheta);
G4double phi = CLHEP::twopi * G4UniformRand();
// Direction pointing downward along negative z towards the phantom
G4double dx = sinTheta * std::cos(phi);
G4double dy = sinTheta * std::sin(phi);
G4double dz = -cosTheta; // negative z
fParticleGun->SetParticleMomentumDirection(G4ThreeVector(dx, dy, dz));
fParticleGun->GeneratePrimaryVertex(anEvent);
}
