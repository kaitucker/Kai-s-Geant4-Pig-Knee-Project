// DetectorConstruction.cc
// Final high-realism pig knee primitive model with patellar tendon detector (mean absorbed dose).
// Place inside Pyrex bowl (open rim facing source). Tendon fascicles run along tendon long axis (X).
#include "DetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4Sphere.hh"
#include "G4Polycone.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4ExtrudedSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
// Scoring
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4PSEnergyDeposit.hh"
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include "G4EmStandardPhysics.hh"

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction() {}
DetectorConstruction::~DetectorConstruction() {}
G4VPhysicalVolume* DetectorConstruction::Construct() {
using namespace CLHEP;
  
// -------------------- TUNABLE PARAMETERS --------------------------------
const G4double tol = 0.01*cm; // geometry tolerance
const G4double tiny_gap = 0.5*mm; // clearance above bowl floor
const G4double shell_margin = 10.0*mm; // soft tissue margin (1 cm)
// Tendon/fascicle parameters (confirmed)
const G4double tendon_len = 61.17*mm; // X axis
const G4double tendon_width = 20.40*mm; // Y axis
const G4double tendon_thick = 11.23*mm; // Z axis
const G4int nFasciclesTarget = 300; // number of fascicles
const G4double fascicle_major = 0.30*mm; // ellipse major axis
const G4double fascicle_minor = 0.20*mm; // ellipse minor axis
const G4int nSegmentsPerFibre = 12; // crimp per fascicle
const G4double crimpAmp = 0.10 * fascicle_major; // 10% amplitude
const G4double crimpLambda = 2.0*mm; // wavelength
const G4double packing_pitch = 0.36*mm; // packing pitch approx
// Bone approximations (from your measurements / defaults)
const G4double femur_condyle_radius = 18.0*mm;
const G4double femoral_separation = 18.0*mm;
const G4double tibia_length = 121.3*mm;
const G4double tibia_shaft_radius = 10.0*mm;
const G4double patella_radius = 10.0*mm;
const G4double fibula_length = 101.8*mm;
const G4double fibula_radius = 6.0*mm;
const G4double cortical_thickness = 3.5*mm;
// Soft tissue composition (75% muscle / 25% fat)
const double muscle_frac = 0.75;
const double fat_frac = 0.25;
// --- Adjustable Collimator Jaw Parameters ---
const G4double x_opening = 12.3*cm; // Total opening width along X
const G4double y_opening = 12.3*cm; // Total opening width along Y
const G4double jaw_thickness = 16.5*cm; // Thickness of jaws along Z
const G4double jaw_length = 20.0*cm; // Length of each jaw block
const G4double jaw_width = 5.0*cm; // Width of each jaw block
// -------------------- Materials ----------------------------------------
G4NistManager* nist = G4NistManager::Instance();
G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
G4Material* lead = nist->FindOrBuildMaterial("G4_Pb");
G4Material* water = nist->FindOrBuildMaterial("G4_WATER");
// Elements
G4Element* elH = new G4Element("Hydrogen","H",1.,1.00794*g/mole);
G4Element* elC = new G4Element("Carbon","C",6.,12.011*g/mole);
G4Element* elN = new G4Element("Nitrogen","N",7.,14.007*g/mole);
G4Element* elO = new G4Element("Oxygen","O",8.,15.999*g/mole);
G4Element* elS = new G4Element("Sulfur","S",16.,32.06*g/mole);
G4Element* elCa = new G4Element("Calcium","Ca",20.,40.078*g/mole);
G4Element* elP = new G4Element("Phosphorus","P",15.,30.973762*g/mole);
G4Element* elB = new G4Element("Boron","B",5.,10.81*g/mole);
G4Element* elNa = new G4Element("Sodium","Na",11.,22.9897*g/mole);
G4Element* elAl = new G4Element("Aluminium","Al",13.,26.9815*g/mole);
G4Element* elSi = new G4Element("Silicon","Si",14.,28.0855*g/mole);
G4Element* elK = new G4Element("Potassium","K",19.,39.0983*g/mole);
G4Element* elSe = new G4Element("Selenium","Se",34.,78.96*g/mole);
// PMMA
G4Material* PMMA = new G4Material("PMMA", 1.19*g/cm3, 3);
PMMA->AddElement(elC, 5);
PMMA->AddElement(elH, 8);
PMMA->AddElement(elO, 2);
// Pyrex glass
G4Material* Pyrex = new G4Material("Pyrex", 2.23*g/cm3, 6);
Pyrex->AddElement(elB, 4.0*perCent);
Pyrex->AddElement(elO, 54.0*perCent);
Pyrex->AddElement(elNa, 2.8*perCent);
Pyrex->AddElement(elAl, 1.1*perCent);
Pyrex->AddElement(elSi,37.7*perCent);
Pyrex->AddElement(elK, 0.3*perCent);
// Amorphous Selenium
G4Material* AmorphSe = new G4Material("AmorphSe", 4.81*g/cm3, 1);
AmorphSe->AddElement(elSe, 1);
// Tendon materials
G4Material* TendonFascicle = new G4Material("TendonFascicle", 1.12*g/cm3, 5);
TendonFascicle->AddElement(elH, 0.064);
TendonFascicle->AddElement(elC, 0.360);
TendonFascicle->AddElement(elN, 0.115);
TendonFascicle->AddElement(elO, 0.454);
TendonFascicle->AddElement(elS, 0.007);
G4Material* TendonMatrix = new G4Material("TendonMatrix", 1.05*g/cm3, 5);
TendonMatrix->AddElement(elH, 0.111);
TendonMatrix->AddElement(elC, 0.104);
TendonMatrix->AddElement(elN, 0.029);
TendonMatrix->AddElement(elO, 0.754); // Corrected from 0.729 to make sum = 1.0
TendonMatrix->AddElement(elS, 0.002);
// Bones and cartilage
G4Material* CorticalBone = new G4Material("CorticalBone", 1.90*g/cm3, 6);
CorticalBone->AddElement(elH, 0.045);
CorticalBone->AddElement(elC, 0.270);
CorticalBone->AddElement(elN, 0.044);
CorticalBone->AddElement(elO, 0.425);
CorticalBone->AddElement(elCa,0.125);
CorticalBone->AddElement(elP, 0.091);
G4Material* TrabecularBone = new G4Material("TrabecularBone", 0.60*g/cm3, 6);
TrabecularBone->AddElement(elH, 0.060);
TrabecularBone->AddElement(elC, 0.330);
TrabecularBone->AddElement(elN, 0.044);
TrabecularBone->AddElement(elO, 0.470);
TrabecularBone->AddElement(elCa,0.040);
TrabecularBone->AddElement(elP, 0.056);
G4Material* Cartilage = new G4Material("Cartilage", 1.05*g/cm3, 5);
Cartilage->AddElement(elH, 0.104);
Cartilage->AddElement(elC, 0.140);
Cartilage->AddElement(elN, 0.033);
Cartilage->AddElement(elO, 0.717);
Cartilage->AddElement(elS, 0.006);
// Muscle & Fat & shell mixture
G4Material* Muscle = new G4Material("Muscle", 1.06*g/cm3, 4);
Muscle->AddElement(elH, 0.102);
Muscle->AddElement(elC, 0.143);
Muscle->AddElement(elN, 0.034);
Muscle->AddElement(elO, 0.721);
G4Material* Fat = new G4Material("Fat", 0.92*g/cm3, 4);
Fat->AddElement(elH, 0.120);
Fat->AddElement(elC, 0.600);
Fat->AddElement(elN, 0.005);
Fat->AddElement(elO, 0.275);
// create shell material as mixture by mass fraction approx using density-weighted approach
G4double rho_shell = muscle_frac * Muscle->GetDensity() + fat_frac * Fat->GetDensity();
G4Material* SoftShell = new G4Material("SoftShell", rho_shell, 2);
SoftShell->AddMaterial(Muscle, muscle_frac);
SoftShell->AddMaterial(Fat, fat_frac);
// -------------------- Scene parameters (unchanged where possible) --------
const G4double boxX = 104.0*cm, boxY = 70.0*cm, boxZ = 196.0*cm;
const G4double detector_from_bottom = 30.0*cm;
const G4double z_bottom = -boxZ/2.0;
const G4double z_detector_plane = z_bottom + detector_from_bottom;
const G4double det_thickness = 1.0*cm;
const G4double z_det_center = z_detector_plane + det_thickness/2.0;
const G4double plexi_in_len = 56.0*cm, plexi_in_dep = 55.5*cm;
const G4double z_top_plexi = 30.0*cm;
const G4double plexi_internal_height = z_top_plexi - z_detector_plane;
const G4double plexi_thickness = 0.5*cm;
const G4double plexi_left_from_boxleft = 39.5*cm;
const G4double x_box_left = -boxX/2.0;
const G4double x_plexi_center = x_box_left + plexi_left_from_boxleft + plexi_in_len/2.0;
const G4double y_plexi_center = 0.0;
const G4double z_plexi_center = (z_top_plexi + z_detector_plane)/2.0;
const G4double top_open_W = 30.0*cm, top_open_L = 30.0*cm;
const G4double holder_size = 40.0*cm, holder_hole_dia = 19.7*cm;
const G4double holder_thickness = 0.5*cm;
const G4double z_holder_top = z_detector_plane + 37.0*cm;
const G4double z_holder_center = z_holder_top - holder_thickness/2.0;
const G4double bowl_total_height = 11.1*cm;
const G4double bowl_rim_outer_dia = 30.0*cm;
const G4double bowl_rim_inner_dia = 27.2*cm;
const G4double bowl_base_dia = 11.1*cm;
const G4double side_wall_thickness = 0.4*cm;
const G4double z_bowl_rim_top = z_holder_top;
const G4double z_bowl_center = z_bowl_rim_top - bowl_total_height/2.0;
const G4double det_X = 47.0*cm, det_Y = 47.0*cm, det_Z = det_thickness;
const G4double det_center_x = x_plexi_center, det_center_y = 0.0, det_center_z = z_det_center;
// Collimator position
const G4double z_coll_exit = z_detector_plane + 106.0*cm;
const G4double z_col_center = z_coll_exit + jaw_thickness/2.0;
const G4double x_col_center = x_plexi_center, y_col_center = 0.0;
// -------------------- World & cavity ----------------------------------
auto solidWorld = new G4Box("WorldSolid", (boxX+20*cm)/2.0, (boxY+20*cm)/2.0, (boxZ+40*cm)/2.0);
auto logicWorld = new G4LogicalVolume(solidWorld, air, "World");
auto physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, false);
auto solidCavity = new G4Box("Cavity", boxX/2.0 - tol, boxY/2.0 - tol, boxZ/2.0 - tol);
auto logicCavity = new G4LogicalVolume(solidCavity, air, "Cavity");
new G4PVPlacement(0, G4ThreeVector(), logicCavity, "Cavity", logicWorld, false, 0, false);
// Plexi shell with hole
auto plexiOuter = new G4Box("PlexiOuter", (plexi_in_len+2*plexi_thickness)/2.0,
(plexi_in_dep+2*plexi_thickness)/2.0,
(plexi_internal_height+2*plexi_thickness)/2.0);
auto plexiInner = new G4Box("PlexiInner", plexi_in_len/2.0 - tol,
plexi_in_dep/2.0 - tol,
plexi_internal_height/2.0 - tol);
auto plexiShell = new G4SubtractionSolid("PlexiShell", plexiOuter, plexiInner, 0, G4ThreeVector(0,0,-plexi_thickness));
auto openBox = new G4Box("TopOpenBox", top_open_L/2.0+0.1*cm, top_open_W/2.0+0.1*cm,
(plexi_internal_height+2*plexi_thickness)/2.0+0.1*cm);
auto plexiShellWithHole = new G4SubtractionSolid("PlexiShellWithHole",
plexiShell, openBox, 0,
G4ThreeVector(0,0,(plexi_internal_height+2*plexi_thickness)/2.0-0.1*cm));
auto plexiLog = new G4LogicalVolume(plexiShellWithHole, PMMA, "PlexiShell");
new G4PVPlacement(0, G4ThreeVector(x_plexi_center, y_plexi_center, z_plexi_center),
plexiLog, "PlexiShell", logicCavity, false, 0, false);
// Holder
auto holderOuter = new G4Box("HolderOuter", holder_size/2.0, holder_size/2.0, holder_thickness/2.0);
auto holderHole = new G4Tubs("HolderHole", 0.0, holder_hole_dia/2.0+0.001*cm,
holder_thickness/2.0+0.001*cm, 0, 360*deg);
auto holderSolid = new G4SubtractionSolid("HolderSolid", holderOuter, holderHole);
auto holderLog = new G4LogicalVolume(holderSolid, PMMA, "Holder");
new G4PVPlacement(0, G4ThreeVector(x_plexi_center,y_plexi_center,z_holder_center),
holderLog, "Holder", logicCavity, false, 0, false);
// Pyrex bowl (polycone), rotate so open rim faces +Z (toward source)
std::vector<G4double> zPlanes_local = {0.0, 2.5*cm, 7.0*cm, bowl_total_height-1.5*cm, bowl_total_height};
std::vector<G4double> rOuter_local = {bowl_rim_outer_dia/2.0,
bowl_rim_outer_dia/2.0-0.8*cm,
bowl_rim_inner_dia/2.0-0.5*cm,
bowl_base_dia/2.0+0.2*cm,
bowl_base_dia/2.0+0.2*cm};
std::vector<G4double> rInner_local = {bowl_rim_inner_dia/2.0,
bowl_rim_inner_dia/2.0-0.8*cm+side_wall_thickness,
bowl_base_dia/2.0+side_wall_thickness,
bowl_base_dia/2.0,
0.0};
int nbPlanes = (int)zPlanes_local.size();
std::vector<G4double> zVec(nbPlanes), rInnerVec(nbPlanes), rOuterVec(nbPlanes);
for(int i=0;i<nbPlanes;++i){
zVec[i] = zPlanes_local[i] - bowl_total_height/2.0;
rOuterVec[i] = rOuter_local[i];
rInnerVec[i] = rInner_local[i];
}
auto bowlOuter = new G4Polycone("BowlOuter",0,360*deg,nbPlanes,zVec.data(),rInnerVec.data(),rOuterVec.data());
std::vector<G4double> zeroInner(nbPlanes,0.0), innerHollow(nbPlanes);
for(int i=0;i<nbPlanes;++i) innerHollow[i]=std::max(0.0*cm,rInnerVec[i]-side_wall_thickness);
auto bowlInner = new G4Polycone("BowlInner",0,360*deg,nbPlanes,zVec.data(),zeroInner.data(),innerHollow.data());
auto bowlShell = new G4SubtractionSolid("BowlShell",bowlOuter,bowlInner);
auto bowlLog = new G4LogicalVolume(bowlShell,Pyrex,"PyrexBowl");
G4RotationMatrix rotBowl;
rotBowl.rotateX(180.*deg);
new G4PVPlacement(G4Transform3D(rotBowl, G4ThreeVector(x_plexi_center,y_plexi_center,z_bowl_center)),
bowlLog,"PyrexBowl",logicCavity,false,0,false);
// -------------------- Build knee primitives (femur condyles, tibia, fibula, patella) ----
// Femoral condyles: two ellipsoids unioned
auto condL = new G4Ellipsoid("CondL", femur_condyle_radius, femur_condyle_radius*0.9, femur_condyle_radius*0.8);
auto condR = new G4Ellipsoid("CondR", femur_condyle_radius, femur_condyle_radius*0.9, femur_condyle_radius*0.8);
G4ThreeVector condLpos(0.0, -femoral_separation/2.0, 0.0);
G4ThreeVector condRpos(0.0, +femoral_separation/2.0, 0.0);
auto condUnion = new G4UnionSolid("CondUnion", condL, condR, 0, condRpos - condLpos);
// Tibia cylinder
auto tibiaCyl = new G4Tubs("TibiaCyl", 0.0, tibia_shaft_radius, tibia_length/2.0, 0.*deg, 360.*deg);
G4ThreeVector tibiaRel(-tibia_length/4.0, 0.0, - (femur_condyle_radius*0.6 + tibia_length/2.0) );
auto kneeWithTibia = new G4UnionSolid("KneeWithTibia", condUnion, tibiaCyl, 0, tibiaRel);
// Fibula
auto fibulaCyl = new G4Tubs("FibulaCyl", 0.0, fibula_radius, fibula_length/2.0, 0.*deg, 360.*deg);
G4ThreeVector fibulaRel(-tibia_length/4.0, femoral_separation/2.0 + 10.0*mm, -tibia_length/4.0);
auto kneeWithFib = new G4UnionSolid("KneeWithFib", kneeWithTibia, fibulaCyl, 0, fibulaRel);
// Patella sphere (small ellipsoid would be better, but sphere acceptable)
auto patella = new G4Sphere("Patella", 0., patella_radius, 0.*deg, 360.*deg, 0.*deg, 180.*deg);
G4ThreeVector patellaRelPos(femur_condyle_radius*0.4, 0.0, femur_condyle_radius*0.9);
auto kneeOuter = new G4UnionSolid("KneeOuter", kneeWithFib, patella, 0, patellaRelPos);
// Inner ellipsoid approximating trabecular core
G4double inner_a = femur_condyle_radius + 10.0*mm;
G4double inner_b = femur_condyle_radius + 6.0*mm;
G4double inner_c = femur_condyle_radius + 12.0*mm;
if (inner_a < 10.0*mm) inner_a = 10.0*mm;
auto innerEll = new G4Ellipsoid("InnerEll", inner_a, inner_b, inner_c);
// Cortical shell and trabecular core (placed at knee center)
auto corticalSolid = new G4SubtractionSolid("CorticalSolid", kneeOuter, innerEll, 0, G4ThreeVector(0.,0.,0.));
auto corticalLog = new G4LogicalVolume(corticalSolid, CorticalBone, "CorticalLog");
auto trabLog = new G4LogicalVolume(innerEll, TrabecularBone, "TrabularLog");
// Place knee in cavity, resting on the bottom of the bowl
const G4double x_pig = x_plexi_center;
const G4double y_pig = 0.0;
// To make the placement realistic, the knee should rest on the bowl floor.
// First, find the lowest point of the bone assembly in its local coordinate system.
// This is determined by the bottom of the tibia cylinder.
const G4double z_lowest_bone_local = tibiaRel.z() - tibia_length/2.0;
// Next, find the Z coordinate of the inner surface of the bowl's bottom.
// The bowl is rotated 180 deg, so its local +Z (the bottom) points to world -Z.
const G4double z_bowl_floor = z_bowl_center - bowl_total_height/2.0;
// Now, set the knee's origin position (z_pig) so that its lowest point
// is sitting tiny_gap above the bowl floor.
const G4double z_pig = z_bowl_floor - z_lowest_bone_local + tiny_gap;
new G4PVPlacement(0, G4ThreeVector(x_pig, y_pig, z_pig), corticalLog, "KneeCortical", logicCavity, false, 0, false);
new G4PVPlacement(0, G4ThreeVector(x_pig, y_pig, z_pig), trabLog, "KneeTrabular", logicCavity, false, 0, false);
// Cartilage & menisci approximations
auto cartRegion = new G4Ellipsoid("CartRegion", inner_a*0.5, inner_b*0.6, inner_c*0.25);
auto cartLog = new G4LogicalVolume(cartRegion, Cartilage, "CartLog");
G4ThreeVector cartWorldPos = G4ThreeVector(x_pig, y_pig, z_pig + inner_c*0.5 - 2.0*mm);
new G4PVPlacement(0, cartWorldPos, cartLog, "Cartilage", logicCavity, false, 0, false);
const G4double men_inner = 6.0*mm, men_outer = 12.0*mm, men_thick = 5.0*mm;
auto menL = new G4Tubs("MenL", men_inner, men_outer, men_thick/2.0, 0.*deg, 180.*deg);
auto menR = new G4Tubs("MenR", men_inner, men_outer, men_thick/2.0, 0.*deg, 180.*deg);
G4ThreeVector menLpos(-8.0*mm, -femoral_separation/2.0, -femur_condyle_radius*0.2);
G4ThreeVector menRpos(-8.0*mm, +femoral_separation/2.0, -femur_condyle_radius*0.2);
auto menLLog = new G4LogicalVolume(menL, Cartilage, "MenLLog");
auto menRLog = new G4LogicalVolume(menR, Cartilage, "MenRLog");
new G4PVPlacement(0, G4ThreeVector(x_pig + menLpos.x(), menLpos.y(), z_pig + menLpos.z()), menLLog, "MenL", logicCavity, false, 0, false);
new G4PVPlacement(0, G4ThreeVector(x_pig + menRpos.x(), menRpos.y(), z_pig + menRpos.z()), menRLog, "MenR", logicCavity, false, 0, false);
// Ligaments (ACL, PCL, MCL, LCL) as cylinders
auto ligRadius = 2.5*mm;
auto ACL = new G4Tubs("ACL", 0., ligRadius, 20.*mm, 0.*deg, 360.*deg);
auto ACLLog = new G4LogicalVolume(ACL, TendonMatrix, "ACLLog");
new G4PVPlacement(0, G4ThreeVector(x_pig - 10.0*mm, 0.0, z_pig - 5.0*mm), ACLLog, "ACL", logicCavity, false, 0, false);
auto PCL = new G4Tubs("PCL", 0., ligRadius, 20.*mm, 0.*deg, 360.*deg);
auto PCLLog = new G4LogicalVolume(PCL, TendonMatrix, "PCLLog");
new G4PVPlacement(0, G4ThreeVector(x_pig + 8.0*mm, 0.0, z_pig - 12.0*mm), PCLLog, "PCL", logicCavity, false, 0, false);
auto MCL = new G4Tubs("MCL", 0., 2.5*mm, 25.*mm, 0.*deg, 360.*deg);
auto MCLLog = new G4LogicalVolume(MCL, TendonMatrix, "MCLLog");
new G4PVPlacement(0, G4ThreeVector(x_pig - 2.0*mm, -femoral_separation/2.0 - 6.0*mm, z_pig), MCLLog, "MCL", logicCavity, false, 0, false);
auto LCL = new G4Tubs("LCL", 0., 2.5*mm, 25.*mm, 0.*deg, 360.*deg);
auto LCLLog = new G4LogicalVolume(LCL, TendonMatrix, "LCLLog");
new G4PVPlacement(0, G4ThreeVector(x_pig - 2.0*mm, +femoral_separation/2.0 + 6.0*mm, z_pig), LCLLog, "LCL", logicCavity, false, 0, false);
// -------------------- Patellar tendon: mother as child of corticalLog (Option A) ----
// Compute tendon center position relative to knee center: anterior face
// anterior Z location ~ z_pig + inner_c*0.5 (approx). Use offset to place tendon slightly outside cortical shell.
G4double anteriorZ = z_pig + inner_c*0.5 + 1.5*mm;
G4ThreeVector tendonWorldPos(x_pig, y_pig, anteriorZ);
// To attach tendon as child of corticalLog, compute tendon position in corticalLog coordinates.
// Because corticalLog was placed at (x_pig, y_pig, z_pig) we can place tendon with translation relative to corticalLog.
G4ThreeVector tendonRelToCortical(0.0, 0.0, anteriorZ - z_pig);
auto tendonBox = new G4Box("TendonBox", tendon_len/2.0, tendon_width/2.0, tendon_thick/2.0);
auto tendonLog = new G4LogicalVolume(tendonBox, TendonMatrix, "PatellarTendonLog");
// Place tendon as daughter of corticalLog so it is anatomically attached (it will override bone material where it overlaps)
new G4PVPlacement(0, tendonRelToCortical, tendonLog, "PatellarTendon", corticalLog, false, 0, false);
// -------------------- Create fascicle segments (extruded elliptical segments) ----
// Helper to create ellipse polygon (Y-Z plane) for G4ExtrudedSolid
auto makeEllipsePolygon = [&](G4double major, G4double minor, G4int nPts=20) {
std::vector<G4TwoVector> poly;
for (G4int k=0;k<nPts;++k) {
double theta = 2.0*M_PI * double(k)/double(nPts);
// After rotY(90), local x -> mother -z, local y -> mother y.
// We want ellipse major axis along mother y, minor along mother z.
double x_local = (minor/2.0) * std::cos(theta); // Corresponds to mother's z-axis extent
double y_local = (major/2.0) * std::sin(theta); // Corresponds to mother's y-axis extent
poly.emplace_back(x_local, y_local);
}
return poly;
};
// Fallback: circular cross-section if G4ExtrudedSolid not available or forced
auto makeFallbackSeg = [&](const G4String& name, G4double segHalfLen, G4double major, G4double minor)->G4VSolid* {
G4double area = M_PI * (major/2.0) * (minor/2.0);
G4double r_eq = std::sqrt(area/M_PI);
return new G4Tubs(name, 0., r_eq, segHalfLen, 0.*deg, 360.*deg);
};
bool forceFallback = false; // set true if your GEANT4 version rejects G4ExtrudedSolid constructor
G4double segHalfLen = tendon_len / (2.0 * nSegmentsPerFibre);
// compute fiber centers in tendon-local coordinates (Y-Z), hex-like packing
// Increased margin from 0.05mm to 0.1mm to prevent splay/crimp from pushing fascicles outside mother volume.
G4double halfY = tendon_width/2.0 - (fascicle_major/2.0) - 0.1*mm;
// Corrected to use fascicle_minor for Z-dimension and increased margin.
G4double halfZ = tendon_thick/2.0 - (fascicle_minor/2.0) - 0.1*mm;
G4double rowSpacing = packing_pitch * std::sqrt(3.0)/2.0;
int maxRows = (int)std::ceil((2.0*halfZ)/rowSpacing) + 4;
int maxCols = (int)std::ceil((2.0*halfY)/packing_pitch) + 4;
std::vector<std::pair<G4double,G4double>> fiberCenters;
for (int rz=-maxRows; rz<=maxRows && (int)fiberCenters.size() < nFasciclesTarget; ++rz) {
G4double zloc = rz * rowSpacing;
bool offset = ((rz % 2) != 0);
for (int cy=-maxCols; cy<=maxCols && (int)fiberCenters.size() < nFasciclesTarget; ++cy) {
G4double yloc = cy * packing_pitch + (offset ? packing_pitch/2.0 : 0.0);
// Use an elliptical boundary condition instead of a rectangular one.
// This prevents splay rotation from pushing fascicles in the "corners" of the grid outside the mother volume.
if ( (yloc/halfY)*(yloc/halfY) + (zloc/halfZ)*(zloc/halfZ) <= 1.0 ) {
fiberCenters.emplace_back(yloc, zloc);
}
}
}
G4cout << "Fiber centers planned: " << fiberCenters.size() << G4endl;
// Place segments as daughters of tendonLog using tendon-local coordinates
// Define a single rotation to align fascicle extrusion axis (local Z) with tendon long axis (mother X)
G4RotationMatrix* rotFascicle = new G4RotationMatrix();
rotFascicle->rotateY(90.*deg);
int globalSegIndex = 0;
for (size_t f=0; f<fiberCenters.size(); ++f) {
G4double y0 = fiberCenters[f].first;
G4double z0 = fiberCenters[f].second;
for (int seg=0; seg<nSegmentsPerFibre; ++seg) {
// center coordinate along tendon (local X)
G4double xLocal = -tendon_len/2.0 + (seg + 0.5) * (tendon_len / nSegmentsPerFibre);
// splay: more at ends
G4double nd = std::abs(xLocal) / (tendon_len/2.0);
G4double splayAngle = 10.0*deg * std::pow(nd, 1.5); // radians
double theta = (y0 >= 0.0) ? splayAngle : -splayAngle;
G4double y_splayed = y0 * std::cos(theta) - z0 * std::sin(theta);
G4double z_splayed = y0 * std::sin(theta) + z0 * std::cos(theta);
// crimp sinusoid added in Y
G4double crimp = crimpAmp * std::sin(2.0*M_PI * xLocal / crimpLambda + double(f % 7));
y_splayed += crimp;
// create segment solid
std::ostringstream segName;
segName << "FibreSeg_f" << f << "_s" << seg;
G4VSolid* segSolid = nullptr;
if (!forceFallback) {
std::vector<G4TwoVector> poly = makeEllipsePolygon(fascicle_major, fascicle_minor, 20);
std::vector<G4ExtrudedSolid::ZSection> sections;
sections.push_back(G4ExtrudedSolid::ZSection(-segHalfLen, G4TwoVector(0.,0.), 1.0));
sections.push_back(G4ExtrudedSolid::ZSection( segHalfLen, G4TwoVector(0.,0.), 1.0));
segSolid = new G4ExtrudedSolid(segName.str(), poly, sections);
} else {
segSolid = makeFallbackSeg(segName.str(), segHalfLen, fascicle_major, fascicle_minor);
}
// logical volume
std::ostringstream lvName;
lvName << "FibreLV_f" << f << "_s" << seg;
G4LogicalVolume* segLog = new G4LogicalVolume(segSolid, TendonFascicle, lvName.str());
// place relative to tendon mother (daughter placement)
G4ThreeVector segPosInTendon(xLocal, y_splayed, z_splayed);
std::ostringstream pvName;
pvName << "FibrePV_f" << f << "_s" << seg;
new G4PVPlacement(rotFascicle, segPosInTendon, segLog, pvName.str(), tendonLog, false, globalSegIndex, false);
++globalSegIndex;
}
}
G4cout << "Placed fascicle segments: " << globalSegIndex << " (~" << globalSegIndex / std::max(1,nSegmentsPerFibre) << " fibers)" << G4endl;
// -------------------- Soft tissue shell (non-spherical envelope) -----------
// A more realistic envelope is formed by the union of four ellipsoids,
// representing the main joint, quadriceps, calf/hamstrings, and tissue along the tibia.
// This creates a more anatomically correct, non-uniform shape.
// 1. Main ellipsoid for the central joint area
G4double shell_f_a = inner_a + shell_margin + 5.0*mm;
G4double shell_f_b = inner_b + shell_margin + 5.0*mm;
G4double shell_f_c = inner_c + shell_margin + 8.0*mm;
auto shellFem = new G4Ellipsoid("ShellFem", shell_f_a, shell_f_b, shell_f_c);
// 2. Tibia-side ellipsoid (shifted down)
auto shellTib = new G4Ellipsoid("ShellTib", shell_f_a*0.8, shell_f_b*0.8, shell_f_c*0.9);
G4ThreeVector shellTibPos(-tibia_length/6.0, 0.0, -tibia_length/4.0);
auto unionShell1 = new G4UnionSolid("UnionShell1", shellFem, shellTib, 0, shellTibPos);
// 3. Quadriceps muscle ellipsoid (superior and anterior)
auto quadShell = new G4Ellipsoid("QuadShell", shell_f_a*0.5, shell_f_b*0.9, shell_f_c*1.0);
G4ThreeVector quadPos(patellaRelPos.x(), 0.0, patellaRelPos.z() + 3.0*cm);
auto unionShell2 = new G4UnionSolid("UnionShell2", unionShell1, quadShell, 0, quadPos);
// 4. Calf/Hamstring muscle ellipsoid (posterior)
auto calfShell = new G4Ellipsoid("CalfShell", shell_f_a*0.6, shell_f_b*0.8, shell_f_c*1.2);
G4ThreeVector calfPos(-femur_condyle_radius - 1.0*cm, 0.0, -femur_condyle_radius);
auto unionShellOuter = new G4UnionSolid("UnionShellOuter", unionShell2, calfShell, 0, calfPos);
// To prevent the soft tissue from extending through the holder, we subtract
// a volume representing the space above the holder's top surface. This simulates
// the specimen being trimmed to fit the apparatus.
G4double z_cutoff_local = z_holder_top - z_pig; // Z-plane of the holder top in the knee's local coords
G4double sub_size = 50.0*cm; // A large size to ensure a complete cut
auto subtractionBox = new G4Box("SubtractionBox", sub_size/2.0, sub_size/2.0, sub_size/2.0);
G4ThreeVector subBoxPos(0.0, 0.0, z_cutoff_local + sub_size/2.0); // Position box to cut everything above the plane
auto unionShellCut = new G4SubtractionSolid("UnionShellCut", unionShellOuter, subtractionBox, 0, subBoxPos);
// Subtract inner bone/cartilage region to create a hollow shell
auto shellInner = new G4Ellipsoid("ShellInner", inner_a, inner_b, inner_c);
auto softShellSolid = new G4SubtractionSolid("SoftShell", unionShellCut, shellInner, 0, G4ThreeVector(0.,0.,0.));
auto softShellLog = new G4LogicalVolume(softShellSolid, SoftShell, "SoftShellLog");
// Place soft shell at knee center (world position x_pig, y_pig, z_pig)
new G4PVPlacement(0, G4ThreeVector(x_pig, y_pig, z_pig), softShellLog, "SoftShell", logicCavity, false, 0, false);
// -------------------- Scoring: Single MFD attached to tendon mother ----------------
// Proper detector registration for modern Geant4

G4SDManager* sdman = G4SDManager::GetSDMpointer();

// Create a MultiFunctionalDetector with energy deposit primitive
G4MultiFunctionalDetector* tendonMFD = new G4MultiFunctionalDetector("TendonMFD");
G4PSEnergyDeposit* edep = new G4PSEnergyDeposit("Edep");
tendonMFD->RegisterPrimitive(edep);

// Register the detector with the SD manager FIRST
sdman->AddNewDetector(tendonMFD);

// Then attach it to the tendon logical volume
SetSensitiveDetector(tendonLog, tendonMFD);

G4cout << "Tendon MFD registered and attached to PatellarTendonLog." << G4endl;
// --- NEW: Adjustable 4-jaw collimator ---
auto jaw_solid = new G4Box("JawSolid", jaw_width/2.0, jaw_length/2.0, jaw_thickness/2.0);
auto jaw_log = new G4LogicalVolume(jaw_solid, lead, "JawLog");
// Position the X jaws
G4double x_jaw_offset = x_opening/2.0 + jaw_width/2.0;
G4ThreeVector jawX1_pos(x_col_center - x_jaw_offset, y_col_center, z_col_center);
G4ThreeVector jawX2_pos(x_col_center + x_jaw_offset, y_col_center, z_col_center);
new G4PVPlacement(0, jawX1_pos, jaw_log, "JawX1", logicCavity, false, 0, false);
new G4PVPlacement(0, jawX2_pos, jaw_log, "JawX2", logicCavity, false, 1, false);
// Position the Y jaws (rotated 90 degrees)
G4RotationMatrix* rot_jaw = new G4RotationMatrix();
rot_jaw->rotateZ(90.*deg);
G4double y_jaw_offset = y_opening/2.0 + jaw_width/2.0;
G4ThreeVector jawY1_pos(x_col_center, y_col_center - y_jaw_offset, z_col_center);
G4ThreeVector jawY2_pos(x_col_center, y_col_center + y_jaw_offset, z_col_center);
new G4PVPlacement(rot_jaw, jawY1_pos, jaw_log, "JawY1", logicCavity, false, 2, false);
new G4PVPlacement(rot_jaw, jawY2_pos, jaw_log, "JawY2", logicCavity, false, 3, false);
// -------------------- Visualization attributes --------------------------
plexiLog->SetVisAttributes(new G4VisAttributes(G4Colour(0.85,0.85,0.65)));
holderLog->SetVisAttributes(new G4VisAttributes(G4Colour(0.85,0.85,0.65)));
bowlLog->SetVisAttributes(new G4VisAttributes(G4Colour(0.6,0.6,0.9,0.35)));
corticalLog->SetVisAttributes(new G4VisAttributes(G4Colour(0.95,0.9,0.8,0.6)));
trabLog->SetVisAttributes(new G4VisAttributes(G4Colour(0.9,0.6,0.6,0.45)));
cartLog->SetVisAttributes(new G4VisAttributes(G4Colour(0.7,0.9,1.0,0.6)));
menLLog->SetVisAttributes(new G4VisAttributes(G4Colour(0.6,0.9,0.9,0.7)));
menRLog->SetVisAttributes(new G4VisAttributes(G4Colour(0.6,0.9,0.9,0.7)));
ACLLog->SetVisAttributes(new G4VisAttributes(G4Colour(1.0,0.6,0.2)));
PCLLog->SetVisAttributes(new G4VisAttributes(G4Colour(1.0,0.6,0.2)));
MCLLog->SetVisAttributes(new G4VisAttributes(G4Colour(1.0,0.6,0.2)));
LCLLog->SetVisAttributes(new G4VisAttributes(G4Colour(1.0,0.6,0.2)));
tendonLog->SetVisAttributes(new G4VisAttributes(G4Colour(0.2,0.8,0.2,0.75)));
softShellLog->SetVisAttributes(new G4VisAttributes(G4Colour(0.95,0.85,0.8,0.45)));
jaw_log->SetVisAttributes(new G4VisAttributes(G4Colour(0.15,0.15,0.15))); // vis for new jaws
logicCavity->SetVisAttributes(new G4VisAttributes(G4Colour(0.95,0.95,0.95,0.0)));
G4cout << "DetectorConstruction: built anatomically-oriented pig knee with tendon fascicles (detector = tendon mother)." << G4endl;
return physWorld;
}
