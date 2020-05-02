// Component for TsSphTEPC
//
// ********************************************************************
// *                                                                  *
// * This  code implementation is an extension to developments by the *
// * TOPAS collaboration.                                             *
// * This extension  is an freely  available in  accordance  with the *
// * freeBSD license:                                                 *
// * Copyright (c) <2015>, <Harald Paganetti>                         *
// * All rights reserved.                                             *
// * Redistribution    and   use in   source and   binary    forms,   *
// * with or without modification, are permitted provided that the    *
// * following conditions are met:                                    *
// *                                                                  *
// *                                                                  *
// * 1. Redistributions of source code must retain the above          *
// * copyright notice, this                                           *
// * list of conditions and the following disclaimer.                 *
// * 2. Redistributions in binary form must reproduce the above       *
// * copyright notice, this list of conditions and the following      *
// * disclaimer in the documentation and/or other materials provided  *
// * with the distribution.                                           *
// *                                                                  *
// * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
// * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
// * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
// * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
// * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR             *
// * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
// * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT *
// * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF *
// * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED  *
// * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT      *
// * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING   *
// * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF   *
// * THE POSSIBILITY OF SUCH DAMAGE.                                  *
// *                                                                  *
// * The views and conclusions contained in the software and          *
// * documentation are those of the authors and should not be         *
// * interpreted as representing official policies, either expressed  *
// * or implied, of the FreeBSD Project.                              *
// *                                                                  *
// * Contacts: Jan Schuemann, jschuemann@mgh.harvard.edu              *
// *           Harald Paganetti, hpaganetti@mgh.harvard.edu           *
// *                                                                  *
// ********************************************************************
//

#include "TsSphTEPC.hh"

#include "TsParameterManager.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcommand.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

//bool
#include "G4BooleanSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4Sphere.hh"


TsSphTEPC::TsSphTEPC(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			   TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
}


TsSphTEPC::~TsSphTEPC()
{
}


G4VPhysicalVolume* TsSphTEPC::Construct()
{
	BeginConstruction();

	G4NistManager* man = G4NistManager::Instance();



	G4Material* G4MaA150 = man->FindOrBuildMaterial("G4_A-150_TISSUE");
	G4Material* G4MaAl = man->FindOrBuildMaterial("G4_Al");
	G4Material* G4MaTEgas = man->FindOrBuildMaterial("TEgas");

	G4String MaAl = "G4_Al";
	G4String MaA150 = "G4_A-150_TISSUE";
	G4String MaTEgas = "TEgas";

	//-------------------------------------------------------------------------------------------
	//-------------------------- Detector Stats for Region Generation ---------------------------
	//-------------------------------------------------------------------------------------------

	G4double AlShellOuterDiameter = fPm->GetDoubleParameter(GetFullParmName("AlShellOuterDiameter"),"Length");
	G4double AlShellThickness = fPm->GetDoubleParameter(GetFullParmName("AlShellThickness"),"Length");

	G4double A150ShellInnerDiameter = fPm->GetDoubleParameter(GetFullParmName("A150ShellInnerDiameter"),"Length");
	G4double A150ShellThickness = fPm->GetDoubleParameter(GetFullParmName("A150ShellThickness"),"Length");


	//------------------------------------------------------------------------------------
	//------------------------- Spherical TEPC Detector  ---------------------------------
	//------------------------------------------------------------------------------------
	G4double outerRadius = AlShellOuterDiameter / 2.;
	G4double innerRadius = 0. * mm;

	G4double SPhi = 0. * deg;
	G4double DPhi = 360. * deg;
	G4double STheta = 0. * deg;
	G4double DTheta = 180. * deg;

	G4Sphere* Detector = new G4Sphere(fName, innerRadius, outerRadius, SPhi, DPhi, STheta, DTheta);
	fEnvelopeLog = CreateLogicalVolume(Detector);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);


	G4String subComponentName;


	//----------------------------------- AlShellOut, Al ---------------------------------------
	outerRadius = AlShellOuterDiameter / 2.;
	innerRadius = AlShellOuterDiameter / 2. - AlShellThickness;

	subComponentName = "AlShellOut";
	G4cout<<"AlShellOut = "<<innerRadius<<","<<outerRadius<<G4endl;
	G4Sphere* solidAlShellOut = new G4Sphere(subComponentName, innerRadius, outerRadius, SPhi, DPhi, STheta, DTheta);

	G4LogicalVolume* logicAlShellOut    = CreateLogicalVolume(subComponentName, MaAl, solidAlShellOut);
	G4VPhysicalVolume* physicAlShellOut = CreatePhysicalVolume(subComponentName, logicAlShellOut, new G4RotationMatrix(), new G4ThreeVector(0,0,0),fEnvelopePhys);


	//----------------------------------- AlShellIn, TEgas ---------------------------------------
	outerRadius = AlShellOuterDiameter / 2. - AlShellThickness;
	innerRadius = A150ShellInnerDiameter / 2. + A150ShellThickness;

	subComponentName = "AlShellIn";
	G4cout<<"AlShellIn = "<<innerRadius<<","<<outerRadius<<G4endl;
	G4Sphere* solidAlShellIn = new G4Sphere(subComponentName, innerRadius, outerRadius, SPhi, DPhi, STheta, DTheta);

	G4LogicalVolume* logicAlShellIn   = CreateLogicalVolume(subComponentName, MaTEgas, solidAlShellIn);
	G4VPhysicalVolume* physicAlShellIn = CreatePhysicalVolume(subComponentName, logicAlShellIn, new G4RotationMatrix(), new G4ThreeVector(0,0,0),fEnvelopePhys);


	//----------------------------------- A150 ---------------------------------------

	outerRadius = A150ShellInnerDiameter / 2. + A150ShellThickness;
	innerRadius = A150ShellInnerDiameter / 2. + A150ShellThickness / 30000.; // this is a hack, very abrupt change if no hack term... it cannot be too small also
	subComponentName = "A150";
	G4cout<<"A150 = "<<innerRadius<<","<<outerRadius<<G4endl;
	G4Sphere* solidA150 = new G4Sphere(subComponentName, innerRadius, outerRadius, SPhi, DPhi, STheta, DTheta);

	G4LogicalVolume* logicA150   = CreateLogicalVolume(subComponentName, MaA150, solidA150);
	G4VPhysicalVolume* physicA150 = CreatePhysicalVolume(subComponentName, logicA150, new G4RotationMatrix(), new G4ThreeVector(0,0,0),fEnvelopePhys);

	//----------------------------------- TEgasSV ---------------------------------------

	outerRadius = A150ShellInnerDiameter / 2. ;
	innerRadius = 0. * mm;

	subComponentName = "TEgasSV";
	G4cout<<"TEgasSV = "<<innerRadius<<","<<outerRadius<<G4endl;
	G4Sphere* solidTEgasSV = new G4Sphere(subComponentName, innerRadius, outerRadius, SPhi, DPhi, STheta, DTheta);

	G4LogicalVolume* logicTEgasSV   = CreateLogicalVolume(subComponentName, MaTEgas, solidTEgasSV);
	G4VPhysicalVolume* physicTEgasSV = CreatePhysicalVolume(subComponentName, logicTEgasSV, new G4RotationMatrix(), new G4ThreeVector(0,0,0),fEnvelopePhys);


	InstantiateChildren(fEnvelopePhys);

	return fEnvelopePhys;
}
