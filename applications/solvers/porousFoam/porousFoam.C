/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    porousFoam

Description
    Dissolution of a porous matrix in three stages:
    Solves for Brinkman flow in a porous matrix
    Solves for transport (steady-state) of reactant and dissolution rate
    Updates porosity field according to local reaction rate (midpoint)

    Author  Tony Ladd (tladd@che.ufl.edu)
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointPatchField.H"
#include "fvMesh.H"
#include "fvIOoptionList.H"
#include "steadyStateControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "readDicts.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    steadyStateControl simple(mesh);        // SIMPLE: no time update

    Info << "Start" << nl;

    if (QCON)
        Info << "Constant flow set: max Q =     " << Qmax << nl;
    if (DARCY)
        Info << "Darcy flow set:    max perm =  " << Kmax << nl;
    if (DEBUG)
        Info << "Debugging turned on" << nl;
    if (Pmax > 1)
        Info << "Pmax > 1: assuming 2D fracture model" << nl << endl;

//  Set time step: dissolution time scale, 1/(k0*gamma)

    dimensionedScalar deltaT = runTime.deltaT()/k0.value();
    Info << "scaled timestep: " << deltaT << endl;

//  Initialize fields (at t=0 and on restart)

    Info << nl << "Time = " << runTime.timeName() << endl;

    #include "createFields.H"               // Create fields
    #include "createUscale.H"               // U(t=0) + velocity scale
    if (runTime.value() == 0)
    {
        #include "solveC.H"                 // C(t=0)
        runTime.writeNow();
    }


//  Main loop: calculate time evolution of the porosity

    while (runTime.run())
    {
        runTime++;
        Info << nl << "Time = " << runTime.timeName() << endl;

        #include "createTmpFields.H"        // Copy Flist -> FTmpList

        for (int n = 0; n < nF; n++)        // Predictor (Heun)
        {
            label Case = n;
            volScalarField& F = FList[n]();
            #include "porosity.H"
            #include "kinetics.H"
            F += R*deltaT;                  // Update mineral fraction
            F.correctBoundaryConditions();
        } 

        #include "porosity.H"               // Solve for new fields
        #include "permModel.H"
        #include "solveU.H"
        #include "solveC.H"

        for (int n = 0; n < nF; n++)        // Corrector (Heun)
        {
            label Case = n;
            volScalarField& F = FList[n]();
            volScalarField& FTmp = FTmpList[n]();
            #include "porosity.H"
            #include "kinetics.H"
            FTmp += R*deltaT;
            F = 0.5*(F + FTmp);             // Average F
            F.correctBoundaryConditions();
        }

        #include "porosity.H"
        #include "permModel.H"
        #include "solveU.H"
        #include "solveC.H"

        runTime.write();
    }

    Info << "End" << endl;
    return 0;
}


// ************************************************************************* //
