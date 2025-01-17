/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    initDynamicMesh

Description
    In-place mesh adaptation to be alternated with setFields or VOFSetFields
    to get a clean interface

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
    Info<< "Reading field alphaGas and U\n" << endl;

    volScalarField alpha
    (
        IOobject
        (
            "alphaGas",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
        
    volScalarField refinementField
    (
        IOobject
        (
            "refinementField",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mag(fvc::grad(alpha))
    );
    refinementField.internalField() *= pow(mesh.V(),1.0/3.0);

    runTime.setDeltaT(1e-6);
    runTime++;
    
	Info<< "Time = " << runTime.timeName() << nl << endl;
	
	scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();
	mesh.update();
	
	Info<< "Execution time for mesh.update() = "
		<< runTime.elapsedCpuTime() - timeBeforeMeshUpdate
		<< " s" << endl;
	
	runTime.writeNow();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
