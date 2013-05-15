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
    VOFSetFields

Description
    Set volume fraction fields that are not aligned with the mesh

    Future development as it becomes useful:
       * Be able to set multiple planes, circles, and spheres, and specify
          build up more complex shapes
       * Add quadratic planes for curved surfaces

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvCFD.H"
#include "cuttableCell.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    const volVectorField& cellCenters = mesh.C();

    word shape(VOFDict.lookup("shape"));
    label solidCells = 0;
    label gasCells = 0;

    if (shape == "circle")
    {
        vector center(VOFDict.lookup("center"));
        scalar radius = readScalar(VOFDict.lookup("radius"));
        vector Usolid(VOFDict.lookup("Usolid"));
        vector Ugas(VOFDict.lookup("Ugas"));
        word planename(VOFDict.lookup("plane"));
        word inv(VOFDict.lookup("inverse"));

        forAll(alpha, cellI)
        {
            vector r = cellCenters[cellI] - center;
            if (planename=="xy")
            {
                r.z() = 0.0;
            }
            else if (planename=="yz")
            {
                r.x() = 0.0;
            }
            else if (planename=="xz")
            {
                r.y() = 0.0;
            }
            else
            {
                FatalError<< "Invalid circle plane "<<planename<<" specified"
                          << "\n  Use 'xy', 'yz', or 'xz'"<<abort(FatalError);
            }

            vector n = r / mag(r);
            point p = center + n*radius;

            cuttableCell pc(mesh, cellI);
            alpha[cellI] = 1.0 - pc.cut( plane(p, n) );

            if (inv == "yes")
            {
                alpha[cellI] = 1.0 - alpha[cellI];
            }

            if (alpha[cellI] > SMALL)
            {
                U[cellI] = Ugas;
                gasCells++;
            }
            else
            {
                U[cellI] = Usolid;
                solidCells++;
            }
        }

    }
    else if (shape == "plane")
    {
        vector point(VOFDict.lookup("point"));
        vector normal(VOFDict.lookup("normal"));
        vector Usolid(VOFDict.lookup("Usolid"));
        vector Ugas(VOFDict.lookup("Ugas"));

        plane p(point, normal);

        forAll(alpha, cellI)
        {
            cuttableCell pc(mesh, cellI);
            alpha[cellI] = max(1.0 - pc.cut( p ), 0.0);

            if (alpha[cellI] > SMALL)
            {
                U[cellI] = Ugas;
                gasCells++;
            }
            else
            {
                U[cellI] = Usolid;
                solidCells++;
            }
        }

    }
    else
    {
        FatalError<< "Invalid shape " << shape << " specified"
                  << abort(FatalError);
    }

	runTime.writeNow();

	Info<< "Set " << gasCells << " gas and mixed cells and "
	    << solidCells << " solid cells" << endl;

    return 0;
}


// ************************************************************************* //
