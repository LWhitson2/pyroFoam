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

\*---------------------------------------------------------------------------*/

#include "cutPlane.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace shapes
{
    defineTypeNameAndDebug(cutPlane, 0);
    addToRunTimeSelectionTable(shape, cutPlane, components);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shapes::cutPlane::cutPlane
(
    const word& name,
    dictionary shapeDict,
    const fvMesh& mesh
)
:
    shape(typeName, name, shapeDict, mesh),
    p_
    (
        coeffDict_.lookup("point"),
        coeffDict_.lookup("normal")
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::shapes::cutPlane::d() const
{
    const volVectorField& cellCenters = mesh_.C();
    
    tmp<scalarField> d = 
        (p_.normal() & (cellCenters.internalField() - p_.refPoint()));
    
    return d;
}

void Foam::shapes::cutPlane::calculate()
{
    //calculate liquid mask        
    forAll(liquidMask_, cellI)
    {        
        cuttableCell pc(mesh_, cellI);
        
        liquidMask_[cellI] = pc.cut( p_ );
    }

    liquidMask_.correctBoundaryConditions();
    
    //calculate vapor mask
    tmp<scalarField> f = Foam::max(dV_ - d(), 0.0);
    
    vaporMask_.internalField() = f * (1.0 - liquidMask_);
    vaporMask_.correctBoundaryConditions();
}




// ************************************************************************* //
