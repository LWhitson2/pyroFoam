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

#include "shape.H"

namespace Foam
{
    defineTypeNameAndDebug(shape, 0);
    defineRunTimeSelectionTable(shape, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shape::shape
(
    const word& type,
    const word& name,
    dictionary shapeDict,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    name_(name),
    shapeDict_(shapeDict),
    coeffDict_(shapeDict.subDict(type+"Coeffs")),
    dV_(readScalar(shapeDict.lookup("delVapor"))),
    Uinit_(shapeDict.lookup("U")),
    liquidSpecies_(shapeDict.lookup("liquidSpecies")),
    vaporSpecies_(shapeDict.lookup("vaporSpecies")),
    liquidMask_
    (
        IOobject
        (
            "liquid_mask_" + name_,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("liquidMask",dimless,0),
        zeroGradientFvPatchScalarField::typeName
    ),
    vaporMask_
    (
        IOobject
        (
            "vapor_mask_" + name_,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("vaporMask",dimless,0),
        zeroGradientFvPatchScalarField::typeName
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::shape> Foam::shape::clone() const
{
    notImplemented("shape::clone() const");
    return autoPtr<shape>(NULL);
}

Foam::List<Foam::word> Foam::shape::species() const
{
    List<word> species = liquidSpecies_.toc();
    species.append(vaporSpecies_.toc());
    return species;
}

void Foam::shape::set
(
    Foam::volScalarField& alphaLiquid,
    Foam::volVectorField& U,
    PtrList<volScalarField>& species
)
{
    calculate();
    
    alphaLiquid.internalField() += liquidMask_;
    U.internalField() += Uinit_*pos(liquidMask_+vaporMask_ - SMALL);

    forAll(species, i)
    {
        if( liquidSpecies_.found(species[i].name()) )
        {
            species[i].internalField() += liquidMask_*liquidSpecies_[species[i].name()];
        }
        
        if( vaporSpecies_.found(species[i].name()) )
        {
            species[i].internalField() += vaporMask_*vaporSpecies_[species[i].name()];
        }
    }
}



// ************************************************************************* //
