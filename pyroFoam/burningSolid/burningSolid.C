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

#include "burningSolid.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burningSolid::burningSolid
(
    const volVectorField& U,
    surfaceScalarField& phi,
    const hsCombustionThermo& thermo
)
:
    mesh_(U.mesh()),
    U_(U),
    phi_(phi),
    thermo_(thermo),

    alpha_
    (
        IOobject
        (
            "alphaGas",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    alphaf_
    (
        IOobject
        (
            "alphaf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphaf", dimless, 0.0)
    ),

    m_pyro_
    (
        IOobject
        (
            "m_pyro",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("m_pyro", dimDensity/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    
    isBurning_
    (
        IOobject
        (
            "isBurning",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("isBurning", dimless, 0.0)
    )
{
    Foam::Info << "Created burning solid class" << Foam::endl;
    alpha_.oldTime();

    //Clip phi by alphaf
    calcAlphaf();
    phi_ *= alphaf_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::burningSolid> Foam::burningSolid::clone() const
{
    notImplemented("burningSolid::clone() const");
    return autoPtr<burningSolid>(NULL);
}


// correct evolves the solid surface, calculates m_pyro_, and updates alpha and
// alphaf. This also clips phi based on gas surface area so it is correct mass
// flux for the rest of the solvers
void Foam::burningSolid::correct()
{
    Foam::Info << "Correcting solid interface" << Foam::endl;
    
    // Determine the burning faces
    // First term 1 for intermediates, 0 otherwise
    // Second term 1 for alpha == 0 and sum(alphaf) /= 0, 0 otherwise
    isBurning_ = pos(alpha_-SMALL)*pos(1-SMALL-alpha_)
                + neg(alpha_-SMALL)*pos(fvc::surfaceSum(alphaf_)-SMALL);
    
    // Calculate m_pyro_ using A = mag(fvc::grad(alpha_))
    m_pyro_ = isBurning_*dimensionedScalar("a",dimDensity/dimTime,1.0);

    // Update alpha_
    solve(fvm::ddt(alpha_));
/*
    fvVectorMatrix alphaEqn
    (
        fvm::ddt(alpha_)
      + fvm::Sp(fvc::ddt(thermo_.p())/thermo_.p(), alpha_)
      + fvc::div(phi_)/thermo_.rho()
     ==
        m_pyro_/thermo_.rho()
    );

    alphaEqn.solve();
*/

    // Update alphaf_
    calcAlphaf();

    // Include alphaf into phi
    phi_ *= alphaf_;

}


// Calculate the gas area fraction on mesh faces
void Foam::burningSolid::calcAlphaf()
{
    alphaf_ = fvc::interpolate(alpha_);
}


// Return the heat generation rate (W/m3)
Foam::tmp<Foam::volScalarField> Foam::burningSolid::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "tSh",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("tSh", dimPower/dimVolume, 0.0)
        )
    );

    return tSh;
}


// ************************************************************************* //
