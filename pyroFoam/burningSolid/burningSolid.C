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
    pyroDict_
    (
        IOobject
        (
            "pyrolysisProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
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
    ),
    
    a_burn_
    (
        IOobject
        (
            "a_burn",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("a_burn_",dimArea, 0.0)
    ),
    rhoS_(pyroDict_.lookup("rhoS")),
    m0_(pyroDict_.lookup("m0"))

{
    Foam::Info << "Created burning solid class" << Foam::endl;
    alpha_.oldTime();
    alphaf_.oldTime();
    
    //Clip phi by alphaf
    calcAlphaf();
    phi_ = (fvc::interpolate(U_*thermo_.rho()) & mesh_.Sf()) * alphaf_;
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
    isBurning_ = pos(alpha_ - SMALL)*pos(1-SMALL - alpha_)
                + neg(alpha_ - SMALL)*pos(fvc::surfaceSum(alphaf_) - SMALL);
                
    // Calculate burning face area
    calcBurningArea();
    
    // Calculate burning face area
    m_pyro_.internalField() = a_burn_ * m0_ / mesh_.V();
    m_pyro_.correctBoundaryConditions();

    // Update alpha_
    solve(fvm::ddt(alpha_) == m_pyro_/rhoS_);

    // Update alphaf_
    calcAlphaf();

    // Include alphaf into phi
    phi_ = (fvc::interpolate(U_*thermo_.rho()) & mesh_.Sf()) * alphaf_;
}


// Calculate the gas area fraction on mesh faces
// LIMITATIONS: This is only valid for 1D serial cases
void Foam::burningSolid::calcAlphaf()
{
    Info<< "Calculating alphaf" << endl;

    //Normal definition applied first, applicable away from interfaces
    alphaf_ = fvc::interpolate(alpha_);

    //Then adjust at the interface regions   
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbor = mesh_.neighbour();

    //Loop through internal faces (internal to this processor)
    forAll(owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbor[faceI];

        if (alpha_[nei] > 1.0-SMALL || alpha_[own] > 1.0-SMALL)
        {
            alphaf_[faceI] = 1.0;
        }
        else if (alpha_[nei] * alpha_[own] < SMALL)
        {
            alphaf_[faceI] = 0.0;
        }
    }
}


// Return the heat generation rate (W/m3)
// LIMITATIONS: NOT IMPLEMENTED
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

// Set the A matrix diagonal to 1/dt in fully solid cells
Foam::tmp<Foam::volScalarField> Foam::burningSolid::setDiag() const
{
    tmp<volScalarField> tSD
    (
        new volScalarField
        (
            IOobject
            (
                "tSD",
                mesh_.time().timeName(),
                mesh_
            ),
            thermo_.rho()/mesh_.time().deltaT()
        )
    );
    
    tSD() *= neg(alpha_ - SMALL);
    return tSD;
}

// Calculate the burning area
// LIMITATIONS: This is only valid for 1D serial cases
void Foam::burningSolid::calcBurningArea()
{
    // Code to get the cell area in the (0 1 0) direction

    const pointField& pf = mesh_.points(); //list of all points in the mesh
    vector dir(0,1,0);
    a_burn_ = dimensionedScalar("zero",dimArea,0.0);

    forAll(mesh_.cells(), cellI)
    {
        if (isBurning_[cellI])
        {
            // a cell is just a list of faces
            const labelList& faces = mesh_.cells()[cellI];
            forAll(faces, faceI)
            {
                const face& f = mesh_.faces()[faces[faceI]];
               
                //This will catch either the face in direction dir, or -dir, but
                // for a blockMesh the area will be the same either way
                if ((dir & (f.normal(pf)/f.mag(pf))) > 0.9)
                {
                    a_burn_[cellI] = f.mag(pf);
                    break;
                }
            }
        }
    }

    a_burn_.correctBoundaryConditions();
} 

// Calculate the burn gas velocity
// LIMITATIONS: This is only valid for 1D serial cases
Foam::tmp<Foam::volVectorField> Foam::burningSolid::burnU() const
{
    tmp<volVectorField> tBurnU
    (
        new volVectorField
        (
            IOobject
            (
                "tBurnU",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector("tBurnU", dimVelocity, vector::zero)
        )
    );
    
    //TODO: Locate "HMXGas" in species list and get its rho here
    tBurnU() = m0_ / thermo_.rho() * vector(0,1,0) * isBurning_;
    return tBurnU;
}


// ************************************************************************* //
