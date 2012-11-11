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

const dimensionedScalar Ru_ = constant::physicoChemical::R*1000;

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

    a_burn_
    (
        IOobject
        (
            "a_burn",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("a_burn",dimArea, 0.0)
    ),

    burnU_
    (
        IOobject
        (
            "burnU",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("burnU", dimVelocity, vector::zero)
    ),
        
    mU_
    (
        IOobject
        (
            "mU",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("mU", dimDensity*dimVelocity/dimTime, vector::zero)
    ),

    USp_
    (
        IOobject
        (
            "USp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("USp", dimDensity/dimTime, 0.0)
    ),

    USu_
    (
        IOobject
        (
            "USu",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("USu", dimDensity*dimVelocity/dimTime, vector::zero)
    ),

    pSp_
    (
        IOobject
        (
            "pSp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("pSp", dimDensity/dimTime/dimPressure, 0.0)
    ),

    pSu_
    (
        IOobject
        (
            "pSu",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("pSu", dimDensity/dimTime, 0.0)
    ),

    iNormal_
    (
        IOobject
        (
            "iNormal",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("iNormal", dimless, vector::zero)
    ),

    rhoS_(pyroDict_.lookup("rhoS")),
    m0_(pyroDict_.lookup("m0")),
    alphaMin_(pyroDict_.lookup("alphaMin")),
    reconstructTol_(pyroDict_.lookup("reconstructTol"))

{
    Foam::Info << "Created burning solid class" << Foam::endl;
    alpha_.oldTime();

    // Calculate new interface position and update alphaf and a_burn
    calculateNewInterface();

    alphaf_.oldTime();
    iNormal_.oldTime();
    a_burn_.oldTime();

    // Calculate mass flux field
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

    // Step 1: Calculate mass flux (m0) using current P and T
    // *** Assumed constant for now ***

    // Step 2: Calculate mass density source using original cell area
    m_pyro_.internalField() = m0_ * a_burn_.oldTime() / mesh_.V();
    m_pyro_.correctBoundaryConditions();

    // Step 3: Evolve interface using calculated burning rate (alpha_)
    solve(fvm::ddt(alpha_) == m_pyro_/rhoS_);

    // Step 4: Calculate burn gas velocity using original interface orientation
    calcBurnU();

    // Step 5: Calculate momentum source
    mU_ = burnU_ * m_pyro_;

    // Step 6: Calculate new interface position and update alphaf and a_burn
    calculateNewInterface();
    
    // Step 7: Fix small cells by transferring momentum (mU) and mass (m_pyro)
    //         to neighbouring larger cells
    fixSmallCells();

    // Step 8: Calculate mass flux field that includes alphaf
    phi_ = (fvc::interpolate(U_*thermo_.rho()) & mesh_.Sf()) * alphaf_;
}


void Foam::burningSolid::calculateInterfaceNormal
(
    const volScalarField& intermeds
)
{
    dimensionedScalar vsrL("vsrL",dimless/dimLength,VSMALL);
    iNormal_ = fvc::grad(alpha_)/(mag(fvc::grad(alpha_))+vsrL)*intermeds;
        //TODO: Use a smoothing procedure to capture the interface better
}


// 3D interface evolution and geometry calculations (replaces calcAlphaf)
//  Update iNormal_, a_burn_, and alphaf_
// LIMITATIONS: SERIAL RUNS ONLY
void Foam::burningSolid::calculateNewInterface()
{
    // Identify intermediate cells based on reconstructTol_
    volScalarField intermeds = pos(alpha_ - reconstructTol_)
                               *pos(1.0-reconstructTol_ - alpha_);

    // Calculate interface normal in intermediate cells
    calculateInterfaceNormal(intermeds);

    // Calculate the cut plane and cut area in intermediate cells
    volVectorField iPoint
    (
        IOobject
        (
            "iPoint",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector("iPoint", dimless, vector::zero)
    );
    
    forAll(iNormal_, cellI)
    {
        if (intermeds[cellI] > SMALL)
        {
            cuttableCell cc(mesh_, cellI);
            plane p = cc.constructInterface(iNormal_[cellI],1.0-alpha_[cellI]);
            iPoint[cellI] = p.refPoint();
            a_burn_[cellI] = cc.cutArea();
            
            // DEBUGGING ONLY, REMOVE ONCE TESTED
            if (a_burn_[cellI] < SMALL)
            {
                Info<< "WARNING: Cut area is zero" 
                    << " with alpha = " << 1.0-alpha_[cellI] << endl;
                iNormal_[cellI] = vector::zero;
                intermeds[cellI] = 0;
            }
        }
        else
        {
            a_burn_[cellI] = 0.0;
        }
    }


    // Calculate alphaf on all faces
    alphaf_ = fvc::interpolate(alpha_); //valid in all homogeneous regions

    // Correct alphaf near interface
    surfaceScalarField hasIntermeds = fvc::interpolate(intermeds);
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbor = mesh_.neighbour();

    // This loop should be modified for parallel cases. This should be doable
    forAll(alphaf_, faceI)
    {
        label own = owner[faceI];
        label nei = neighbor[faceI];

        // Faces where there is a cut plane in one or both bounding cells
        if (hasIntermeds[faceI] > SMALL)
        {
            cuttableFace cf(mesh_, faceI);

            scalar alphafNei = 1.0;
            scalar alphafOwn = 1.0;

            if (mag(iNormal_[nei]) > SMALL)
            {
                Foam::plane p(iPoint[nei], iNormal_[nei]);
                alphafNei = cf.cut(p);
            }

            if (mag(iNormal_[own]) > SMALL)
            {
                Foam::plane p(iPoint[own], iNormal_[own]);
                alphafOwn = cf.cut(p);
            }
            alphaf_[faceI] = Foam::min(alphafNei, alphafOwn);
        }

        // Catch sharp face-coincident interfaces. In this case, the solid cell
        // is the one that will be burning at the next time step.
        else if (mag(alpha_[own] - alpha_[nei]) > 0.1)
        {
            alphaf_[faceI] = 0.0;

            //set iNormal and a_burn for this case
            label solidcell = (alpha_[own] < 0.1) ? own : nei;

            // Increment a_burn since a solidcell could technically have more
            // than one sharp boundary
            a_burn_[solidcell] += mesh_.magSf()[faceI];
            
            // Calculate face normal
            vector norm = mesh_.Sf()[faceI] / mesh_.magSf()[faceI];
            
            // Flip norm if pointed wrong way
            vector vSF = mesh_.Cf()[faceI] - mesh_.C()[solidcell];
            
            if ((vSF & norm) < 0.0)
            {
                norm *= -1.0;
            }
            
            // Increment iNormal_ for the same reason a_burn is incremented
            iNormal_[solidcell] += norm;
        }
    }
    
    a_burn_.correctBoundaryConditions();
    
    //Re-normalize iNormal (only needed for the cases when it is incremented)
    iNormal_ /= (mag(iNormal_) + VSMALL);
}


// Return the heat generation rate (W/m3)
// LIMITATIONS: NOT IMPLEMENTED
tmp<volScalarField> Foam::burningSolid::Sh() const
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

// Calculate alpha that excludes small cells
Foam::tmp<Foam::volScalarField> Foam::burningSolid::alphaUsed() const
{
    return alpha_ * pos(alpha_ - alphaMin_);
}


// Calculate the burn gas velocity
void Foam::burningSolid::calcBurnU()
{
    //TODO: Locate "HMXGas" in species list and get its rho here
    burnU_ = m0_ / thermo_.rho() * iNormal_.oldTime();
}


// Transfer mass and momentum from small cells to larger neighbour cells and
// fix the values in the small cells to:
//   U = Uburn
//   p = weighted average of p in neighbouring gas cells

void Foam::burningSolid::fixSmallCells()
{
    scalarField m_transferred = m_pyro_*mesh_.V()*(1.0 - thermo_.rho()/rhoS_);

    // Value to force small cells to designated velocity
    dimensionedScalar rhordT
    (
        "rhordT",
        dimDensity/dimTime,
        1.0/mesh_.time().deltaTValue()
    );

    // Value to force small cells to designated pressure
    dimensionedScalar psirdT
    (
        "psirdT",
        dimDensity/dimTime/dimPressure,
        1e-5/mesh_.time().deltaTValue()
    );

    // If negative, alpha is a small cell
    volScalarField alphaShift = alpha_ - alphaMin_;

    // Intra-cell transfer weights
    surfaceScalarField w
    (
        IOobject
        (
            "w",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("w", dimless, 0.0)
    );

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbor = mesh_.neighbour();

    // Calculate transfer weights
    Info<< "Calculating transfer weights" << endl;
    forAll(w, faceI)
    {
        label own = owner[faceI];
        label nei = neighbor[faceI];

        if (alphaShift[own] * alphaShift[nei] * alphaf_[faceI] < 0.0)
        { //one is small, one is not, and they share a gas boundary

            label sc = (alphaShift[own] < 0.0) ? own : nei; // Small Cell

            w[faceI] = mag
            (
                burnU_[sc] & mesh_.Sf()[faceI]
            ) * alphaf_[faceI];
        }
    }

    volScalarField wtot = fvc::surfaceSum(w);

    // Sets diagonal terms for SMALL and SOLID cells
    USp_ = neg(alphaShift)*rhordT;
    USu_ = dimensionedVector("zero",dimVelocity*dimDensity/dimTime,vector::zero);

    // pSp set in SMALL and SOLID cells
    pSp_ = neg(alphaShift)*psirdT;

    // Only set pSu in SOLID cells, not SMALL cells. SMALL cells are set in the
    // loop below
    pSu_ = neg(alphaShift)*dimensionedScalar("ps",dimPressure,1e5)*psirdT;


    // Do mass transfers and "de-activate" small cells
    Info<< "Doing transfers" << endl;
    forAll(w, faceI)
    {
        if (w[faceI] > 0.0)
        {
            label own = owner[faceI];
            label nei = neighbor[faceI];

            label sc = (alphaShift[own] < 0.0) ? own : nei; //Small Cell
            label rc = (alphaShift[own] < 0.0) ? nei : own; //Receiving Cell

            m_pyro_[rc] += w[faceI]/wtot[sc] * m_transferred[sc] / mesh_.V()[sc];
            mU_[rc] += w[faceI]/wtot[sc] * m_transferred[sc] * burnU_[sc] / mesh_.V()[sc];
            alphaf_[faceI] = 0.0;

            m_pyro_[sc] = 0.0;
            mU_[sc] = vector::zero;
            USu_[sc] = burnU_[sc] * rhordT.value();
            //pSu_[sc] += w[faceI]/wtot[sc] * thermo_.p()[sc] * psirdT.value();
        }
    }
}

// ************************************************************************* //
