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

    isBurning_
    (
        IOobject
        (
            "isBurning",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
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

    TsSp_
    (
        IOobject
        (
            "TsSp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("TsSp", dimPower/dimVolume/dimTemperature, 0.0)
    ),

    TsSu_
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
        dimensionedScalar("TsSu", dimPower/dimVolume, 0.0)
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

    Ts_
    (
        IOobject
        (
            "Ts",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Ts",dimTemperature, 300.0)
    ),

    rhoS_(pyroDict_.lookup("rhoS")),
    m0_(pyroDict_.lookup("m0")),
    alphaMin_(pyroDict_.lookup("alphaMin")),
    Ac_(pyroDict_.lookup("Ac")),
    Ec_(pyroDict_.lookup("Ec")),
    Qc_(pyroDict_.lookup("Qc")),
    kc_(pyroDict_.lookup("kc")),
    Cpc_(pyroDict_.lookup("Cpc"))

{
    Foam::Info << "Created burning solid class" << Foam::endl;
    alpha_.oldTime();

    // Calculate new interface position and update alphaf and a_burn
    calcAlphaf();

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
    calcAlphaf();

    // Step 7: Fix small cells by transferring momentum (mU) and mass (m_pyro)
    //         to neighbouring larger cells
    fixSmallCells();

    // Step 8: Calculate mass flux field that includes alphaf
    phi_ = (fvc::interpolate(U_*thermo_.rho()) & mesh_.Sf()) * alphaf_;
}


// 3D interface evolution and geometry calculations (replaces calcAlphaf)
//  Update iNormal_, a_burn_, and alphaf_
void Foam::burningSolid::calculateNewInterface()
{
    // Identify intermediate cells
    volScalarField intermeds = pos(alpha_ - SMALL)*pos(1.0-SMALL - alpha_);

    // Calculate interface normal in intermediate cells
    iNormal_ = fvc::grad(alpha_)/mag(fvc::grad(alpha_)+VSMALL) * intermeds;
        //TODO: Use a smoothing procedure to capture the interface better

    // Set a_burn equal to zero everywhere
    a_burn_ = dimensionedScalar("zero",dimArea,0.0);

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
        /*
            cuttableCell cc(mesh_, cellI);
            plane p = cc.constructInterface(iNormal_[cellI], alpha_[cellI]);
            iPoint[cellI] = p.refPoint();
            a_burn_[cellI] = cc.cutArea();
        */
        }
    }


    // Calculate alphaf on all faces
    alphaf_ = fvc::interpolate(alpha_); //valid in all homogeneous regions

    // Correct alphaf near interface
    surfaceScalarField avgNorm = fvc::interpolate(mag(iNormal_));

    forAll(alphaf_, faceI)
    {
    /*
        // Faces where there is a cut plane in one or both neighbour cells
        if (avgNorm[faceI] > SMALL)
        {

            cuttableFace cf(mesh_, faceI);

            scalar alphafNei = 1.0;
            scalar alphafOwn = 1.0;

            if (mag(iNormal_[nei]) > SMALL)
            {
                alphafNei = cf.cutArea(iPoint_[nei], iNormal_[nei]);
            }

            if (mag(iNormal_[own]) > SMALL)
            {
                alphafOwn = cf.cutArea(iPoint_[own], iNormal_[own]);
            }

            alphaf_[face] = min(alphafNei, alphafOwn);

        }

        // catch sharp face-coincident boundaries
        else if (alpha_[own]*alpha_[nei] < SMALL && mag(alpha_[own]+alpha_[nei]-1)<SMALL)
        {

            alphaf_[faceI] = 0.0;

            //set iNormal and a_burn for this case


        }
    */
    }


}


// Calculate the gas area fraction on mesh faces
// LIMITATIONS: This is only valid for 1D serial cases
void Foam::burningSolid::calcAlphaf()
{
    calcBurningArea();

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

    isBurning_ = pos(alpha_ - SMALL)*pos(1-SMALL - alpha_)
                + neg(alpha_ - SMALL)*pos(fvc::surfaceSum(alphaf_) - SMALL);
    iNormal_ = dimensionedVector("n",dimless,vector(0,1,0)) * isBurning_;
}


// Return the heat generation rate (W/m3)
// LIMITATIONS: NOT IMPLEMENTED
volScalarField Foam::burningSolid::Sh()
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

    tSh = Qc_*rhoS_*Ac_*exp(-Ec_/(Ts_*Ru_));

    return tSh;
}

// Calculate alpha that excludes small cells
Foam::tmp<Foam::volScalarField> Foam::burningSolid::alphaUsed() const
{
    return alpha_ * pos(alpha_ - alphaMin_);
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
    //scalarField m_generated = m_pyro_*mesh_.V(); //a_burn_ * m0_;     // (kg_gas/s)
    //scalarField m_stored = fvc::ddt(alpha_)*thermo_.rho()*mesh_.V();  // (kg_gas/s)

    // ddt(alpha_) = m_pyro_/rhoS, which can be substituted into the expression above

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

    // Momentum being generated in current cell
    mU_ = m_pyro_ * burnU_;

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

void Foam::burningSolid::solveTs()
{
    // Copy gas temperature field for full gas cells
    Ts_ = Ts_ + (thermo_.T() - Ts_)*pos(alpha_.value() - 1. + SMALL);

    // Calculate explicit source term
    TsSu_ = (1. - alpha_)*Sh();

    // Calculate implicit source term

    fvScalarMatrix TsEqn
    (
        fvm::ddt((1.-alpha_)*rhoS_*Cpc_, Ts_)
     ==
        fvm::laplacian(kc_, Ts_)
      + TsSu_ + fvm::Sp(TsSp_,Ts_)
    );

    TsEqn.relax();
    TsEqn.solve();

//    thermo.correct(); //we may need to create a 'specie' for the solid so that this works

    Info<< "T solid min/max   = " << min(Ts_).value() << ", "
        << max(Ts_).value() << endl;
}

// ************************************************************************* //
