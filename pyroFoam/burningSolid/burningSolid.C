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
    const fvMesh& mesh,
    immersedBoundary& ib,
    const hsCombustionThermo& gasThermo
)
:
    mesh_(mesh),
    
    solidThermo_
    (
        basicSolidThermo::New( mesh )
    ),
    
    ib_(ib),
    
    gasThermo_(gasThermo),
    
    pyroDict_
    (
        IOobject
        (
            "pyrolysisProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
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
    
    m0_(pyroDict_.lookup("m0")),
    
    solidName_(pyroDict_.lookup("solidName")),
    
    gasName_(pyroDict_.lookup("gasName")),
    
    EMg_(NULL)
{
    Info<< "Created burningSolid" << endl;
    
    // Rename T in solidThermo to Ts to avoid name conflicts with the T
    //  in gasThermo.
    solidThermo_->T().rename("Ts");
    
    // Find gasName in species
    const multiComponentMixture<gasThermoPhysics>& composition =
        dynamic_cast<const multiComponentMixture<gasThermoPhysics>& >
        ( gasThermo_.composition() );
        
    const PtrList<gasThermoPhysics>& speciesData = composition.speciesData();
    
    forAll(speciesData, s)
    {
        if (speciesData[s].name() == gasName_)
        {
            Info<< "  Linked burningSolid with gas specie" << endl;
            EMg_ = speciesData(s);
            break;
        }
    }
    
    if (EMg_ == NULL)
    {
        FatalErrorIn
        (
            "burningSolid::burningSolid()"
        )   << "Could not find gas " << gasName_
            << " in loaded gas species"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::burningSolid::fixSmallCells()
{
    volScalarField m_transferred = m_pyro_*
        (1.0 - gasThermo_.rho()/solidThermo_->rho());
        
    volVectorField mU_transferred = m_transferred*burnU_;
    
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

    // Values to set in the solid region
    dimensionedScalar pSolid("ps",dimPressure,1e5);
    dimensionedVector USolid("Us",dimVelocity,vector::zero);
    
    tmp<surfaceScalarField> tw = ib_.scTransferWeights();
    const surfaceScalarField& w = tw();

    // Transfer mass and momentum out of small cells
    ib_.transfer<scalar>(w, m_transferred, m_pyro_, 0.0);
    ib_.transfer<vector>(w, mU_transferred, mU_, vector::zero);
    
    // Set source terms in small and solid cells to specify value
    ib_.setScValue<vector>(w, USu_, USp_, burnU_, USolid, rhordT, "fix");
    ib_.setScValue<scalar>(w, pSu_, pSp_, gasThermo_.p(), 
                           pSolid, psirdT, "avg");
}

// Calculate the burn gas velocity
void Foam::burningSolid::calcBurnU()
{
    //Get the gas constant for the pyrolysis gas
    dimensionedScalar Rg
    (
        "Rg",
        dimVelocity*dimVelocity/dimTemperature,
        EMg_->R()
    );
    
    burnU_ = m0_*gasThermo_.T()*Rg/gasThermo_.p()*ib_.normal().oldTime();
}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

void Foam::burningSolid::correct
(
    const volVectorField& U,
    surfaceScalarField& phi
)
{
    Foam::Info << "Correcting burningSolid" << Foam::endl;

    // Step 1: Calculate mass flux (m0) using current P and T
    // *** Assumed constant for now ***

    // Step 2: Calculate mass density source using original cell area
    m_pyro_.internalField() = m0_ * ib_.area().oldTime() / mesh_.V();
    m_pyro_.correctBoundaryConditions();
    
    // Step 3: Calculate burn gas velocity using current interface orientation
    calcBurnU();

    // Step 4: Calculate momentum source
    mU_ = burnU_ * m_pyro_;
    
    // Step 5: Evolve interface using calculated burning rate (vol frac/s)
    ib_.moveInterface( m_pyro_ / solidThermo_->rho() );

    // Step 6: Fix small cells by transferring momentum (mU) and mass (m_pyro)
    //         to neighbouring larger cells
    fixSmallCells();

    // Step 7: Calculate mass flux field that includes alphaf
    phi = (ib_.interpolate(U*gasThermo_.rho()) & mesh_.Sf())*ib_.alphaf();
}


tmp<volScalarField> Foam::burningSolid::YSu(const word& Yname) const
{
    dimensionedScalar onerDt
    (
        "onerDt",
        dimless/dimTime,
        1.0/mesh_.time().deltaTValue()
    );
    
    if (Yname == gasName_)
    {
        return ib_.gasCells()*m_pyro_
               - ib_.smallCells()*onerDt*gasThermo_.rho();
    }
    else
    {
        return ib_.noCells()*onerDt*gasThermo_.rho();
    }
}


tmp<volScalarField> Foam::burningSolid::YSp() const
{
    return ib_.smallAndSolidCells()*gasThermo_.rho()
            *dimensionedScalar
            (
                "onerDt",
                dimless/dimTime,
                1.0/mesh_.time().deltaTValue()
            );
}



// ************************************************************************* //
