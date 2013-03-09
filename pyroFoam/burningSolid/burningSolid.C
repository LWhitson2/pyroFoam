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

    mCM_
    (
        dynamic_cast<const multiComponentMixture<gasThermoPhysics>& >
        ( gasThermo_.composition() )
    ),

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

    pyroModel_(pyrolysisModel::New(pyroDict_)),

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

    Ti_
    (
        IOobject
        (
            "Ti",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Ti", dimTemperature, 298.0)
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
        mesh_
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

    hsSp_
    (
        IOobject
        (
            "hsSp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("hsSp", dimDensity/dimTime, 0.0)
    ),

    hsSu_
    (
        IOobject
        (
            "hsSu",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("hsSu", dimDensity*dimPower/dimMass, 0.0)
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
            "TsSu",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("TsSu", dimPower/dimVolume, 0.0)
    ),

    m0_(pyroDict_.lookup("m0")),

    solidName_(pyroDict_.lookup("solidName")),

    gasName_(pyroDict_.lookup("gasName")),

    EMg_(NULL),

    Qt_g_
    (
        IOobject
        (
            "Qt_g_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qtg", dimPower/dimVolume, 0.0)
    ),

    Qt_s_
    (
        IOobject
        (
            "Qt_s_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qts", dimPower/dimVolume, 0.0)
    ),

    surfStress_
    (
        IOobject
        (
            "surfStress_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
//         dimensionedVector("surfStress", dimless/(dimLength*dimTime), vector::zero)
        dimensionedScalar("surfStress", dimless/dimArea, 0.0)
    )
{
    Info<< "Created burningSolid" << endl;

    // Rename T in solidThermo to Ts to avoid name conflicts with the T
    //  in gasThermo.
    solidThermo_->T().rename("Ts");
    solidThermo_->rho().rename("rhos");

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

    volScalarField Qtg_transferred = Qt_g_;
    volScalarField Qts_transferred = Qt_s_;

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

    // Value to force small cells to designated enthalpy
    dimensionedScalar hsrdT
    (
        "hsrdT",
        dimDensity/dimTime,
        1.0/mesh_.time().deltaTValue()
    );

    // Value to force gas cells to designated temperature
    dimensionedScalar TsrdT
    (
        "TsrdT",
        dimPower/dimVolume/dimTemperature,
        1.0/mesh_.time().deltaTValue()
    );

    // Values to set in the solid region
    dimensionedScalar pSolid("ps",dimPressure,1e5);
    dimensionedVector USolid("Us",dimVelocity,vector::zero);
    dimensionedScalar hsSolid("hsSolid",dimEnergy/dimMass,0.0);
    dimensionedScalar TsGas("TsGas",dimTemperature,298.0);

    // Calculate transfer weights
    tmp<surfaceScalarField> tw = ib_.scTransferWeights("gas");
    tmp<surfaceScalarField> tws = ib_.scTransferWeights("solid");
    const surfaceScalarField& w = tw();
    const surfaceScalarField& ws = tws();
    Info << "Min/Max Gas Transfer Weights: " << min(w).value() << ", " << max(w).value() << endl;
    Info << "Min/Max Solid Transfer Weights: " << min(ws).value() << ", " << max(ws).value() << endl;

    // Transfer mass, momentum, and energy out of small cells
    ib_.transfer<scalar>(w, m_transferred, m_pyro_, 0.0, "gas");
    ib_.transfer<vector>(w, mU_transferred, mU_, vector::zero, "gas");
    ib_.transfer<scalar>(w, Qtg_transferred, Qt_g_, 0.0, "gas");
    ib_.transfer<scalar>(ws, Qts_transferred, Qt_s_, 0.0, "solid");

    // Set source terms in small and solid cells to specify value
    ib_.setScValue<vector>(w, USu_, USp_, burnU_,
                           USolid, rhordT, "fix", "gas");
    ib_.setScValue<scalar>(w, pSu_, pSp_, gasThermo_.p(),
                           pSolid, psirdT, "avg", "gas");
    ib_.setScValue<scalar>(w, hsSu_, hsSp_, gasThermo_.hs(),
                           hsSolid, hsrdT, "avg", "gas");
    ib_.setScValue<scalar>(ws, TsSu_, TsSp_, Ts_,
                           TsGas, TsrdT, "avg", "solid");
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
    surfaceScalarField& phi,
    bool allowHtx
)
{
    Foam::Info << "Correcting burningSolid" << Foam::endl;

    // Step 1: Calculate mass flux (m0) using current P and Ts
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

    // Step 4b: Calculate heat transfer sources
    if (allowHtx)
    {
        calcHeatTransfer();
    }

    // Step 6: Fix small cells by transferring momentum (mU) and mass (m_pyro)
    //         to neighbouring larger cells
    fixSmallCells();

    // Step 7: Calculate mass flux field that includes alphaf
    phi = (ib_.interpolate(U*gasThermo_.rho()) & mesh_.Sf())*ib_.alphaf();

    // Step 8: Calculate surface stress in interface cells
    calcSurfaceStress();
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

void Foam::burningSolid::calcHeatTransfer()
{
    // Conduction coefficients
    volScalarField Ks = solidThermo_->K();
    volScalarField Kg = gasThermo_.alpha()*gasThermo_.Cp();

    // Conduction lengths
    const volScalarField& Lg = ib_.gasL();//.oldTime();
    const volScalarField& Ls = ib_.solidL();//.oldTime();

    // Interface area and volume
    const volScalarField& Ai = ib_.area();//.oldTime();
    const volScalarField::DimensionedInternalField& Vc = mesh_.V();

    // Initialize cells to no heat transfer
    Qt_g_ = dimensionedScalar("zero", dimPower/dimVolume, 0.0);
    Qt_s_ = dimensionedScalar("zero", dimPower/dimVolume, 0.0);

    // Cell identification
    volScalarField normalCell = ib_.mixedCells();
    volScalarField solidCell = ib_.solidCells();
    volScalarField fullCell = solidCell*pos(Ai
                            - dimensionedScalar("tmp", dimArea, SMALL));

    forAll(Ts_,cellI)
    {
        // Normal mixed cell conduction transfer
        if (normalCell[cellI])
        {
            // Calculate Thermal Resistances
            scalar Rg = Lg[cellI]/(Kg[cellI]*Ai[cellI]);
            scalar Rs = Ls[cellI]/(Ks[cellI]*Ai[cellI]);
            scalar Req = Rg + Rs;

            // Calculate transfer to gas from constant temperature solid
            Qt_g_[cellI] = (Ts_[cellI] - gasThermo_.T()[cellI])/(Req*Vc[cellI]);
            Qt_s_[cellI] = -Qt_g_[cellI];

        }
    }

    // Full solid cell to face neighbor conduction transfer
    forAll(mesh_.magSf(), faceI)
    {
        label own = mesh_.owner()[faceI];
        label nei = mesh_.neighbour()[faceI];
        label sc = (solidCell[own]) ? own:nei;
        label mc = (sc == own) ? nei:own;

        if (fullCell[sc] && !solidCell[mc])
        {
            scalar tmpA = ib_.alphafs()[faceI]
                                    * mesh_.magSf()[faceI];
            scalar tmpL = mag((mesh_.Cf()[faceI]
                        - mesh_.C()[mc]) & mesh_.Sf()[faceI])
                        / mesh_.magSf()[faceI];
            Qt_g_[mc] += Kg[mc]*(Ts_[sc] - gasThermo_.T()[mc])
                        * tmpA/(tmpL*Vc[mc]);
            Qt_s_[sc] -= Qt_g_[mc]*Vc[mc]/Vc[sc];
        }
    }

    // Full solid cell to parallel neighbor conduction transfer
    const fvPatchList& patches = mesh_.boundary();

    // Get boundary fields for required values
    const volScalarField::GeometricBoundaryField& KgBf = Kg.boundaryField();
    const volScalarField::GeometricBoundaryField& fullCellBf = fullCell.boundaryField();
    const volScalarField::GeometricBoundaryField& solidCellBf = solidCell.boundaryField();
    const volVectorField::GeometricBoundaryField& meshCBf = mesh_.C().boundaryField();
    const volScalarField::GeometricBoundaryField& TBf = gasThermo_.T().boundaryField();
    const volScalarField::GeometricBoundaryField& TsBf = Ts_.boundaryField();

    forAll(patches, patchI)
    {
        const fvPatch& curPatch = patches[patchI];
        const labelList& pFaceCells = patches[patchI].faceCells();

        // Get values for current patch
        const fvPatchScalarField& KgPf = KgBf[patchI];
        const fvPatchScalarField& fullCellPf = fullCellBf[patchI];
        const fvPatchScalarField& solidCellPf = solidCellBf[patchI];
        const fvPatchVectorField& meshCPf = meshCBf[patchI];
        const fvPatchScalarField& TPf = TBf[patchI];
        const fvPatchScalarField& TsPf = TsBf[patchI];

        if (curPatch.coupled()) //returns true for parallel and cyclic patches
        {
            // Get values across parallel patch
            const scalarField KgPNf(KgPf.patchNeighbourField());
            const scalarField fullCellPNf(fullCellPf.patchNeighbourField());
            const scalarField solidCellPNf(solidCellPf.patchNeighbourField());
            const vectorField meshCPNf(meshCPf.patchNeighbourField());
            const scalarField TPNf(TPf.patchNeighbourField());
            const scalarField TsPNf(TsPf.patchNeighbourField());

            forAll(curPatch, pFaceI)
            {
                label pfCellI = pFaceCells[pFaceI];

                scalar tmpA = ib_.alphafs()[pFaceI]
                                    * mesh_.magSf()[pFaceI];

                // Boundary cell solid, neighbor cell mixed
                if (fullCell[pfCellI] && !solidCellPNf[pFaceI])
                {

                    //Info << "Parallel Solid Cell" << endl;
                    scalar tmpL = mag((mesh_.Cf()[pFaceI]
                                - meshCPNf[pFaceI]) & mesh_.Sf()[pFaceI])
                                / mesh_.magSf()[pFaceI];


                    Qt_s_[pfCellI] -= KgPNf[pFaceI]*(Ts_[pfCellI]
                                    - TPNf[pFaceI])
                                    * tmpA/(tmpL*Vc[pfCellI]);
                }
                // Boundary cell mixed, neighbor cell solid
                else if (fullCellPNf[pFaceI] && !solidCell[pfCellI])
                {
                    //Info << "Parallel Mixed Cell" << endl;
                    scalar tmpL = mag((mesh_.Cf()[pFaceI]
                                - mesh_.C()[pfCellI]) & mesh_.Sf()[pFaceI])
                                / mesh_.magSf()[pFaceI];

                    Qt_g_[pfCellI] += Kg[pfCellI]*(TsPNf[pFaceI]
                                    - gasThermo_.T()[pfCellI])
                                    * tmpA/(tmpL*Vc[pfCellI]);
                }
            }
        }
    }
}

void Foam::burningSolid::calcSurfaceStress()
{
    // Centroid lengths
    const volScalarField& Lg = ib_.gasL();

    // Interface area and volume
    const volScalarField& Ai = ib_.area();
    const volScalarField::DimensionedInternalField& Vc = mesh_.V();

    // Initialize cells to no heat transfer
    surfStress_ = dimensionedScalar("tmpSurfStress", dimless/dimArea, 0.0);

    // Cell identification
    volScalarField normalCell = ib_.mixedCells();
    volScalarField solidCell = ib_.solidCells();
    volScalarField fullCell = solidCell*pos(Ai
                            - dimensionedScalar("tmp", dimArea, SMALL));

    
    forAll(Ts_,cellI)
    {
        // Normal mixed cell surface stress
        if (normalCell[cellI])
        {
            // Calculate transfer to gas from constant temperature solid
            surfStress_[cellI] = Ai[cellI]/(Lg[cellI]*Vc[cellI]);
        }
    }
}

// ************************************************************************* //
