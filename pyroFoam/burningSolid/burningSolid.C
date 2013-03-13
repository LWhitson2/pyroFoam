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

    QgSp_
    (
        IOobject
        (
            "QgSp_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("QgSp", dimDensity/dimTime, 0.0)
    ),

    QsSp_
    (
        IOobject
        (
            "QsSp_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("QsSp", dimPower/dimTemperature/dimVolume, 0.0)
    ),

    QgSu_
    (
        IOobject
        (
            "QgSu_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("QgSu", dimPower/dimVolume, 0.0)
    ),

    QsSu_
    (
        IOobject
        (
            "QsSu_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("QsSu", dimPower/dimVolume, 0.0)
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

    // Set source terms in small and solid cells to specify value
    ib_.setScValue<vector>(w, USu_, USp_, burnU_,
                           USolid, rhordT, "fix", "gas");
    ib_.setScValue<scalar>(w, pSu_, pSp_, gasThermo_.p(),
                           pSolid, psirdT, "avg", "gas");
    ib_.setScValue<scalar>(w, hsSu_, hsSp_, gasThermo_.hs(),
                           hsSolid, hsrdT, "avg", "gas");
    ib_.setScValue<scalar>(ws, TsSu_, TsSp_, Ts_,
                           TsGas, TsrdT, "avg", "solid");

    TsSp_ = TsSp_*(1. - ib_.smallSolidCells());
    TsSu_ = TsSu_*(1. - ib_.smallSolidCells());
    hsSp_ = hsSp_*(1. - ib_.smallCells());
    hsSu_ = hsSu_*(1. - ib_.smallCells());

    // Transfer mass and momentum out of small cells
    ib_.transfer<scalar>(w, m_transferred, m_pyro_, 0.0, "gas");
    ib_.transfer<vector>(w, mU_transferred, mU_, vector::zero, "gas");
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
    tmp<volScalarField> Cpg = gasThermo_.Cp();
    volScalarField alphag = gasThermo_.alpha();
    volScalarField Kg = alphag*Cpg();

    // Conduction lengths
    const volScalarField& Lg = ib_.gasL();
    const volScalarField& Ls = ib_.solidL();

    // Interface area and volume
    const volScalarField& Ai = ib_.area();//.oldTime();
    const volScalarField::DimensionedInternalField& Vc = mesh_.V();

    // Initialize cells to no heat transfer
    QgSp_ = dimensionedScalar("zero", dimDensity/dimTime, 0.0);
    QsSp_ = dimensionedScalar("zero", dimPower/dimTemperature/dimVolume, 0.0);
    QgSu_ = dimensionedScalar("zero", dimPower/dimVolume, 0.0);
    QsSu_ = dimensionedScalar("zero", dimPower/dimVolume, 0.0);

    // Initialize Interface Temperature to Solid Temperature
    Ti_ = Ts_;

    // Cell identification
    volScalarField normalCell = ib_.mixedCells();
    volScalarField solidCell = ib_.solidCells();
    volScalarField smallGasCell = ib_.smallCells();
    volScalarField smallSolidCell = ib_.smallSolidCells();
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

            // Calculate Interface Temperature
            Ti_[cellI] = (Rg*Ts_[cellI] + Rs*gasThermo_.T()[cellI])/Req;

            // Gas source terms
            QgSp_[cellI] = 1./(Rg*Cpg()[cellI]*Vc[cellI]);
            QgSu_[cellI] = mCM_.cellMixture(cellI).Hs(Ti_[cellI])
                         / (Rg*Cpg()[cellI]*Vc[cellI]);

            // Solid source terms
            QsSp_[cellI] = 1./(Rs*Vc[cellI]);
            QsSu_[cellI] = Ti_[cellI]/(Rs*Vc[cellI]);
        }
    }

    volScalarField scReq = Ti_*0.;
    forAll(mesh_.magSf(), faceI)
    {
        label own = mesh_.owner()[faceI];
        label nei = mesh_.neighbour()[faceI];
        label sc = (solidCell[own]) ? own:nei;
        label mc = (sc == own) ? nei:own;

        // Full solid cell to face neighbor conduction transfer
        if (fullCell[sc] && !solidCell[mc])
        {
            scalar tmpA = ib_.alphafU()[faceI]*mesh_.magSf()[faceI];

            if (tmpA > 0.)
            {
                scalar tmpLg = mag((mesh_.Cf()[faceI] - ib_.gasC()[mc])
                             & mesh_.Sf()[faceI])/mesh_.magSf()[faceI];
                scalar tmpLs = mag((mesh_.Cf()[faceI] - mesh_.C()[sc])
                             & mesh_.Sf()[faceI])/mesh_.magSf()[faceI];

                // Calculate thermal resistance
                scalar Rg = tmpLg/(Kg[mc]*tmpA);
                scalar Rs = tmpLs/(Ks[sc]*tmpA);
                scalar Req = Rg + Rs;

                // Calculate current interface temperature
                scalar tmpTi = (Rg*Ts_[sc] + Rs*gasThermo_.T()[mc])/Req;

                // Gas source terms
                QgSp_[mc] += 1./(Rg*Cpg()[mc]*Vc[mc]);
                QgSu_[mc] += mCM_.cellMixture(mc).Hs(tmpTi)
                           / (Rg*Cpg()[mc]*Vc[mc]);

                // Solid source terms
                QsSp_[sc] += 1./(Rs*Vc[sc]);
                QsSu_[sc] += tmpTi/(Rs*Vc[sc]);

                // Partial calculation for Ti of solid cell
                scReq[sc] += 1./Req;
            }
        }

        // Set small cell value based on steady state conduction
        if (smallGasCell[own] + smallGasCell[nei])
        {
            scalar tmpA = ib_.alphafU()[faceI]*mesh_.magSf()[faceI];

            if (tmpA > 0.)
            {
                // Calculate conduction lengths
                scalar Lown = mag((mesh_.Cf()[faceI] - mesh_.C()[own])
                            & mesh_.Sf()[faceI])/mesh_.magSf()[faceI];
                scalar Lnei = mag((mesh_.Cf()[faceI] - mesh_.C()[nei])
                            & mesh_.Sf()[faceI])/mesh_.magSf()[faceI];

                // Calculate thermal resistance
                scalar Rown = Lown/(tmpA*alphag[own]);
                scalar Rnei = Lnei/(tmpA*alphag[nei]);
                scalar Req = Rown + Rnei;

                // Add thermal resistance to system
                if (smallGasCell[own])
                {
                    QgSp_[own] += 1./(Req*Vc[own]);
                    QgSu_[own] += gasThermo_.hs()[nei]/(Req*Vc[own]);
                }
                if (smallGasCell[nei])
                {
                    QgSp_[nei] += 1./(Req*Vc[nei]);
                    QgSu_[nei] += gasThermo_.hs()[own]/(Req*Vc[nei]);
                }
            }
        }
        if (smallSolidCell[own] + smallSolidCell[nei])
        {
            scalar tmpA = ib_.alphafsU()[faceI]*mesh_.magSf()[faceI];

            if (tmpA > 0.)
            {
                // Calculate conduction lengths
                scalar Lown = mag((mesh_.Cf()[faceI] - mesh_.C()[own])
                            & mesh_.Sf()[faceI])/mesh_.magSf()[faceI];
                scalar Lnei = mag((mesh_.Cf()[faceI] - mesh_.C()[nei])
                            & mesh_.Sf()[faceI])/mesh_.magSf()[faceI];

                // Calculate thermal resistance
                scalar Rown = Lown/(tmpA*Ks[own]);
                scalar Rnei = Lnei/(tmpA*Ks[nei]);
                scalar Req = Rown + Rnei;

                // Add thermal resistance to system
                if (smallSolidCell[own])
                {
                    QsSp_[own] += 1./(Req*Vc[own]);
                    QsSu_[own] += Ts_[nei]/(Req*Vc[own]);
                }
                if (smallSolidCell[nei])
                {
                    QsSp_[nei] += 1./(Req*Vc[nei]);
                    QsSu_[nei] += Ts_[own]/(Req*Vc[nei]);
                }
            }
        }
    }

    // Full solid cell to parallel neighbor conduction transfer
    const fvPatchList& patches = mesh_.boundary();

    // Get boundary fields for required values
    const volScalarField::GeometricBoundaryField& KgBf = Kg.boundaryField();
    const volScalarField::GeometricBoundaryField& KsBf = Ks.boundaryField();
    const volScalarField::GeometricBoundaryField& alphagBf = alphag.boundaryField();
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
        const fvPatchScalarField& KsPf = KsBf[patchI];
        const fvPatchScalarField& alphagPf = alphagBf[patchI];
        const fvPatchScalarField& fullCellPf = fullCellBf[patchI];
        const fvPatchScalarField& solidCellPf = solidCellBf[patchI];
        const fvPatchVectorField& meshCPf = meshCBf[patchI];
        const fvPatchScalarField& TPf = TBf[patchI];
        const fvPatchScalarField& TsPf = TsBf[patchI];

        if (curPatch.coupled()) //returns true for parallel and cyclic patches
        {
            // Get values across parallel patch
            const scalarField KgPNf(KgPf.patchNeighbourField());
            const scalarField KsPNf(KsPf.patchNeighbourField());
            const scalarField alphagPNf(alphagPf.patchNeighbourField());
            const scalarField fullCellPNf(fullCellPf.patchNeighbourField());
            const scalarField solidCellPNf(solidCellPf.patchNeighbourField());
            const vectorField meshCPNf(meshCPf.patchNeighbourField());
            const scalarField TPNf(TPf.patchNeighbourField());
            const scalarField TsPNf(TsPf.patchNeighbourField());

            forAll(curPatch, pFaceI)
            {
                label pfCellI = pFaceCells[pFaceI];

                // Boundary cell solid, neighbor cell mixed
                if (fullCell[pfCellI] && !solidCellPNf[pFaceI])
                {
                    scalar tmpA = ib_.alphafU()[pFaceI]
                                    * mesh_.magSf()[pFaceI];

                    if (tmpA > 0.)
                    {
                        scalar tmpLg = mag((mesh_.Cf()[pFaceI]
                                    - meshCPNf[pFaceI]) & mesh_.Sf()[pFaceI])
                                    / mesh_.magSf()[pFaceI];
                        scalar tmpLs = mag((mesh_.Cf()[pFaceI]
                                    - mesh_.C()[pfCellI]) & mesh_.Sf()[pFaceI])
                                    / mesh_.magSf()[pFaceI];

                        // Calculate thermal resistance
                        scalar Rg = tmpLg/(tmpA*KgPNf[pFaceI]);
                        scalar Rs = tmpLs/(tmpA*Ks[pfCellI]);
                        scalar Req = Rg + Rs;

                        // Solid source terms
                        QsSp_[pfCellI] += 1./(Req*Vc[pfCellI]);
                        QsSu_[pfCellI] += TPNf[pFaceI]/(Req*Vc[pfCellI]);

                        // Partial calculation for Ti of solid cell
                        scReq[pfCellI] += 1./Req;
                    }
                }

                // Boundary cell mixed, neighbor cell solid
                else if (fullCellPNf[pFaceI] && !solidCell[pfCellI])
                {
                    scalar tmpA = ib_.alphafU()[pFaceI]
                                    * mesh_.magSf()[pFaceI];
                    if (tmpA > 0.)
                    {
                        scalar tmpLg = mag((mesh_.Cf()[pFaceI]
                                     - mesh_.C()[pfCellI]) & mesh_.Sf()[pFaceI])
                                     / mesh_.magSf()[pFaceI];
                        scalar tmpLs = mag((mesh_.Cf()[pFaceI]
                                    - meshCPNf[pFaceI]) & mesh_.Sf()[pFaceI])
                                    / mesh_.magSf()[pFaceI];

                        // Calculate thermal resistance
                        scalar Rg = tmpLg/(Kg[pfCellI]*tmpA);
                        scalar Rs = tmpLs/(KsPNf[pFaceI]*tmpA);
                        scalar Req = Rg + Rs;

                        // Calculate interface temperature for source terms
                        scalar tmpTi = (Rg*TsPNf[pFaceI]
                                     + Rs*gasThermo_.T()[pfCellI])/Req;

                        // Gas source terms
                        QgSp_[pfCellI] += 1./(Rg*Cpg()[pfCellI]*Vc[pfCellI]);
                        QgSu_[pfCellI] += mCM_.cellMixture(pfCellI).Hs(tmpTi)
                                        / (Rg*Cpg()[pfCellI]*Vc[pfCellI]);
                    }
                }

                // Set small cell value based on steady state conduction
                if (smallGasCell[pfCellI])
                {
                    scalar tmpA = ib_.alphafU()[pFaceI] * mesh_.magSf()[pFaceI];
                    if (tmpA > 0.)
                    {
                        scalar Lown = mag((mesh_.Cf()[pFaceI]
                                    - ib_.gasC()[pfCellI]) & mesh_.Sf()[pFaceI])
                                    / mesh_.magSf()[pFaceI];
                        scalar Lnei = mag((mesh_.Cf()[pFaceI] - meshCPNf[pFaceI])
                                    & mesh_.Sf()[pFaceI])/mesh_.magSf()[pFaceI];

                        // Calculate Thermal Resistances
                        scalar Rown = Lown/(tmpA*alphag[pfCellI]);
                        scalar Rnei = Lnei/(tmpA*alphagPNf[pFaceI]);
                        scalar Req = Rown + Rnei;

                        // Calculate interface temperature
                        scalar tmpTi = (Rown*gasThermo_.T()[pfCellI]
                                     + Rnei*TPNf[pFaceI])/Req;

                        // Add thermal resistance to system
                        QgSp_[pfCellI] += 1./(Rown*Vc[pfCellI]);
                        QgSu_[pfCellI] += mCM_.cellMixture(pfCellI).Hs(tmpTi)
                                        / (Rown*Vc[pfCellI]);
                    }
                }
                if (smallSolidCell[pfCellI])
                {
                    scalar tmpA = ib_.alphafsU()[pFaceI] * mesh_.magSf()[pFaceI];
                    if (tmpA > 0.)
                    {
                        scalar Lown = mag((mesh_.Cf()[pFaceI]
                                    - ib_.solidC()[pfCellI]) & mesh_.Sf()[pFaceI])
                                    / mesh_.magSf()[pFaceI];
                        scalar Lnei = mag((mesh_.Cf()[pFaceI] - meshCPNf[pFaceI])
                                    & mesh_.Sf()[pFaceI])/mesh_.magSf()[pFaceI];

                        // Caclulate Thermal Resistance
                        scalar Rown = Lown/(tmpA*Ks[pfCellI]);
                        scalar Rnei = Lnei/(tmpA*KsPNf[pFaceI]);
                        scalar Req = Rown + Rnei;

                        // Add thermal resistance to system
                        QsSp_[pfCellI] += 1./(Req*Vc[pfCellI]);
                        QsSu_[pfCellI] += TsPNf[pFaceI]/(Req*Vc[pfCellI]);
                    }
                }
            }
        }
    }

    dimensionedScalar hsrdT
    (
        "hsrdT",
        dimDensity/dimTime,
        1.0/mesh_.time().deltaTValue()
    );
    dimensionedScalar TsrdT
    (
        "TsrdT",
        dimPower/dimVolume/dimTemperature,
        1.0/mesh_.time().deltaTValue()
    );

    // Solve conduction for small cells only
    volScalarField rhog = gasThermo_.rho();
    volScalarField rhos = solidThermo_->rho();
    volScalarField Cps = solidThermo_->Cp();
    forAll(QgSu_, cellI)
    {
        if (smallGasCell[cellI])
        {
            scalar qcond = (QgSu_[cellI] - QgSp_[cellI]*gasThermo_.hs()[cellI])
                         * mesh_.time().deltaTValue();
            scalar hsGas = QgSu_[cellI]/QgSp_[cellI];
            scalar hsTime = (rhog.oldTime()[cellI]*gasThermo_.hs().oldTime()[cellI]
                          + qcond)/rhog[cellI];
            if (qcond > 0.) hsGas = min(hsGas, hsTime);
            else hsGas = max(hsGas, hsTime);
            hsSp_[cellI] = hsrdT.value();
            hsSu_[cellI] = hsGas*hsrdT.value();
            QsSu_[cellI] += (gasThermo_.hs().oldTime()[cellI] - hsGas)*rhog[cellI]
                          / mesh_.time().deltaTValue();
            QgSu_[cellI] = 0.0;
            QgSp_[cellI] = 0.0;
        }

        if (smallSolidCell[cellI])
        {
            scalar qcond = (QsSu_[cellI] - QsSp_[cellI]*Ts_[cellI])
                         * mesh_.time().deltaTValue();
            scalar Tsolid = QsSu_[cellI]/QsSp_[cellI];
            scalar TsTime = Ts_.oldTime()[cellI]
                          + qcond/(rhos[cellI]*Cps[cellI]);
            if (qcond > 0.) Tsolid = min(Tsolid, TsTime);
            else Tsolid = max(Tsolid, TsTime);
            TsSp_[cellI] = TsrdT.value();
            TsSu_[cellI] = Tsolid*TsrdT.value();
            QgSu_[cellI] += (Ts_.oldTime()[cellI] - Tsolid)
                          * rhos[cellI]*Cps[cellI]/mesh_.time().deltaTValue();
            QsSu_[cellI] = 0.0;
            QsSp_[cellI] = 0.0;
        }

        // Calculate interface temperature in fullCells
//         if (scReq[cellI] > 0.0)
//         {
//             Ti_[cellI] = Ts_[cellI] + (QsSu_[cellI] - QsSp_[cellI]*Ts_[cellI])
//                        / scReq[cellI];
//         }
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

// void Foam::burningSolid::calcSurfaceStress2()
// {
//     // Centroid lengths
//     volVectorField Lg = ib_.gasL()*ib_.iNormal();
//
//     // Interface area and volume
//     volScalarField& Ai = ib_.area();
//     const volScalarField::DimensionedInternalField& Vc = mesh_.V();
//
//     // Initialize cells to no heat transfer
//     surfStress_ = dimensionedScalar("tmpSurfStress", dimless/dimArea, 0.0);
//
//     // Cell identification
//     volScalarField normalCell = ib_.mixedCells();
//     volScalarField solidCell = ib_.solidCells();
//     volScalarField fullCell = solidCell*pos(Ai
//                             - dimensionedScalar("tmp", dimArea, SMALL));
//
//     forAll(Ts_,cellI)
//     {
//         // Normal mixed cell surface stress
//         if (normalCell[cellI])
//         {
//             // Calculate transfer to gas from constant temperature solid
//             surfStress_[cellI] = Ai[cellI]/(Lg[cellI]*Vc[cellI]);
//         }
//     }
// }

void Foam::burningSolid::calcHeatTransferOld()
{
    // Conduction coefficients
    volScalarField Ks = solidThermo_->K();
    tmp<volScalarField> Cpg = gasThermo_.Cp();
    volScalarField Kg = gasThermo_.alpha()*Cpg();

    // Conduction lengths
    const volScalarField& Lg = ib_.gasL();//.oldTime();
    const volScalarField& Ls = ib_.solidL();//.oldTime();

    // Interface area and volume
    const volScalarField& Ai = ib_.area();//.oldTime();
    const volScalarField::DimensionedInternalField& Vc = mesh_.V();

    // Initialize cells to no heat transfer
    QgSp_ = dimensionedScalar("zero", dimDensity/dimTime, 0.0);
    QsSp_ = dimensionedScalar("zero", dimPower/dimTemperature/dimVolume, 0.0);
    QgSu_ = dimensionedScalar("zero", dimPower/dimVolume, 0.0);
    QsSu_ = dimensionedScalar("zero", dimPower/dimVolume, 0.0);

    // Cell identification
    volScalarField normalCell = ib_.mixedCells();
    volScalarField solidCell = ib_.solidCells();
    volScalarField smallGasCell = ib_.smallCells();
    volScalarField smallSolidCell = ib_.smallSolidCells();
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

            // Calculate source terms
            QgSp_[cellI] = 1./(Req*Vc[cellI]);
            QsSp_[cellI] = 1./(Req*Vc[cellI]);
//             QgSu_[cellI] = mCM_.cellMixture(cellI).Hs(Ts_[cellI])
//                          / (Cpg()[cellI]*Req*Vc[cellI]);
            QgSu_[cellI] = Ts_.oldTime()[cellI]/(Req*Vc[cellI]);
            QsSu_[cellI] = gasThermo_.T().oldTime()[cellI]/(Req*Vc[cellI]);
        }
    }


    forAll(mesh_.magSf(), faceI)
    {
        label own = mesh_.owner()[faceI];
        label nei = mesh_.neighbour()[faceI];
        label sc = (solidCell[own]) ? own:nei;
        label mc = (sc == own) ? nei:own;

        // Full solid cell to face neighbor conduction transfer
        if (fullCell[sc] && !solidCell[mc])
        {
            scalar tmpA = ib_.alphafU()[faceI]*mesh_.magSf()[faceI];
            scalar tmpL = mag((mesh_.Cf()[faceI] - ib_.gasC()[mc])
                        & mesh_.Sf()[faceI])/mesh_.magSf()[faceI];

//             Qt_g_[mc] += Kg[mc]*(Ts_[sc] - gasThermo_.T()[mc])
//                         * tmpA/(tmpL*Vc[mc]);
//             Qt_s_[sc] -= Qt_g_[mc]*Vc[mc]/Vc[sc];

            QgSp_[mc] += Kg[mc]*tmpA/(tmpL*Vc[mc]);
            QsSp_[sc] += Kg[mc]*tmpA/(tmpL*Vc[sc]);
            QgSu_[mc] += Kg[mc]*Ts_.oldTime()[sc]*tmpA/(tmpL*Vc[mc]);
//             QgSu_[mc] += mCM_.cellMixture(mc).Hs(Ts_[sc])*Kg[mc]*tmpA
//                        / (Cpg()[mc]*tmpL*Vc[mc]);
            QsSu_[sc] += Kg[mc]*gasThermo_.T().oldTime()[mc]*tmpA/(tmpL*Vc[sc]);
        }

        // Set small cell value based on steady state conduction
        if (smallGasCell[own] + smallGasCell[nei])
        {
            scalar tmpA = ib_.alphafU()[faceI]*mesh_.magSf()[faceI];

            // Calculate thermal resistance
            scalar Lown = mag((mesh_.Cf()[faceI] - mesh_.C()[own])
                        & mesh_.Sf()[faceI])/mesh_.magSf()[faceI];
            scalar Lnei = mag((mesh_.Cf()[faceI] - mesh_.C()[nei])
                        & mesh_.Sf()[faceI])/mesh_.magSf()[faceI];
            scalar Rown = Lown/Kg[own];
            scalar Rnei = Lnei/Kg[nei];
            scalar Req = Rown + Rnei;

            // Add thermal resistance to system
            if (smallGasCell[own])
            {
                QgSp_[own] += tmpA/(Req*Vc[own]);
//                 QgSu_[own] = gasThermo_.hs()[nei]
//                            / (Cpg()[own]*Req*Vc[own]);
                QgSu_[own] += gasThermo_.T().oldTime()[nei]*tmpA/(Req*Vc[own]);
            }
            if (smallGasCell[nei])
            {
                QgSp_[nei] += tmpA/(Req*Vc[nei]);
//                 QgSu_[nei] = gasThermo_.hs()[own]
//                            / (Cpg()[nei]*Req*Vc[nei]);
                QgSu_[nei] += gasThermo_.T().oldTime()[own]*tmpA/(Req*Vc[nei]);
            }
        }
        if (smallSolidCell[own] + smallSolidCell[nei])
        {
            scalar tmpA = ib_.alphafsU()[faceI]*mesh_.magSf()[faceI];

            // Calculate thermal resistance
            scalar Lown = mag((mesh_.Cf()[faceI] - mesh_.C()[own])
                        & mesh_.Sf()[faceI])/mesh_.magSf()[faceI];
            scalar Lnei = mag((mesh_.Cf()[faceI] - mesh_.C()[nei])
                        & mesh_.Sf()[faceI])/mesh_.magSf()[faceI];
            scalar Rown = Lown/Ks[own];
            scalar Rnei = Lnei/Ks[nei];
            scalar Req = Rown + Rnei;

            // Add thermal resistance to system
            if (smallSolidCell[own])
            {
                QsSp_[own] += tmpA/(Req*Vc[own]);
                QsSu_[own] += Ts_.oldTime()[nei]*tmpA/(Req*Vc[own]);
            }
            if (smallSolidCell[nei])
            {
                QsSp_[nei] += tmpA/(Req*Vc[nei]);
                QsSu_[nei] += Ts_.oldTime()[own]*tmpA/(Req*Vc[nei]);
            }
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

                // Boundary cell solid, neighbor cell mixed
                if (fullCell[pfCellI] && !solidCellPNf[pFaceI])
                {
                    scalar tmpA = ib_.alphafU()[pFaceI]
                                    * mesh_.magSf()[pFaceI];
                    //Info << "Parallel Solid Cell" << endl;
                    scalar tmpL = mag((mesh_.Cf()[pFaceI]
                                - meshCPNf[pFaceI]) & mesh_.Sf()[pFaceI])
                                / mesh_.magSf()[pFaceI];


//                     Qt_s_[pfCellI] -= KgPNf[pFaceI]*(Ts_[pfCellI]
//                                     - TPNf[pFaceI])
//                                     * tmpA/(tmpL*Vc[pfCellI]);
                    QsSp_[pFaceI] += KgPNf[pFaceI]*tmpA/(tmpL*Vc[pfCellI]);
                    QsSu_[pFaceI] += KgPNf[pFaceI]*TPNf[pFaceI]*tmpA
                                   / (tmpL*Vc[pfCellI]);
                }
                // Boundary cell mixed, neighbor cell solid
                else if (fullCellPNf[pFaceI] && !solidCell[pfCellI])
                {
                    scalar tmpA = ib_.alphafU()[pFaceI]
                                    * mesh_.magSf()[pFaceI];
                    //Info << "Parallel Mixed Cell" << endl;
                    scalar tmpL = mag((mesh_.Cf()[pFaceI]
                                - mesh_.C()[pfCellI]) & mesh_.Sf()[pFaceI])
                                / mesh_.magSf()[pFaceI];

//                     Qt_g_[pfCellI] += Kg[pfCellI]*(TsPNf[pFaceI]
//                                     - gasThermo_.T()[pfCellI])
//                                     * tmpA/(tmpL*Vc[pfCellI]);
                    QgSp_[pFaceI] += Kg[pFaceI]*tmpA/(tmpL*Vc[pfCellI]);
                    QgSu_[pFaceI] += TsPNf[pFaceI]* Kg[pFaceI]*tmpA
                                   / (tmpL*Vc[pfCellI]);
//                     QgSu_[pFaceI] += mCM_.cellMixture(pfCellI).Hs(TsPNf[pFaceI])
//                                    * Kg[pFaceI]*tmpA
//                                    / (Cpg()[pfCellI]*tmpL*Vc[pfCellI]);
                }

                // Set small cell value based on steady state conduction
                if (smallGasCell[pfCellI])
                {
                    scalar tmpA = ib_.alphafU()[pFaceI] * mesh_.magSf()[pFaceI];

                    // Calculate thermal resistance
                    scalar Lown = mag((mesh_.Cf()[pFaceI]
                                - mesh_.C()[pfCellI]) & mesh_.Sf()[pFaceI])
                                / mesh_.magSf()[pFaceI];
//                     scalar Lnei = mag((mesh_.Cf()[faceI] - mesh_.C()[nei])
//                                 & mesh_.Sf()[faceI])/mesh_.magSf()[faceI];
                    scalar Rown = Lown/Kg[pfCellI];
//                     scalar Rnei = Lnei/Kg[nei];
                    scalar Req = Rown; // + Rnei;

                    // Add thermal resistance to system
                    QgSp_[pfCellI] += tmpA/(Req*Vc[pfCellI]);
//                    QgSu_[own] = gasThermo_.hs()[nei]
//                              / (Cpg()[own]*Req*Vc[own]);
                    QgSu_[pfCellI] += TPNf[pFaceI]*tmpA/(Req*Vc[pfCellI]);
                }
                if (smallSolidCell[pfCellI])
                {
                    scalar tmpA = ib_.alphafsU()[pFaceI] * mesh_.magSf()[pFaceI];

                    // Calculate thermal resistance
                    scalar Lown = mag((mesh_.Cf()[pFaceI]
                                - mesh_.C()[pfCellI]) & mesh_.Sf()[pFaceI])
                                / mesh_.magSf()[pFaceI];
//                     scalar Lnei = mag((mesh_.Cf()[faceI] - mesh_.C()[nei])
//                                 & mesh_.Sf()[faceI])/mesh_.magSf()[faceI];
                    scalar Rown = Lown/Ks[pFaceI];
//                     scalar Rnei = Lnei/Ks[nei];
                    scalar Req = Rown; //+ Rnei;

                    // Add thermal resistance to system
                    QsSp_[pFaceI] += tmpA/(Req*Vc[pFaceI]);
                    QsSu_[pFaceI] += TsPNf[pFaceI]*tmpA/(Req*Vc[pFaceI]);
                }
            }
        }
    }

    // Emmulate explicit behaviour
//     QgSu_ = QgSu_ - QgSp_*gasThermo_.T();
//     QgSp_ = QgSp_*0.0;
//     QgSu_[cellI] = QgSu_[cellI] - ib_.gasCells()()[cellI]*QgSp_[cellI]*gasThermo_.T()[cellI];
//     QgSp_[cellI] = ib_.smallCells()()[cellI]*QgSp_[cellI];
//     QsSu_ = QsSu_ - QsSp_*Ts_;
//     QsSp_ = QsSp_*0.0;

    dimensionedScalar hsrdT
    (
        "hsrdT",
        dimDensity/dimTime,
        1.0/mesh_.time().deltaTValue()
    );
    dimensionedScalar TsrdT
    (
        "TsrdT",
        dimPower/dimVolume/dimTemperature,
        1.0/mesh_.time().deltaTValue()
    );

    // Solve conduction for small cells only
    volScalarField rhog = gasThermo_.rho();
    volScalarField rhos = solidThermo_->rho();
    volScalarField Cps = solidThermo_->Cp();
    volScalarField gasCell = ib_.gasCells();
    forAll(QgSu_, cellI)
    {

        QgSu_[cellI] = QgSu_[cellI] - gasCell[cellI]*QgSp_[cellI]*gasThermo_.T()[cellI];
        QgSp_[cellI] = smallGasCell[cellI]*QgSp_[cellI];

        if (smallGasCell[cellI])
        {
            scalar Tgas = QgSu_[cellI]/QgSp_[cellI];
            scalar hsGas = mCM_.cellMixture(cellI).Hs(Tgas);
            hsSp_[cellI] = hsrdT.value();
            hsSu_[cellI] = hsGas*hsrdT.value();
            QsSu_[cellI] += (gasThermo_.hs().oldTime()[cellI] - hsGas)
                         *  rhog[cellI]/mesh_.time().deltaTValue();
            QgSu_[cellI] = 0.0;
            QgSp_[cellI] = 0.0;
        }

//         if (smallSolidCell[cellI])
//         {
//             scalar Tsolid = QsSu_[cellI]/QsSp_[cellI];
//             TsSp_[cellI] = TsrdT.value();
//             TsSu_[cellI] = Tsolid*TsrdT.value();
//             QgSu_[cellI] += (Ts_[cellI] - Tsolid)
//                          *  rhos[cellI]*Cps[cellI]/mesh_.time().deltaTValue();
//             QsSu_[cellI] = 0.0;
//             QsSp_[cellI] = 0.0;
//         }

    }
}

// ************************************************************************* //
