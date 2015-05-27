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

//#include "burningSolid.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class GasThermoType, class ReactionThermoType>
Foam::burningSolid<GasThermoType,ReactionThermoType>::burningSolid
(
    const fvMesh& mesh,
    immersedBoundary& ib,
    const ReactionThermoType& gasThermo
)
:
    mesh_(mesh),

    solidDict_
    (
        IOobject
        (
            "solidThermophysicalProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    solidThermo_
    (
        solidThermo::New( mesh , solidDict_)
    ),

    ib_(ib),

    gasThermo_(gasThermo),

    mCM_
    (
        dynamic_cast<const multiComponentMixture<GasThermoType>& >
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

    pyroModel_(pyrolysisModel::New(pyroDict_, solidThermo_.operator->())),

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

    mgen_
    (
        IOobject
        (
            "mgen",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("mgen", dimDensity/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Ti_
    (
        IOobject
        (
            "Ti",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
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

    rhoSp_
    (
        IOobject
        (
            "rhoSp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("rhoSp", dimless/dimTime, 0.0)
    ),

    rhoSu_
    (
        IOobject
        (
            "rhoSu",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("rhoSu", dimDensity/dimTime, 0.0)
    ),

    mflux_
    (
        IOobject
        (
            "mflux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("mflux", dimMass/dimArea/dimTime, 0.0)
    ),

    qflux_
    (
        IOobject
        (
            "qflux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qflux", dimPower/dimArea, 0.0)
    ),

    qgeng_
    (
        IOobject
        (
            "qgeng",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qgeng", dimPower/dimVolume, 0.0)
    ),

    qgens_
    (
        IOobject
        (
            "qgens",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qgens", dimPower/dimVolume, 0.0)
    ),

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
    ),

    ignFlux_(pyroDict_.lookup("ignFlux")),

    ignTime_(pyroDict_.lookup("ignTime")),

    ignRelax_(pyroDict_.lookup("ignRelax")),

    testPyro_(pyroDict_.lookup("testPyro"))

{
    Info<< "Created burningSolid" << endl;

    // Rename T in solidThermo to Ts to avoid name conflicts with the T
    //  in gasThermo.
    solidThermo_->T().rename("Ts");
    solidThermo_->rho().rename("rhos");

    // Find gasName in species
    const multiComponentMixture<GasThermoType>& composition =
        dynamic_cast<const multiComponentMixture<GasThermoType>& >
        ( gasThermo_.composition() );

    const PtrList<GasThermoType>& speciesData = composition.speciesData();

    forAll(speciesData, s)
    {
        const specie& specieS = dynamic_cast<const specie&>(speciesData[s]);
        
        if (specieS.name() == gasName_ )
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
template<class GasThermoType, class ReactionThermoType>
void Foam::burningSolid<GasThermoType,ReactionThermoType>::fixSmallCells()
{
    tmp<volScalarField> TsGas = gasThermo_.T();
    
    tmp<volScalarField> rhoSolid = solidThermo_->rho();

    volVectorField mU_transferred = mU_;

    // Small cell density value
    dimensionedScalar Rg
    (
        "Rg",
        dimVelocity*dimVelocity/dimTemperature,
        EMg_->R()
    );

    volScalarField rhog = gasThermo_.p()/(Ti_*Rg);

    volScalarField mgen_transferred = mgen_*(1. - rhog/rhoSolid);
    volScalarField qgeng_transferred = qgeng_*(1. - rhog/rhoSolid);

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

    // Value to force gas cells to designated density
    dimensionedScalar rdT
    (
        "rdT",
        dimless/dimTime,
        1.0/mesh_.time().deltaTValue()
    );

    // Values to set in the solid region
    tmp<volScalarField> pSolid = (gasThermo_.p() - gasThermo_.p())
                               + gasThermo_.p().weightedAverage(ib_.area());
                               
    //TODO: Fix these to do a proper creation rather than this X-X
    tmp<volVectorField> USolid = (burnU_ - burnU_)
                               + dimensionedVector("Us",dimVelocity,vector::zero);
                               
    tmp<volScalarField> heSolid = (gasThermo_.he() - gasThermo_.he());
    
    forAll(heSolid(),cellI)
    {
        heSolid()[cellI] = mCM_.cellMixture(cellI).Hs
        (
            pSolid()[cellI],
            Ti_[cellI]
        );
    }

    // Calculate transfer weights
    tmp<surfaceScalarField> tw = ib_.scTransferWeights("gas");
    tmp<surfaceScalarField> tws = ib_.scTransferWeights("solid");
    const surfaceScalarField& w = tw();
    const surfaceScalarField& ws = tws();
    Info << "Min/Max Gas Transfer Weights: " << min(w).value() << ", " << max(w).value() << endl;
    Info << "Min/Max Solid Transfer Weights: " << min(ws).value() << ", " << max(ws).value() << endl;

    // Set source terms in small and solid cells to specify value
    Info << "Setting source terms" << endl;
    ib_.setScValue<vector>(w, USu_, USp_, burnU_,
                           USolid, rhordT, "fix", "gas");
    ib_.setScValue<scalar>(w, pSu_, pSp_, gasThermo_.p(),
                           pSolid, psirdT, "avg", "gas");
    ib_.setScValue<scalar>(w, hsSu_, hsSp_, gasThermo_.he(),
                           heSolid, hsrdT, "avg", "gas");
    ib_.setScValue<scalar>(ws, TsSu_, TsSp_, Ts_,
                           TsGas, TsrdT, "avg", "solid");
    ib_.setScValue<scalar>(w, rhoSu_, rhoSp_, rhog,
                           rhoSolid, rdT, "avg", "gas");

    // Transfer mass and momentum out of small cells
    Info << "Transferring out of small cells" << endl;
    ib_.transfer<scalar>(w, mgen_transferred, mgen_, 0.0, "gas");
    ib_.transfer<vector>(w, mU_transferred, mU_, vector::zero, "gas");
    ib_.transfer<scalar>(w, qgeng_transferred, qgeng_, 0.0, "gas");

    // Solve conduction for small cells only
    forAll(QgSu_, cellI)
    {
        if (ib_.smallCells()()[cellI])
        {
            hsSp_[cellI] = 0.0;
            hsSu_[cellI] = 0.0;
        }

        if (ib_.smallSolidCells()()[cellI])
        {
            TsSp_[cellI] = 0.0;
            TsSu_[cellI] = 0.0;
        }
    }
}

// Calculate the burn gas velocity
template<class GasThermoType, class ReactionThermoType>
void Foam::burningSolid<GasThermoType,ReactionThermoType>::calcBurnU()
{
    //Get the gas constant for the pyrolysis gas
    dimensionedScalar Rg
    (
        "Rg",
        dimVelocity*dimVelocity/dimTemperature,
        EMg_->R()
    );

    burnU_ = (Ti_*Rg/gasThermo_.p() - 1./solidThermo_->rho())
           * mflux_*ib_.normal().oldTime(); //should use oldTime??
}

template<class GasThermoType, class ReactionThermoType>
void Foam::burningSolid<GasThermoType,ReactionThermoType>::calcInterfaceFlux()
{
    mflux_ = dimensionedScalar("zero", dimMass/dimArea/dimTime, 0.);
    qflux_ = dimensionedScalar("zero", dimPower/dimArea, 0.);
    qgens_ = dimensionedScalar("zero", dimPower/dimVolume, 0.);
    qgeng_ = dimensionedScalar("zero", dimPower/dimVolume, 0.);

    // Surface pyrolysis fluxes
    forAll(Ts_, cellI)
    {
        if (ib_.area().oldTime()[cellI] > SMALL)
        {
            mflux_[cellI] = pyroModel_->mass_burning_rate(
                Ti_[cellI], gasThermo_.p()[cellI], cellI).value();

            qflux_[cellI] = pyroModel_->energy_generation(
                Ti_[cellI], gasThermo_.p()[cellI], cellI).value();
        }
    }

    // Laser ignition flux
    volScalarField Cps = solidThermo_->Cp();
    if (mesh_.time().timeOutputValue() < ignTime_.value())
    {
        qflux_ += ignFlux_;
    }
    else if (mesh_.time().timeOutputValue() < (ignTime_ + ignRelax_).value())
    {
        scalar time = mesh_.time().timeOutputValue();
        scalar maxTime = (ignTime_.value() + ignRelax_.value());
        dimensionedScalar ignFlux = ignFlux_*(maxTime - time)/ignRelax_.value();
        qflux_ += ignFlux;
    }
}

template<class GasThermoType, class ReactionThermoType>
void Foam::burningSolid<GasThermoType,ReactionThermoType>::calcSurfaceEnergy()
{
    // Convert surface flux to source
    // TODO add an AbyV or ArV function to ib
    // Handled in the surface balance for Ti.
    // qgens_.internalField() += qflux_*ib_.area().oldTime()/mesh_.V();
    // if (testPyro_ == "enthalpy")
    // {
    //     qgeng_.internalField() += qflux_*ib_.area().oldTime()/mesh_.V();
    // }

    // Interface area and volume
    const volScalarField& Ai = ib_.area();
    const volScalarField::DimensionedInternalField& Vc = mesh_.V();

    // Calculate energy transferred between solid and gas
    tmp<volScalarField> Cps = solidThermo_->Cp();
    forAll(mflux_, cellI)
    {
        // Energy lost by solid
        qgens_.internalField()[cellI] -=
              m_pyro_[cellI]*solidThermo_->Cp()()[cellI]
            * (Ti_[cellI] - Ts_.oldTime()[cellI]);

        // Energy gained by gas
        qgeng_.internalField()[cellI] += m_pyro_[cellI]
            * mCM_.cellMixture(cellI).Hs(gasThermo_.p()[cellI], Ti_[cellI]);
            // + pyroModel_->energy_generation(
            //     Ts_[cellI], gasThermo_.p()[cellI], cellI).value()*Ai[cellI]/Vc[cellI];
    }
}

template<class GasThermoType, class ReactionThermoType>
void Foam::burningSolid<GasThermoType,ReactionThermoType>::calcSurfaceMomentum(const volVectorField& U)
{
    // Convert surface momentum flux to source
    // TODO add an AbyV or ArV function to ib
    mU_ = burnU_*m_pyro_;
}

template<class GasThermoType, class ReactionThermoType>
void Foam::burningSolid<GasThermoType,ReactionThermoType>::calcSurfaceMass()
{
    // Convert surface mass flux to source
    // TODO add an AbyV or ArV function to ib
    mgen_ = m_pyro_;
}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class GasThermoType, class ReactionThermoType>
void Foam::burningSolid<GasThermoType,ReactionThermoType>::correct
(
    const volVectorField& U,
    surfaceScalarField& phi,
    bool allowHtx
)
{
    Foam::Info << "Correcting burningSolid" << Foam::endl;

    // Calculate mass and energy flux (mflux, qflux) using current P and Ts
    Foam::Info << "Calculate mass/energy flux at surface" << Foam::endl;
    calcInterfaceFlux();

    // Calculate mass and energy source using original cell area
    Foam::Info << "Calculate mass burning rate" << Foam::endl;
    m_pyro_.internalField() = mflux_ * ib_.area().oldTime() / mesh_.V();
    Info << "Interface Area: " << max(ib_.area().oldTime()).value() << ", " << max(ib_.area()).value() << endl;
    m_pyro_.correctBoundaryConditions();

    // Calculate burn gas velocity using current interface orientation
    Foam::Info << "Calculate gas velocity" << Foam::endl;
    calcBurnU();

    // Evolve interface using calculated burning rate (vol frac/s)
    Foam::Info << "Move Interface" << Foam::endl;
    ib_.moveInterface( m_pyro_ / solidThermo_->rho() );

    // // Calculate mass source
    Foam::Info << "Calculate mass source" << Foam::endl;
    calcSurfaceMass();

    // Calculate momentum source
    Foam::Info << "Calculate momentum source" << Foam::endl;
    calcSurfaceMomentum(U);

    // Calculate energy source
    Foam::Info << "Calculate energy source" << Foam::endl;
    calcSurfaceEnergy();

    Foam::Info << "Calculate gas/solid interface" << Foam::endl;
    calcInterfaceTransfer();

    // Fix small cells by transferring momentum (mU) and mass (m_pyro)
    // to neighbouring larger cells
    Foam::Info << "Fix small cells" << Foam::endl;
    fixSmallCells();

    // Calculate mass flux field that includes alphaf
    phi = (ib_.interpolate(U*gasThermo_.rho(), "gas") & mesh_.Sf())*ib_.alphaf();

    // Calculate surface stress in interface cells
    Foam::Info << "Calculate surface stress" << Foam::endl;
    // calcSurfaceStress();
}


template<class GasThermoType, class ReactionThermoType>
tmp<volScalarField> Foam::burningSolid<GasThermoType,ReactionThermoType>::YSu
(
    const word& Yname
) const
{
    dimensionedScalar onerDt
    (
        "onerDt",
        dimless/dimTime,
        1.0/mesh_.time().deltaTValue()
    );

    if (Yname == gasName_)
    {
        return mgen_ + ib_.smallCells()*onerDt*gasThermo_.rho();
    }
    else
    {
        return ib_.noCells()*onerDt*gasThermo_.rho();
    }
}

template<class GasThermoType, class ReactionThermoType>
tmp<volScalarField> Foam::burningSolid<GasThermoType,ReactionThermoType>::YSp() const
{
    return ib_.smallAndSolidCells()*gasThermo_.rho()
            *dimensionedScalar
            (
                "onerDt",
                dimless/dimTime,
                1.0/mesh_.time().deltaTValue()
            );
}

// template<class GasThermoType, class ReactionThermoType>
// void Foam::burningSolid<GasThermoType,ReactionThermoType>::calcSurfaceStress()
// {
//     // Centroid lengths
//     const volScalarField& Lg = ib_.gasL();

//     // Interface area and volume
//     const volScalarField& Ai = ib_.area();
//     const volScalarField::DimensionedInternalField& Vc = mesh_.V();

//     // Initialize cells to no heat transfer
//     surfStress_ = dimensionedScalar("tmpSurfStress", dimless/dimArea, 0.0);

//     // Cell identification
//     volScalarField normalCell = ib_.mixedCells();
//     volScalarField solidCell = ib_.solidCells();
//     volScalarField fullCell = solidCell*pos(Ai
//                             - dimensionedScalar("tmp", dimArea, SMALL));

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

template<class GasThermoType, class ReactionThermoType>
void Foam::burningSolid<GasThermoType,ReactionThermoType>::calcInterfaceTransfer()
{
    // Solid Properties
    volScalarField rhos = solidThermo_->rho();
    volScalarField Cps = solidThermo_->Cp();
    volScalarField Ks = solidThermo_->kappa();
    // if (testPyro_ == "enthalpy") Ks = Ks*0.;

    // Gas Properties
    volScalarField rhog = gasThermo_.rho();
    tmp<volScalarField> Cpg = gasThermo_.Cp();
    tmp<volScalarField> Kg = gasThermo_.kappa();
    if (testPyro_ == "solid") Kg = Kg*0.;
    const volScalarField& Tg = gasThermo_.T();

    const dimensionedScalar vsL = dimensionedScalar("vsL", dimLength, SMALL);

//     Info << "Kg min/max: " << min(Kg) << ", " << max(Kg) << endl;

    // Conduction lengths
    const volScalarField& Lg = ib_.gasL();
    const volScalarField& Ls = ib_.solidL();

    // Interface area and volume
    const volScalarField& Ai = ib_.area();
    const volScalarField::DimensionedInternalField& Vc = mesh_.V();

    // Initialize cells to no heat transfer
    QgSp_ = dimensionedScalar("zero", dimDensity/dimTime, 0.0);
    QsSp_ = dimensionedScalar("zero", dimPower/dimTemperature/dimVolume, 0.0);
    QgSu_ = dimensionedScalar("zero", dimPower/dimVolume, 0.0);
    QsSu_ = dimensionedScalar("zero", dimPower/dimVolume, 0.0);

    // Cell identification
    volScalarField normalCell = ib_.reconstructedCells()*pos(Lg - vsL)*pos(Ls - vsL);
    volScalarField solidCell = ib_.solidCells();
    volScalarField smallGasCell = ib_.smallCells();
    volScalarField smallSolidCell = ib_.smallSolidCells();
    volScalarField fullCell = pos(Ai - dimensionedScalar("tmp", dimArea, SMALL))
                            * neg(ib_.reconstructedCells() - SMALL);

    // Initialize Ti to Ts
    Ti_ = Ts_;

    Info << "Normal Cell Conduction" << endl;
    forAll(Ts_,cellI)
    {
        // Normal mixed cell conduction transfer
        if (normalCell[cellI])
        {
            // Calculate thermal constants
            // Info << ib_.alpha()[cellI] << endl;
            scalar Cg = Kg()[cellI]/Lg[cellI];
            scalar Cs = Ks[cellI]/Ls[cellI];
            
            // Info << "Cell " << cellI <<" is normal" << endl;
            // Info << "Lg/Ls: " << Lg[cellI] << ", " << Ls[cellI] << endl;
            // Info << "Kg/Ks: " << Kg[cellI] << ", " << Ks[cellI] << endl;
            // Info << "Cg/Cs: " << Cg << ", " << Cs << endl;

            // Info << Cg << " " << Cs << endl;

            // Calculate interface temperature
            Ti_[cellI] = (Cs*Ts_[cellI] + Cg*Tg[cellI] + qflux_[cellI])
                       / (Cs + Cg);

            // Set Ti explicitly
            // Ti_[cellI] = 714.64;

            if ((testPyro_ != "solid"))
            {
                // Calculate thermal resistance
                scalar Req = 0.;
                if (Cg > SMALL) Req += 1./Cg;
                if (Cs > SMALL) Req += 1./Cs;

                // Solid source terms (dirchlet boundary at Ti)
                // QsSp_[cellI] = Ai[cellI]/(Req*Vc[cellI]);
                // QsSu_[cellI] = Tg[cellI]
                //             * Ai[cellI]/(Req*Vc[cellI]);
                QsSp_[cellI] = Ai[cellI]*Cs/Vc[cellI];
                QsSu_[cellI] = Ti_[cellI]
                            * Ai[cellI]*Cs/Vc[cellI];

                // Gas source terms (neumann boundary)
                QgSp_[cellI] = Ai[cellI]*Cg/(Cpg()[cellI]*Vc[cellI]);
                QgSu_[cellI] = mCM_.cellMixture(cellI).Hs(gasThermo_.p()[cellI], Ti_[cellI])
                            * Ai[cellI]*Cg/(Cpg()[cellI]*Vc[cellI]);
                // QgSp_[cellI] = 0.;
                // QgSu_[cellI] = (Ti_[cellI] - Tg[cellI])
                //             * Ai[cellI]*Cg/Vc[cellI];

                

                // Transfer to Qgen for small cells
                // if (smallGasCell[cellI])
                // {
                //     Info << "QgSu/QgSp: " << QgSu_[cellI] << ", " << QgSp_[cellI] << endl;
                //     qgeng_[cellI] += QgSu_[cellI]
                //                   - QgSp_[cellI]*mCM_.cellMixture(cellI).Hs(gasThermo_.p()[cellI], Tg[cellI]);
                //     QgSu_[cellI] = 0.;
                //     QgSp_[cellI] = 0.;
                // }
            }
        }
    }

    Info << "Full Cell Conduction" << endl;
    forAll(mesh_.magSf(), faceI)
    {
        label own = mesh_.owner()[faceI];
        label nei = mesh_.neighbour()[faceI];
        label sc = (fullCell[own]) ? own:nei;
        label mc = (sc == own) ? nei:own;

        // Full solid cell to face neighbor conduction transfer
        if (fullCell[sc] && !solidCell[mc] && !fullCell[mc])
        {
            scalar tmpA = mesh_.magSf()[faceI];

            if (tmpA > SMALL)
            {
                scalar tmpLg = (mag((mesh_.Cf()[faceI] - ib_.gasC()[mc])
                             & mesh_.Sf()[faceI])/mesh_.magSf()[faceI]);
                scalar tmpLs = (mag((mesh_.Cf()[faceI] - mesh_.C()[sc])
                             & mesh_.Sf()[faceI])/mesh_.magSf()[faceI]);

                // Calculate thermal resistance
                scalar Cg = Kg()[mc]/tmpLg;
                scalar Cs = Ks[sc]/tmpLs;

                // Info << "Cell " << sc <<" is full" << endl;
                // Info << "Lg/Ls: " << tmpLg << ", " << tmpLs << endl;
                // Info << "Kg/Ks: " << Kg[mc] << ", " << Ks[sc] << endl;
                // Info << "Cg/Cs: " << Cg << ", " << Cs << endl;

                // Set Ti explicitly
                // Ti_[sc] = 714.64;
                // tmpTi = 714.64;

                if (testPyro_ != "solid")
                {
                    // Calculate interface temperature
                    scalar tmpTi = (Cs*Ts_[sc] + Cg*Tg[mc] + qflux_[sc])
                           / (Cs + Cg);
                    Ti_[sc] = tmpTi; // TODO Add function to calculate average Ti

                    // Calculate thermal resistance
                    scalar Req = 0.;
                    if (Cg > SMALL) Req += 1./Cg;
                    if (Cs > SMALL) Req += 1./Cs;

                    // // Gas source terms
                    QgSp_[mc] += tmpA*Cg/(Cpg()[mc]*Vc[mc]);
                    QgSu_[mc] += mCM_.cellMixture(mc).Hs(gasThermo_.p()[mc], tmpTi)
                               * tmpA*Cg/(Cpg()[mc]*Vc[mc]);
                    // QgSp_[mc] += 0.;
                    // QgSu_[mc] += (tmpTi - Tg[mc])
                    //            * tmpA*Cg/Vc[mc];

                    // Solid source terms
                    // QsSp_[sc] += tmpA/(Req*Vc[sc]);
                    // QsSu_[sc] += Tg[mc]
                    //            * tmpA/(Req*Vc[sc]);
                    QsSp_[sc] += tmpA*Cs/Vc[sc];
                    QsSu_[sc] += tmpTi
                               * tmpA*Cs/Vc[sc];

                    // Solid side gas source terms
                    QgSu_[sc] += mCM_.cellMixture(mc).Hs(gasThermo_.p()[sc], Ti_[sc])
                               * Cs*tmpA/(Cpg()[sc]*Vc[sc])
                               - qflux_[sc]*tmpA/Vc[sc];
                    QgSp_[sc] += tmpA*Cs/(Cpg()[sc]*Vc[sc]);
                }
                else
                {
                    // Calculate interface temperature
                    scalar tmpTi = (Cs*Ts_[sc] + qflux_[sc]) / Cs;
                    Ti_[sc] = tmpTi;

                    // Solid Source Terms
                    
                }
            }
        }
    }

    // TODO Update parallel routine for calculating interface transfer
    /*// Full solid cell to parallel neighbor conduction transfer
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

                        // Calculate current interface temperature
                        scalar tmpTi = (Rg*Ts_[pfCellI] + Rs*TPNf[pFaceI])/Req;

                        // Solid source terms
                        QsSp_[pfCellI] += 1./(Req*Vc[pfCellI]);
                        QsSu_[pfCellI] += TPNf[pFaceI]/(Req*Vc[pfCellI]);

                        // Partial calculation for Ti of solid cell
                        sumAiTi[pfCellI] += tmpTi*tmpA;
                        sumAi[pfCellI] += tmpA;
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
    }*/
}

// ************************************************************************* //
