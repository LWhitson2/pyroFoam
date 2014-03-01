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

Class
    Foam::burningSolid

Description


SourceFiles
    burningSolid.C

\*---------------------------------------------------------------------------*/

#ifndef burningSolid_H
#define burningSolid_H

#include "volFields.H"
#include "dictionaryEntry.H"
#include "fvCFD.H"
#include "dimensionedScalar.H"
#include "dimensionedScalarFwd.H"
#include "cuttableCell.H"
#include "cuttableFace.H"
#include "specie.H"
#include "multiComponentMixture.H"
#include "solidThermo.H"
#include "immersedBoundary.H"
#include "pyrolysisModel.H"
//#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class phase Declaration
\*---------------------------------------------------------------------------*/
template<class GasThermoType, class ReactionThermoType>
class burningSolid
{
    // Private data
        //- Mesh
        const fvMesh& mesh_;

        //- Solid thermo pointer
        autoPtr<solidThermo> solidThermo_;

        //- Immersed solid boundary
        immersedBoundary& ib_;

        //- Gas thermo
        const ReactionThermoType& gasThermo_;

        //- Multicomponent Mixture
        const multiComponentMixture<GasThermoType>& mCM_;

        //- Input dictionary
        IOdictionary pyroDict_;

        //- Pyrolysis model
        autoPtr<pyrolysisModel> pyroModel_;

        //- Solid burning rate (kg/m^3-s)
        volScalarField m_pyro_;

        //- Gas mass generation from burning (kg/m^3-s)
        volScalarField mgen_;

        //- Interface temperature (K)
        volScalarField Ti_;

        //- Burning gas velocity (m/s)
        volVectorField burnU_;

        //- Burning gas momentum source
        volVectorField mU_;

        //- solid Temperature field
        volScalarField Ts_;

        //- Momentum implicit source term
        volScalarField USp_;

        //- Momentum explicit source term
        volVectorField USu_;

        //- Pressure implicit source term
        volScalarField pSp_;

        //- Pressure explicit source term
        volScalarField pSu_;

        //- Gas Enthalpy implicit source term
        volScalarField hsSp_;

        //- Gas Enthalpy explicit source term
        volScalarField hsSu_;

        //- Solid Temperature implicit source term
        volScalarField TsSp_;

        //- Solid Temperature explicit source term
        volScalarField TsSu_;

        //- Density implicit source term
        volScalarField rhoSp_;

        //- Density explicit source term
        volScalarField rhoSu_;

        //- Solid burning flux (kg/m2/s)
        volScalarField mflux_;

        //- Solid burning heat flux (J/m2/s)
        volScalarField qflux_;

        //- Gas burning heat generation (J/m3/s)
        volScalarField qgeng_;

        //- Solid burning heat generation (J/m3/s)
        volScalarField qgens_;

        //- Name of the solid specie
        word solidName_;

        //- Name of the pyrolized gas specie
        word gasName_;

        //- Pointer to the thermo for the pyro gas
        const GasThermoType* EMg_;

        //- Gas heat transfer implicit source
        volScalarField QgSp_;

        //- Solid heat transfer implicit source
        volScalarField QsSp_;

        //- Gas heat transfer explicit source
        volScalarField QgSu_;

        //- Solid heat transfer explicit source
        volScalarField QsSu_;

        //- Solid stress term
        volScalarField surfStress_;

        //- Igniter Flux
        dimensionedScalar ignFlux_;

        //- Igniter Time
        dimensionedScalar ignTime_;

        //- Igniter Relax Time
        dimensionedScalar ignRelax_;

        //- Flag to test various parts of code
        word testPyro_;

    // Private functions
        void fixSmallCells();
        void calcBurnU();
        void calcSurfaceStress();
        void calcInterfaceFlux();
        void calcSurfaceMass();
        void calcSurfaceMomentum(const volVectorField& U);
        void calcSurfaceEnergy();

    public:

        // Constructor
        burningSolid
        (
            const fvMesh& mesh,
            immersedBoundary& ib,
            const ReactionThermoType& gasThermo
        );

        //Functions
        const word& gasName()
        {
            return gasName_;
        }

        const word& solidName()
        {
            return solidName_;
        }

        const volVectorField& mU()
        {
            return mU_;
        }

        const volVectorField& USu()
        {
            return USu_;
        }

        const volScalarField& USp()
        {
            return USp_;
        }

        const volScalarField& pSp()
        {
            return pSp_;
        }

        const volScalarField& pSu()
        {
            return pSu_;
        }

        const volScalarField& m_pyro()
        {
            return m_pyro_;
        }

        const volScalarField& mgen()
        {
            return mgen_;
        }

        void correct
        (
            const volVectorField& U,
            surfaceScalarField& phi,
            bool allowHtx = true
        );

        tmp<volScalarField> YSu(const word& Yname) const;

        tmp<volScalarField> YSp() const;

        const volScalarField& hsSu()
        {
            return hsSu_;
        }

        const volScalarField& hsSp()
        {
            return hsSp_;
        }

        const volScalarField& TsSu()
        {
            return TsSu_;
        }

        const volScalarField& TsSp()
        {
            return TsSp_;
        }

        const volScalarField& rhoSu()
        {
            return rhoSu_;
        }

        const volScalarField& rhoSp()
        {
            return rhoSp_;
        }

        volScalarField& Ts()
        {
            return Ts_;
        }

        const volScalarField& rhos()
        {
            return solidThermo_->rho();
        }

        const tmp<volScalarField> Ks()
        {
            return solidThermo_->kappa();
        }

        const tmp<volScalarField> Cps()
        {
            return solidThermo_->Cp();
        }

        solidThermo& thermo()
        {
            return solidThermo_();
        }

        const volScalarField& QgSp()
        {
            return QgSp_;
        }

        const volScalarField& QsSp()
        {
            return QsSp_;
        }

        const volScalarField& QgSu()
        {
            return QgSu_;
        }

        const volScalarField& QsSu()
        {
            return QsSu_;
        }

        const volScalarField& Stress()
        {
            return surfStress_;
        }

        const volScalarField& Qgens()
        {
            return qgens_;
        }

        const volScalarField& Qgeng()
        {
            return qgeng_;
        }

        const volScalarField& Ti()
        {
            return Ti_;
        }

        const word testPyro()
        {
            return testPyro_;
        }

        void calcInterfaceTemp();

        void calcHeatTransferOld();

        void calcInterfaceTransfer();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "burningSolid.C"
#endif

#endif

// ************************************************************************* //
