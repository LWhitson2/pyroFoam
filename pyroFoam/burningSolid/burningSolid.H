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
#include "hsCombustionThermo.H"
#include "cuttableCell.H"
#include "cuttableFace.H"
#include "thermoPhysicsTypes.H"
#include "multiComponentMixture.H"
#include "basicSolidThermo.H"
#include "immersedBoundary.H"
#include "pyrolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class phase Declaration
\*---------------------------------------------------------------------------*/

class burningSolid
{
    // Private data
        //- Mesh
        const fvMesh& mesh_;

        //- Solid thermo pointer
        autoPtr<basicSolidThermo> solidThermo_;

        //- Immersed solid boundary
        immersedBoundary& ib_;

        //- Gas thermo
        const hsCombustionThermo& gasThermo_;

        //- Multicomponent Mixture
        const multiComponentMixture<gasThermoPhysics>& mCM_;

        //- Input dictionary
        IOdictionary pyroDict_;

        //- Pyrolysis model
        autoPtr<pyrolysisModel> pyroModel_;

        //- Solid burning rate (kg/m^3-s)
        volScalarField m_pyro_;

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

        //- Constant burning rate (for now, kg/m2/s)
        dimensionedScalar m0_;

        //- Name of the solid specie
        word solidName_;

        //- Name of the pyrolized gas specie
        word gasName_;

        //- Pointer to the thermo for the pyro gas
        const gasThermoPhysics* EMg_;

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

    // Private functions
        void fixSmallCells();
        void calcBurnU();
//         void calcHeatTransfer();
        void calcSurfaceStress();

public:

    // Constructor
    burningSolid
    (
        const fvMesh& mesh,
        immersedBoundary& ib,
        const hsCombustionThermo& gasThermo
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

    volScalarField& Ts()
    {
        return Ts_;
    }

    const volScalarField& rhos()
    {
        return solidThermo_->rho();
    }

    const volScalarField& Ks()
    {
        return solidThermo_->K();
    }

    const tmp<volScalarField> Cps()
    {
        return solidThermo_->Cp();
    }

    basicSolidThermo& solidThermo()
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

    void calcHeatTransfer();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
