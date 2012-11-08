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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class phase Declaration
\*---------------------------------------------------------------------------*/

class burningSolid
{
    // Private data

        // Defining mesh and dictionary
        const fvMesh& mesh_;                // Read-Only mesh
        IOdictionary pyroDict_;             // Input dictionary

        // Field references to the outside world

        const volVectorField& U_;           // Read-Only velocity vector
        surfaceScalarField& phi_;           // Read-Write face flux
        const hsCombustionThermo& thermo_;  // Read-Only thermo properties

        // Fields owned by this class
        volScalarField alpha_;        // Volume fraction of GAS
        surfaceScalarField alphaf_;   // Face area fraction of gas
        volScalarField m_pyro_;       // Solid burning rate
        volScalarField isBurning_;    // Flag for burning cells
        volScalarField a_burn_;       // Burning area
        volVectorField burnU_;        // Burning gas velocity

        volVectorField a_burn_norm_;  // Burning area normal vector
        volVectorField mU_;           // Burning gas small cell momentum
        volVectorField USp_;
        volVectorField USu_;
        volScalarField pSp_;
        volScalarField pSu_;

        dimensionedScalar rhoS_;      // Density of energetic material
        dimensionedScalar m0_;        // Constant burning rate (for now, kg/m2/s)

    // Private member functions
        void calcAlphaf();            // Calculate gas face fraction
        void calcBurningArea();       // Calculate burning face areas
        void calcBurnU();             // Calculate burning gas velocity

public:

    // Constructors

        //- Construct from components
        burningSolid
        (
            const volVectorField& U,
            surfaceScalarField& phi,
            const hsCombustionThermo& thermo
        );

        //- Return clone
        autoPtr<burningSolid> clone() const;



    // Public Member Functions

        const surfaceScalarField& alphaf() const
        {
            return alphaf_;
        }

        const volScalarField& alpha() const
        {
            return alpha_;
        }

        const volScalarField& m_pyro() const
        {
            return m_pyro_;
        }

        const volScalarField& isBurning() const
        {
            return isBurning_;
        }

        const dimensionedScalar& rhoS() const
        {
            return rhoS_;
        }

        const volScalarField& a_burn() const
        {
            return a_burn_;
        }

        const volVectorField& burnU() const
        {
            return burnU_;
        }

        const volVectorField& a_burn_norm() const
        {
            return a_burn_norm_;
        }

        tmp<volScalarField> Sh() const;

        tmp<volScalarField> setDiag() const;

        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
