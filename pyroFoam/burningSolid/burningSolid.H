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

        // Field references to the outside world
        const fvMesh& mesh_;                // Read-Only mesh
        const volVectorField& U_;           // Read-Only velocity vector
        surfaceScalarField& phi_;           // Read-Write face flux
        const hsCombustionThermo& thermo_;  // Read-Only thermo properties

        // Fields owned by this class
        volScalarField alpha_;        // Volume fraction of GAS
        surfaceScalarField alphaf_;   // Face fraction of gas
        volScalarField m_pyro_;       // Solid burning rate
        volScalarField isBurning_;    // Flag for burning cells
        dimensionedScalar rhoS_;      // Density of energetic material

    // Private member functions
        void calcAlphaf();            // Calculate gas face fraction

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

        tmp<volScalarField> Sh() const;

        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
