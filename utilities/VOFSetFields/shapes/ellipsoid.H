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
    Foam::phase

Description


SourceFiles
    ellipsoid.C

\*---------------------------------------------------------------------------*/

#ifndef ellipsoid_H
#define ellipsoid_H

#include "volFields.H"
#include "PtrDictionary.H"
#include "dictionaryEntry.H"
#include "fvCFD.H"
#include "dimensionedScalar.H"
#include "dimensionedScalarFwd.H"
#include "cuttableCell.H"
#include "shape.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace shapes
{
/*---------------------------------------------------------------------------*\
                         Class ellipsoid Declaration
\*---------------------------------------------------------------------------*/

class ellipsoid
:
    public shape
{
    // Private data
        vector center_;
        vector radius_;

        tmp<scalarField> r() const;
public:

    //- Runtime type information
        TypeName("ellipsoid");
    
    // Constructors

        //- Construct from components
        ellipsoid
        (
            const word& name,
            dictionary dropletDict,
            const fvMesh& mesh
        );


    // Member Functions

        void calculate();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace shapes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
