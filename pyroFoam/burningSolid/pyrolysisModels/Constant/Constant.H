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
    Foam::Constant

Description
    Constant pyrolysis model

SourceFiles
    Constant.C

\*---------------------------------------------------------------------------*/

#ifndef Constant_H
#define Constant_H

#include "pyrolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pyrolysisModels
{
/*---------------------------------------------------------------------------*\
                      Class evaporationModel Declaration
\*---------------------------------------------------------------------------*/

class Constant
:
    public pyrolysisModel
{

    //model-specific stuff

        dimensionedScalar mb_;  //constant burn rate (kg/m2/s)

public:

    //- Runtime type information
    TypeName("Constant");

    // Constructors

        //- construct from components
        Constant
        (
            dictionary pyroDict,
            solidThermo const* solidThermo_
        );

    //- Destructor
        virtual ~Constant() {}


    // Member Functions
        virtual dimensionedScalar mass_burning_rate
        (
            const dimensionedScalar& T,
            const dimensionedScalar& p,
            const label& cell
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace pyrolysisModel

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
