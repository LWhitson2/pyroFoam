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
    Foam::Arrhenius

Description
    Arrhenius pyrolysis model

SourceFiles
    Arrhenius.C

\*---------------------------------------------------------------------------*/

#ifndef Arrhenius_H
#define Arrhenius_H

#include "pyrolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pyrolysisModels
{
/*---------------------------------------------------------------------------*\
                      Class evaporationModel Declaration
\*---------------------------------------------------------------------------*/

class Arrhenius
:
    public pyrolysisModel
{

    //model-specific stuff

        dimensionedScalar A_;  //should be in kg/(m2-s)
        dimensionedScalar Ea_; //should be in J/kmol

public:

    //- Runtime type information
    TypeName("Arrhenius");

    // Constructors

        //- construct from components
        Arrhenius
        (
            dictionary pyroDict,
            basicSolidThermo const* solidThermo_
        );

    //- Destructor
        virtual ~Arrhenius() {}


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
