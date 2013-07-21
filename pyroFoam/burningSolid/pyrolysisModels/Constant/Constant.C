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

#include "Constant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pyrolysisModels
{
    defineTypeNameAndDebug(Constant, 0);
    addToRunTimeSelectionTable(pyrolysisModel, Constant, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::pyrolysisModels::Constant::Constant
(
    dictionary pyroDict,
    basicSolidThermo const* solidThermo
)
:
    pyrolysisModel(typeName, pyroDict, solidThermo),
    mb_(pyroModelDict_.lookup("mb"))
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::dimensionedScalar Foam::pyrolysisModels::Constant::mass_burning_rate
(
    const dimensionedScalar& T,
    const dimensionedScalar& p,
    const label& cell
)
{
    return mb_;
}

// ************************************************************************* //
