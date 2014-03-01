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

#include "Arrhenius.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pyrolysisModels
{
    defineTypeNameAndDebug(Arrhenius, 0);
    addToRunTimeSelectionTable(pyrolysisModel, Arrhenius, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::pyrolysisModels::Arrhenius::Arrhenius
(
    dictionary pyroDict,
    solidThermo const* solidThermo
)
:
    pyrolysisModel(typeName, pyroDict, solidThermo),
    A_(pyroModelDict_.lookup("A")),
    Ea_(pyroModelDict_.lookup("Ea"))
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::dimensionedScalar Foam::pyrolysisModels::Arrhenius::mass_burning_rate
(
    const dimensionedScalar& T,
    const dimensionedScalar& p,
    const label& cell
)
{
    dimensionedScalar Ru("Ru",dimEnergy/dimMoles/dimTemperature,8314);
    return A_*Foam::exp(-Ea_/(Ru*T));
}

// ************************************************************************* //
