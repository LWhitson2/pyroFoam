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

#include "WSB.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pyrolysisModels
{
    defineTypeNameAndDebug(WSB, 0);
    addToRunTimeSelectionTable(pyrolysisModel, WSB, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::pyrolysisModels::WSB::WSB
(
    dictionary pyroDict,
    basicSolidThermo const* solidThermo_
)
:
    pyrolysisModel(typeName, pyroDict, solidThermo_),
    Ac_(pyroModelDict_.lookup("Ac")),
    Ec_(pyroModelDict_.lookup("Ec")),
    T0_(pyroModelDict_.lookup("T0"))
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::dimensionedScalar Foam::pyrolysisModels::WSB::mass_burning_rate
(
    const dimensionedScalar& T,
    const dimensionedScalar& p,
    const label& cell
)
{
    dimensionedScalar Ru("Ru",dimEnergy/dimMoles/dimTemperature,8314.);
    scalar mb = 0.;

    scalar rho = solidThermo_->rho()[cell];
    scalar k = solidThermo_->K()()[cell];
    scalar cp = solidThermo_->Cp()()[cell];
    scalar denom = Ec_.value()*(cp*(T.value() - T0_.value()) - Qc_.value()/2.);
    if (denom > 0.0)
    {
        mb = Foam::sqrt(Ac_.value()*Ru.value()*Foam::sqr(T.value())*k*rho
           * Foam::exp(-Ec_.value()/(Ru.value()*T.value()))
           / denom);
    }
    else mb = 0.0;
    return dimensionedScalar("mb", dimMass/dimArea/dimTime, mb);
}

// ************************************************************************* //
