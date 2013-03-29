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

#include "pyrolysisModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pyrolysisModel, 0);
    defineRunTimeSelectionTable(pyrolysisModel, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pyrolysisModel::pyrolysisModel
(
    const word& type,
    dictionary pyroDict,
    basicSolidThermo const* solidThermo
)
:   //allow an empty dictionary for NoPyrolysis model, all others will crash
    pyroModelDict_(pyroDict.subOrEmptyDict(type + "Coeffs")),
    solidThermo_(solidThermo),
    Qc_(dimensionedScalar("zero", dimEnergy/dimMass, 0.))
{
    if(!pyroModelDict_.empty()) Qc_ = pyroModelDict_.lookup("Qc");
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

 //Define any general pyrolysisModel functions here. These will be things
 // that are the same for every model.
Foam::dimensionedScalar Foam::pyrolysisModel::energy_generation
    (
        const dimensionedScalar& T,
        const dimensionedScalar& p,
        const label& cell
    )
{
    return mass_burning_rate(T, p, cell)*Qc_;
}

// ************************************************************************* //

