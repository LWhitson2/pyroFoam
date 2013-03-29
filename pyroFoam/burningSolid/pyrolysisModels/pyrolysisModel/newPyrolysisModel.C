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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::pyrolysisModel>
Foam::pyrolysisModel::New
(
    dictionary pyroDict,
    basicSolidThermo const* solidThermo
)
{
    // Look in pyroDict to see which model is specified
    word pyrolysisModelTypeName
    (
        pyroDict.lookup("pyrolysisModel")
    );

    Info<< "Selecting pyrolysis model "
        << pyrolysisModelTypeName << endl;

    // Attempt to load the specified model from the list of valid models
    componentsConstructorTable::iterator cstrIter =
        componentsConstructorTablePtr_
            ->find(pyrolysisModelTypeName);

    // Die if the specified model cannot be found
    if (cstrIter == componentsConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "pyrolysisModel::New"
        )   << "Unknown pyrolysisModel type "
            << pyrolysisModelTypeName << endl << endl
            << "Valid  pyrolysisModel are : " << endl
            << componentsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    // Return a pointer to the newly created pyrolysis model using macro magic
    return autoPtr<pyrolysisModel>(cstrIter()(pyroDict, solidThermo));
}


// ************************************************************************* //
