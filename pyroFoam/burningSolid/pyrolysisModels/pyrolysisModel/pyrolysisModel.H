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
    Foam::pyrolysisModel

Description
    Abstract pyrolysis base class

SourceFiles
    pyrolysisModel.C

\*---------------------------------------------------------------------------*/

#ifndef pyrolysisModel_H
#define pyrolysisModel_H

#include "volFields.H"
#include "PtrDictionary.H"
#include "dictionaryEntry.H"
#include "fvCFD.H"
#include "dimensionedScalar.H"
#include "dimensionedScalarFwd.H"
#include "solidThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class pyrolysisModel Declaration
\*---------------------------------------------------------------------------*/

class pyrolysisModel
{
    // Protected data - all pyrolysis models inherit this stuff
    // Possible things to add include:
    //    pyrolysis gas molar mass
    protected:

        dictionary pyroModelDict_;

        pyrolysisModel(const pyrolysisModel&);

        void operator=(const pyrolysisModel&);

        solidThermo const* solidThermo_;

        dimensionedScalar Qc_; // J/kg

    public:

        //- Runtime type information
        TypeName("pyrolysisModel");

        // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            pyrolysisModel,
            components,
            (
                dictionary pyroDict,
                solidThermo const* solidThermo
            ),
            (
                pyroDict,
                solidThermo
            )
        );


        // Selectors

            //- Return a reference to the selected phaseChange model
            static autoPtr<pyrolysisModel> New
            (
                dictionary pyroDict,
                solidThermo const* solidThermo
            );

        // Constructors

            pyrolysisModel
            (
                const word& type,
                dictionary pyroDict,
                solidThermo const* solidThermo
            );

            virtual ~pyrolysisModel() {}


        // Member Functions

            virtual dimensionedScalar mass_burning_rate
            (
                const dimensionedScalar& T,
                const dimensionedScalar& p,
                const label& cell
            ) = 0;

            virtual dimensionedScalar energy_generation
            (
                const dimensionedScalar& T,
                const dimensionedScalar& p,
                const label& cell
            );

            virtual dimensionedScalar dmdT
            (
                const dimensionedScalar& T,
                const dimensionedScalar& p,
                const label& cell
            );

            dimensionedScalar Qc()
            {
                return Qc_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
