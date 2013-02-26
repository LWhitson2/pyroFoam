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
                dictionary pyroDict
            ),
            (pyroDict)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<pyrolysisModel> New
        (
            dictionary pyroDict
        );
        
    // Constructors
    
        pyrolysisModel
        (
            const word& type,
            dictionary pyroDict
        );
        
        virtual ~pyrolysisModel() {}


    // Member Functions
    
        virtual dimensionedScalar mass_burning_rate
        (
            const dimensionedScalar& T,
            const dimensionedScalar& p
        ) = 0;
        
        // Other possible functions include:
        //  Derivatives to help solving, dmdT, dmdP
        //  Heat generation rates
        //
        //

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
