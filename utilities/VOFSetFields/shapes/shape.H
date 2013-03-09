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
    Abstract shape class

SourceFiles
    shape.C

\*---------------------------------------------------------------------------*/

#ifndef shape_H
#define shape_H

#include "volFields.H"
#include "PtrDictionary.H"
#include "dictionaryEntry.H"
#include "fvCFD.H"
#include "dimensionedScalar.H"
#include "dimensionedScalarFwd.H"
#include "cuttableCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class droplet Declaration
\*---------------------------------------------------------------------------*/

class shape
{
    // Protected data
protected:

        const fvMesh& mesh_;

        //- Shape name
        word name_;
        
        //- Shape dictionary
        dictionary shapeDict_;
        dictionary coeffDict_;

        scalar dV_;
        vector Uinit_;
        
        HashTable<scalar> liquidSpecies_;
        HashTable<scalar> vaporSpecies_;
        
        volScalarField liquidMask_; //TODO: Change to scalarField
        volScalarField vaporMask_; //TODO: Change to scalarField
        
        
        shape(const shape&);
        
        void operator=(const shape&);
        
public:
    
        TypeName("shape");
    
        declareRunTimeSelectionTable
        (
            autoPtr,
            shape,
            components,
            (
                const word& name,
                dictionary shapeDict,
                const fvMesh& mesh
            ),
            (name, shapeDict, mesh)
        );
        
    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<shape> New
        (
            const word& name,
            dictionary shapeDict,
            const fvMesh& mesh
        );
        
    // Constructors

        //- Construct from components
        shape
        (
            const word& type,
            const word& name,
            dictionary shapeDict,
            const fvMesh& mesh
        );

        virtual ~shape() {}
        
        //- Return clone
        autoPtr<shape> clone() const;

        //- Return a pointer to a new shape created on freestore
        //  from Istream
        class iNew
        {
            const fvMesh& mesh_;

        public:

            iNew
            (
                const fvMesh& mesh
            )
            :
                mesh_(mesh)
            {}

            autoPtr<shape> operator()(Istream& is) const
            {
                dictionaryEntry ent(dictionary::null, is);
                return autoPtr<shape>(shape::New(ent.keyword(), ent, mesh_));
            }
        };


    // Member Functions

        const dictionary& dict() const
        {
            return shapeDict_;
        }

        const word& name() const
        {
            return name_;
        }

        const word& keyword() const
        {
            return name();
        }
        
        const volScalarField& liquidMask() const
        {
            return liquidMask_;
        }
        
        const volScalarField& vaporMask() const
        {
            return vaporMask_;
        }
        
        virtual void set
        (
            volScalarField& alphaLiquid,
            volVectorField& U,
            PtrList<volScalarField>& species
        );
        
        void calcDropMask();
        
        List<word> species() const;

        virtual void calculate() =0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
