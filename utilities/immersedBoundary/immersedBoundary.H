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
    Foam::immersedBoundary

Description


SourceFiles
    immersedBoundary.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundary_H
#define immersedBoundary_H

#include "volFields.H"
#include "dictionaryEntry.H"
#include "fvCFD.H"
#include "dimensionedScalar.H"
#include "dimensionedScalarFwd.H"
#include "cuttableCell.H"
#include "cuttableFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class phase Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundary
{
    // Private data

        //- Host mesh
        const fvMesh& mesh_;

        //- Parameter dictionary
        IOdictionary ibDict_;

        //- Volume fraction of gas
        volScalarField alpha_;

        volScalarField alphas_;

        volScalarField alphaCorr_;

        volScalarField alphasCorr_;

        //- Face area fraction of gas
        surfaceScalarField alphaf_;

        surfaceScalarField alphafs_;

        surfaceScalarField alphafCorr_;

        surfaceScalarField alphafsCorr_;


        volScalarField sumalphaf_;
        volScalarField sumalphafs_;

        //- Interface area
        volScalarField iArea_;

        //- Interface regression velocity
        //volVectorField iU_;

        //- Interface normal vector
        volVectorField iNormal_;

        //- Interface plane point
        volVectorField iPoint_;

        //- Centroid of gas portion of a cell (zero in solid)
        volVectorField gasC_;

        //- Centroid of solid portion of a cell (zero in gas)
        volVectorField solidC_;

        //- Distance from gas centroid to interface centroid (zero otherwise)
        volScalarField gasL_;

        //- Distance from solid centroid to interface centroid (zero otherwise)
        volScalarField solidL_;

        //- Mixed cells that are solid and gas
        volScalarField intermedsOut_;

        //- Alpha cutoff for small cell treatment
        dimensionedScalar alphaMin_;

        //- Alpha interface reconstruction cutoff
        dimensionedScalar reconstructTol_;

    // Private member functions
        //- Calculate the interface normals
        void calculateInterfaceNormal(const volScalarField& intermeds);

        //- Get the outward facing normal vector from faceI relative to cellI
        vector outwardNormal(label faceI, label cellI) const;

public:

    // Constructors

        //- Construct on host mesh
        immersedBoundary
        (
            const fvMesh& mesh
        );

        //- Return clone
        autoPtr<immersedBoundary> clone() const;

    // Public Member Functions

        const surfaceScalarField& alphaf() const
        {
            return alphaf_;
        }

        const surfaceScalarField& alphafs() const
        {
            return alphafs_;
        }

        const surfaceScalarField& alphafU() const
        {
            return alphafCorr_;
        }

        const surfaceScalarField& alphafsU() const
        {
            return alphafsCorr_;
        }

        const volScalarField& alpha() const
        {
            return alpha_;
        }

//         const volScalarField& alphas() const
//         {
//             return alphas_;
//         }
//
//         const volScalarField& alphaCorr() const
//         {
//             return alphaCorr_;
//         }
//
//         const volScalarField& alphasCorr() const
//         {
//             return alphasCorr_;
//         }
        const tmp<volScalarField> alphas() const
        {
            return (1. - alpha_);
        }

        const volScalarField& area() const
        {
            return iArea_;
        }

        const volVectorField& normal() const
        {
            return iNormal_;
        }

        const volVectorField& point() const
        {
            return iPoint_;
        }

        const volScalarField& gasL() const
        {
            return gasL_;
        }

        const volScalarField& solidL() const
        {
            return solidL_;
        }

        const volVectorField& gasC() const
        {
            return gasC_;
        }

        const volVectorField& solidC() const
        {
            return solidC_;
        }

        //- Get the "used" alpha (zero in small cells)
        tmp<volScalarField> alphaCorr() const;
        tmp<volScalarField> alphasCorr() const;

        //- Evolve the interface given a solid consumption rate
        void moveInterface(const volScalarField& ddtalpha);

        //- Calculate the interface parameters given alpha
        void correct();

        //- Update the interface after a mesh update
        void update();

        tmp<surfaceScalarField> lapCorr
        (
            const word& input
        );

        //- Get the mesh refinement criteria field
        tmp<volScalarField> getRefinementField
        (
            const volVectorField& U
        ) const;

        //- Calculate the small cell transfer weights
        tmp<surfaceScalarField> scTransferWeights(const word& input);

        //- Transfer a field from small cells to neighbour cells
        template<class Type>
        void transfer
        (
            const surfaceScalarField& w,
            const GeometricField<Type, fvPatchField, volMesh>& src,
            GeometricField<Type, fvPatchField, volMesh>& dest,
            Type zero,
            const word& input
        );

        //- Set source terms in small cells to set values in governing eqns
        template<class Type>
        void setScValue
        (
            const surfaceScalarField& w,
            GeometricField<Type, fvPatchField, volMesh>& Su,
            GeometricField<scalar, fvPatchField, volMesh>& Sp,
            const GeometricField<Type, fvPatchField, volMesh>& src,
//             const dimensioned<Type>& solidValue,
            const GeometricField<Type, fvPatchField, volMesh>& solidValue,
            const dimensionedScalar& diag,
            const word& mode,
            const word& input
        );

        //- Return a mask for small gas cells
        tmp<volScalarField> smallCells() const;

        //- Return a mask for full gas cells
        tmp<volScalarField> fullGasCells() const;

        //- Return a mask for small and solid cells
        tmp<volScalarField> smallAndSolidCells() const;

        //- Return a mask for small solid and gas cells
        tmp<volScalarField> smallSolidAndGasCells() const;

        //- Return a mask for solid cells
        tmp<volScalarField> solidCells() const;

        //- Return a mask for small solid cells
        tmp<volScalarField> smallSolidCells() const;

        //- Return a mask for gas cells
        tmp<volScalarField> gasCells() const;

        //- Return a mask for mixed cells
        tmp<volScalarField> mixedCells() const;

        //- Return a mask for no cells
        tmp<volScalarField> noCells() const;

        //- Return a mask for large solid cells
        tmp<volScalarField> largeSolidCells() const;

        //- IB-aware interpolation
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > ibInterpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;


        template<class Type>
        tmp
        <
            GeometricField
            <
                typename outerProduct<vector,Type>::type, fvPatchField, volMesh
            >
        >
        reconstruct
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
        ) const;

        template<class Type>
        tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > reconstruct
        (
            const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >&
        ) const;

        // Overloading

        inline tmp<surfaceScalarField> interpolate
        (
            const volScalarField& vf
        ) const
        {
            return ibInterpolate<scalar>(vf);
        }

        inline tmp<surfaceScalarField> interpolate
        (
            tmp<volScalarField> tvf
        ) const
        {
            return ibInterpolate<scalar>(tvf());
        }

        inline tmp<surfaceVectorField> interpolate
        (
            const volVectorField& vf
        ) const
        {
            return ibInterpolate<vector>(vf);
        }

        inline tmp<surfaceVectorField> interpolate
        (
            tmp<volVectorField> tvf
        ) const
        {
            return ibInterpolate<vector>(tvf());
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "immersedBoundaryTemplates.H"

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
