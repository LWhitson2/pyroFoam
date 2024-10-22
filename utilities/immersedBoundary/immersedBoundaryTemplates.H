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

#ifndef immersedBoundaryTemplates_H
#define immersedBoundaryTemplates_H

template<class Type>
void immersedBoundary::transfer
(
    const surfaceScalarField& w,
    const GeometricField<Type, fvPatchField, volMesh>& src,
    GeometricField<Type, fvPatchField, volMesh>& dest,
    Type zero,
    const word& input
)
{
    const tmp<volScalarField> alphastmp = 1.0 - alpha_;
    const volScalarField& alpha = (input == "gas") ? alpha_:alphastmp();
    const volScalarField::DimensionedInternalField& Vc = mesh_.V();

    volScalarField wtot = fvc::surfaceSum(w);


    const labelUList& owner = mesh_.owner();
    const labelUList& neighbor = mesh_.neighbour();

    forAll(w, faceI)
    {
        if (w[faceI] > 0.0)
        {
            label own = owner[faceI];
            label nei = neighbor[faceI];

            label sc = (alpha[own] < alphaMin_.value()) ? own : nei;
            label rc = (sc == own) ? nei : own;

            dest[rc] += w[faceI]/wtot[sc] * src[sc] * Vc[sc]/Vc[rc];
            dest[sc] = zero;
        }
    }

    // Now do transfer at parallel patches
    // TODO: Add V(src)/V(dest) to the final calculation.
    const typename GeometricField<Type, fvPatchField, volMesh>::GeometricBoundaryField&
        srcBf = src.boundaryField();

    const volScalarField::GeometricBoundaryField&
        wtotBf = wtot.boundaryField();

    const surfaceScalarField::GeometricBoundaryField&
        wBf = w.boundaryField();

    forAll(srcBf, patchI)
    {
        const fvPatchField<Type>& srcPf = srcBf[patchI];
        const fvPatchScalarField& wtotPf = wtotBf[patchI];
        const scalarField& wPf = wBf[patchI];

        const labelList& pFaceCells = mesh_.boundary()[patchI].faceCells();

        if (srcPf.coupled())
        {
            // Get values across parallel patch
            const Field<Type> srcPNf(srcPf.patchNeighbourField());
            const scalarField wtotPNf(wtotPf.patchNeighbourField());

            forAll(srcPf, pFaceI)
            {
                label pfCellI = pFaceCells[pFaceI];

                if ( wPf[pFaceI] > 0.0)
                {
                    if (alpha[pfCellI] < alphaMin_.value())
                    {//I own the small cell
                        dest[pfCellI] = zero;
                    }
                    else
                    {//I own the receiving cell
                        dest[pfCellI] += wPf[pFaceI]/wtotPNf[pFaceI]
                                          * srcPNf[pFaceI];
                    }
                }
            }
        }
    }
}


template<class Type>
void immersedBoundary::setScValue
(
    const surfaceScalarField& w,
    GeometricField<Type, fvPatchField, volMesh>& Su,
    GeometricField<scalar, fvPatchField, volMesh>& Sp,
    const GeometricField<Type, fvPatchField, volMesh>& src,
    const GeometricField<Type, fvPatchField, volMesh>& solidValue,
    const dimensionedScalar& diag,
    const word& mode,
    const word& input
)
{
    tmp<volScalarField> alphastmp = 1.0 - alpha_;
    const volScalarField& alpha = (input == "gas") ? alpha_:alphastmp();

    Su = ((input == "gas") ? solidCells():fullGasCells())*solidValue*diag;
    Sp = ((input == "gas") ? smallAndSolidCells():smallSolidAndGasCells())*diag;

    if (mode == "fix")
    {
        Su += ((input == "gas") ? smallCells():smallSolidCells())*src*diag;
    }
    else if (mode == "avg")
    {
        volScalarField wtot = fvc::surfaceSum(w);

        const labelUList& owner = mesh_.owner();
        const labelUList& neighbor = mesh_.neighbour();

        forAll(w, faceI)
        {
            if (w[faceI] > SMALL)
            {
                label own = owner[faceI];
                label nei = neighbor[faceI];

                label sc = (alpha[own] < alphaMin_.value()) ? own : nei;
                label rc = (sc == own) ? nei : own;

                if (alpha[sc] > SMALL)
                {
                    Su[sc] += w[faceI]/wtot[sc] * src[rc] * diag.value();
                }
            }
        }

        // Now do transfer at parallel patches
        // const typename GeometricField
        //     <Type, fvPatchField, volMesh>::GeometricBoundaryField&
        //     srcBf = src.boundaryField();

        // const surfaceScalarField::GeometricBoundaryField&
        //     wBf = w.boundaryField();

        // forAll(srcBf, patchI)
        // {
        //     const fvPatchField<Type>& srcPf = srcBf[patchI];
        //     const scalarField& wPf = wBf[patchI];

        //     const labelList& pFaceCells = mesh_.boundary()[patchI].faceCells();

        //     if (srcPf.coupled())
        //     {
        //         // Get values across parallel patch
        //         const Field<Type> srcPNf(srcPf.patchNeighbourField());

        //         forAll(srcPf, pFaceI)
        //         {
        //             label pfCellI = pFaceCells[pFaceI];

        //             if ( wPf[pFaceI] > 0.0)
        //             {
        //                 if (alpha[pfCellI] < alphaMin_.value())
        //                 {//I own the small cell
        //                     Su[pfCellI] += wPf[pFaceI]/wtot[pfCellI]
        //                                   * srcPNf[pFaceI] * diag.value();
        //                 }
        //             }
        //         }
        //     }
        // }
    }
    else
    {
        FatalErrorIn
        (
            "immersedBoundary::setScValue()"
        )   << "Setting mode  '" << mode
            << "' not recognized"
            << exit(FatalError);
    }
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
immersedBoundary::ibInterpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& input
) const
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
    (
        fvc::interpolate(vf)
    );
    
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf();
       
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbor = mesh_.neighbour();

    // Internal faces first
    forAll(alphaf_, faceI)
    {
        label own = owner[faceI];
        label nei = neighbor[faceI];

        // Check if either cell has a cut plane in it
        if( mag(iPoint_[own]) > SMALL || mag(iPoint_[nei]) > SMALL )
        {
            // get the proper centroid from each cell
            vector ownC = (mag(iPoint_[own]) > SMALL) ?
                        ( (input == "gas") ? gasC_[own] : solidC_[own] )
                      : ( mesh_.C()[own] );


            vector neiC = (mag(iPoint_[nei]) > SMALL) ?
                        ( (input == "gas") ? gasC_[nei] : solidC_[nei] )
                      : ( mesh_.C()[nei] );
            
            scalar deltaOwn = mag((ownC - mesh_.Cf()[faceI]) & mesh_.Sf()[faceI]) / mesh_.magSf()[faceI];
            scalar deltaNei = mag((neiC - mesh_.Cf()[faceI]) & mesh_.Sf()[faceI]) / mesh_.magSf()[faceI];
            
            sf[faceI] = (deltaNei*vf[own] + deltaOwn*vf[nei])/(deltaOwn+deltaNei);
        }
    }

    //Repeat for parallel patches (and boundary patches?)
    
    return tsf;
}

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, fvPatchField, volMesh
    >
>
immersedBoundary::reconstruct
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<GeometricField<GradType, fvPatchField, volMesh> > treconField
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                "volIntegrate("+ssf.name()+')',
                ssf.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            inv(fvc::surfaceSum(sqr(mesh_.Sf()*(alphaf_ + SMALL))/(mesh_.magSf())))
          & fvc::surfaceSum((mesh_.Sf()/mesh_.magSf())*ssf),
            zeroGradientFvPatchField<GradType>::typeName
        )
    );

    treconField().correctBoundaryConditions();

    return treconField;
}

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
immersedBoundary::reconstruct
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tssf
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    tmp<GeometricField<GradType, fvPatchField, volMesh> > tvf
    (
        reconstruct(tssf())
    );
    tssf.clear();
    return tvf;
}

#endif
