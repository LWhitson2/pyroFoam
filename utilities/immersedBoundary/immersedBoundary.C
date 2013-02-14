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

#include "immersedBoundary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundary::immersedBoundary
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    ibDict_
    (
        IOobject
        (
            "immersedBoundaryDict",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    alpha_
    (
        IOobject
        (
            "alphaGas",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    alphaf_
    (
        IOobject
        (
            "alphaf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphaf", dimless, 0.0)
    ),

    alphafs_
    (
        IOobject
        (
            "alphafs",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphafs", dimless, 0.0)
    ),

    sumalphaf_
    (
        IOobject
        (
            "sumalphaf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("sumalphaf",dimless,0.0)
    ),

    sumalphafs_
    (
        IOobject
        (
            "sumalphafs",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("sumalphafs",dimless,0.0)
    ),

    iArea_
    (
        IOobject
        (
            "iArea",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("iArea",dimArea, 0.0)
    ),

    iNormal_
    (
        IOobject
        (
            "iNormal",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("iNormal", dimless, vector::zero)
    ),

    iPoint_
    (
        IOobject
        (
            "iPoint",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("iPoint", dimless, vector::zero)
    ),

    gasC_
    (
        IOobject
        (
            "gasC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("gasC", dimless, vector::zero)
    ),

    solidC_
    (
        IOobject
        (
            "solidC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("solidC", dimless, vector::zero)
    ),

    intermedsOut_
    (
        IOobject
        (
            "intermedsOut",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("iNormal", dimless, 0)
    ),

    alphaMin_(ibDict_.lookup("alphaMin")),
    reconstructTol_(ibDict_.lookup("reconstructTol"))

{
    Foam::Info << "Created immersed boundary" << Foam::endl;

    // Calculate new interface position and fields
    correct();

    alpha_.oldTime();
    alphaf_.oldTime();
    iNormal_.oldTime();
    iPoint_.oldTime();
    iArea_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::immersedBoundary> Foam::immersedBoundary::clone() const
{
    notImplemented("immersedBoundary::clone() const");
    return autoPtr<immersedBoundary>(NULL);
}

// Perform a correct update of the interface after the mesh is adapted
void Foam::immersedBoundary::update()
{
    //Mesh update will simply apply alpha of the parent cell to its child cells
    // which is incorrect for the sharp interface. Fortunately, it will also
    // write iNormal and iPoint to each child cell so that alpha can be
    // recalculated easily.

    // When coarsening, this gets more complicated. To alleviate this, we will
    // require the boundary to always be refined to the maximum level. It can
    // only unrefine when the boundary has moved away.

    forAll(alpha_, cellI)
    {
        if (mag(iNormal_[cellI]) > SMALL && mag(iPoint_[cellI]) > SMALL)
        {
            cuttableCell cc(mesh_, cellI);
            plane p(iPoint_[cellI],iNormal_[cellI]);
            alpha_[cellI] = 1.0 - cc.cut(p);
        }
    }
    alpha_.correctBoundaryConditions();

    correct();
}


void Foam::immersedBoundary::calculateInterfaceNormal
(
    const volScalarField& intermeds
)
{
    //Use a smoothing procedure to capture the interface better

    // Get gradient
    iNormal_ = fvc::grad(alpha_)*dimensionedScalar("one",dimLength,1.0);

    surfaceVectorField iNormalf = interpolate(iNormal_);
    iNormal_ = 0.7*fvc::average(iNormalf) + 0.3*iNormal_;

    for (label i = 0; i < 1; ++i)
    {
        iNormalf = interpolate(iNormal_);
        iNormal_ = 0.7*fvc::average(iNormalf) + 0.3*iNormal_;
    }

    /*
    const labelListList& cellCells = mesh_.cellCells();
    const scalarField& cellVolumes = mesh_.V();

    // Do spatial smoothing of gradient by volume-weighted local averaging
    forAll(mesh_.cells(), cellI)
    {
        vector normSum = iNormal_[cellI]*cellVolumes[cellI];
        scalar Vtot = cellVolumes[cellI];

        const labelList& nbList = cellCells[cellI];

        forAll(nbList, nbI)
        {
            normSum += iNormal_[ nbList[nbI] ]*cellVolumes[ nbList[nbI] ];
            Vtot += cellVolumes[ nbList[nbI] ];
        }

        iNormal_[cellI] = normSum / Vtot;
    }
    */

    // Normalize and limit iNormal only to intermediate cells
    iNormal_ *= intermeds / (mag(iNormal_) + SMALL);
    iNormal_.correctBoundaryConditions();
}


// Get the outward facing normal vector on faceI relative to cellI
Foam::vector Foam::immersedBoundary::outwardNormal
(
    label faceI,
    label cellI
) const
{
    const face& f = mesh_.faces()[faceI];
    const pointField& points = mesh_.points();

    // Calculate face normal
    vector norm = f.normal(points);
    norm /= mag(norm);

    // Flip norm if pointed wrong way
    vector vSF = f.centre(points) - mesh_.C()[cellI];

    if ((vSF & norm) < 0.0)
    {
        norm *= -1.0;
    }

    return norm;
}


// Given alpha, calculate derived interface fields
// LIMITATIONS: SERIAL RUNS ONLY
void Foam::immersedBoundary::correct()
{
    // Step 1: Identify intermediate cells based on alpha_ and reconstructTol_
    volScalarField intermeds = pos(alpha_ - reconstructTol_)
                               *pos(1.0 - reconstructTol_ - alpha_);

    // Step 2: Calculate interface normal in intermediate cells
    calculateInterfaceNormal(intermeds);
    intermedsOut_ = intermeds;

    // Step 3: Calculate the cut plane and cut area in intermediate cells
    //         setting a_burn to zero in all non-intermediate cells
    forAll(iNormal_, cellI)
    {
        if (intermeds[cellI] > SMALL)
        {
            cuttableCell cc(mesh_, cellI);
            plane p = cc.constructInterface(iNormal_[cellI],1.0-alpha_[cellI]);
            iPoint_[cellI] = p.refPoint();
            iArea_[cellI] = cc.cutArea();

            // Save gas and solid portion centroids
            gasC_[cellI] = cc.lostCentroid();
            solidC_[cellI] = cc.cutCentroid();

            // DEBUGGING PURPOSES
            if (iArea_[cellI] < SMALL)
            {
                Info<< "WARNING: Cut area is " << iArea_[cellI]
                    << " with alphaSolid = " << 1.0-alpha_[cellI] << endl;
            }
        }
        else
        {
            iArea_[cellI] = 0.0;
            iPoint_[cellI] = vector::zero;
            if (alpha_[cellI] > 0.5)
            { //gas cell
                gasC_[cellI] = mesh_.C()[cellI];
                solidC_[cellI] = vector::zero;
            }
            else
            { //solid cell
                gasC_[cellI] = vector::zero;
                solidC_[cellI] = mesh_.C()[cellI];
            }
        }
    }

    iPoint_.correctBoundaryConditions();
    iNormal_.correctBoundaryConditions();


    // Step 4: Calculate alphaf on all faces, valid only in the homogeneous
    //         regions away from the interface
    alphaf_ = fvc::interpolate(alpha_);

    // Step 5: Use the planes in intermediate cells to correct alphaf
    //         near the interface. Also catch solid cells that have a partial
    //         burning face as calculated from an intermediate cell cut plane.
    surfaceScalarField hasIntermeds = fvc::interpolate(intermeds);

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbor = mesh_.neighbour();

    // This loop should be modified for parallel cases. This should be doable
    forAll(alphaf_, faceI)
    {
        label own = owner[faceI];
        label nei = neighbor[faceI];

        // Faces where there is a cut plane in one or both bounding cells
        if (hasIntermeds[faceI] > SMALL)
        {
            cuttableFace cf(mesh_, faceI);

            scalar alphafNei = 1.0;
            scalar alphafOwn = 1.0;

            if (mag(iNormal_[nei]) > SMALL)
            {
                Foam::plane p(iPoint_[nei], iNormal_[nei]);
                alphafNei = cf.cut(p);
            }

            if (mag(iNormal_[own]) > SMALL)
            {
                Foam::plane p(iPoint_[own], iNormal_[own]);
                alphafOwn = cf.cut(p);
            }

            alphaf_[faceI] = Foam::min(alphafNei, alphafOwn);

            // Catch the cases where:
            //
            //   +---------+
            //   |\        |
            //   | \   g   |
            //   |s \      |
            //   +---======+  <- Face in question
            //   |         |
            //   |    s    |
            //   |         |
            //   +---------+
            //
            //  The value of alphaf (===) for the face is calculated correctly,
            //  but needs to be added to the interface area of the solid cell

            if (alpha_[own] < reconstructTol_.value()
                && alphaf_[faceI] > SMALL)
            {
                //own is solid
                iArea_[own] += mesh_.magSf()[faceI] * alphaf_[faceI];
                iNormal_[own] += outwardNormal(faceI, own)*alphaf_[faceI];
            }
            else if (alpha_[nei] < reconstructTol_.value()
                     && alphaf_[faceI] > SMALL)
            {
                //nei is solid
                iArea_[nei] += mesh_.magSf()[faceI] * alphaf_[faceI];
                iNormal_[nei] += outwardNormal(faceI, nei)*alphaf_[faceI];
            }
        }

        // Catch sharp face-coincident interfaces. In this case, the solid cell
        // is the one that will be burning at the next time step.
        else if (mag(alpha_[own] - alpha_[nei]) > 0.1)
        {
            alphaf_[faceI] = 0.0;

            //set iNormal and a_burn for this case
            label solidcell = (alpha_[own] < 0.1) ? own : nei;

            // Increment a_burn since a solidcell could technically have more
            // than one sharp boundary
            iArea_[solidcell] += mesh_.magSf()[faceI];

            // Increment iNormal_
            iNormal_[solidcell] += outwardNormal(faceI, solidcell);
        }
    }

    // Now set alphaf on parallel patches
    const volScalarField::GeometricBoundaryField& alphaBf =
        alpha_.boundaryField();
    const volVectorField::GeometricBoundaryField& iPointBf =
        iPoint_.boundaryField();
    volVectorField::GeometricBoundaryField& iNormalBf =
        iNormal_.boundaryField();
    volScalarField::GeometricBoundaryField& iAreaBf =
        iArea_.boundaryField();
    const surfaceScalarField::GeometricBoundaryField& hasIntermedsBf =
        hasIntermeds.boundaryField();

    surfaceScalarField::GeometricBoundaryField& alphafBf =
        alphaf_.boundaryField();

    forAll(alphafBf, patchI)
    {
        const fvPatchScalarField& alphaPf = alphaBf[patchI];
        const fvPatchVectorField& iPointPf = iPointBf[patchI];

        fvPatchVectorField& iNormalPf = iNormalBf[patchI];
        fvPatchScalarField& iAreaPf = iAreaBf[patchI];

        const scalarField& hasIntermedsPf = hasIntermedsBf[patchI];
        scalarField& alphafPf = alphafBf[patchI];

        const labelList& pFaceCells = mesh_.boundary()[patchI].faceCells();

        if (alphaPf.coupled()) //returns true for parallel and cyclic patches
        {
            // Get values across parallel patch
            const vectorField iPointPNf(iPointPf.patchNeighbourField());
            const vectorField iNormalPNf(iNormalPf.patchNeighbourField());
            const scalarField alphaPNf(alphaPf.patchNeighbourField());
            const scalarField iAreaPNf(iAreaPf.patchNeighbourField());

            //patch face starting IDs in mesh.faces()
            label patchFs = alphaPf.patch().start();
            const fvPatch& meshPf = mesh_.boundary()[patchI];

            forAll(alphafPf, pFaceI)
            {
                label pfCellI = pFaceCells[pFaceI];

                if (hasIntermedsPf[pFaceI] > SMALL)
                {
                    cuttableFace cf(mesh_, patchFs+pFaceI);

                    scalar alphafNei = 1.0;
                    scalar alphafOwn = 1.0;

                    if (mag(iNormalPNf[pFaceI]) > SMALL)
                    {
                        Foam::plane p(iPointPNf[pFaceI], iNormalPNf[pFaceI]);
                        alphafNei = cf.cut(p);
                    }

                    if (mag(iNormal_[pfCellI]) > SMALL)
                    {
                        Foam::plane p(iPoint_[pfCellI], iNormal_[pfCellI]);
                        alphafOwn = cf.cut(p);
                    }

                    alphafPf[pFaceI] = Foam::min(alphafNei, alphafOwn);

                    // now catch sharp edges
                    if (alpha_[pfCellI] < reconstructTol_.value()
                        && alphafPf[pFaceI] > SMALL)
                    {
                        //own is solid
                        iArea_[pfCellI] += meshPf.magSf()[pFaceI]
                                           *alphafPf[pFaceI];

                        iNormal_[pfCellI] +=
                                outwardNormal(patchFs+pFaceI, pfCellI)
                                *alphafPf[pFaceI];
                    }

                }
                else if (mag(alpha_[pfCellI] - alphaPNf[pFaceI]) > 0.1)
                {
                    alphafPf[pFaceI] = 0.0;

                    //set iNormal and a_burn for this case
                    if (alpha_[pfCellI] < 0.1)
                    {
                        // than one sharp boundary
                        iArea_[pfCellI] += meshPf.magSf()[pFaceI];

                        // Increment iNormal_
                        iNormal_[pfCellI] +=
                            outwardNormal(patchFs+pFaceI, pfCellI);
                    }
                }
            }
        }
    }

    //Re-normalize iNormal (only needed for the cases when it is incremented)
    iNormal_ /= (mag(iNormal_) + VSMALL);

    sumalphaf_ = fvc::surfaceSum(alphaf_);
    alphafs_ = alphaf_;

    iArea_.correctBoundaryConditions();
    iNormal_.correctBoundaryConditions();
    gasC_.correctBoundaryConditions();
    solidC_.correctBoundaryConditions();
}



// Calculate alpha that excludes small cells
Foam::tmp<Foam::volScalarField> Foam::immersedBoundary::alphaCorr() const
{
    return alpha_ * pos(alpha_ - alphaMin_);
}

Foam::tmp<Foam::volScalarField> Foam::immersedBoundary::alphasCorr() const
{
    return alphas() * pos(alphas() - alphaMin_);
}


//If no dir is available, this could use iNormal_ instead
Foam::tmp<Foam::surfaceScalarField>
Foam::immersedBoundary::scTransferWeights(const word& input)
{
    tmp<surfaceScalarField> tw
    (
        new surfaceScalarField
        (
            IOobject
            (
                "tw",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("tw", dimless, 0.0)
        )
    );

    surfaceScalarField& w = tw();

    const tmp<volScalarField> alphastmp = 1.0 - alpha_;
    const volScalarField& alpha = (input == "gas") ? alpha_:alphastmp();
    surfaceScalarField& alphaf = (input == "gas") ? alphaf_:alphafs_;
    volScalarField& sumalphaf = (input == "gas") ? sumalphaf_:sumalphafs_;
    scalar iFlip = (input=="gas") ? 1 : -1;

    // If negative, alpha is a small cell
    volScalarField alphaShift = alpha - alphaMin_;

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbor = mesh_.neighbour();

    // Calculate transfer weights
    forAll(w, faceI)
    {
        label own = owner[faceI];
        label nei = neighbor[faceI];

        if (alphaShift[own] * alphaShift[nei] * alphaf[faceI] < 0.0)
        { //one is small, one is not, and they share a gas boundary

            label sc = (alphaShift[own] < 0.0) ? own : nei;

            w[faceI] = mag
            (
                iNormal_[sc] & mesh_.Sf()[faceI]
            ) * alphaf[faceI] * iFlip;

            alphaf[faceI] = 0.0;
        }

        // if both cells are small, zero out alphaf between then
        if (alphaShift[own] < 0.0 && alphaShift[nei] < 0.0)
        {
            alphaf[faceI] = 0.0;
        }
    }



    // Now set alphaf on parallel patches
    const volScalarField::GeometricBoundaryField& alphaShiftBf =
        alphaShift.boundaryField();
    const volVectorField::GeometricBoundaryField& iNormalBf =
        iNormal_.boundaryField();


    surfaceScalarField::GeometricBoundaryField& wBf =
        w.boundaryField();
    surfaceScalarField::GeometricBoundaryField& alphafBf =
        alphaf.boundaryField();

    forAll(alphafBf, patchI)
    {
        const fvPatchScalarField& alphaShiftPf = alphaShiftBf[patchI];
        const fvPatchVectorField& iNormalPf = iNormalBf[patchI];

        scalarField& alphafPf = alphafBf[patchI];
        scalarField& wPf = wBf[patchI];

        const labelList& pFaceCells = mesh_.boundary()[patchI].faceCells();

        if (alphaShiftPf.coupled())
        {
            // Get values across parallel patch
            const scalarField alphaShiftPNf(alphaShiftPf.patchNeighbourField());
            const vectorField iNormalPNf(iNormalPf.patchNeighbourField()*iFlip);

            const fvPatch& meshPf = mesh_.boundary()[patchI];

            forAll(alphafPf, pFaceI)
            {
                label pfCellI = pFaceCells[pFaceI];

                //Calculate weight on small-large cell faces and close
                // small cell faces
                if ( alphaShift[pfCellI] * alphaShiftPNf[pFaceI]
                     * alphafPf[pFaceI] < 0.0)
                {

                    vector scNorm = (alphaShift[pfCellI] < 0.0)
                                    ? iFlip*iNormal_[pfCellI]
                                    : iNormalPNf[pFaceI];

                    wPf[pFaceI] = mag
                    (
                        scNorm & meshPf.Sf()[pFaceI]
                    ) * alphafPf[pFaceI];

                    alphafPf[pFaceI] = 0.0;
                }

                //Close faces between small cells
                if (alphaShift[pfCellI] < 0.0 && alphaShiftPNf[pFaceI] < 0.0)
                {
                    alphafPf[pFaceI] = 0.0;
                }
            }
        }
    }

    sumalphaf = fvc::surfaceSum(alphaf);

    return tw;
}




Foam::tmp<Foam::volScalarField>
Foam::immersedBoundary::getRefinementField
(
    const volVectorField& U
) const
{
    // Force all cells near the interface to refine to the maximum level
    dimensionedScalar C("C",dimLength,1e6);
    tmp<volScalarField> tRefinementField = C*mag(fvc::grad(alpha_));

    //Include curl criteria from Popinet (Gerris), scaled by 0.5, to also
    // refine key fluid flow regions
    tRefinementField().internalField() = max
    (
        tRefinementField().internalField(),
        Foam::mag(fvc::curl(U)) * Foam::pow(mesh_.V(),1.0/3.0) * 0.5
    );

    return tRefinementField;
}


void Foam::immersedBoundary::moveInterface(const volScalarField& ddtalpha)
{
    //Update alpha
    solve(fvm::ddt(alpha_) == ddtalpha);

    //Correct the interface parameters using the new alpha
    correct();
}



Foam::tmp<Foam::volScalarField>
Foam::immersedBoundary::smallAndSolidCells() const
{
    return neg(alpha_ - alphaMin_);
}

Foam::tmp<Foam::volScalarField>
Foam::immersedBoundary::smallSolidAndGasCells() const
{
    return neg(alphas() - alphaMin_);
}

Foam::tmp<Foam::volScalarField> Foam::immersedBoundary::fullGasCells() const
{
    return neg(alphas() - SMALL);
}

Foam::tmp<Foam::volScalarField> Foam::immersedBoundary::smallCells() const
{
    return neg(alpha_ - alphaMin_)*pos(alpha_ - SMALL);
}


Foam::tmp<Foam::volScalarField> Foam::immersedBoundary::solidCells() const
{
    return neg(alpha_ - SMALL);
}

Foam::tmp<Foam::volScalarField> Foam::immersedBoundary::smallSolidCells() const
{
    return neg(alphas() - alphaMin_)*pos(alphas() - SMALL);
}

Foam::tmp<Foam::volScalarField> Foam::immersedBoundary::gasCells() const
{
    return pos(alpha_ - alphaMin_);
}

Foam::tmp<Foam::volScalarField> Foam::immersedBoundary::noCells() const
{
    return neg(alpha_ + 100.0);
}


// ************************************************************************* //
