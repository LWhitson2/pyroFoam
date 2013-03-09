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

#include "ellipsoid.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace shapes
{
    defineTypeNameAndDebug(ellipsoid, 0);
    addToRunTimeSelectionTable(shape, ellipsoid, components);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shapes::ellipsoid::ellipsoid
(
    const word& name,
    dictionary shapeDict,
    const fvMesh& mesh
)
:
    shape(typeName, name, shapeDict, mesh),
    center_(coeffDict_.lookup("center")),
    radius_(coeffDict_.lookup("radius"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::scalarField> Foam::shapes::ellipsoid::r() const
{
    vector sr = radius_;
        
    const volVectorField& cellCenters = mesh_.C();
    
    //vector Isr(1/(sr.x()+SMALL), 1/(sr.y()+SMALL), 1/(sr.z()+SMALL));
    
    if (sr.x() < SMALL) { sr.x() = GREAT; }
    if (sr.y() < SMALL) { sr.y() = GREAT; }
    if (sr.z() < SMALL) { sr.z() = GREAT; }
    
    tensor Isr(1/sr.x(), 0, 0, 
               0, 1/sr.y(), 0,
               0, 0, 1/sr.z());
    
    tmp<scalarField> r = Foam::mag(Isr & (cellCenters.internalField() - center_));
    
    return r;
}

void Foam::shapes::ellipsoid::calculate()
{
    //calculate liquid mask
    vector sr = radius_;
    const volVectorField& cellCenters = mesh_.C();
    if (sr.x() < SMALL) { sr.x() = GREAT; }
    if (sr.y() < SMALL) { sr.y() = GREAT; }
    if (sr.z() < SMALL) { sr.z() = GREAT; }
    
    tmp<scalarField> tr = r();
        
    forAll(liquidMask_, cellI)
    {
        vector cp = cellCenters[cellI] - center_;

        //scalar d = Foam::sqrt(1.0/( cp.x()*cp.x()/sr.x()/sr.x() + 
        //    cp.y()*cp.y()/sr.y()/sr.y() + cp.z()*cp.z()/sr.z()/sr.z() ));
            
        scalar d = 1.0/tr()[cellI];
            
        point p = center_ + d*cp;
        
        vector n(cp.x()/sr.x()/sr.x(),
                 cp.y()/sr.y()/sr.y(),
                 cp.z()/sr.z()/sr.z());
                 
        n /= mag(n);
        
        cuttableCell pc(mesh_, cellI);
        
        liquidMask_[cellI] = pc.cut( plane(p, n) );
    }

    liquidMask_.correctBoundaryConditions();
    
    //calculate vapor mask
    scalar d_layer = Foam::mag(radius_)/dV_;

    tmp<scalarField> f = Foam::max(1.0 - d_layer*(r() - 1.0), 0.0);
    
    vaporMask_.internalField() = f * (1.0 - liquidMask_);
    vaporMask_.correctBoundaryConditions();
}




// ************************************************************************* //
