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

Global
    rhoEqn

Description
    Solve the continuity for density.

\*---------------------------------------------------------------------------*/

{

    Info<< "Solving rhoEqn" << endl;

    solve
    (
        fvm::ddt(ib.alphaCorr(), rho)
      + fvc::div(phi)
     ==
        solid.rhoSu()
      - fvm::Sp(solid.rhoSp(), rho)
      + solid.mgen()
    );

    // Info<< rho[37] << endl;
    // Info<< ib.alphaCorr().oldTime()[37] << endl;
    // Info<< ib.alphaCorr()[37] << endl;
    // Info << rho.oldTime()[37] << endl;
    // Info << fvc::div(phi)()[37] << endl;
    // Info<< solid.mgen()[37] << endl;
    // Info<< "Min,max rho = " << Foam::min(rho).value()
    //     << ", " << Foam::max(rho).value() << endl;
    // Info<< "Max alphaDiff = " << Foam::max(Foam::mag(
    //   ib.alphaCorr().oldTime() - ib.alpha().oldTime())).value() << endl;
    // Info<< "Min,max alphaCorr.oldTime = " << Foam::min(ib.alphaCorr().oldTime()).value()
    //     << ", " << Foam::max(ib.alphaCorr().oldTime()).value() << endl;
    // Info<< "Min,max ddtalpha = " << Foam::min(fvc::ddt(ib.alpha())).value()
    //     << ", " << Foam::max(fvc::ddt(ib.alpha())).value() << endl;
    // Info<< "Min,max ddtalphaCorr = " << Foam::min(fvc::ddt(ib.alphaCorr())).value()
    //      << ", " << Foam::max(fvc::ddt(ib.alphaCorr())).value() << endl;
    // Info<< "Min,max ddtrho = " << Foam::min(fvc::ddt(rho)).value()
    //     << ", " << Foam::max(fvc::ddt(rho)).value() << endl;
    // Info<< "Min,max ddtrhoalpha = " << Foam::min(fvc::ddt(ib.alpha(), rho)).value()
    //     << ", " << Foam::max(fvc::ddt(ib.alpha(), rho)).value() << endl;
    // Info<< "Min,max ddtrhoalphaCorr = " << Foam::min(fvc::ddt(ib.alphaCorr(), rho)).value()
    //     << ", " << Foam::max(fvc::ddt(ib.alphaCorr(), rho)).value() << endl;
}

// ************************************************************************* //
