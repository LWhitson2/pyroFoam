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

Application
    reactingFoam

Description
    Solver for combustion with chemical reactions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "turbulenceModel.H"
#include "psiCombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "burningSolid.H"
#include "OFstream.H"

#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"
//     #include "diffusionNo.H"
    #include "setInitialDeltaT.H"

    //Time for fluid flow to run before enableing heat transfer to solid
    scalar flowRelaxTime =
        runTime.controlDict().lookupOrDefault<scalar>("flowRelaxTime", -1.0);

    // Overall conservation fields
    scalar coMass = 0.;
    vector coMom = vector(0., 0., 0.);
    scalar coEnergy = 0.;

    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    dimensionedScalar TsAvg = dimensionedScalar("TsAvg", dimTemperature, 300.0);

    while (runTime.run())
    {
        #include "readControls.H"
        #include "compressibleCourantNo.H"
        #include "setPyroDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        {
            // Store divU from the previous mesh for correctPhi
            volScalarField divU(fvc::div(phi));

            //Identify regions near ANY interface or reaction zone
            refinementField = ib.getRefinementField(U);

            // Do any mesh changes
            mesh.update();

            if (mesh.changing() && correctPhi)
            {
                ib.update();
                #include "correctPhi.H"
            }

            if (mesh.changing() && checkMeshCourantNo)
            {
                #include "meshCourantNo.H"
            }
        }

        // Evolve the solid surface
        solid.correct(U, phi, runTime.timeOutputValue() > flowRelaxTime);
        if (solid.testPyro() != "solid")
        {
            #include "rhoEqn.H"
        }

        while (pimple.loop())
        {
            if (solid.testPyro() == "solid")
            {
                #include "TsEqn.H"
            }
            else if (solid.testPyro() == "momentum")
            {
                #include "UEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }
            }
            else if (solid.testPyro() == "species")
            {
                #include "UEqn.H"
                #include "YEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }
            }
            else if (solid.testPyro() == "enthalpy")
            {
                #include "UEqn.H"
                #include "hsEqn.H"
                #include "TsEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }
            }
            else
            {
                #include "UEqn.H"
                #include "YEqn.H"
                #include "hsEqn.H"
                #include "TsEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }
            }

            // Info << nl << rho << nl <<endl;

            if (solid.testPyro() != "solid")
            {
                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }
            }
        }

        Tall = Ts*ib.alphas() + T*ib.alpha();

        heOut = thermo.he();

        runTime.write();

        // Record average Interface Temperature
        {
            dimensionedScalar smallA("dA", dimArea, SMALL);
            forAll(Y, i)
            {
                // Info << Y[i].name() << endl;
                if (Y[i].name() == "EMg")
                {
                    volScalarField& Yi = Y[i];
                    logTi << runTime.value() << token::TAB
                          << ib.alpha().weightedAverage(ib.area()).value() << token::TAB
                          << Ts.weightedAverage(ib.area()).value() << token::TAB
                          << solid.Ti().weightedAverage(ib.area()).value() << token::TAB
                          << T.weightedAverage(ib.area()).value() << token::TAB
                          << p.weightedAverage(ib.area()).value() << token::TAB
                          << Yi.weightedAverage(ib.area()).value() << token::TAB
                          << rho.weightedAverage(ib.area()).value() << token::TAB
                          << U.weightedAverage(ib.area()).component(1).value() << endl;
                }
            }
        }

        // #include "conservationCheck.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
