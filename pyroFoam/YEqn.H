{
    reaction->correct();
    dQ = reaction->dQ();
    volScalarField Yt(0.0*Y[0]);

    // Calculate corrected species diffusivity
    surfaceScalarField DEff = ib.interpolate(turbulence->alphaEff(), "gas")
                               *ib.alphaf()*ib.lapCorr("gas"); // unity Lewis No.

    forAll(Y, i)
    {
        if (Y[i].name() == solid.solidName())
        {
            Y[i] = ib.solidCells();
            Y[i].max(0.0);
            Yt += Y[i];
        }
        else
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rhoAlpha, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(DEff, Yi)
              ==
                ib.alphaCorr()*reaction->R(Yi)
//                 combustion->R(Yi)
              + solid.YSu( Yi.name() )
              - fvm::Sp(solid.YSp(), Yi)
            );

            YiEqn.relax();
            YiEqn.solve(mesh.solver("Yi"));
            Yi.max(0.0);

            Yt += Yi;
        }
    }

    // Added to get sum(Yi) = 1
    //Yt += Y[inertIndex];

    forAll(Y, i)
    {
        Y[i] = Y[i] / Yt;

        Info<< "Min,max for " << Y[i].name() << " = "
            << Foam::min(Y[i]).value() << ", "
            << Foam::max(Y[i]).value() << endl;
    }

    Info<< "Min,max for Total Y = "
            << Foam::min(Yt).value() << ", "
            << Foam::max(Yt).value() << endl;

    thermo.correct();
}
