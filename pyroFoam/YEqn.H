tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    combustion->correct();
    dQ = combustion->dQ();
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (Y[i].name() == "EMs")
        {
            Y[i] = neg(solid.alpha() - SMALL);
            Y[i].max(0.0);
            Yt += Y[i];
        }
        else if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
                (
                    fvm::ddt(solid.alphaUsed()*rho, Yi)
                +   mvConvection->fvmDiv(phi, Yi)
                // - fvm::laplacian(turbulence->muEff(), Yi)
                ==
                    solid.alphaUsed()*combustion->R(Yi)
                +   solid.YSu(Yi)
                +   fvm::Sp(solid.YSp(Yi), Yi)
                );

            YiEqn.relax();
            YiEqn.solve(mesh.solver("Yi"));
            Yi.max(0.0);
            Yt += Yi;
            solid.YSuC() = solid.YSu(Yi);
            solid.YSpC() = solid.YSp(Yi);
        }
        else
        {
            inertIndex = i;
        }
    }

    Y[inertIndex] = scalar(1.0) - Yt;
    Y[inertIndex].max(0.0);

    // Added to get sum(Yi) = 1
    /*Yt += Y[inertIndex];

    forAll(Y, i)
    {
        Y[i] = Y[i] / Yt;
    }*/
}
