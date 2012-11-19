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
    volScalarField Yt(0.0*Y[0]);
    
    //Unity Schmidt number assumed (mu ~ rho*D)
    surfaceScalarField DEff = ib.interpolate(turbulence->muEff())
                               *ib.alphaf(); //unity Schmidt No.

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
                fvm::ddt(ib.alphaCorr()*rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(DEff, Yi)
              ==
                ib.alphaCorr()*combustion->R(Yi)
              + solid.YSu( Yi.name() )
              + fvm::Sp(solid.YSp(), Yi)
            );

            YiEqn.relax();
            YiEqn.solve(mesh.solver("Yi"));
            Yi.max(0.0);
                           
            Yt += Yi;
            solid.YSuC() = solid.YSu(Yi);
            solid.YSpC() = solid.YSp(Yi);
        }
    }
    
    thermo.correct();
    
    // Added to get sum(Yi) = 1
    //Yt += Y[inertIndex];

    forAll(Y, i)
    {
        Y[i] = Y[i] / Yt;
        
        Info<< "Min,max for " << Y[i].name() << " = "
            << Foam::min(Y[i]).value() << ", "
            << Foam::max(Y[i]).value() << endl;
    }
}
