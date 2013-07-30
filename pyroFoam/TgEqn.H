{
    // Calculate corrected thermal conductivity
    surfaceScalarField KbyCp = ib.interpolate(Ks)
                               *ib.alphafs(); //TODO: edit

    fvScalarMatrix TgEqn
    (
        fvm::ddt(ib.alphaCorr()*rho, Tg)
      + mvConvection->fvmDiv(phi, Tg)
      - fvm::Sp(fvc::ddt(ib.alphaCorr()*rho) + fvc::div(phi), Tg)
      - fvm::laplacian(KbyCp, Tg)
     ==
        solid.rCp()*
        (
            ib.alphaCorr()*dpdt
          - fvc::ddt(ib.alphaCorr()*rho, K)
          - fvc::div(phi, K)
          + (fvc::ddt(ib.alphaCorr()*rho) + fvc::div(phi))*K
          + ib.alphaCorr()*combustion->Sh()
        )
      + solid.TgSu()
      - fvm::Sp(solid.TgSp()*solid.rCp(), Ts)
    );

    TgEqn.relax();
    TgEqn.solve();



}
