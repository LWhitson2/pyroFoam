{
    // Calculate corrected thermal conductivity
    const tmp<volScalarField> Cps = solid.Cps();
    const tmp<volScalarField> alphas = ib.alphasCorr();
    const tmp<volScalarField> Ks = solid.Ks();
    surfaceScalarField KEff = ib.interpolate(Ks(), "solid")
                            * ib.alphafsU()*ib.lapCorr("solid");

    fvScalarMatrix TsEqn
    (
        alphas*rhos*Cps*fvm::ddt(Ts)
      - fvm::laplacian(KEff, Ts)
     ==
        solid.TsSu()
      - fvm::Sp(solid.TsSp(), Ts)
      + solid.QsSu()
      - fvm::Sp(solid.QsSp(), Ts)
      + solid.Qgens()
    );

    TsEqn.relax();
    TsEqn.solve();

    // Info << "ddtTs: " << fvc::ddt(Ts)()[99] << endl;
    // Info << "ddtAlphaTs: " << fvc::ddt(ib.alphaCorr(), Ts)()[99] << endl;
    // Info << "laplacian: " << -fvc::laplacian(KEff, Ts)()[99] << endl;
    // Info << "Qs Terms: " << solid.QsSu()[99] - fvc::Sp(solid.QsSp(), Ts)()[99] << endl;
    // Info << "Ts Terms: " << solid.TsSu()[99] - fvc::Sp(solid.TsSp(), Ts)()[99] << endl;

    Info<< "T solid min/max   = " <<
        min(ib.fullGasCells()*dimensionedScalar("tmpT", dimTemperature, 100000.0)
        + Ts).value() << ", "
        << max((ib.solidCells()+ib.mixedCells())*Ts).value() << endl;

    // Update solidThermo
    solid.thermo().T() = Ts;
    solid.thermo().correct();
}
