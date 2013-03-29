{
    //solid.calcHeatTransfer();
    // Calculate corrected thermal conductivity
    tmp<volScalarField> Cps = solid.Cps();
    surfaceScalarField KEff = ib.interpolate(Ks)
                               *ib.alphafs();

    fvScalarMatrix TsEqn
    (
        fvm::ddt(ib.alphasCorr()*rhos*Cps, Ts)
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

    Info<< "T solid min/max   = " <<
        min(ib.fullGasCells()*dimensionedScalar("tmpT", dimTemperature, 100000.0)
        + Ts).value() << ", " << max(ib.solidCells()*Ts).value()
        << endl;

    // Update solidThermo
    solid.solidThermo().T() = Ts;
    solid.solidThermo().correct();
}
