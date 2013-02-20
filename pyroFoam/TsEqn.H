{
    // Calculate corrected thermal conductivity
    tmp<volScalarField> Cps = solid.Cps();
    surfaceScalarField KEff = ib.interpolate(Ks)
                               *(1. - ib.alphaf());

    fvScalarMatrix TsEqn
    (
        fvm::ddt(ib.alphasCorr()*rhos*Cps, Ts)
      - fvm::laplacian(KEff, Ts)
     ==
        solid.TsSu()
      - fvm::Sp(solid.TsSp(), Ts)
    );

    TsEqn.relax();
    TsEqn.solve();

    Info<< "T solid min/max   = " <<
        min(ib.gasCells()*dimensionedScalar("tmpT", dimTemperature, 1000.0)
        + ib.solidCells()*Ts).value() << ", " << max(ib.solidCells()*Ts).value()
        << endl;

    // Update solidThermo
    solid.solidThermo().T() = Ts;
    solid.solidThermo().correct();
}
