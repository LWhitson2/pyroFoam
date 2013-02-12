{
    // Calculate corrected thermal conductivity
    surfaceScalarField KEff = ib.interpolate(Ks)
                               *(1. - ib.alphaf());

    fvScalarMatrix TsEqn
    (
        fvm::ddt((1. - ib.alpha())*rhos*Cps, Ts)
      - fvm::laplacian(KEff, Ts)
     ==
        solid.TsSu()
      - fvm::Sp(solid.TsSp(), Ts)
    );

    TsEqn.relax();
    TsEqn.solve();

    thermo.correct(); //we may need to create a 'specie' for the solid so that this works

    Info<< "T solid min/max   = " << min(Ts).value() << ", "
        << max(Ts).value() << endl;
}
