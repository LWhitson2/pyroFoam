{
    // Update solid temperature field
    solid.Ts() = Ts;
    tmp<volScalarField> Cps = solid.Cps();

    // Calculate corrected thermal conductivity
    Info << "Calculating KEff" << endl;
    surfaceScalarField KEff = ib.interpolate(Ks)
                               *(1. - ib.alphaf());

    Info << "Setting up TsEqn" << endl;
    fvScalarMatrix TsEqn
    (
        fvm::ddt((1. - ib.alpha())*rhos*Cps, Ts)
      - fvm::laplacian(KEff, Ts)
     ==
        solid.TsSu()
      - fvm::Sp(solid.TsSp(), Ts)
    );

    Info << "Relaxing TsEqn" << endl;
    TsEqn.relax();

    Info << "Solving TsEqn" << endl;
    TsEqn.solve();

    Info<< "T solid min/max   = " << min(Ts).value() << ", "
        << max(Ts).value() << endl;
}
