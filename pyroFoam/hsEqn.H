{
    fvScalarMatrix hsEqn
    (
        fvm::ddt(rho, hs)
      + mvConvection->fvmDiv(phi, hs)
      - fvm::laplacian(turbulence->alphaEff(), hs) //needs to be adjusted
     ==
        dpdt
      - (fvc::ddt(rho, K) + fvc::div(phi, K))
      + combustion->Sh()
      + solid.Sh()   //<--
    );

    hsEqn.relax();
    hsEqn.solve();

    thermo.correct(); //we may need to create a 'specie' for the solid so that this works

    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}
