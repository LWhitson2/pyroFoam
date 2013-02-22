{
    // Calculate corrected thermal diffusivity
    surfaceScalarField alphaEff = ib.interpolate(turbulence->alphaEff())
                               *ib.alphaf();

    fvScalarMatrix hsEqn
    (
        fvm::ddt(ib.alphaCorr()*rho, hs)
      + mvConvection->fvmDiv(phi, hs)
      - fvm::laplacian(alphaEff, hs)
     ==
        ib.alphaCorr()*dpdt
      - (fvc::ddt(ib.alphaCorr()*rho, K) + fvc::div(phi, K))
      //+ ib.alphaCorr()*combustion->Sh()
      //+ solid.Sh()   //<--
      + solid.hsSu()
      - fvm::Sp(solid.hsSp(), hs)
      + solid.Qtg()
    );

    hsEqn.relax();
    hsEqn.solve();

    thermo.correct(); //we may need to create a 'specie' for the solid so that this works

    Info<< "T gas min/max   = " <<
        min(ib.solidCells()*dimensionedScalar("tmpT", dimTemperature, 1000.0)
        + ib.gasCells()*T).value() << ", " << max(ib.gasCells()*T).value()
        << endl;
}
