{
//     solid.calcInterfaceTransfer();

    // Calculate corrected thermal diffusivity
    surfaceScalarField alphaEff = ib.interpolate(turbulence->alphaEff())
                               *ib.alphafU()*ib.lapCorr("gas");

//     forAll(hs, cellI)
//     {
//         if (ib.smallCells()()[cellI]) Info << solid.QgSu()[cellI] << ", " << solid.QgSp()[cellI] << endl;
//     }

    fvScalarMatrix heEqn
    (
        fvm::ddt(ib.alphaCorr()*rho, he)
      + mvConvection->fvmDiv(phi, he)
      - fvm::laplacian(alphaEff, he)
     ==
        ib.alphaCorr()*dpdt
      - (fvc::ddt(ib.alphaCorr()*rho, K) + fvc::div(phi, K))
      + ib.alphaCorr()*reaction->Sh()
//       + combustion->Sh()
      + solid.hsSu()
      - fvm::Sp(solid.hsSp(), he)
      + solid.QgSu()
      - fvm::Sp(solid.QgSp(), he)
      + solid.Qgeng()
    );

    heEqn.relax();
    heEqn.solve();

    thermo.correct(); //we may need to create a 'specie' for the solid so that this works

    Info<< "T gas min/max   = " <<
        min(ib.solidCells()*dimensionedScalar("tmpT", dimTemperature, 100000.0)
        + T).value() << ", " << max(T).value()
        << endl;
}
