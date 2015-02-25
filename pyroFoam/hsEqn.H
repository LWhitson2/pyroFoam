{
    // solid.calcInterfaceTransfer();

    // Calculate corrected thermal diffusivity
    surfaceScalarField alphaEff = ib.interpolate(turbulence->alphaEff(),"gas")
                               *ib.alphaf()*ib.lapCorr("gas");

//     forAll(hs, cellI)
//     {
//         if (ib.smallCells()()[cellI]) Info << solid.QgSu()[cellI] << ", " << solid.QgSp()[cellI] << endl;
//     }

    // Info << "Setting up gas enthalpy equation" << endl;
    fvScalarMatrix heEqn
    (
        fvm::ddt(rhoAlphaCorr, he)
      + mvConvection->fvmDiv(phi, he)
      - fvm::laplacian(alphaEff, he)
     ==
      //   ib.alphaCorr()*dpdt
      // - (fvc::ddt(rhoAlphaCorr, K) + fvc::div(phi, K))
       ib.alphaCorr()*reaction->Sh()
      + solid.hsSu()
      - fvm::Sp(solid.hsSp(), he)
      + solid.QgSu()
      - fvm::Sp(solid.QgSp(), he)
      + solid.Qgeng()
    );

    // Info << "Relaxing gas enthalpy equation" << endl;
    heEqn.relax();
    // Info << "Solving gas enthalpy equation" << endl;
    heEqn.solve();

    // Info<< "hs Curr: " << he[68] << endl;
    // Info << "hs Old: " << he.oldTime()[68] << endl;
    // Info<< "rhoAlph Curr: " << rhoAlpha.oldTime()[68] << endl;
    // Info<< "rhoAlph Old: " << rhoAlpha[68] << endl;
    // Info<< "rhoAlphCorr Curr: " << rhoAlphaCorr.oldTime()[68] << endl;
    // Info<< "rhoAlphCorr Old: " << rhoAlphaCorr[68] << endl;
    // Info << "Sum phi hs: " << mvConvection->fvcDiv(phi,he)()[68] << endl;
    // Info << "Laplacian: " << -fvc::laplacian(alphaEff, he)()[68] << endl;
    // Info << "dpdt: " << dpdt[68] << endl;
    // Info << "Kinetic Energy: " << -(fvc::ddt(rhoAlphaCorr, K) + fvc::div(phi, K))()[68] << endl;
    // Info << "Combustion: " << ib.alphaCorr()[68]*reaction->Sh()()[68] << endl;
    // Info<< "hs Gen: " << solid.Qgeng()[68] << endl;

    Info << "Correcting gas thermo" << endl;
    thermo.correct();

    Info<< "T gas min/max   = " <<
        min(ib.solidCells()*dimensionedScalar("tmpT", dimTemperature, 100000.0)
        + T).value() << ", " << max(T).value()
        << endl;
}
