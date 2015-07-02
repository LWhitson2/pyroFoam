{
    // solid.calcInterfaceTransfer();

    // Calculate corrected thermal diffusivity
    surfaceScalarField alphaEff = ib.interpolate(turbulence->alphaEff(),"gas")
                               *ib.alphafCorr()*ib.lapCorr("gas");

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

    scalar ci = 99;
    Info << "hs Curr: " << he[ci] << endl;
    Info << "hs Old: " << he.oldTime()[ci] << endl;
    Info << "rhoAlph Curr: " << rhoAlpha.oldTime()[ci] << endl;
    Info << "rhoAlph Old: " << rhoAlpha[ci] << endl;
    Info << "rhoAlphCorr Curr: " << rhoAlphaCorr.oldTime()[ci] << endl;
    Info << "rhoAlphCorr Old: " << rhoAlphaCorr[ci] << endl;
    Info << "Sum phi hs: " << mvConvection->fvcDiv(phi,he)()[ci] << endl;
    Info << "Laplacian: " << -fvc::laplacian(alphaEff, he)()[ci] << endl;
    Info << "dpdt: " << dpdt[ci] << endl;
    Info << "Kinetic Energy: " << -(fvc::ddt(rhoAlphaCorr, K) + fvc::div(phi, K))()[ci] << endl;
    Info << "Combustion: " << ib.alphaCorr()[ci]*reaction->Sh()()[ci] << endl;
    Info << "hs Gen: " << solid.Qgeng()[ci] << endl;
    Info << "hsSu: " << solid.hsSu()[ci] << endl;
    Info << "hsSp: " << solid.hsSp()[ci] << endl;

    Info << "Correcting gas thermo" << endl;
    thermo.correct();

    Info<< "T gas min/max   = " <<
        min(ib.solidCells()*dimensionedScalar("tmpT", dimTemperature, 100000.0)
        + T).value() << ", " << max(T).value()
        << endl;
}
