{
//     solid.calcInterfaceTransfer();

    // Calculate corrected thermal diffusivity
    surfaceScalarField alphaEff = ib.interpolate(turbulence->alphaEff(),"solid")
                               *ib.alphaf()*ib.lapCorr("gas");

//     forAll(hs, cellI)
//     {
//         if (ib.smallCells()()[cellI]) Info << solid.QgSu()[cellI] << ", " << solid.QgSp()[cellI] << endl;
//     }

    fvScalarMatrix heEqn
    (
        fvm::ddt(rhoAlpha, hs)
      + mvConvection->fvmDiv(phi, he)
      - fvm::laplacian(alphaEff, he)
     ==
        ib.alphaCorr()*dpdt
      - (fvc::ddt(rhoAlpha, K) + fvc::div(phi, K))
      + ib.alphaCorr()*reaction->Sh()
      + solid.hsSu()
      - fvm::Sp(solid.hsSp(), he)
      + solid.QgSu()
      - fvm::Sp(solid.QgSp(), he)
      + solid.Qgeng()
    );

    heEqn.relax();
    heEqn.solve();

    Info<< "hs Curr: " << hs[99] << endl;
    Info << "hs Old: " << hs.oldTime()[99] << endl;
    Info<< "rhoAlph Curr: " << rhoAlpha.oldTime()[99] << endl;
    Info<< "rhoAlph Old: " << rhoAlpha[99] << endl;
    Info << "Sum phi hs: " << mvConvection->fvcDiv(phi,hs)()[99] << endl;
    Info << "Laplacian: " << -fvc::laplacian(alphaEff, hs)()[99] << endl;
    Info << "dpdt: " << dpdt[99] << endl;
    Info << "Kinetic Energy: " << -(fvc::ddt(rhoAlpha, K) + fvc::div(phi, K))()[99] << endl;
    Info << "Combustion: " << ib.alphaCorr()[99]*combustion->Sh()()[99] << endl;
    Info<< "hs Gen: " << solid.Qgeng()[99] << endl;

    thermo.correct(); //we may need to create a 'specie' for the solid so that this works

    Info<< "T gas min/max   = " <<
        min(ib.solidCells()*dimensionedScalar("tmpT", dimTemperature, 100000.0)
        + T).value() << ", " << max(T).value()
        << endl;
}
