tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

    Info<< "Solving UEqn" << endl;

    // surfaceScalarField muEff = ib.interpolate(turbulence->muEff(), "gas")
    //                            *ib.alphaf();
    // surfaceVectorField gradUeff = fvc::interpolate(dev2(Foam::T(fvc::grad(U))))
    //                             & mesh.Sf();
    volScalarField rhoAlphaCorr = rho*ib.alphaCorr();
    rhoAlphaCorr.oldTime() = rho.oldTime() * ib.alphaCorr().oldTime();
    
    volScalarField rhoAlpha = rho*ib.alpha();
    rhoAlpha.oldTime() = rho.oldTime() * ib.alpha().oldTime();

    fvVectorMatrix UEqn
    (
        fvm::ddt(rhoAlphaCorr, U)
      + fvm::div(phi, U)
      // - fvm::laplacian(muEff*ib.lapCorr("gas"), U)
      // - fvc::div(muEff*gradUeff)
      // - rho*ib.alphaCorr()*g
     ==
        solid.USu()
      - fvm::Sp(solid.USp(), U)
      + solid.mU()
      // - rho.oldTime()*U.oldTime()*fvc::ddt(ib.alphaCorr())
      // + fvm::Sp(solid.mgen(), U)

      //+ fvm::Sp(turbulence->muEff()*solid.Stress(), U)
    );

    UEqn.relax();

    // Info << "Ueqn Coefficients: rho*U*ddt(alpha), div(phi,U)" << endl;
    // Info << "Balance: " << (fvc::ddt(rhoAlpha, U) + fvc::div(phi, U))()[99] << endl;
    Info << "Cell 99: " << (rhoAlphaCorr.oldTime()*fvc::ddt(U))()[99] << ", " << (U.oldTime()*fvc::ddt(rhoAlphaCorr))()[99] << ", " << fvc::div(phi, U)()[99] << endl;
    Info << "Cell 98: " << (rhoAlphaCorr.oldTime()*fvc::ddt(U))()[98] << ", " << (U.oldTime()*fvc::ddt(rhoAlphaCorr))()[98] << ", " << fvc::div(phi, U)()[98] << endl;
    // Info << "Cell 6: " << (ib.alphaCorr().oldTime()*rho*fvc::ddt(U))()[6] << ", "  << (rho.oldTime()*U.oldTime()*fvc::ddt(ib.alphaCorr()))()[6] << ", " << fvc::div(phi, U)()[6] << endl;
    // Info << "Cell 20: " << fvc::div(phi, U)()[20] << endl;

    if (pimple.momentumPredictor())
    {
        solve
        (
           UEqn == ib.reconstruct(-fvc::snGrad(p)*mesh.magSf()*ib.alphafCorr())
        );
        K = 0.5*magSqr(U); // scale this by ib.alpha() ?
    }

    Info<< "Min,max U = " << Foam::min(Foam::mag(U)).value()
        << ", " << Foam::max(Foam::mag(U)).value() << endl;
