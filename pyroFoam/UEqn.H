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

    surfaceScalarField muEff = ib.interpolate(turbulence->muEff(), "gas")
                               *ib.alphaf();
    surfaceVectorField gradUeff = fvc::interpolate(dev2(Foam::T(fvc::grad(U))))
                                & mesh.Sf();
    volScalarField rhoAlpha = rho*ib.alphaCorr();
    rhoAlpha.oldTime() = rho.oldTime() * ib.alphaCorr().oldTime();

    fvVectorMatrix UEqn
    (
        fvm::ddt(rhoAlpha, U)
      + fvm::div(phi, U)
      // - fvm::laplacian(muEff*ib.lapCorr("gas"), U)
      // - fvc::div(muEff*gradUeff)
      // - rho*ib.alphaCorr()*g
     ==
        solid.USu()
      - fvm::Sp(solid.USp(), U)
      // + solid.mU()
      // - rho.oldTime()*U.oldTime()*fvc::ddt(ib.alphaCorr())
      // + fvm::Sp(solid.mgen(), U)

      //+ fvm::Sp(turbulence->muEff()*solid.Stress(), U)
    );

    UEqn.relax();

    // Info << "Ueqn Coefficients: rho*U*ddt(alpha), div(phi,U)" << endl;
    // Info << "Balance: " << (fvc::ddt(rhoAlpha, U) + fvc::div(phi, U))()[4] << endl;
    // Info << "Cell 4: " << (ib.alphaCorr().oldTime()*rho*fvc::ddt(U))()[4] << ", " << (rho.oldTime()*U.oldTime()*fvc::ddt(ib.alphaCorr()))()[4] << ", " << fvc::div(phi, U)()[4] << endl;
    // Info << "Cell 5: " << (ib.alphaCorr().oldTime()*rho*fvc::ddt(U))()[5] << ", "  << (rho.oldTime()*U.oldTime()*fvc::ddt(ib.alphaCorr()))()[5] << ", " << fvc::div(phi, U)()[5] << endl;
    // Info << "Cell 6: " << (ib.alphaCorr().oldTime()*rho*fvc::ddt(U))()[6] << ", "  << (rho.oldTime()*U.oldTime()*fvc::ddt(ib.alphaCorr()))()[6] << ", " << fvc::div(phi, U)()[6] << endl;
    // Info << "Cell 20: " << fvc::div(phi, U)()[20] << endl;

    if (pimple.momentumPredictor())
    {
        solve
        (
           UEqn == ib.reconstruct(-fvc::snGrad(p)*mesh.magSf()*ib.alphaf())
        );
        K = 0.5*magSqr(U); // scale this by ib.alpha() ?
    }

    Info<< "Min,max U = " << Foam::min(Foam::mag(U)).value()
        << ", " << Foam::max(Foam::mag(U)).value() << endl;
