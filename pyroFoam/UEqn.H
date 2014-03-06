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
      - fvm::laplacian(muEff*ib.lapCorr("gas"), U)
      - fvc::div(muEff*gradUeff)
      - rho*ib.alphaCorr()*g
     ==
        solid.USu()
      - fvm::Sp(solid.USp(), U)
//       + solid.Ugen()
      + solid.mU()

      //+ fvm::Sp(turbulence->muEff()*solid.Stress(), U)
    );

    UEqn.relax();

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
