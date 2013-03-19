
    Info<< "Solving UEqn" << endl;

    surfaceScalarField muEff = ib.interpolate(turbulence->muEff())
                               *ib.alphaf();
    surfaceVectorField gradUeff = fvc::interpolate(dev2(Foam::T(fvc::grad(U))))
                                & mesh.Sf();

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho*ib.alphaCorr(), U)
      + fvm::div(phi, U)
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*gradUeff)
      - solid.mU()
      - rho*ib.alpha()*g
      + fvm::Sp(solid.USp(), U)
      - solid.USu()
      + fvm::Sp(turbulence->muEff()*solid.Stress(), U)
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
