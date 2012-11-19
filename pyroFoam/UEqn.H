
    Info<< "Solving UEqn" << endl;
    
    surfaceScalarField muEff = ib.interpolate(turbulence->muEff())
                               *ib.alphaf();

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho*ib.alphaCorr(), U)
      + fvm::div(phi, U)
      - fvm::laplacian(muEff, U)
      //- fvc::div(muEff*dev2(T(fvc::grad(U))))
      - solid.mU()
      - rho*ib.alpha()*g
      + fvm::Sp(solid.USp(), U)
      - solid.USu()
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
