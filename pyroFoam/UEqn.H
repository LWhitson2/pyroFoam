
    Info<< "Solving UEqn" << endl;
    
    surfaceScalarField muEff = fvc::interpolate(turbulence->muEff())*solid.alphaf();

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho*solid.alphaUsed(), U)
      + fvm::div(phi, U)
      //+ turbulence->divDevRhoReff(U)
      // - fvm::laplacian(muEff, U) - fvc::div(muEff*dev2(T(fvc::grad(U))))
      - fvm::laplacian(muEff, U)
      - solid.mU()
      - rho*solid.alpha()*g
      + fvm::Sp(solid.USp(), U) //switch to set U = 0 in the solid and Uburn in small cells
      - solid.USu()
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
           UEqn == fvc::reconstruct(-fvc::snGrad(p)*mesh.magSf()*solid.alphaf())
        );
        K = 0.5*magSqr(U); // scale this by solid.alpha() ?
    }

    Info<< "Min,max U = " << Foam::min(Foam::mag(U)).value()
        << ", " << Foam::max(Foam::mag(U)).value() << endl;
