
    Info<< "Solving UEqn" << endl;

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho*solid.alphaUsed(), U)
      + fvm::div(phi, U)
      //+ turbulence->divDevRhoReff(U)
      - solid.mU()
      - rho*solid.alpha()*g
      + fvm::Sp(solid.USp(), U) //switch to set U = 0 in the solid
      + solid.USu()
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
