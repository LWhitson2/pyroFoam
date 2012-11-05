
    Info<< "Solving UEqn" << endl;
    
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho*solid.alpha(), U)
      + fvm::div(phi, U)
      //+ turbulence->divDevRhoReff(U)
      - solid.m_pyro()*solid.burnU()
      - rho*solid.alpha()*g
      + fvm::Sp(solid.setDiag(), U) //switch to set U = 0 in the solid
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
