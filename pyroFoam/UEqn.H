
    Info<< "Solving UEqn" << endl;
    
    // The continuous forcing approach is used for the immersed boundary
    // through the solid permeability model.
    // This method becomes instable for high Re, so a discrete forcing
    // approach should be considered for future development.
       
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho*solid.alpha(), U)                   //<--
      + fvm::div(phi, U)
      //+ turbulence->divDevRhoReff(U)                   //<--
     ==
        rho*solid.alpha()*g                              //<--
      - fvm::Sp(turbulence->muEff()/solid.kappa(), U)
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));
        K = 0.5*magSqr(U); // scale this by solid.alpha() ?
    }

    Info<< "Min,max U = " << Foam::min(Foam::mag(U)).value() 
        << ", " << Foam::max(Foam::mag(U)).value() << endl;
