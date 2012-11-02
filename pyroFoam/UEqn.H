
    Info<< "Solving UEqn" << endl;
    
    // Construct a "momentum density field"
    // In cells with more than 1e-5 gas, the equation considers the gas phase
    // only. In cells which are more then 99.99% solid, the solid can act as a
    // momentum sink with its large mass. Momentum flux, however, always uses
    // the gas density.
    volScalarField rhoM = rho*solid.alpha() + (solid.rhoS() - rho)*pos(1e-5 - solid.alpha());

    fvVectorMatrix UEqn
    (
        fvm::ddt(rhoM, U)                   //<--
      + fvm::div(phi, U)
      //+ turbulence->divDevRhoReff(U)                   //<--
     ==
        rho*solid.alpha()*g                              //<--
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));
        K = 0.5*magSqr(U); // scale this by solid.alpha() ?
    }

    Info<< "Min,max U = " << Foam::min(Foam::mag(U)).value() 
        << ", " << Foam::max(Foam::mag(U)).value() << endl;
