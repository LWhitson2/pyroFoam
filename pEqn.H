
Info<< "Solving pEqn" << endl;

rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
U = rAU*UEqn.H();

Info<< "Min,max Uest = " << Foam::min(Foam::mag(U)).value() 
    << ", " << Foam::max(Foam::mag(U)).value() << endl;

/*if (pimple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rAU, rho, U, phi)
        )
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rho*rAU, p)
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi == pEqn.flux();
        }
    }
}
else*/
{
    surfaceScalarField& phiEst = phi; //for clarity, let's rename phi here

    phiEst =
        fvc::interpolate(rho)*solid.alphaf() //force phiEst to zero inside the solid
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rAU, rho, U, phi)
        );

    // treat the solid as "slightly" compressible 
    //  (rho changes 0.001 kg/m3 per bar)
    dimensionedScalar psiS("psiS",dimDensity/dimPressure,1e-8); 

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psiS+solid.alpha()*(psi-psiS), p)
          + fvc::div(phiEst)
          - fvm::laplacian(rho*rAU*solid.alpha(), p)
         == 
            solid.m_pyro()
          - fvm::Sp(psi*fvc::ddt(solid.alpha()), p)
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi += pEqn.flux(); // * solid.alphaf(); // ?
        }
    }
}

Info<< "Min,max p = " << Foam::min(p).value() 
    << ", " << Foam::max(p).value() << endl;

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U -= rAU*fvc::grad(p)*solid.alpha();
U.correctBoundaryConditions();

Info<< "Min,max U = " << Foam::min(Foam::mag(U)).value() 
    << ", " << Foam::max(Foam::mag(U)).value() << endl;

K = 0.5*magSqr(U);

dpdt = fvc::ddt(p);
