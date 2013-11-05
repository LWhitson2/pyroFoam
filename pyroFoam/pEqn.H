
Info<< "Solving pEqn" << endl;

rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf(rAU.name() + 'f', ib.interpolate(rho*rAU));
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
            (fvc::interpolate(U) & mesh.Sf())*solid.alphaf()
         // + fvc::ddtPhiCorr(rAU, rho, U, phi)
        )
    );

    dimensionedScalar psiS("psiS",dimDensity/dimPressure,SMALL);
    volScalarField psiEff = psi*solid.alpha() + psiS*(1-solid.alpha());

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psiEff, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rhorAUf*solid.alphaf(), p)
         ==
            solid.m_pyro()
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi == pEqn.flux();
            U = rAU*fvc::reconstruct(pEqn.flux()/rhorAUf);
            U.correctBoundaryConditions();
        }
    }
}
else*/
{
    surfaceScalarField phiEst
    (
        ib.interpolate(rho)*
        (
            (ib.interpolate(U) & mesh.Sf())*ib.alphaf()
          //+ ib.ddtPhiCorr(rAU, rho, U, phi)
        )
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi*ib.alphaCorr(), p)
          + fvc::div(phiEst)
          - fvm::laplacian(rhorAUf*ib.alphaf()*ib.lapCorr("gas"), p)
          + fvm::Sp(solid.pSp(), p)
         ==
            solid.m_pyro()
          + solid.pSu()
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiEst + pEqn.flux();
            U += rAU*ib.reconstruct(pEqn.flux()/rhorAUf);
            U.correctBoundaryConditions();
        }
    }
}

Info<< "Min,max p = " << Foam::min(p).value()
    << ", " << Foam::max(p).value() << endl;

// #include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

Info<< "Min,max U = " << Foam::min(Foam::mag(U)).value()
    << ", " << Foam::max(Foam::mag(U)).value() << endl;

K = 0.5*magSqr(U);

dpdt = fvc::ddt(p);
