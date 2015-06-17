
Info<< "Solving pEqn" << endl;

rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf(rAU.name() + 'f', ib.interpolate(rho*rAU, "gas"));
U = rAU*UEqn.H();

UestOut = U;

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
        ib.interpolate(rho, "gas")*
        (
            (ib.interpolate(U, "gas") & mesh.Sf())*ib.alphafCorr()
            // + fvc::ddtPhiCorr(rAU, rho, U, phi)
            // + ib.ddtPhiCorr(rAU, rho, U, phi)
        )
    );

    while (pimple.correctNonOrthogonal())
    {
        volScalarField psiAlphaCorr = psi*ib.alphaCorr();
        psiAlphaCorr.oldTime() = psi.oldTime() * ib.alphaCorr().oldTime();
        fvScalarMatrix pEqn
        (
            fvm::ddt(psiAlphaCorr, p)
          + fvc::div(phiEst)
          - fvm::laplacian(rhorAUf*ib.alphafCorr()*ib.lapCorr("gas"), p)
         ==
            solid.pSu()
          - fvm::Sp(solid.pSp(), p)
          + solid.mgen()
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

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

Info<< "Min,max U = " << Foam::min(Foam::mag(U)).value()
    << ", " << Foam::max(Foam::mag(U)).value() << endl;

// Info << "Balance: " << (fvc::ddt(rhoAlpha, U) + fvc::div(phi, U))()[4] << endl;

K = 0.5*magSqr(U);

dpdt = fvc::ddt(ib.alphaCorr(), p);
