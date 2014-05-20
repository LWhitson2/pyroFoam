// Find outlet patch at top
word patchName = "top";
label patchID = mesh.boundaryMesh().findPatchID(patchName);
tmp<volScalarField> Cps = solid.Cps();
tmp<volScalarField> alphas = ib.alphasCorr();
volScalarField rhoAlpha = rho*ib.alphaCorr();
rhoAlpha.oldTime() = rho.oldTime() * ib.alphaCorr().oldTime();

if( patchID >= 0 )
{
    Info << "Massflow" << endl;
    scalar massflow = Foam::sum(rho.boundaryField()[patchID]*phi.boundaryField()[patchID]);

    Info << "energyflow" << endl;
    scalar energyFlow = Foam::sum(rho.boundaryField()[patchID]
                      *phi.boundaryField()[patchID]
                      *hsOut.boundaryField()[patchID]);

    Info << "energychange" << endl;
    scalar energyChange = Foam::fvc::domainIntegrate(
        Foam::fvc::ddt(rhoAlpha, hs) + alphas*rhos*Cps*fvc::ddt(Ts)
        ).value();

    Info << "totalenergy" << endl;
    scalar totalEnergy = Foam::fvc::domainIntegrate(
        rhoAlpha*hs
        + ib.alphasCorr()*solid.rhos()*solid.Cps()*(Ts - dimensionedScalar("T0", dimTemperature, 273.15)
        )).value();

    Info << "reduce" << endl;
    Foam::reduce(massflow, sumOp<scalar>());
    Info << "reduce" << endl;
    Foam::reduce(energyFlow, sumOp<scalar>());
    Info << "reduce" << endl;
    Foam::reduce(energyChange, sumOp<scalar>());
    Info << "coe" << endl;
    coe = coe + (energyChange + energyFlow)*runTime.deltaTValue();
    Info << "print" << endl;
    Info<< runTime.value() << token::TAB 
          << energyFlow << token::TAB
          << energyChange << token::TAB
          << totalEnergy << token::TAB
          << coe/totalEnergy << endl;

}