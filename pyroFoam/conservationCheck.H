// Find outlet patch at top
word patchName = "top";
label patchID = mesh.boundaryMesh().findPatchID(patchName);
tmp<volScalarField> Cps = solid.Cps();

// Solid Fields
volScalarField rhosAlphas = rhos*(1. - ib.alpha());
rhosAlphas.oldTime() = rhos*(1. - ib.alpha().oldTime());

// Gas Fields
volScalarField rhoAlpha = rho*ib.alpha();
rhoAlpha.oldTime() = rho.oldTime() * ib.alpha().oldTime();

if( patchID >= 0 )
{
    // Mass conservation
    scalar solidMassChange = Foam::fvc::domainIntegrate(
        Foam::fvc::ddt(rhosAlphas)).value();
    scalar gasMassChange = Foam::fvc::domainIntegrate(
        Foam::fvc::ddt(rhoAlpha)).value();
    scalar massFlow = Foam::sum(phi.boundaryField());
    scalar massConserve = solidMassChange + gasMassChange + massFlow;

    // scalar energyFlow = Foam::sum(phi.boundaryField()[patchID]
    //                   *heOut.boundaryField()[patchID])
    //                   +Foam::sum(solid.qflux().internalField()*ib.area().internalField());

    // scalar energyChange = Foam::fvc::domainIntegrate(
    //     Foam::fvc::ddt(rhoAlpha, he) + alphas*rhos*Cps*fvc::ddt(Ts)).value();

    // scalar totalEnergy = Foam::fvc::domainIntegrate(
    //     rhoAlpha*he
    //     + ib.alphasCorr()*solid.rhos()*solid.Cps()*(Ts - dimensionedScalar("T0", dimTemperature, 273.15)
    //     )).value();

    // Foam::reduce(energyFlow, sumOp<scalar>());

    // Foam::reduce(energyChange, sumOp<scalar>());

    // coe = coe + (energyChange + energyFlow)*runTime.deltaTValue();

    logConserve << runTime.value() << token::TAB 
                << massConserve << endl;
}
