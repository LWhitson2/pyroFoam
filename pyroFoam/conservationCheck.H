// Find outlet patch at top
word patchName = "top";
label patchID = mesh.boundaryMesh().findPatchID(patchName);

// Solid Fields
tmp<volScalarField> Cps = solid.Cps();
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
    scalar massFlow = Foam::sum(phi.boundaryField()[patchID]);
    scalar massConserve = solidMassChange + gasMassChange + massFlow;
    coMass = coMass + massConserve*runTime.deltaTValue();

    // Momentum conservation
    vector gasMomChange = Foam::fvc::domainIntegrate(
        Foam::fvc::ddt(rhoAlpha, U)).value();
    vector momFlow = Foam::sum(phi.boundaryField()[patchID]
                             * U.boundaryField()[patchID]);
    vector momConserve = gasMomChange + momFlow;
    coMom = coMom + momConserve*runTime.deltaTValue();

    // Energy conservation
    scalar solidCombustion = Foam::sum(ib.area()*solid.qflux()).value();
    scalar gasCombustion = Foam::fvc::domainIntegrate(ib.alphaCorr()*reaction->Sh()).value();

    scalar solidEnergyChange = Foam::fvc::domainIntegrate(
        Cps*Foam::fvc::ddt(rhosAlphas, Ts)).value();
    scalar gasEnthalpyChange = Foam::fvc::domainIntegrate(
        Foam::fvc::ddt(rhoAlpha, he)).value();
    scalar gasKineticChange = Foam::fvc::domainIntegrate(
        Foam::fvc::ddt(rhoAlpha, K)).value();
    scalar gasPressureChange = Foam::fvc::domainIntegrate(
        Foam::fvc::ddt(ib.alphaCorr(), p)).value();
    scalar energyFlow = Foam::sum(phi.boundaryField()[patchID]
                                * he.boundaryField()[patchID]);
    scalar energyConserve = solidEnergyChange + gasEnthalpyChange
                          // + gasKineticChange + gasPressureChange
                          + energyFlow - solidCombustion - gasCombustion;
    coEnergy = coEnergy + energyConserve*runTime.deltaTValue();

    // Debug outputs
    // tmp<volScalarField> Is = solid.Is();
    // Info << endl << endl << max(Is()).value() << ", " << max(Is().oldTime()).value() << endl << endl; 

    // Output values to logs
    logConserve << runTime.value() << token::TAB 
                << massConserve << token::TAB
                << momConserve << token::TAB
                << energyConserve << endl;

    logConserve2 << runTime.value() << token::TAB
                 << coMass << token::TAB
                 << coMom << token::TAB
                 << coEnergy << endl;
}
