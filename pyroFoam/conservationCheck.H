
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
    // Info << "Mass Conservation Check" << endl;
    scalar totalMass = Foam::fvc::domainIntegrate(rhosAlphas).value()
                     + Foam::fvc::domainIntegrate(rhoAlpha).value();
    scalar solidMassChange = Foam::fvc::domainIntegrate(
        Foam::fvc::ddt(rhosAlphas)).value();
    scalar gasMassChange = Foam::fvc::domainIntegrate(
        Foam::fvc::ddt(rhoAlpha)).value();
    scalar massFlow = Foam::sum(phi.boundaryField()[patchID]);
    scalar massConserve = solidMassChange + gasMassChange + massFlow;
    coMass = coMass + massConserve*runTime.deltaTValue();

    // Momentum conservation
    // Info << "Momentum Conservation Check" << endl;
    vector totalMomentum = Foam::fvc::domainIntegrate(rhoAlpha*U).value();
    vector gasMomChange = Foam::fvc::domainIntegrate(
        Foam::fvc::ddt(rhoAlpha, U)).value();
    vector momFlow = Foam::sum(phi.boundaryField()[patchID]
                             * U.boundaryField()[patchID]);
    vector momConserve = gasMomChange + momFlow;
    coMom = coMom + momConserve*runTime.deltaTValue();

    // Energy conservation
    // Info << "Energy Conservation Check" << endl;
    scalar totalEnergy = Foam::fvc::domainIntegrate(Cps()*rhosAlphas*Ts).value();
    totalEnergy += Foam::fvc::domainIntegrate(rhosAlphas*(400.e3+3018.e3)).value();
    totalEnergy += Foam::fvc::domainIntegrate(rhoAlpha*he).value(); 
    forAll(Y, i)
    {
        if (Y[i].name() == "EMg")
        {
            volScalarField& Yi = Y[i];
            totalEnergy += Foam::fvc::domainIntegrate(rhoAlpha*3018.e3*Yi).value();
        }
    }

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
    // Info << "Output Conservation" << endl;
    logConserve << runTime.value() << token::TAB 
                << massConserve << token::TAB
                << momConserve << token::TAB
                << energyConserve << endl;

    logConserve2 << runTime.value() << token::TAB
                 << coMass/totalMass << token::TAB
                 << coMom[0]/totalMomentum[0] << token::TAB
                 << coMom[1]/totalMomentum[1] << token::TAB
                 << coMom[2]/totalMomentum[2] << token::TAB
                 << coEnergy/totalEnergy << endl;
}
