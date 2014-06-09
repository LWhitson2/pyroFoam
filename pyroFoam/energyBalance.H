// Find outlet patch at top
word patchName = "top";
label patchID = mesh.boundaryMesh().findPatchID(patchName);
tmp<volScalarField> Cps = solid.Cps();
tmp<volScalarField> alphas = ib.alphasCorr();
volScalarField rhoAlpha = rho*ib.alphaCorr();
rhoAlpha.oldTime() = rho.oldTime() * ib.alphaCorr().oldTime();

if( patchID >= 0 )
{
    scalar energyFlow = Foam::sum(phi.boundaryField()[patchID]
                      *heOut.boundaryField()[patchID])
                      +Foam::sum(solid.qflux().internalField()*ib.area().internalField());

    scalar energyChange = Foam::fvc::domainIntegrate(
        Foam::fvc::ddt(rhoAlpha, he) + alphas*rhos*Cps*fvc::ddt(Ts)).value();

    scalar totalEnergy = Foam::fvc::domainIntegrate(
        rhoAlpha*he
        + ib.alphasCorr()*solid.rhos()*solid.Cps()*(Ts - dimensionedScalar("T0", dimTemperature, 273.15)
        )).value();

    Foam::reduce(energyFlow, sumOp<scalar>());

    Foam::reduce(energyChange, sumOp<scalar>());

    coe = coe + (energyChange + energyFlow)*runTime.deltaTValue();

    Info<< runTime.value() << token::TAB 
          << energyFlow << token::TAB
          << energyChange << token::TAB
          << totalEnergy << token::TAB
          << coe/totalEnergy << endl;
}
