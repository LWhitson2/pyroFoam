Info<< "Creating combustion model\n" << endl;

autoPtr<combustionModels::psiChemistryCombustionModel> combustion
(
    combustionModels::psiChemistryCombustionModel::New
    (
        mesh
    )
);

psiChemistryModel& chemistry = combustion->pChemistry();

hsCombustionThermo& thermo = chemistry.thermo();

basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

//word inertSpecie(thermo.lookup("inertSpecie"));

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField refinementField
(
    IOobject
    (
        "refinementField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("refinementField",dimless,0)
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();
volScalarField& hs = thermo.hs();
const volScalarField& T = thermo.T();

#include "compressibleCreatePhi.H"

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

// Set the turbulence into the combustion model
combustion->setTurbulence(turbulence());


Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt("dpdt", fvc::ddt(p));

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(hs);

volScalarField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
);




// Create an immersed boundary
immersedBoundary ib(mesh);

// Create the burning solid object defined by the immersed boundary
burningSolid solid(mesh, ib, thermo);

// Create Solid Fields
volScalarField& Ts = solid.Ts();
solid.solidThermo().T() = Ts;
solid.solidThermo().correct();


volScalarField Tall
(
    IOobject
    (
        "Tall",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Ts*ib.alphas() + T*ib.alpha()
);


const volScalarField& rhos = solid.rhos();
const tmp<volScalarField> Cps = solid.Cps();
const volScalarField& Ks = solid.Ks();

OFstream log(args.path()/"TsAvg.out");

if(Pstream::master())
{
    log<< "Time" << token::TAB << "Avg. Solid Temp. [K]" << token::TAB
        << endl;
}
