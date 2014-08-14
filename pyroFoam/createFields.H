Info<< "Creating reaction model\n" << endl;

autoPtr<combustionModels::psiCombustionModel> reaction
(
    combustionModels::psiCombustionModel::New(mesh)
);


psiReactionThermo& thermo = reaction->thermo();
thermo.validate(args.executable(), "h", "e");


basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();


word inertSpecie(thermo.lookup("inertSpecie"));

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
volScalarField& he = thermo.he();
const volScalarField& T = thermo.T();

volScalarField heOut
(
    IOobject
    (
        "heOut",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.he()
);

volVectorField UestOut
(
    IOobject
    (
        "UestOut",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U
);

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
reaction->setTurbulence(turbulence());


Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt("dpdt", fvc::ddt(p));

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(he);

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
//  For now, just set thermo type here to match what is in dictionary
burningSolid<constGasHThermoPhysics,psiReactionThermo> solid(mesh, ib, thermo);

// Create Solid Fields
volScalarField& Ts = solid.Ts();
solid.thermo().T() = Ts;
solid.thermo().correct();

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
//const tmp<volScalarField> Cps = solid.Cps();
//const tmp<volScalarField> Ks = solid.Ks();

OFstream log(args.path()/"TsAvg.out");
OFstream logTi(args.path()/"TiAvg.out");

log << "Time" << token::TAB << "Avg. Solid Temp. [K]" << token::TAB
    << endl;

logTi << "Time" << token::TAB
      << "alpha" << token::TAB
      << "Ts" << token::TAB
      << "Ti" << token::TAB
      << "Tg" << token::TAB
      << "p" << token::TAB
      << "EMg" << token::TAB
      << "Lg" << token::TAB
      << "Ls" << endl;

// Conservation Check log
OFstream logConserve(args.path()/"TimeStepConserve.out");
logConserve << "Time" << token::TAB
            << "Mass" << token::TAB
            << "Momentum" << token::TAB
            << "Energy" << endl;

// Overall Conservation log
OFstream logConserve2(args.path()/"OverallConserve.out");
logConserve << "Time" << token::TAB
            << "Mass" << token::TAB
            << "Momentum" << token::TAB
            << "Energy" << endl;