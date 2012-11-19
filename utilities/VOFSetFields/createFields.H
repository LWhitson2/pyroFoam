// * * * * * ALPHA FIELDS * * * * * * * * * * * * * * * * * * * * * * * * //
    IOdictionary VOFDict
    (
        IOobject
        (
            "VOFSetFieldsDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    
    
    word alphaName(VOFDict.lookup("alphaName"));



    volScalarField alpha
    (
        IOobject
        (
            alphaName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );



// * * * * * VELOCITY FIELDS * * * * * * * * * * * * * * * * * * * * * * * * //
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


