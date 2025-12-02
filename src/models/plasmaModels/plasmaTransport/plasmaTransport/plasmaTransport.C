/*---------------------------------------------------------------------------*\
  File: plasmaTransport.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaTransport.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "ScharfetterGummel.H"
#include "plasmaTransport.H"
#include "plasmaTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(plasmaTransport, 0);

// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * *  //

void plasmaTransport::constructModels()
{
    // Loop over species and create a transport model for each one
    for (label i = 0; i < species_.nSpecies(); ++i)
    {
        const word& sName = species_.speciesNames()[i];
        const dictionary& sDict = species_.speciesDict(sName);

        if (!sDict.found("transportModel"))
        {
            FatalIOErrorInFunction(sDict)
                << "Species '" << sName
                << "' is missing required entry 'transportModel' in "
                << species_.dictName() << nl
                << exit(FatalIOError);
        }

        word modelName;
        sDict.lookup("transportModel") >> modelName;

        // Construct the model using the runtime selection system
        transportModels_.set
        (
            i,
            plasmaTransportModel::New(modelName, sDict, mesh_, species_, i, E_)
        );

        Info<< "plasmaTransport: Created transport model '" << modelName
            << "' for species '" << sName << "'\n";
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaTransport::plasmaTransport
(
    plasmaSpecies& species,
    const fvMesh& mesh,
    const volVectorField& E
)
:
    mesh_(mesh),
    species_(species),
    E_(E),
    transportModels_(species.nSpecies()),
    flux_(),
    surfaceFlux_()
{
    constructModels();

    flux_.setSize(species.nSpecies());
    surfaceFlux_.setSize(species.nSpecies());
    for (label i = 0; i < species.nSpecies(); ++i)
    {
        flux_.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "Gamma_" + species.speciesNames()[i], 
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector
                (
                    "zero",
                    dimensionSet(0, -2, -1, 0, 0, 0, 0),
                    vector::zero
                )
            )
        );

        surfaceFlux_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phi_" + species.speciesNames()[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "zero",
                    dimensionSet(0, 3, -1, 0, 0, 0, 0),
                    0.0)
            )
        );
    }
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void plasmaTransport::correct()
{
    const label nSpecies = species_.nSpecies();

    for (label i = 0; i < nSpecies; ++i)
    {
        transportModels_[i].correct();

        const volVectorField& U = transportModels_[i].driftVelocity();

        surfaceFlux_[i] = fvc::flux(U);

        volScalarField& n = species_.numberDensity(i);

        // fvScalarMatrix eqn
        // (
        //     fvm::ddt(n) + fvm::div(surfaceFlux_[i], n)
        // );
const volScalarField* D = transportModels_[i].diffusivity();
fvScalarMatrix eqn
(
    // You must dereference the pointer D to pass the actual field object
    fvm::ddt(n) + fvm::ScharfetterGummel(n, surfaceFlux_[i], *D) // FIX HERE
);

        // if (const volScalarField* D = transportModels_[i].diffusivity())
        // {
        //     eqn -= fvm::laplacian(*D, n);
        // }

        eqn.solve();
    }
}

tmp<volScalarField> plasmaTransport::elecConductionCoeff() const
{
    tmp<volScalarField> elecConductionCoeff =
        tmp<volScalarField>::New
        (
            volScalarField
            (
                IOobject
                (
                    "elecConductionCoeff",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "zeroConduction",
                    dimensionSet(-1,-3,3,0,0,2,0), 
                    0.0
                )
            )
        );

    const label nSpecies = species_.nSpecies();
    for (label i = 0; i < nSpecies; ++i)
    {
        elecConductionCoeff.ref() += transportModels_[i].elecConductionCoeff();
    }

    return elecConductionCoeff;
}

tmp<volScalarField> plasmaTransport::elecDiffusionCharge() const
{
    tmp<volScalarField> elecDiffusionCharge =
        tmp<volScalarField>::New
        (
            volScalarField
            (
                IOobject
                (
                    "elecDiffusionCharge",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "zeroDiffusionCharge",
                    dimensionSet(0,-3,0,0,0,1,0), 
                    0.0
                )
            )
        );

    const label nSpecies = species_.nSpecies();
    for (label i = 0; i < nSpecies; ++i)
    {
        elecDiffusionCharge.ref() += transportModels_[i].elecDiffusionCharge();
    }

    return elecDiffusionCharge;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //












