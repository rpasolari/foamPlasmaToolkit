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
    regIOobject
    (
        IOobject
        (
            "plasmaTransport",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    species_(species),
    E_(E),
    transportModels_(species.nSpecies()),
    particleFlux_(),
    Gamma_()
{
    constructModels();

    particleFlux_.setSize(species.nSpecies());
    Gamma_.setSize(species.nSpecies());

    for (label i = 0; i < species.nSpecies(); ++i)
    {
        // Surface Scalar Field
        particleFlux_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "particleFlux_" + species.speciesNames()[i], 
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "zero",
                    dimensionSet(0, 0, -1, 0, 0, 0, 0), 
                    0.0
                )
            )
        );

        // Volumetric Vector Field
        Gamma_.set
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
    }
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void plasmaTransport::correct()
{
    const label nSpecies = species_.nSpecies();

    for (label i = 0; i < nSpecies; ++i)
    {
        // Update coefficients (mobility, diffusivity, driftVelocity etc.)
        transportModels_[i].correct();

        // Build the matrix
        tmp<fvScalarMatrix> tEqn = transportModels_[i].nEqn();
        fvScalarMatrix& nEqn = tEqn.ref();

        // Solve the transport equation
        nEqn.solve();

        // Update fluxes
        particleFlux_[i] = transportModels_[i].particleFlux();
        Gamma_[i] = fvc::reconstruct(particleFlux_[i]);
    }
}

tmp<volScalarField> plasmaTransport::electricalConductivity() const
{
    tmp<volScalarField> tSigma
    (
        new volScalarField
        (
            IOobject
            (
                "electricalConductivity",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "zero", 
                dimensionSet(-1, -3, 3, 0, 0, 2, 0), 
                0.0
            )
        )
    );

    // Get reference
    volScalarField& sigma = tSigma.ref();

    const label nSpecies = species_.nSpecies();
    for (label i = 0; i < nSpecies; ++i)
    {
        sigma.ref() += transportModels_[i].electricalConductivity();
    }

    return tSigma;
}

tmp<volScalarField> plasmaTransport::diffusiveChargeSource() const
{
    tmp<volScalarField> tRhoDiff
    (
        new volScalarField
        (
            IOobject
            (
                "diffusiveChargeSource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "zero", 
                dimensionSet(0, -3, 0, 0, 0, 1, 0), 
                0.0
            )
        )
    );

    // Get reference
    volScalarField& rhoDiff = tRhoDiff.ref();

    const label nSpecies = species_.nSpecies();
    for (label i = 0; i < nSpecies; ++i)
    {
        rhoDiff.ref() += transportModels_[i].diffusiveChargeSource();
    }

    return tRhoDiff;
}

bool plasmaTransport::writeData(Ostream& os) const
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //












