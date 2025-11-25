/*---------------------------------------------------------------------------*\
  File: constantMobility.C
  Part of: foamPlasmaToolkit
\*---------------------------------------------------------------------------*/

#include "constantMobility.H"
#include "volFields.H"

namespace Foam
{

// * * * * * * * * * Runtime Type Information * * * * * * * * * //

defineTypeNameAndDebug(constantMobility, 0);
addToRunTimeSelectionTable(mobilityModel, constantMobility, dictionary);


// * * * * * * * * * * Constructors * * * * * * * * * * //

constantMobility::constantMobility
(
    const dictionary& coeffs,
    const fvMesh& mesh
)
:
    mobilityModel(coeffs, mesh),
    mu_(coeffs.lookupOrDefault<scalar>("mu", 0.0))
{}


// * * * * * * * * * * Member Functions * * * * * * * * * * //

tmp<volScalarField> constantMobility::mu() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "mu",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("mu", dimless/dimTime, mu_)
        )
    );
}

} // End namespace Foam
// ************************************************************************* //
