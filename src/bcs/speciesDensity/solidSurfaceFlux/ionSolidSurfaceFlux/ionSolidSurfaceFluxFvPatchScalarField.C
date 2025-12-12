/*---------------------------------------------------------------------------*\
  File: ionSolidSurfaceFluxFvPatchScalarField.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::ionSolidSurfaceFluxFvPatchScalarField.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

#include "ionSolidSurfaceFluxFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Register the class to the runtime selection table
makePatchTypeField
(
    fvPatchScalarField,
    ionSolidSurfaceFluxFvPatchScalarField
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> ionSolidSurfaceFluxFvPatchScalarField::calcWallVelocity
(
    const dimensionedScalar& m,
    const dimensionedScalar& T,
    const scalarField& uDriftNormal
) const
{
    // Calculate the base thermal velocity
    dimensionedScalar uTh = calcThermalVelocity(m, T);

    // u_wall =  u_th/4
    scalar uWallRaw = 0.25 * uTh.value();

    // Add drift velocity only if it is positive
    scalarField uDriftSurf = max(uDriftNormal, scalar(0.0));

    // Return as a field (uniform value across the patch)
    return uWallRaw + uDriftSurf;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Standard Constructor
ionSolidSurfaceFluxFvPatchScalarField::
ionSolidSurfaceFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    solidSurfaceFluxFvPatchScalarField(p, iF)
{}

// Dictionary Constructor
ionSolidSurfaceFluxFvPatchScalarField::
ionSolidSurfaceFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    solidSurfaceFluxFvPatchScalarField(p, iF, dict)
{}

// Mapping Constructor
ionSolidSurfaceFluxFvPatchScalarField::
ionSolidSurfaceFluxFvPatchScalarField
(
    const ionSolidSurfaceFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidSurfaceFluxFvPatchScalarField(ptf, p, iF, mapper)
{}

// Copy Constructor
ionSolidSurfaceFluxFvPatchScalarField::
ionSolidSurfaceFluxFvPatchScalarField
(
    const ionSolidSurfaceFluxFvPatchScalarField& ptf
)
:
    solidSurfaceFluxFvPatchScalarField(ptf)
{}

// Copy Constructor (with new internal field)
ionSolidSurfaceFluxFvPatchScalarField::
ionSolidSurfaceFluxFvPatchScalarField
(
    const ionSolidSurfaceFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    solidSurfaceFluxFvPatchScalarField(ptf, iF)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ionSolidSurfaceFluxFvPatchScalarField::write(Ostream& os) const
{
    solidSurfaceFluxFvPatchScalarField::write(os);   
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
