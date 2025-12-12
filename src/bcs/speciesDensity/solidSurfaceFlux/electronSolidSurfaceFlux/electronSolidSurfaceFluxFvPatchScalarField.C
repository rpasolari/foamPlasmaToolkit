/*---------------------------------------------------------------------------*\
  File: electronSolidSurfaceFluxFvPatchScalarField.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::electronSolidSurfaceFluxFvPatchScalarField.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

#include "electronSolidSurfaceFluxFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Register the class to the runtime selection table
makePatchTypeField
(
    fvPatchScalarField,
    electronSolidSurfaceFluxFvPatchScalarField
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> electronSolidSurfaceFluxFvPatchScalarField::calcWallVelocity
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

    if(includeDriftFlux_)
    {
        return uWallRaw + max(uDriftNormal, scalar(0.0));
    }
    else
    {
        // Return as a field (uniform value across the patch)
        return tmp<scalarField>::New(this->size(), uWallRaw);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Standard Constructor
electronSolidSurfaceFluxFvPatchScalarField::
electronSolidSurfaceFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    solidSurfaceFluxFvPatchScalarField(p, iF),
    includeDriftFlux_(false)
{}

// Dictionary Constructor
electronSolidSurfaceFluxFvPatchScalarField::
electronSolidSurfaceFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    solidSurfaceFluxFvPatchScalarField(p, iF, dict),
    includeDriftFlux_(dict.lookupOrDefault<bool>("includeDriftFlux", false))
{}

// Mapping Constructor
electronSolidSurfaceFluxFvPatchScalarField::
electronSolidSurfaceFluxFvPatchScalarField
(
    const electronSolidSurfaceFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidSurfaceFluxFvPatchScalarField(ptf, p, iF, mapper),
    includeDriftFlux_(ptf.includeDriftFlux_)
{}

// Copy Constructor
electronSolidSurfaceFluxFvPatchScalarField::
electronSolidSurfaceFluxFvPatchScalarField
(
    const electronSolidSurfaceFluxFvPatchScalarField& ptf
)
:
    solidSurfaceFluxFvPatchScalarField(ptf),
    includeDriftFlux_(ptf.includeDriftFlux_)
{}

// Copy Constructor (with new internal field)
electronSolidSurfaceFluxFvPatchScalarField::
electronSolidSurfaceFluxFvPatchScalarField
(
    const electronSolidSurfaceFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    solidSurfaceFluxFvPatchScalarField(ptf, iF),
    includeDriftFlux_(ptf.includeDriftFlux_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void electronSolidSurfaceFluxFvPatchScalarField::write(Ostream& os) const
{
    solidSurfaceFluxFvPatchScalarField::write(os);   

    os.writeEntry("includeDriftFlux", includeDriftFlux_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
