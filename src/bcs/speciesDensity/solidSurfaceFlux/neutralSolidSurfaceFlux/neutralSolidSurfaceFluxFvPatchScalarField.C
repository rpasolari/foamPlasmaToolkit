/*---------------------------------------------------------------------------*\
  File: neutralSolidSurfaceFluxFvPatchScalarField.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::neutralSolidSurfaceFluxFvPatchScalarField.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

#include "neutralSolidSurfaceFluxFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Register the class to the runtime selection table
makePatchTypeField
(
    fvPatchScalarField,
    neutralSolidSurfaceFluxFvPatchScalarField
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> neutralSolidSurfaceFluxFvPatchScalarField::calcWallVelocity
(
    const dimensionedScalar& m,
    const dimensionedScalar& T,
    const scalarField& uDriftNormal
) const
{
    // Calculate the base thermal velocity
    dimensionedScalar uTh = calcThermalVelocity(m, T);

    // u =  u_th/4
    scalar uWallRaw = 0.25 * uTh.value();

    // Return as a field (uniform value across the patch)
    return tmp<scalarField>::New(this->size(), uWallRaw);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Standard Constructor
neutralSolidSurfaceFluxFvPatchScalarField::
neutralSolidSurfaceFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    solidSurfaceFluxFvPatchScalarField(p, iF)
{
    if (this->enableSurfaceCharging_)
    {
        WarningInFunction
            << "Surface charging was enabled in dictionary for patch '" 
            << p.name() << "', but species is neutral." << nl
            << "    Forcing enableSurfaceCharging to FALSE."
            << endl;

        this->enableSurfaceCharging_ = false;
    }
}

// Dictionary Constructor
neutralSolidSurfaceFluxFvPatchScalarField::
neutralSolidSurfaceFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    solidSurfaceFluxFvPatchScalarField(p, iF, dict)
{
    if (this->enableSurfaceCharging_)
    {
        WarningInFunction
            << "Surface charging was enabled in dictionary for patch '" 
            << p.name() << "', but species is neutral." << nl
            << "    Forcing enableSurfaceCharging to FALSE."
            << endl;

        this->enableSurfaceCharging_ = false;
    }
}

// Mapping Constructor
neutralSolidSurfaceFluxFvPatchScalarField::
neutralSolidSurfaceFluxFvPatchScalarField
(
    const neutralSolidSurfaceFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidSurfaceFluxFvPatchScalarField(ptf, p, iF, mapper)
{
    if (this->enableSurfaceCharging_)
    {
        WarningInFunction
            << "Surface charging was enabled in dictionary for patch '" 
            << p.name() << "', but species is neutral." << nl
            << "    Forcing enableSurfaceCharging to FALSE."
            << endl;

        this->enableSurfaceCharging_ = false;
    }
}

// Copy Constructor
neutralSolidSurfaceFluxFvPatchScalarField::
neutralSolidSurfaceFluxFvPatchScalarField
(
    const neutralSolidSurfaceFluxFvPatchScalarField& ptf
)
:
    solidSurfaceFluxFvPatchScalarField(ptf)
{
    this->enableSurfaceCharging_ = false;
}

// Copy Constructor (with new internal field)
neutralSolidSurfaceFluxFvPatchScalarField::
neutralSolidSurfaceFluxFvPatchScalarField
(
    const neutralSolidSurfaceFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    solidSurfaceFluxFvPatchScalarField(ptf, iF)
{
    this->enableSurfaceCharging_ = false;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void neutralSolidSurfaceFluxFvPatchScalarField::write(Ostream& os) const
{
    solidSurfaceFluxFvPatchScalarField::write(os);   
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
