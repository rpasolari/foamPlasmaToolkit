/*---------------------------------------------------------------------------*\
License
    This file is part of the foamPlasmaToolkit.

    The foamPlasmaToolkit is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2025 Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::solidSurfaceFluxFvPatchScalarField

Description
    Abstract base class for scalar flux boundary conditions at solid surfaces,
    derived from mixedFvPatchScalarField.

    This class serves as a foundation for specific particle or charge flux 
    conditions (e.g., ion or electron fluxes) where particles impact a wall.

    This is an abstract class. It cannot be used directly in a case file
    as the 'type'. Derived classes must implement the pure virtual function.

Usage
    This class is not instantiated directly. Derived classes are instantiated.
    
SourceFiles
    solidSurfaceFluxFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_solidSurfaceFluxFvPatchScalarField_H
#define Foam_solidSurfaceFluxFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of plasmaTransport class
class plasmaTransport;

/*---------------------------------------------------------------------------*\
         Class solidSurfaceFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class solidSurfaceFluxFvPatchScalarField
:
    public mixedFvPatchScalarField
{
protected:

    // Protected Data

        //- Constant Temperature Value
        dimensionedScalar T_;

        //- Switch to enable surface charging
        Switch enableSurfaceCharging_;

    // Protected Member Functions

        //- Calculate thermal velocity
        dimensionedScalar calcThermalVelocity
        (
            const dimensionedScalar& m,
            const dimensionedScalar& T
        ) const;

        //- Calculate the effective wall absorption velocity
        virtual tmp<scalarField> calcWallVelocity
        (
            const dimensionedScalar& m,
            const dimensionedScalar& T,
            const scalarField& uDriftNormal
        ) const = 0;

public:

    //- Runtime type information
    TypeName("solidSurfaceFlux");

    // Constructors

        //- Construct from patch and internal field
        solidSurfaceFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidSurfaceFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given BC onto a new patch
        solidSurfaceFluxFvPatchScalarField
        (
            const solidSurfaceFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidSurfaceFluxFvPatchScalarField
        (
            const solidSurfaceFluxFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        solidSurfaceFluxFvPatchScalarField
        (
            const solidSurfaceFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const = 0;

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const = 0;

        // Member Functions

            //- Update the coefficients
            virtual void updateCoeffs();

            //- Get the switch enableSurfaceCharging
            bool enableSurfaceCharging() const;

            //- Write dictionary entries
            virtual void write(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
