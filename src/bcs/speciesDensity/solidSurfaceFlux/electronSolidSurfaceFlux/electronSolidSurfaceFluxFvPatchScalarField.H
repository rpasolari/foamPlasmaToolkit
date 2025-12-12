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
    Foam::electronSolidSurfaceFluxFvPatchScalarField

Description
    Concrete implementation of solidSurfaceFluxFvPatchScalarField for 
    electron species.

    This boundary condition calculates the electron flux to a solid surface 
    based on thermal motion.

Usage
    Example boundary condition entry:
    \verbatim
    <patchName>
    {
        type            electronSolidSurfaceFlux;

        // Electron temperature [K]
        T               Te; 

        // Optional: Enable surface charge accumulation
        enableSurfaceCharging_ false;
    }
    \endverbatim
    
SourceFiles
    electronSolidSurfaceFluxFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef electronSolidSurfaceFluxFvPatchScalarField_H
#define electronSolidSurfaceFluxFvPatchScalarField_H

#include "solidSurfaceFluxFvPatchScalarField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class electronSolidSurfaceFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class electronSolidSurfaceFluxFvPatchScalarField
:
    public solidSurfaceFluxFvPatchScalarField
{
protected:

    // Protected data

        //- Switch to include drift velocity in wall flux calculation
        Switch includeDriftFlux_;

    // Protected Member Functions

        //- Calculate the electron-specific wall absorption velocity
        virtual tmp<scalarField> calcWallVelocity
        (
            const dimensionedScalar& m,
            const dimensionedScalar& T,
            const scalarField& uDriftNormal
        ) const override;

public:

    //- Runtime type information
    TypeName("electronSolidSurfaceFlux");

    // Constructors

        //- Construct from patch and internal field
        electronSolidSurfaceFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        electronSolidSurfaceFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given BC onto a new patch
        electronSolidSurfaceFluxFvPatchScalarField
        (
            const electronSolidSurfaceFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        electronSolidSurfaceFluxFvPatchScalarField
        (
            const electronSolidSurfaceFluxFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        electronSolidSurfaceFluxFvPatchScalarField
        (
            const electronSolidSurfaceFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new electronSolidSurfaceFluxFvPatchScalarField(*this)
            );
        }

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new electronSolidSurfaceFluxFvPatchScalarField(*this, iF)
            );
        }

        // Member Functions

            //- Write dictionary entries
            virtual void write(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
