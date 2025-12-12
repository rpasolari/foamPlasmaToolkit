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
    Foam::ionSolidSurfaceFluxFvPatchScalarField

Description
    Concrete implementation of solidSurfaceFluxFvPatchScalarField for 
    ion species.

    This boundary condition calculates the ion flux to a solid surface 
    based on thermal motion.

Usage
    Example boundary condition entry:
    \verbatim
    <patchName>
    {
        type            ionSolidSurfaceFlux;

        // Electron temperature [K]
        T               Tion; 

        // Optional: Enable surface charge accumulation
        enableSurfaceCharging_ false;
    }
    \endverbatim
    
SourceFiles
    ionSolidSurfaceFluxFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef ionSolidSurfaceFluxFvPatchScalarField_H
#define ionSolidSurfaceFluxFvPatchScalarField_H

#include "solidSurfaceFluxFvPatchScalarField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class ionSolidSurfaceFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class ionSolidSurfaceFluxFvPatchScalarField
:
    public solidSurfaceFluxFvPatchScalarField
{
protected:

    // Protected Member Functions

        //- Calculate the ion-specific wall absorption velocity
        virtual tmp<scalarField> calcWallVelocity
        (
            const dimensionedScalar& m,
            const dimensionedScalar& T,
            const scalarField& uDriftNormal
        ) const override;

public:

    //- Runtime type information
    TypeName("ionSolidSurfaceFlux");

    // Constructors

        //- Construct from patch and internal field
        ionSolidSurfaceFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ionSolidSurfaceFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given BC onto a new patch
        ionSolidSurfaceFluxFvPatchScalarField
        (
            const ionSolidSurfaceFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        ionSolidSurfaceFluxFvPatchScalarField
        (
            const ionSolidSurfaceFluxFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        ionSolidSurfaceFluxFvPatchScalarField
        (
            const ionSolidSurfaceFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new ionSolidSurfaceFluxFvPatchScalarField(*this)
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
                new ionSolidSurfaceFluxFvPatchScalarField(*this, iF)
            );
        }

        // Member Functions

            //- Write dictionary entries
            virtual void write(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
