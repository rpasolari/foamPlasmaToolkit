/*---------------------------------------------------------------------------*\
License
    This file is part of the foamPlasmaToolkit.

    The foamPlasmaToolkit is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2025
        Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::mobilityModel

Description
    Abstract base class for plasma species mobility modelling.

    This class defines the common interface for all mobility models used
    by plasma transport implementations such as drift–diffusion.

    Derived models provide the actual evaluation of the mobility μ, which
    may be constant, computed from Einstein-type relations, or based on
    more complex physics. Models are selected at runtime via the standard
    OpenFOAM run-time selection mechanism using the 'mobilityModel' entry
    in the species dictionary.

SourceFiles
    mobilityModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_mobilityModel_H
#define Foam_mobilityModel_H

#include "fvMesh.H"
#include "dictionary.H"
#include "volFields.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class mobilityModel Declaration
\*---------------------------------------------------------------------------*/

class mobilityModel
{
protected:

    const fvMesh& mesh_;
    dictionary coeffs_;

    mobilityModel
    (
        const dictionary& coeffs,
        const fvMesh& mesh
    )
    :
        mesh_(mesh),
        coeffs_(coeffs)
    {}

public:

    TypeName("mobilityModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        mobilityModel,
        dictionary,
        (const dictionary& coeffs, const fvMesh& mesh),
        (coeffs, mesh)
    );

    static autoPtr<mobilityModel> New
    (
        const word& modelName,
        const dictionary& coeffs,
        const fvMesh& mesh
    );

    virtual ~mobilityModel() {}

    virtual tmp<volScalarField> mu() const = 0;

    inline const fvMesh& mesh() const { return mesh_; }
    inline const dictionary& coeffs() const { return coeffs_; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
