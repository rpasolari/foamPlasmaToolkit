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
    Foam::plasmaMobilityModel

Description
    Abstract base class for plasma species mobility modelling.

    This class defines the common interface for all mobility models used
    by plasma transport implementations such as drift–diffusion.

    Derived models provide the actual evaluation of the mobility μ, which
    may be constant, based on lookupTables, etc. Models are selected at runtime
    via the standard OpenFOAM run-time selection mechanism using the 
    'mobilityModel' entry in the species dictionary.

SourceFiles
    plasmaMobilityModelI.H
    plasmaMobilityModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_plasmaMobilityModel_H
#define Foam_plasmaMobilityModel_H

#include "fvMesh.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"

#include "plasmaSpecies.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class plasmaMobilityModel Declaration
\*---------------------------------------------------------------------------*/

class plasmaMobilityModel
{
protected:

    //- Model name (e.g. "constant", "lookupTable", etc.)
    const word modelName_;

    //- Mesh reference
    const fvMesh& mesh_;

    //- Reference to species database
    const plasmaSpecies& species_;

    //- Index of the species this model applies to
    const label specieIndex_;

    //- Mobility model dictionary
    const dictionary dict_;

public:

    //- Runtime type information
    TypeName("plasmaMobilityModel");

    //- Runtime selection table declaration
    declareRunTimeSelectionTable
    (
        autoPtr,
        plasmaMobilityModel,
        dictionary,
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex
        ),
        (modelName, dict, mesh, species, specieIndex)
    );

    // Constructors

        //- Construct base model
        plasmaMobilityModel
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex
        );

    //- Selector
    static autoPtr<plasmaMobilityModel> New
    (
        const word& modelName,
        const dictionary& dict,
        const fvMesh& mesh,
        const plasmaSpecies& species,
        const label specieIndex
    );

    //- Destructor
    virtual ~plasmaMobilityModel() {}

    // Member Functions

        //- Correct mobility
        virtual void correct(volScalarField& mu) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
