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
    Foam::plasmaDiffusivityModel

Description
    Abstract base class for plasma species diffusivity modelling.

    This class defines the common interface for all diffusivity models used
    by plasma transport implementations such as drift–diffusion.

    Derived models provide the actual evaluation of the diffusivity D, which
    may be constant, based on lookupTables, etc. Models are selected at runtime
    via the standard OpenFOAM run-time selection mechanism using the 
    'diffusivityModel' entry in the species dictionary.

SourceFiles
    plasmaDiffusivityModelI.H
    plasmaDiffusivityModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_plasmaDiffusivityModel_H
#define Foam_plasmaDiffusivityModel_H

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
                    Class plasmaDiffusivityModel Declaration
\*---------------------------------------------------------------------------*/

class plasmaDiffusivityModel
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

    //- Diffusivity model dictionary
    const dictionary dict_;

public:

    //- Runtime type information
    TypeName("plasmaDiffusivityModel");

    //- Runtime selection table declaration
    declareRunTimeSelectionTable
    (
        autoPtr,
        plasmaDiffusivityModel,
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
        plasmaDiffusivityModel
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex
        );

    //- Selector
    static autoPtr<plasmaDiffusivityModel> New
    (
        const word& modelName,
        const dictionary& dict,
        const fvMesh& mesh,
        const plasmaSpecies& species,
        const label specieIndex
    );

    //- Destructor
    virtual ~plasmaDiffusivityModel() {}

    // Member Functions

        //- Correct diffusivity
        virtual void correct(volScalarField& D) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
