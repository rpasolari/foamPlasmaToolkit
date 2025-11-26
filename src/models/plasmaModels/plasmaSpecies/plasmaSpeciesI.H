/*---------------------------------------------------------------------------*\
  File: plasmaSpeciesI.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Inline member function definitions for the plasmaSpecies class.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef Foam_plasmaSpeciesI_H
#define Foam_plasmaSpeciesI_H

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Mesh and species metadata access
inline const Foam::fvMesh&
Foam::plasmaSpecies::mesh() const
{
    return mesh_;
}

inline const Foam::wordList&
Foam::plasmaSpecies::speciesNames() const
{
    return speciesNames_;
}

inline Foam::label
Foam::plasmaSpecies::nSpecies() const
{
    return nSpecies_;
}

inline const Foam::List<Foam::scalar>&
Foam::plasmaSpecies::speciesCharge() const
{
    return speciesCharge_;
}


//- Number density access
inline const Foam::volScalarField&
Foam::plasmaSpecies::numberDensity(Foam::label i) const
{
    return numberDensity_[i];
}

inline Foam::volScalarField&
Foam::plasmaSpecies::numberDensity(Foam::label i)
{
    return numberDensity_[i];
}

inline const Foam::PtrList<Foam::volScalarField>&
Foam::plasmaSpecies::numberDensity() const
{
    return numberDensity_;
}

inline Foam::PtrList<Foam::volScalarField>&
Foam::plasmaSpecies::numberDensity()
{
    return numberDensity_;
}


//- Flux access
inline const Foam::volVectorField&
Foam::plasmaSpecies::flux(Foam::label i) const
{
    return flux_[i];
}


inline Foam::volVectorField&
Foam::plasmaSpecies::flux(Foam::label i)
{
    return flux_[i];
}


inline const Foam::PtrList<Foam::volVectorField>&
Foam::plasmaSpecies::flux() const
{
    return flux_;
}


inline Foam::PtrList<Foam::volVectorField>&
Foam::plasmaSpecies::flux()
{
    return flux_;
}


//- Species dictionaries access
inline const Foam::dictionary&
Foam::plasmaSpecies::speciesDict(Foam::label i) const
{
    return speciesDicts_[speciesNames_[i]];
}

inline const Foam::dictionary&
Foam::plasmaSpecies::speciesDict(const Foam::word& name) const
{
    return speciesDicts_[name];
}

inline const Foam::dictionary&
Foam::plasmaSpecies::defaultDict() const
{
    return defaultDict_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
