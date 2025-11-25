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

inline const Foam::dictionary& 
Foam::plasmaSpecies::speciesDict(label i) const
{
    return speciesDicts_[speciesNames_[i]];
}

inline const Foam::dictionary& 
Foam::plasmaSpecies::speciesDict(const word& name) const
{
    return speciesDicts_[name];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
