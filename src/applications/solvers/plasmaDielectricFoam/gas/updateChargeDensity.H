/*---------------------------------------------------------------------------*\
  File: updateChargeDensity.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

Description:
    Updates the charge density field by summing the contributions from the
    number densities of all participating species.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

chargeDensity = dimensionedScalar("zero", chargeDensity.dimensions(),0);

for (label i = 0; i < nSpecies; ++i)
{
    chargeDensity += 
        speciesCharge[i]
      * constant::plasma::eCharge
      * numberDensity[i]
}
