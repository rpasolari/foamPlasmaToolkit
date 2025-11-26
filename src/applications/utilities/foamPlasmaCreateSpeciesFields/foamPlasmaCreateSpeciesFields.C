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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Application
    foamPlasmaCreateSpeciesFields

Description
    Utility to generate template number-density fields (n_speciesName) in the
    0/ directory for all species listed in the plasmaSpecies dictionary.

    The tool:
      - detects and uses the gas region (if present in a multi-region case),
      - reads species from the plasmaSpecies dictionary,
      - loads the mesh of the region,
      - creates volScalarField files with default values and zeroGradient
        boundary conditions for all patches.

    These generated fields are **only templates**.
    The user must review and modify the internalField values, boundary
    conditions, and any other field settings according to the needs of the
    specific plasma simulation.

Usage
    \b plasmaDielectricFoam [OPTIONS]

    Example:
        foamPlasmaCreateSpeciesFields -case plasmaCase

Author
    Rention Pasolari
    Contact: r.pasolari@gmail.com
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "regionProperties.H"
#include "fvMesh.H"
#include "IOdictionary.H"
#include <fstream>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    argList::addNote
    (
        "Create species fields in 0 folder"
    );

    argList::noBanner();
    argList::noParallel();
    argList::noJobInfo();
    argList::noFunctionObjects();

    argList args(argc, argv);

    // Create minimal Time
    Time runTime
    (
        Time::controlDictName,
        args.rootPath(),
        args.caseName()
    );

    word gasRegionName;

    fileName regionPropsFile = runTime.constant()/"regionProperties";

    // Find the gas region, or default region in single region cases
    if (isFile(regionPropsFile))
    {
        Info << "Found regionProperties: "
            << regionPropsFile << nl;

        regionProperties rp(runTime);

        if (rp.found("gas"))
        {
            const wordList& names = rp["gas"];

            if (names.size() == 0)
            {
                FatalErrorInFunction
                    << "A 'gas' region type is defined in regionProperties, "
                    << "but no region names are listed under it." << nl
                    << exit(FatalError);
            }

            if (names.size() > 1)
            {
                FatalErrorInFunction
                    << "Multiple gas regions detected in "
                    << "constant/regionProperties:" << nl
                    << "  gas: " << names << nl
                    << "This utility supports ONLY ONE gas region." << nl
                    << exit(FatalError);
            }

            // Exactly one gas region
            gasRegionName = names[0];

            Info << "Detected gas region: "
                << gasRegionName << nl;
        }
    }

    // Report region used
    if (gasRegionName.empty())
    {
        Info << "No gas region found → using constant/ as search path." << nl;
    }
    else
    {
        Info << "Using region: " << gasRegionName << nl;
    }

    fileName plasmaSpeciesPath;

    // Find the plasmaSpecies file in the case
    if (!gasRegionName.empty())
    {
        // Multi-region path
        fileName regionDir = runTime.constant()/gasRegionName;

        if (isFile(regionDir/"plasmaSpecies"))
        {
            plasmaSpeciesPath = regionDir/"plasmaSpecies";
        }
        else
        {
            FatalErrorInFunction
                << "Region '" << gasRegionName << "' exists but has no "
                << "plasmaSpecies file:" << nl
                << "  " << regionDir/"plasmaSpecies" << nl
                << exit(FatalError);
        }
    }
    else
    {
        // Single-region path
        if (isFile(runTime.constant()/"plasmaSpecies"))
        {
            plasmaSpeciesPath = runTime.constant()/"plasmaSpecies";
        }
        else
        {
            FatalErrorInFunction
                << "No 'constant/plasmaSpecies' file found." << nl 
                << exit(FatalError);
        }
    }

    // Read plasmaSpecies dict
    Info << "Reading plasmaSpecies from: " << plasmaSpeciesPath << nl;

    IOdictionary speciesDict
    (
        IOobject
        (
            "plasmaSpecies",
            plasmaSpeciesPath.path(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (!speciesDict.found("species"))
    {
        FatalErrorInFunction
            << "The plasmaSpecies dictionary has no 'species' entry!"
            << exit(FatalError);
    }

    wordList species(speciesDict.lookup("species"));

    Info << "Found " << species.size() << " species:" << nl;

    for (const word& s : species)
    {
        Info << "  - " << s << nl;
    }

    // Determine region path for mesh creation
    word regionToLoad = gasRegionName;

    // Create the mesh
    fvMesh mesh
    (
        IOobject
        (
            regionToLoad,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info << "Boundaries: " << nl;

    // Loop over patches
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    for (const polyPatch& p : patches)
    {
        Info << "  - " << p.name() << nl;
    }

    // Determine correct "0" directory
    fileName zeroDir;

    if (gasRegionName.empty())
    {
        zeroDir = runTime.path()/"0";
    }
    else
    {
        zeroDir = runTime.path()/"0"/gasRegionName;
    }

    // Ensure directory exists
    mkDir(zeroDir);

    Info << nl << "Creating species number-density fields in: "
         << zeroDir << nl;

    // Loop over species
    for (const word& s : species)
    {
        const word fieldName = "n_" + s;
        const fileName fieldPath = zeroDir/fieldName;

        Info << "  - Creating field: " << fieldName << nl;

        // Construct volScalarField
        IOdictionary header
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                gasRegionName,  // empty for default region
                runTime,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            )
        );

        // Write to file manually using standard C++ ofstream
        std::ofstream os(fieldPath.c_str());

        if (!os)
        {
            FatalErrorInFunction
                << "Cannot open file for writing: " << fieldPath << nl
                << exit(FatalError);
        }

        os << "/*--------------------------------*- C++ -*----------------------------------*\\\n"
              "| =========                 |                                                 |\n"
              "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
              "|  \\\\    /   O peration     | Version:  v2412                                 |\n"
              "|   \\\\  /    A nd           | Website:  www.openfoam.com                      |\n"
              "|    \\\\/     M anipulation  |                                                 |\n"
              "\\*---------------------------------------------------------------------------*/\n";

        os << "FoamFile\n"
              "{\n"
              "    version     2.0;\n"
              "    format      ascii;\n"
              "    class       volScalarField;\n"
              "    object      " << fieldName << ";\n"
              "}\n"
              "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

        os << "dimensions      [0 -3 0 0 0 0 0];\n\n";
        os << "internalField   uniform 0;\n\n";

        os << "boundaryField\n{\n";

        // Write each boundary patch with zeroGradient
        const polyBoundaryMesh& patches2 = mesh.boundaryMesh();
        for (const polyPatch& p : patches2)
        {
            os << "    " << p.name() << "\n"
               << "    {\n"
               << "        type    zeroGradient;\n"
               << "    }\n\n";
        }

        os << "}\n\n";
        os << "// ************************************************************************* //\n";

        os.close();
    }

    return 0;
}
