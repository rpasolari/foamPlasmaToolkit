#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "fvMesh.H"

using namespace Foam;

int main(int argc, char* argv[])
{
    argList::noParallel();
    argList args(argc, argv);

    // Create minimal Time object (needed for reading constant/)
    Time runTime
    (
        Time::controlDictName,
        args.rootPath(),
        args.caseName()
    );

    fileName regionPropsFile = runTime.constant()/"regionProperties";
    word gasRegionName;

    if (isFile(regionPropsFile))
    {
        Info << "Found regionProperties: " << regionPropsFile << nl;

        IOdictionary regionProps
        (
            IOobject
            (
                "regionProperties",
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        if (regionProps.found("regions"))
        {
            const dictionary& regions = regionProps.subDict("regions");

            for (const entry& e : regions)
            {
                if (!e.isDict()) continue;

                const word& regionName = e.keyword();
                const dictionary& sub = e.dict();

                word type = sub.lookupOrDefault<word>("type", "unknown");

                if (type == "gas" || type == "fluid")
                {
                    gasRegionName = regionName;
                    Info << "Detected gas/fluid region: " << gasRegionName << nl;
                    break;
                }
            }
        }
    }

    // If regionProperties NOT found OR no gas region defined → default
    if (gasRegionName.empty())
    {
        Info << "No gas region detected. Assuming single-region case." << nl;
        gasRegionName = fvMesh::defaultRegion;
    }

    Info << "Using region: " << gasRegionName << nl;


    // ------------------------------------------------------------
    // 2. Locate plasmaSpecies dictionary
    // ------------------------------------------------------------

    fileName plasmaSpeciesPath;

    // Try multi-region location first: constant/<region>/plasmaSpecies
    fileName regionDir = runTime.constant()/gasRegionName;

    if (isFile(regionDir/"plasmaSpecies"))
    {
        plasmaSpeciesPath = regionDir/"plasmaSpecies";
        Info << "Found plasmaSpecies at: " << plasmaSpeciesPath << nl;
    }
    // Try single-region location: constant/plasmaSpecies
    else if (isFile(runTime.constant()/"plasmaSpecies"))
    {
        plasmaSpeciesPath = runTime.constant()/"plasmaSpecies";
        Info << "Found plasmaSpecies at: " << plasmaSpeciesPath << nl;
    }
    else
    {
        Warning
            << "No plasmaSpecies file found in:\n"
            << " - " << regionDir/"plasmaSpecies" << "\n"
            << " - " << runTime.constant()/"plasmaSpecies" << "\n"
            << "Continuing without species information.\n" << nl;
    }

    // ------------------------------------------------------------
    // 3. If plasmaSpecies exists → load it
    // ------------------------------------------------------------

    if (!plasmaSpeciesPath.empty())
    {
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

        Info << "Loaded species dictionary with "
             << speciesDict.size() << " entries." << nl;
    }

    return 0;
}
