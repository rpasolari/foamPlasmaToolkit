const wordList fluidNames(rp["fluid"]);

PtrList<dynamicFvMesh> fluidRegions(fluidNames.size());

forAll(fluidNames, i)
{
    Info << "Create fluid mesh for region " << fluidNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

    fluidRegions.set
    (
        i,
        dynamicFvMesh::New
        (
            IOobject
            (
                fluidNames[i],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).ptr()
    );

    // Handle Dry-Run flags matching the gasMesh implementation
    if (args.dryRun() || args.found("dry-run-write"))
    {
        Info<< "Operating in 'dry-run' mode for region " 
            << fluidNames[i] << endl;
        
        FieldBase::allowConstructFromLargerSize = true;

        if (args.found("dry-run-write"))
        {
            // Using saWriteNow triggers function objects execute(), write()
            runTime.stopAt(Time::saWriteNow);

            // Ensure mesh gets output to the current time step
            fluidRegions[i].setInstance(runTime.timeName());
        }
        else
        {
            // Using saNoWriteNow triggers function objects execute(), but not write()
            runTime.stopAt(Time::saNoWriteNow);
        }

        functionObject::outputPrefix = "postProcessing-dry-run";
    }
}
