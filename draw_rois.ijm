 // draw_rois.ijm
//
// Semi-automated ROI drawing helper.
// First checks that ALL supported images are multi-channel.
// If any image has only one channel, the macro stops before
// any ROI drawing begins.

TRACE_CHANNEL_PREFIX = "C1-";   // <-- set to your membrane/red channel prefix

// ================================================================
// SELECT FOLDERS
// ================================================================

imageDir = getDirectory("Choose the folder containing images");
roiDir   = getDirectory("Choose (or create) the folder to save ROIs into");

list = getFileList(imageDir);

setTool("rectangle");


// ================================================================
// STEP 1: CHECK ALL IMAGES ARE MULTI-CHANNEL
// ================================================================

print("========================================");
print("Checking all images for multiple channels");
print("========================================");

singleChannelCount = 0;
checkedImages = 0;

for (i = 0; i < list.length; i++) {

    filename = list[i];

    if (!(endsWith(filename, ".dv") ||
          endsWith(filename, ".tif") ||
          endsWith(filename, ".tiff"))) {
        continue;
    }

    checkedImages++;

    path = imageDir + filename;

    print("Checking: " + filename);

    // Open image temporarily
    run("Bio-Formats Importer",
        "open=[" + path + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");

    setMinAndMax(0, 65535);

    // Get hyperstack dimensions
    Stack.getDimensions(width, height, channels, slices, frames);

    print("  Channels: " + channels);

    if (channels <= 1) {
        singleChannelCount++;

        print("  *** WARNING: SINGLE-CHANNEL IMAGE ***");
    }

    // Close the image before checking the next one
    close();
}


// ================================================================
// STOP IF ANY SINGLE-CHANNEL IMAGES WERE FOUND
// ================================================================

if (singleChannelCount > 0) {

    print("");
    print("========================================");
    print("ERROR: SINGLE-CHANNEL IMAGES FOUND");
    print("========================================");
    print("Number of single-channel images: " +
        singleChannelCount);
    print("");
    print("No ROIs have been processed.");
    print("All images must contain multiple channels.");
    print("========================================");

    showMessage(
        "Multi-channel check failed",
        "Found " + singleChannelCount +
        " single-channel image(s).\n\n" +
        "No ROIs have been processed.\n\n" +
        "See the Log window for the filenames."
    );

    exit();
}


// ================================================================
// ALL IMAGES PASSED
// ================================================================

print("");
print("========================================");
print("Multi-channel check PASSED");
print("Images checked: " + checkedImages);
print("Starting ROI drawing...");
print("========================================");
print("");

skipped = newArray(0);


// ================================================================
// STEP 2: PROCESS IMAGES
// ================================================================

for (i = 0; i < list.length; i++) {

    filename = list[i];

    if (!(endsWith(filename, ".dv") ||
          endsWith(filename, ".tif") ||
          endsWith(filename, ".czi") ||
          endsWith(filename, ".tiff"))) {
        continue;
    }

    path = imageDir + filename;

    run("Bio-Formats Importer",
        "open=[" + path + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");

    title = getTitle();

    run("Split Channels");

    print("Open windows after split:");
    windowList = getList("image.titles");
    for (w = 0; w < windowList.length; w++) {
        print("  " + windowList[w]);
    }

    traceWindow = TRACE_CHANNEL_PREFIX + title;

    if (!isOpen(traceWindow)) {

        print("WARNING: expected window '" +
            traceWindow +
            "' not found - check TRACE_CHANNEL_PREFIX.");

        print("Skipping: " + filename);

        run("Close All");

        continue;
    }


    // ------------------------------------------------------------
    // SHOW TRACE CHANNEL
    // ------------------------------------------------------------

    selectWindow(traceWindow);
    setBatchMode(false);


    // ------------------------------------------------------------
    // WAIT FOR USER TO DRAW ROI
    // ------------------------------------------------------------

    waitForUser(
        "Draw ROI: " + filename,
        "Draw a rectangle ROI on the image now.\n" +
        "Click OK when done.\n\n" +
        "To skip this image, click OK WITHOUT drawing a selection."
    );


    // ------------------------------------------------------------
    // CHECK FOR ROI
    // ------------------------------------------------------------

    if (selectionType() == -1) {

        print("No selection drawn for " +
            filename +
            " - skipping.");

        skipped = Array.concat(skipped, filename);

        run("Close All");

        continue;
    }


    // ------------------------------------------------------------
    // SAVE ROI
    // ------------------------------------------------------------

    roiManager("reset");

    roiManager("Add");

    base = filename;

    base = replace(base, ".dv", "");
    base = replace(base, ".tif", "");
    base = replace(base, ".tiff", "");

    savePath = roiDir + base + "_RoiSet.zip";

    roiManager("Select", 0);

    roiManager("Save", savePath);

    print("Saved ROI for " +
        filename +
        " -> " +
        savePath);

    run("Close All");
}


// ================================================================
// FINISHED
// ================================================================

print("");
print("========================================");
print("DONE");
print("========================================");

if (skipped.length > 0) {

    print("Skipped " +
        skipped.length +
        " image(s) because no ROI was drawn:");

    for (i = 0; i < skipped.length; i++) {
        print("  " + skipped[i]);
    }
}