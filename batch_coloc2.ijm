// batch_coloc2.ijm
//
// ROI-driven batch Coloc2 analysis.
//
// The ROI directory is the source of truth:
//   - Only images that have a matching *_RoiSet.zip are processed.
//   - Images in imageDir without an ROI are ignored.
//   - ROIs without a corresponding image are reported and skipped.
//
// Runs Coloc2 on every pairwise combination of channels 1, 2, 3
// (channel 4 is excluded entirely): C1 vs C2, C1 vs C3, C2 vs C3.
// One row is written per (image, channel pair) to coloc_results.csv.
//
// Expected naming:
//   image: example.dv
//   ROI:   example_RoiSet.zip
//
// SETUP NOTES:
//   - CHANNEL_PREFIXES must match the channel naming Fiji assigns
//     after "Split Channels" for channels 1, 2, 3 (channel 4 is
//     deliberately not included here).
//   - Display range is left at Fiji's default (resetMinAndMax()) on
//     import rather than being forced to a fixed range - this keeps
//     the "color balance" at whatever Fiji's default reset state is
//     for that bit depth, rather than the previous behaviour of
//     forcing every image to a hardcoded 0-65535 display range.
//   - "manders'_correlation" is explicitly included in the Coloc 2
//     run call below - without it, Coloc2 does not compute Manders'
//     coefficients at all, and they simply won't appear in the log
//     (this was previously causing M1/M2 to always come out as "NA").
//   - Coloc2's log output is COMMA-separated ("Label, value"), not
//     colon-separated - parseColocValue() below splits on the last
//     comma in the matching line. If a value still comes out "NA",
//     print(logText) on one test image/pair and check the exact
//     label wording matches what parseColocValue() is searching for.
//   - BACKGROUND_SUBTRACT / ROLLING_BALL_RADIUS / USE_SLIDING_PARABOLOID:
//     applies background subtraction to each of channels 1-3 (once
//     per channel, before any pairwise Coloc2 run) - this removes a
//     constant/uneven background offset, which otherwise can push
//     Costes' auto-threshold up near the channel mean and trigger
//     "y-intercept far from zero" / "threshold too high" warnings.
//     Set BACKGROUND_SUBTRACT = false to disable entirely.
//   - APPLY_HARD_THRESHOLD / HARD_THRESHOLD_METHOD: applied AFTER
//     background subtraction, BEFORE the ROI/Coloc2 - runs Fiji's
//     built-in auto-threshold on each of channels 1-3 and clears
//     everything below that cutoff to 0. Set APPLY_HARD_THRESHOLD =
//     false to skip this step.
//   - COSTES_RANDOMISATIONS: number of Costes shuffle iterations
//     Coloc2 runs for its P-value - set to 100 (recommended minimum
//     for a reasonably stable P-value estimate).

CHANNEL_PREFIXES = newArray("C1-", "C3-", "C4-");   // channels 1, 2, 3 only - channel 4 excluded
BACKGROUND_SUBTRACT  = true;
ROLLING_BALL_RADIUS  = 50;
USE_SLIDING_PARABOLOID = true;   // tracks uneven backgrounds more tightly than a fixed-radius ball
APPLY_HARD_THRESHOLD = true;     // clears anything below the auto-threshold to 0, before ROI/Coloc2
HARD_THRESHOLD_METHOD = "Default";   // Fiji auto-threshold method (e.g. "Default", "Otsu", "Triangle")
COSTES_RANDOMISATIONS = 100;
DEBUG_PRINT_FIRST_LOG = true;   // prints the full raw Coloc2 log for the FIRST pair only, then turns itself off

imageDir  = getDirectory("Select image directory");
roiDir    = getDirectory("Select ROI directory");
outputDir = getDirectory("Select results output directory");

outputFile = outputDir + "coloc_results.csv";


// ================================================================
// WRITE CSV HEADER
// ================================================================

f = File.open(outputFile);
print(f, "filename,channel_1,channel_2,M1,M2,pearson,pearson_above_threshold,pearson_below_threshold," +
         "mask_type_used," +
         "ch1_mean_raw,ch1_min_raw,ch1_max_raw,ch2_mean_raw,ch2_min_raw,ch2_max_raw," +
         "ch1_mean,ch2_mean,ch1_max_threshold,ch2_max_threshold," +
         "percent_zero_zero,percent_saturated_ch1,percent_saturated_ch2," +
         "slope,y_intercept,y_intercept_to_mean_ratio," +
         "costes_pvalue,costes_ratio_rand_ge_actual," +
         "warning_y_intercept_far_from_zero,warning_ch1_threshold_too_high,warning_ch2_threshold_too_high");
File.close(f);


// ================================================================
// GET ROI FILES
// ================================================================

roiList = getFileList(roiDir);

processed = 0;
missingImages = 0;


// ================================================================
// LOOP THROUGH ROI FILES
// ================================================================

for (i = 0; i < roiList.length; i++) {

    roiFilename = roiList[i];

    // Only process ROI zip files with the expected naming convention
    if (!endsWith(roiFilename, "_RoiSet.zip")) {
        continue;
    }


    // ------------------------------------------------------------
    // Derive image filename from ROI filename
    // ------------------------------------------------------------

    base = replace(roiFilename, "_RoiSet.zip", "");

    imageFilename = "";

    if (File.exists(imageDir + base)) {
        imageFilename = base;
    }


    // ------------------------------------------------------------
    // ROI exists but matching image does not
    // ------------------------------------------------------------

    if (imageFilename == "") {

        print("WARNING: ROI found but matching image not found:");
        print("  ROI: " + roiFilename);

        missingImages++;

        continue;
    }


    // ------------------------------------------------------------
    // PROCESS THIS IMAGE
    // ------------------------------------------------------------

    print("");
    print("Processing: " + imageFilename);
    print("Using ROI:  " + roiFilename);

    path = imageDir + imageFilename;

    run("Bio-Formats Importer",
        "open=[" + path + "] color_mode=Default view=Hyperstack stack_order=XYCZT");

    title = getTitle();


    // ------------------------------------------------------------
    // SPLIT CHANNELS
    // ------------------------------------------------------------

    run("Split Channels");


    // ------------------------------------------------------------
    // BUILD WINDOW NAMES FOR CHANNELS 1-3 AND CHECK THEY EXIST
    // ------------------------------------------------------------

    channelWindows = newArray(CHANNEL_PREFIXES.length);
    allChannelsFound = true;

    for (c = 0; c < CHANNEL_PREFIXES.length; c++) {
        channelWindows[c] = CHANNEL_PREFIXES[c] + title;

        if (!isOpen(channelWindows[c])) {
            print("WARNING: channel not found:");
            print("  Expected: " + channelWindows[c]);
            allChannelsFound = false;
        }
    }

    if (!allChannelsFound) {
        print("Skipping " + imageFilename);
        run("Close All");
        continue;
    }


    // ------------------------------------------------------------
    // DEFAULT DISPLAY RANGE (color balance) FOR EACH CHANNEL
    //
    // resetMinAndMax() sets the display range to Fiji's normal
    // default for the image's bit depth (equivalent to clicking
    // "Reset" in the Brightness/Contrast dialog) - this replaces
    // any prior explicit/forced display range.
    // ------------------------------------------------------------

    for (c = 0; c < channelWindows.length; c++) {
        selectWindow(channelWindows[c]);
        resetMinAndMax();
    }


    // ------------------------------------------------------------
    // MEASURE RAW (PRE-SUBTRACTION) CHANNEL STATS
    //
    // Captured before any processing, over the WHOLE image (not yet
    // ROI-restricted) - ground truth for whether an image was
    // genuinely acquired dim versus becoming dim only after
    // background subtraction.
    // ------------------------------------------------------------

    ch_mean_raw = newArray(channelWindows.length);
    ch_min_raw  = newArray(channelWindows.length);
    ch_max_raw  = newArray(channelWindows.length);

    for (c = 0; c < channelWindows.length; c++) {
        selectWindow(channelWindows[c]);
        getStatistics(area, mean_c, min_c, max_c);
        ch_mean_raw[c] = mean_c;
        ch_min_raw[c]  = min_c;
        ch_max_raw[c]  = max_c;
    }


    // ------------------------------------------------------------
    // BACKGROUND SUBTRACTION (applied to whole image, before ROI,
    // once per channel)
    // ------------------------------------------------------------

    if (BACKGROUND_SUBTRACT) {
        bgOptions = "rolling=" + ROLLING_BALL_RADIUS;
        if (USE_SLIDING_PARABOLOID) {
            bgOptions = bgOptions + " sliding";
        }

        for (c = 0; c < channelWindows.length; c++) {
            selectWindow(channelWindows[c]);
            run("Subtract Background...", bgOptions);
        }
    }

    if (APPLY_HARD_THRESHOLD) {
        for (c = 0; c < channelWindows.length; c++) {
            applyHardThreshold(channelWindows[c], HARD_THRESHOLD_METHOD);
        }
    }


    // ------------------------------------------------------------
    // LOAD ROI
    // ------------------------------------------------------------

    roiPath = roiDir + roiFilename;

    roiManager("reset");
    roiManager("Open", roiPath);

    roiCount = roiManager("count");

    if (roiCount == 0) {

        print("WARNING: ROI file contains no ROIs:");
        print("  " + roiFilename);

        run("Close All");

        continue;
    }


    // ------------------------------------------------------------
    // APPLY ROI TO EACH CHANNEL
    // ------------------------------------------------------------

    for (c = 0; c < channelWindows.length; c++) {
        selectWindow(channelWindows[c]);
        roiManager("Select", 0);
        if (selectionType() == -1) {
            print("WARNING: ROI did not apply to channel " + channelWindows[c] + " for " + imageFilename);
        }
    }


    // ------------------------------------------------------------
    // RUN COLOC2 ON EVERY PAIRWISE COMBINATION OF CHANNELS 1-3
    // ------------------------------------------------------------

    for (a = 0; a < channelWindows.length; a++) {
        for (b = a + 1; b < channelWindows.length; b++) {

            chan1 = channelWindows[a];
            chan2 = channelWindows[b];

            print("");
            print("  Coloc2: " + chan1 + " vs " + chan2);

            // CLEAR LOG BEFORE COLOC2
            print("\\Clear");

            run("Coloc 2",
			    "channel_1=[" + chan1 + "] channel_2=[" + chan2 + "] " +
			    "roi_or_mask=[ROI Manager] " +
			    "threshold_regression=Costes " +
			    "display_images_in_result " +
			    "manders'_correlation " +
			    "costes'_significance_test " +
			    "psf=3 costes_randomisations=" + COSTES_RANDOMISATIONS);

            // ------------------------------------------------------------
            // READ COLOC2 LOG
            // ------------------------------------------------------------

            logText = getInfo("log");

            if (DEBUG_PRINT_FIRST_LOG) {
                print("========== RAW COLOC2 LOG (debug) ==========");
                print(logText);
                print("========== END RAW LOG ==========");
                DEBUG_PRINT_FIRST_LOG = false;  // only do this once
            }

            m1      = parseColocValue(logText, "Manders' tM1");
            m2      = parseColocValue(logText, "Manders' tM2");
            mask_type_used        = parseColocValue(logText, "Mask Type Used");
            pearson               = parseColocValue(logText, "Pearson's R value (no threshold)");
            pearson_above         = parseColocValue(logText, "Pearson's R value (above threshold)");
            pearson_below         = parseColocValue(logText, "Pearson's R value (below threshold)");

            ch1_mean              = parseColocValue(logText, "Channel 1 Mean");
            ch2_mean              = parseColocValue(logText, "Channel 2 Mean");
            ch1_max_threshold     = parseColocValue(logText, "Ch1 Max Threshold");
            ch2_max_threshold     = parseColocValue(logText, "Ch2 Max Threshold");

            percent_zero_zero     = parseColocValue(logText, "% zero-zero pixels");
            percent_saturated_ch1 = parseColocValue(logText, "% saturated ch1 pixels");
            percent_saturated_ch2 = parseColocValue(logText, "% saturated ch2 pixels");

            slope                 = parseColocValue(logText, "m (slope)");
            y_intercept           = parseColocValue(logText, "b (y-intercept)");
            y_intercept_ratio     = parseColocValue(logText, "b to y-mean ratio");

            costes_pvalue         = parseColocValue(logText, "Costes P-Value");
            costes_ratio          = parseColocValue(logText, "Ratio of rand. Pearsons >= actual Pearsons value");

            warning_intercept = "FALSE";
            if (indexOf(logText, "y-intercept far from zero") != -1) {
                warning_intercept = "TRUE";
            }

            warning_ch1_thresh = "FALSE";
            if (indexOf(logText, "Threshold of ch. 1 too high") != -1) {
                warning_ch1_thresh = "TRUE";
            }

            warning_ch2_thresh = "FALSE";
            if (indexOf(logText, "Threshold of ch. 2 too high") != -1) {
                warning_ch2_thresh = "TRUE";
            }

            if (mask_type_used != "ROI") {
                print("WARNING: Mask Type Used = '" + mask_type_used + "' (expected 'ROI') for " + imageFilename + " (" + chan1 + " vs " + chan2 + ")");
            }

            File.append(
                imageFilename + "," +
                chan1 + "," +
                chan2 + "," +
                m1 + "," +
                m2 + "," +
                pearson + "," +
                pearson_above + "," +
                pearson_below + "," +
                mask_type_used + "," +
                ch_mean_raw[a] + "," +
                ch_min_raw[a] + "," +
                ch_max_raw[a] + "," +
                ch_mean_raw[b] + "," +
                ch_min_raw[b] + "," +
                ch_max_raw[b] + "," +
                ch1_mean + "," +
                ch2_mean + "," +
                ch1_max_threshold + "," +
                ch2_max_threshold + "," +
                percent_zero_zero + "," +
                percent_saturated_ch1 + "," +
                percent_saturated_ch2 + "," +
                slope + "," +
                y_intercept + "," +
                y_intercept_ratio + "," +
                costes_pvalue + "," +
                costes_ratio + "," +
                warning_intercept + "," +
                warning_ch1_thresh + "," +
                warning_ch2_thresh,
                outputFile
            );

            processed++;

            // close Coloc2's own results/warnings frame between pairs
            // so it doesn't accumulate or get confused with the next
            // pair's run
            closeFramesContaining("coloc");
            closeFramesContaining("warning");
        }
    }


    // ------------------------------------------------------------
    // CLEAN UP
    // ------------------------------------------------------------

    if (nImages() > 0) {
        run("Close All");
    }

    winTitles = getList("window.titles");
    for (w = 0; w < winTitles.length; w++) {
        wTitle = winTitles[w];
        if (wTitle == "ROI Manager") {
            continue;
        }
        if (isOpen(wTitle)) {
            close(wTitle);
        }
    }

    closeFramesContaining("coloc");
    closeFramesContaining("warning");
}


// ================================================================
// FINISHED
// ================================================================

print("");
print("========================================");
print("Done.");
print("========================================");
print("Image/channel-pair rows processed: " + processed);
print("ROIs with no matching image: " + missingImages);
print("Results written to:");
print(outputFile);
print("========================================");


function applyHardThreshold(windowTitle, method) {
    selectWindow(windowTitle);

    setAutoThreshold(method + " dark");
    getThreshold(lower, upper);

    setThreshold(lower, upper);
    run("Create Selection");
    run("Clear Outside");
    resetThreshold();
    run("Select None");

    print("  Applied hard threshold (" + method + ") to '" + windowTitle + "': cutoff = " + lower);
}

function closeFramesContaining(substr) {
    substrLower = toLowerCase(substr);
    code =
        "importClass(Packages.java.awt.Frame);\n" +
        "var frames = Frame.getFrames();\n" +
        "for (var i = 0; i < frames.length; i++) {\n" +
        "    var t = frames[i].getTitle();\n" +
        "    if (t != null && t.toLowerCase().indexOf('" + substrLower + "') >= 0) {\n" +
        "        frames[i].dispose();\n" +
        "    }\n" +
        "}\n";
    eval("script", code);
}

function parseColocValue(text, label) {

    idx = indexOf(text, label);

    if (idx == -1)
        return "NA";

    sub = substring(text, idx);

    lineEnd = indexOf(sub, "\n");

    if (lineEnd == -1)
        lineEnd = lengthOf(sub);

    line = substring(sub, 0, lineEnd);

    commaIdx = lastIndexOf(line, ",");

    if (commaIdx == -1)
        return "NA";

    val = String.trim(substring(line, commaIdx + 1));

    return val;
}