#include <cairo/cairo.h>
#include <iostream>
#include <cmath>
#include <string>

#include "Karyogram.h"
#include "SDDutilities.h"


// Bookkeeping struct - remembers each drawn chromosome's pixel
// position and height so the DSB/SSB-drawing code further down can
// look it back up by chromosome/chromatid number.
struct ChromosomeGeometry
{
    int chromosomeNumber;						// Chromosome Number (i.e. ID).
    int chromatidNumber;						// Chromatid Number (1 for unduplicated, 1 or 2 for duplicated)
    double x;								// X pixel position
    double y;								// Y pixel position
    double height;							// Chromosome height in pixels
};



bool Karyogram::generateKaryogram(                                         // Function that generates the overall Karyogram structure
    const std::vector<double>& chromosomeSizes,                            // Requires chromosome sizes in the SDD file header for scaling on the image
    const std::vector<double>& cellCyclePhase,                             // Requires cell cycle phase in the SDD file header for determining presence of duplicated chromosomes.
    const std::vector<double>& doseOrFluence,                              // Accesses the SDD dose or fluence listed in the header
    const std::vector<int>& incidentParticles,                             // Accesses the SDD incident particles in the header
    const std::vector<Exposure>& exposures,                                // Accesses the damage locations associated with each chromosome in a given exposure.
    bool humanGenome,                                                      // Boolean to determine if karyogram should be drawn with human centromere ranges or generic centromere locations for non-human genomes.
    const std::string& outputFilename)                                     // Outputs the image to the outputFilename.
{

    if (chromosomeSizes.empty())                                           // Make sure SDD header 'Chromosome sizes' Is non-empty.
    {
        std::cerr << "ERROR: No chromosome sizes provided.\n";
        return false;
    }

    std::vector<DamageLocation> dsbVec = getDoubleStrandBreaks(exposures); // Use this vector for plotting DSBs onto karyogram

    std::vector<DamageLocation> ssbVec = getSingleStrandBreaks(exposures);

    std::vector<ChromosomeGeometry> chromosomeGeometry;                    // Store the geometry of the chromosomes in x and y pixel coordinates to draw the DSBs.

    const int chromosomeCount = static_cast<int>(chromosomeSizes[0]);      // The first entry specifies the number of chromosomes.

    //const int homologousSetSize = (chromosomeCount - 2) / 2;             // Store number of homologs for DSB drawing

    // Generic metacentric chromosomes similar to MEDRAS-MC
    const double genericCentromereStart = 0.49;                            // In case non-human chromosomes specified, draw generic centromere regions
    const double genericCentromereEnd = 0.51;                              // Will be replaced later if human genome is specified with actual centromere positions.

    if (chromosomeSizes.size() !=                                          // Check that the vector contains the expected number of chromosome sizes.
        static_cast<size_t>(chromosomeCount + 1))
    {
        std::cerr << "ERROR: Chromosome count does not match "
                  << "the number of chromosome sizes provided.\n";
        return false;
    }

    bool doubleChromatid = false;                                          // Check cell cycle phase to determine number of chromatids per chromosomes to draw.

    if (!cellCyclePhase.empty())                                           // Cell cycle phase of 3, 4, or 5 indicates post-replicated DNA, two chromatids per chromosome.
    {
        const int phase = static_cast<int>(cellCyclePhase[0]);

        if (phase == 3 || phase == 4 || phase == 5)
        {
            doubleChromatid = true;
        }
    }

    // ------------------------------------------------------------
    // Determine chromosome layout
    // ------------------------------------------------------------

    const ChromosomeLayout chromosomeLayout = determineChromosomeLayout(chromosomeSizes); // Determine which layout the chromosome sizes are listed in for plotting logic.

    const bool hasHomologs = chromosomeLayout != ChromosomeLayout::NON_HOMOLOGOUS;  // Default layout is non-homologous 1,...,N,Y,X

    const int homologousPairs = hasHomologs ? (chromosomeCount - 2) / 2 : 0;   // Number of homologous pairs if homologs present.

    const int drawableGroups = hasHomologs ? homologousPairs : chromosomeCount - 2; // Number of chromosome groups that will be drawn before the Y and X

    const int totalColorGroups = drawableGroups + 2;                       // Use totalColorGroups to have evenly distrbuted colors across all chromosomes.

    double maxBP = 0.0;                                                    // Find largest chromosome, scale all chromosomes with respect to the largest

    for (size_t i = 1 ; i < chromosomeSizes.size(); i++)
    {
        if (chromosomeSizes[i] > maxBP)
        {
            maxBP = chromosomeSizes[i];
        }
    }

    if (maxBP <= 0.0)
    {
        std::cerr << "ERROR: Chromosome sizes must contain positive values.\n";
        return false;
    }

    // ----------------------------------------------------
    // Karyogram dimensions
    // ----------------------------------------------------

    const int imgWidth = 1000;                                             // Image width in number of pixels
    const int columns = 4;                                                 // Desired number of columns in the karyogram
    const double colWidth = static_cast<double>(imgWidth) / columns;       // Column width also in number of pixels
    const double rowHeight = 250.0;                                        // Row height in number of pixels
    const double startY = 188.0;                                           // Choosing the starting Y position for the first row of chromosomes
    const int rows = (drawableGroups + columns - 1) / columns;             // Determine number of rows based on number of chromosomes passed.
    const double legendHeight = 70.0;
    const double legendBottomMargin = 25.0;
    const int imgHeight = static_cast<int>(startY + rows * rowHeight + legendHeight + legendBottomMargin);          // Adjust image height based on the number of rows, also in pixels.
    const double maxRenderHeight = 180.0;                                  // Limit the maximum height of the chromosome, also in pixels
    const double chromosomeWidth = 20.0;                                   // Individual chromosomes are 20 pixels wide.
    const double chromatidGap = 1.5;                                       // The gaps between homologs in a replicated chromosome are 1.5 pixels.


    // ----------------------------------------------------
    // Generate the background Karyogram template
    // ----------------------------------------------------
    cairo_surface_t* surface =
        cairo_image_surface_create(
      	    CAIRO_FORMAT_ARGB32,
            imgWidth,
            imgHeight
        );

    cairo_t* cr = cairo_create(surface);

    cairo_set_antialias(cr, CAIRO_ANTIALIAS_BEST);

    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    cairo_paint(cr);

    // ----------------------------------------------------
    // Draw SDD summary
    // ----------------------------------------------------

    drawSDDsummary(cr, cellCyclePhase, exposures, doseOrFluence, incidentParticles);


    // -----------------------------------------------------------
    // Draw chromosomes depending on chromosome size layout given
    // -----------------------------------------------------------

    for (int i = 0; i < drawableGroups; ++i)
    {
        size_t firstIndex = 0;                                             // Declare index of first homolog to be determined depending on chromosome sizes layout.
        size_t secondIndex = 0;                                            // Second homolog begins after the first

        // -----------------------------------------------------
        // Determine indexing based off chromosome sizes layout
        // -----------------------------------------------------

        if (chromosomeLayout == ChromosomeLayout::SPLIT_HOMOLOGS)          // Chromosomes listed [count,1,2,3,..., 1,2,3,...,Y,X]
        {
            firstIndex = static_cast<size_t>(i + 1);
            secondIndex = static_cast<size_t>(i + 1 + homologousPairs);
        }
        else if (chromosomeLayout == ChromosomeLayout::ADJACENT_HOMOLOGS)  // Chromosomes listed [count,1,1,2,2,...,Y,X]
        {
            firstIndex = static_cast<size_t>(2 * i + 1);
            secondIndex = static_cast<size_t>(2 * i + 2);
        }
        else                                                               // Chromosomes listed [count,1,2,3,4...,N-2,Y,X]
        {
            firstIndex = static_cast<size_t>(i + 1);
        }


        // Determine position of chromosome group
        int col = i % columns;
        int row = i / columns;

        double groupCenterX = (col * colWidth) + (colWidth / 2.0);
        double posY = startY + (row * rowHeight);

        RGB chromosomeColor = generateChromosomeColor(i, totalColorGroups);// Homologs receive same color, non-homologous chromosomes receive their own color.


        // -----------------------------------------------------------------
        // NON-HOMOLOGOUS CHROMOSOME DRAWING
        // -----------------------------------------------------------------

        if (chromosomeLayout == ChromosomeLayout::NON_HOMOLOGOUS)
        {
            const double chromosomeSize = chromosomeSizes[firstIndex];     // Only care about firstIndex because there are no homologs.

            double chromosomeHeight = (chromosomeSize / maxBP) * maxRenderHeight; // Normalize chromosome heights relative to the largest chromosoem maxBP

            if (chromosomeHeight < 30.0)                                   // Prevent too small chromosome sizes.
            {
                chromosomeHeight = 30.0;
            }

            // If non-human genome passed, use generic centromere start and end locations
            double centromereStart = genericCentromereStart;
            double centromereEnd = genericCentromereEnd;

            // Adjust centromere locations to the ones listed at the start if a humanGenome specified.
            if (humanGenome)
            {
                int humanChromosomeID = static_cast<int>(firstIndex);

                if (chromosomeCount == 24)
                {
                    if (firstIndex == 23)                                  // When chromosome sizes index 23 is reached, associate that with the Y chromosome
                    {
                        humanChromosomeID = 45;   // Y
                    }
                    else if (firstIndex == 24)                             // When chromosome sizes index 24 is reached, associate that with the X chromosome.
                    {
                        humanChromosomeID = 46;   // X
                    }
                }


                const CentromerePosition* centromere = getHumanCentromere(humanChromosomeID);       // Get corresponding centromere positions per chromosome ID

                if (centromere != nullptr)
                {
                    // Set corresponding centromere ranges and convert to fraction of the chromosome length for conversion to pixels later.
                    centromereStart = static_cast<double>(centromere->start) / (chromosomeSize * 1000000.0);        // chromosome sizes listed in Mbp, convert to bp.
                    centromereEnd = static_cast<double>(centromere->end)/ (chromosomeSize*1000000.0);
                }
            }


            // ------------------------------------------------
            // Single-chromatid drawing
            // ------------------------------------------------
            if (!doubleChromatid)
            {
                const double posX = groupCenterX - chromosomeWidth / 2.0;  // Set the horizontal position of the single chromatid drawing

                // Draw single-chromatid chromosome and store the ID, chromatid number, position and size in chromosomeGeometry.
                drawChromosome(cr, posX, posY, chromosomeHeight, chromosomeWidth, chromosomeColor, centromereStart, centromereEnd);
                chromosomeGeometry.push_back({static_cast<int>(firstIndex), 1, posX, posY, chromosomeHeight});

            }

            // ------------------------------------------------
            // Two chromatids drawing
            // ------------------------------------------------
            else
            {
                // Determine horizontal X positions in the karyogram of the two chromatids for the corresponding chromosome.
                const double posX = groupCenterX - chromosomeWidth - chromatidGap / 2.0;
                const double posX2 = groupCenterX + chromosomeWidth + chromatidGap;

                // Draw both chromatids of the chromosome in their respective X positions.
                drawChromosome(cr, posX, posY, chromosomeHeight, chromosomeWidth, chromosomeColor, centromereStart, centromereEnd);
                drawChromosome(cr, posX2, posY, chromosomeHeight, chromosomeWidth, chromosomeColor, centromereStart, centromereEnd);

                // Both chromatids in the non-homologous case are assigned the same chromosome ID, it is just in the post-replicated cell phase, so there are two.
                chromosomeGeometry.push_back({static_cast<int>(firstIndex), 1, posX, posY, chromosomeHeight});
                chromosomeGeometry.push_back({static_cast<int>(firstIndex), 2, posX2, posY, chromosomeHeight});
            }

            // ------------------------------------------------
            // Label chromosome
            // ------------------------------------------------
            cairo_set_source_rgb(
                cr,
                0.1,
                0.1,
                0.1
            );

            cairo_select_font_face(
                cr,
                "Sans",
                CAIRO_FONT_SLANT_NORMAL,
                CAIRO_FONT_WEIGHT_BOLD
            );

            cairo_set_font_size(cr, 16.0);

            std::string label =
                "Chr " + std::to_string(i + 1);

            cairo_move_to(
                cr,
                groupCenterX - 20.0,
                posY + chromosomeHeight + 25.0
            );

            cairo_show_text(
                cr,
                label.c_str()
            );
        }

        // -----------------------------------------------------
        // HOMOLOGOUS CHROMOSOMES DRAWING
        // -----------------------------------------------------

        else
        {
            // Initialize centromere positions per chromosome ID for each homolog
            double firstCentromereStart = genericCentromereStart;

            double firstCentromereEnd = genericCentromereEnd;

            double secondCentromereStart = genericCentromereStart;

            double secondCentromereEnd = genericCentromereEnd;

            // If human centromeres, get the positions corresponding to the correct ID.
            if (humanGenome)
            {
                const CentromerePosition* centromere = getHumanCentromere(i + 1);

                if(centromere != nullptr)
                {
                // Obtain centromre start and end positions for both homologs from the CentromerePosition struct.
                    firstCentromereStart = static_cast<double>(centromere->start) / (chromosomeSizes[firstIndex] * 1000000.0);
                    firstCentromereEnd = static_cast<double>(centromere->end) / (chromosomeSizes[firstIndex] * 1000000.0);
                    secondCentromereStart = static_cast<double>(centromere->start) / (chromosomeSizes[secondIndex] * 1000000.0);
                    secondCentromereEnd = static_cast<double>(centromere->end) / (chromosomeSizes[secondIndex] * 1000000.0);
                }

            }


            // ----------------------------------------------
            // Scale chromosome heights
            // ----------------------------------------------

            double firstHeight = (chromosomeSizes[firstIndex] / maxBP) * maxRenderHeight;   // Normalize height of first homolog to the largest chromosome height

            double secondHeight = (chromosomeSizes[secondIndex] / maxBP) * maxRenderHeight; // Normalize height of second homolog to the largest chromosome height.

	    // Clamp minimum chromosome height to 30.0 pixels
            if (firstHeight < 30.0)
            {
                firstHeight = 30.0;
            }

            if (secondHeight < 30.0)
            {
                secondHeight = 30.0;
            }


            // ----------------------------------------------
            // Single-chromatid homologs
            // ----------------------------------------------

            if (!doubleChromatid)
            {
                const double homologGap = 18.0;                            			// Define spacing between homologs in pixels.

                const double leftChromPosX = groupCenterX - chromosomeWidth - homologGap / 2.0;	// Determine the X positions of the left homolog

                const double rightChromPosX = groupCenterX + homologGap / 2.0;              	// Determine the X positions of the right homolog

                // Draw first homolog at the desired location with the desired dimensions
                drawChromosome(cr, leftChromPosX, posY, firstHeight, chromosomeWidth, chromosomeColor, firstCentromereStart, firstCentromereEnd);

                // Draw second homolog at the desired location with the desired dimensions
                drawChromosome(cr, rightChromPosX, posY, secondHeight, chromosomeWidth, chromosomeColor, secondCentromereStart, secondCentromereEnd);

                // Store actual SDD chromosome IDs and their corresponding dimensions and locations --> for damage marker drawings later on.
                chromosomeGeometry.push_back({static_cast<int>(firstIndex), 1, leftChromPosX, posY, firstHeight});
                chromosomeGeometry.push_back({static_cast<int>(secondIndex), 1, rightChromPosX, posY, secondHeight});
            }


            // ----------------------------------------------
            // Double-chromatid homologs
            // ----------------------------------------------

            else
            {
                const double homologGap = 40.0;                            			// Increase homolog gap to account for second chromatids

                // Determine homolog 1 left and right chromatid positions
                const double leftChromatid1PosX = groupCenterX - chromosomeWidth - homologGap / 2.0 - chromatidGap / 2.0;
                const double leftChromatid2PosX = leftChromatid1PosX + chromosomeWidth + chromatidGap;

                // Determine homolog 2 left and right chromatid positions
                const double rightChromatid1PosX = groupCenterX + homologGap / 2.0;
                const double rightChromatid2PosX = rightChromatid1PosX + chromosomeWidth + chromatidGap;


                // Draw homolog 1 chromatids at the desired locations with the desired dimensions
                drawChromosome(cr, leftChromatid1PosX, posY, firstHeight, chromosomeWidth, chromosomeColor, firstCentromereStart, firstCentromereEnd);
                drawChromosome(cr, leftChromatid2PosX, posY, firstHeight, chromosomeWidth, chromosomeColor, firstCentromereStart, firstCentromereEnd);


                // Draw homolog 2 chromatids at the desired locations with the desired dimensions
                drawChromosome(cr, rightChromatid1PosX, posY, secondHeight, chromosomeWidth, chromosomeColor, secondCentromereStart, secondCentromereEnd);
                drawChromosome(cr, rightChromatid2PosX, posY, secondHeight, chromosomeWidth, chromosomeColor, secondCentromereStart, secondCentromereEnd);

                // Store actual SDD chromosome IDs with their corresponding dimensions and locations --> for damage marker drawings later on
                chromosomeGeometry.push_back({static_cast<int>(firstIndex), 1, leftChromatid1PosX, posY, firstHeight});
                chromosomeGeometry.push_back({static_cast<int>(firstIndex), 2, leftChromatid2PosX, posY, firstHeight});
                chromosomeGeometry.push_back({static_cast<int>(secondIndex), 1, rightChromatid1PosX, posY, secondHeight});
                chromosomeGeometry.push_back({static_cast<int>(secondIndex), 2, rightChromatid2PosX, posY, secondHeight});
            }


            // ----------------------------------------------
            // Label homologous pair
            // ----------------------------------------------

            const double labelHeight = std::max(firstHeight, secondHeight);

            cairo_set_source_rgb(
                cr,
                0.1,
                0.1,
                0.1
            );

            cairo_select_font_face(
                cr,
                "Sans",
                CAIRO_FONT_SLANT_NORMAL,
                CAIRO_FONT_WEIGHT_BOLD
            );

            cairo_set_font_size(cr, 16.0);

            std::string label =
                "Chr " + std::to_string(i + 1);

            cairo_move_to(
                cr,
                groupCenterX - 20.0,
                posY + labelHeight + 25.0
            );

            cairo_show_text(
                cr,
                label.c_str()
            );
        }
    }


    // --------------------------------------------------
    // X and Y chromosome drawing
    // --------------------------------------------------

    // Always assume the user gives the Y chromosome then the X chromosome as the last two listed chromosome sizes.
    const size_t xIndex = static_cast<size_t>(chromosomeCount);            // Determining the index of the X chromosome
    const size_t yIndex = static_cast<size_t>(chromosomeCount - 1);        // Determining the index of the Y chromosome

    // Set the centromere locations to the generic positions, replaced if user passed '--karyogram human' option.
    double xCentromereStart = genericCentromereStart;
    double xCentromereEnd = genericCentromereEnd;
    double yCentromereStart = genericCentromereStart;
    double yCentromereEnd = genericCentromereEnd;

    if (humanGenome)                                                       // If user specified a human genome '--karyogram human'
    {
        int xHumanID = 46;
        int yHumanID = 45;

        // Get the base pair locations of the centromeres for the X and Y chromosomes
        const CentromerePosition* xCentromere = getHumanCentromere(xHumanID);
        const CentromerePosition* yCentromere = getHumanCentromere(yHumanID);

        if (xCentromere != nullptr)
        {
            // Scale centromere locations from Mbp to bp
            xCentromereStart = static_cast<double>(xCentromere->start) / (chromosomeSizes[xIndex] * 1000000.0);
            xCentromereEnd = static_cast<double>(xCentromere->end) / (chromosomeSizes[xIndex] * 1000000.0);
        }

        if (yCentromere != nullptr)
        {
            // Scale centromere locations from Mbp to bp
            yCentromereStart = static_cast<double>(yCentromere->start) / (chromosomeSizes[yIndex] * 1000000.0);
            yCentromereEnd = static_cast<double>(yCentromere->end) / (chromosomeSizes[yIndex] * 1000000.0);
        }
    }

    // Scale the X and Y chromosome heights relative to the largest chromosome.
    double xHeight = (chromosomeSizes[xIndex] / maxBP) * maxRenderHeight;
    double yHeight = (chromosomeSizes[yIndex] / maxBP) * maxRenderHeight;

    // Prevent X and Y chromosomes from being smaller than 30 pixels
    if (xHeight < 30.0)
    {
        xHeight = 30.0;
    }
    if (yHeight < 30.0)
    {
        yHeight = 30.0;
    }

    // -------------------------------------------------------------
    // X and Y chromosome pixel dimensions
    // -------------------------------------------------------------
    // X and Y occupy the remaining positions in the final row.
    const int sexChromosomeRow = (drawableGroups > 0) ? (drawableGroups - 1) / columns : 0;
    // Plot X first then Y next.
    const int xColumn = 2;
    const int yColumn = 3;
    const double xGroupCenterX = (xColumn * colWidth) + (colWidth / 2.0);
    const double yGroupCenterX = (yColumn * colWidth) + (colWidth / 2.0);
    const double sexChromPosY = startY + (sexChromosomeRow * rowHeight);
    const double xChromPosX = xGroupCenterX - chromosomeWidth / 2.0;
    const double yChromPosX = yGroupCenterX - chromosomeWidth / 2.0;


    // Set X and Y chromosome Colors
    RGB yColor = generateChromosomeColor(drawableGroups, totalColorGroups);
    RGB xColor = generateChromosomeColor(drawableGroups + 1, totalColorGroups);

    if (!doubleChromatid)                                                  // Depending on if SDD header 'Cell cycle phase' specified a post-replicated chromosome.
    {
        // Draw X
        drawChromosome(cr, xChromPosX, sexChromPosY, xHeight, chromosomeWidth, xColor, xCentromereStart, xCentromereEnd);
        // Store geometry of X chromosome for damage drawing
        chromosomeGeometry.push_back({chromosomeCount, 1, xChromPosX, sexChromPosY, xHeight});

        // Draw Y
        drawChromosome(cr, yChromPosX, sexChromPosY, yHeight, chromosomeWidth, yColor, yCentromereStart, yCentromereEnd);
        // Store geometry of Y chromosome for damage drawing
        chromosomeGeometry.push_back({chromosomeCount - 1, 1, yChromPosX, sexChromPosY, yHeight});

    }
    else
    {
        // Draw X chromatids 1 and 2
        drawChromosome(cr, xChromPosX, sexChromPosY, xHeight, chromosomeWidth, xColor, xCentromereStart, xCentromereEnd);
        drawChromosome(cr, xChromPosX + chromosomeWidth + chromatidGap, sexChromPosY, xHeight, chromosomeWidth, xColor, xCentromereStart, xCentromereEnd);
        // Store X chromatids 1 and 2 geometries for damage drawing
        chromosomeGeometry.push_back({chromosomeCount, 1, xChromPosX, sexChromPosY, xHeight});
        chromosomeGeometry.push_back({chromosomeCount, 2, xChromPosX + chromosomeWidth + chromatidGap, sexChromPosY, xHeight});

        // Draw Y chromatids 1 and 2
        drawChromosome(cr, yChromPosX, sexChromPosY, yHeight, chromosomeWidth, yColor, yCentromereStart, yCentromereEnd);
        drawChromosome(cr, yChromPosX + chromosomeWidth + chromatidGap, sexChromPosY, yHeight, chromosomeWidth, yColor, yCentromereStart, yCentromereEnd);

        // Store Y chromatids 1 and 2 geometries for damage drawing
        chromosomeGeometry.push_back({chromosomeCount - 1, 1, yChromPosX, sexChromPosY, yHeight});
        chromosomeGeometry.push_back({chromosomeCount - 1, 2, yChromPosX + chromosomeWidth + chromatidGap, sexChromPosY, yHeight});

    }


    // --------------------------------------------------
    // Draw double-strand break markers
    // --------------------------------------------------

    for (const auto& dsb : dsbVec)                                      // Loop through the double-strand break vector containing exposure data (damages per chromosome per exposure)
    {
        const int chromosomeID = dsb.chromosomeID.chromosomeNumber;     // Chromosome ID 1-46 for humans

        const int chromatidID = dsb.chromosomeID.chromatidNumber;       // Can be either 1 for unduplicated chromosomes, or 1 or 2 for duplicated chromosomes.

        for (const auto& geometry : chromosomeGeometry)                 // Get the associated pixel coordinates for the given chromosome ID and chromatid to draw the double-strand breaks.
        {
            if (geometry.chromosomeNumber != chromosomeID)
            {
                continue;
            }

            if (geometry.chromatidNumber != chromatidID)
            {
                continue;
            }


            // ------------------------------------------
            // Determine chromosome size for the given ID
            // ------------------------------------------
            const size_t sizeIndex = static_cast<size_t>(chromosomeID);

            if (sizeIndex >= chromosomeSizes.size())
            {
                continue;
            }

            const double chromosomeSize = chromosomeSizes[sizeIndex];


            // ------------------------------------------
            // Convert damage position to fraction of chromosome length
            // ------------------------------------------
            double damageFraction = getDamageFraction(dsb, chromosomeSize);

            // If fraction appears outside of the chromosome, clamp the damage to appear at the end of the chromosome.
            damageFraction = std::clamp(damageFraction, 0.0, 1.0);


            // ------------------------------------------
            // Convert fraction to image pixel coordinates
            // ------------------------------------------
            const double markerX = geometry.x + chromosomeWidth / 2.0;
            const double markerY = geometry.y + damageFraction * geometry.height;


            // ------------------------------------------
            // Draw DSB marker
            // ------------------------------------------
            drawDoubleStrandBreakMarker(cr, markerX, markerY, chromosomeWidth);				// Takes damage marker x location, y location, and width of the chromosome to scale with width of the line.

            break;
        }
    }

    // --------------------------------------------------
    // Draw single-strand break markers
    // --------------------------------------------------

    for (const auto& ssb : ssbVec)
    {
        const int chromosomeID = ssb.chromosomeID.chromosomeNumber;

        const int chromatidID = ssb.chromosomeID.chromatidNumber;

        // Find the chromosome/chromatid geometry
        for (const auto& geometry : chromosomeGeometry)
        {
            if (geometry.chromosomeNumber != chromosomeID)
            {
                continue;
            }

            if (geometry.chromatidNumber != chromatidID)
            {
                continue;
            }


            // ------------------------------------------
            // Determine chromosome size
            // ------------------------------------------
            const size_t sizeIndex = static_cast<size_t>(chromosomeID);

            if (sizeIndex >= chromosomeSizes.size())
            {
                continue;
            }

            const double chromosomeSize = chromosomeSizes[sizeIndex];


            // ------------------------------------------
            // Convert SSB position to fraction
            // ------------------------------------------
            double damageFraction = getDamageFraction(ssb, chromosomeSize);

            damageFraction = std::clamp(damageFraction, 0.0, 1.0);				// damageFraction is the location of the damage along the length of the chromosome, clamped to 1 for safety.


            // ------------------------------------------
            // Convert fraction to image coordinates
            // ------------------------------------------
            const double markerX = geometry.x + chromosomeWidth / 2.0;
            const double markerY = geometry.y + damageFraction * geometry.height;


            // ------------------------------------------
            // Determine SSB marker length
            // ------------------------------------------
            const double ssbBaseLength = 4.0;                              // Scaling single-strand break line length based on number of damages in the damage site
            const double ssbLengthPerBreak = 4.0;                          // Marker length per individual single-strand break to be summed over the whole site
            const double maxSSBLength = chromosomeWidth;                   // Maximum single-strand break marker length is the width of the chromosome

            double markerLength = ssbBaseLength + ssbLengthPerBreak * (ssb.numSingleStrandBreaks - 1);
            markerLength = std::min(markerLength, maxSSBLength);		// Scale damage marker length depending on the number of SSBs in the damage site.


            // ------------------------------------------
            // Draw SSB marker
            // ------------------------------------------
            drawSingleStrandBreakMarker(cr, markerX, markerY, markerLength);			// Take damage marker X location, Y location, and scaled markerLength depending on the number of SSBs in the damage site.

            break;
        }
    }


    // X label
    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);

    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);

    cairo_set_font_size(cr, 16.0);

    cairo_move_to(cr, xGroupCenterX - 4.0, sexChromPosY + xHeight + 25.0);

    cairo_show_text(cr, "X");

    // Y label
    cairo_move_to(cr, yGroupCenterX - 4.0, sexChromPosY + yHeight + 25.0);

    cairo_show_text(cr, "Y");

    // ----------------------------------------------------
    // Draw Karyogram legend at bottom
    // ----------------------------------------------------
    const double legendY = static_cast<double>(imgHeight) - 100.0;
    drawLegend(cr, legendY);

    // Write the karyogram to the desired output filename "'SDDfileName'_karyogram.png"
    cairo_status_t status = cairo_surface_write_to_png(surface, outputFilename.c_str());
    cairo_destroy(cr);
    cairo_surface_destroy(surface);

    // Check if Karyogram was generated successfully
    if (status != CAIRO_STATUS_SUCCESS)
    {
        std::cerr << "ERROR: Failed to write karyogram image.\n";
        return false;
    }

    return true;
}








void Karyogram::drawChromosome(
    cairo_t* cr,
    double x,                                           // Pixel X coordinate
    double y,                                           // Pixel Y coordinate
    double height,                                      // Scaled chromosome pixel height
    double width,                                       // Chromosome width in pixels
    RGB color,                                          // Associated chromosome color
    double centromereStart,                             // Centromere start location as fraction of the chromosome length from p arm to q arm
    double centromereEnd)                               // Centromere end location as fraction of the chromosome length from p arm to q arm.
{

    // --------------------------------------
    // Centromere position and dimensions
    // --------------------------------------
    const double capRadius = width / 2.0;
    const double centromereCenter = (centromereStart + centromereEnd) / 2.0;
    const double centromereY = y + (height * centromereCenter);
    const double constrictionHeight = 7.0;
    const double constrictionAmount = 2.0;
    const double constrictionTop = centromereY - constrictionHeight / 2.0;
    const double constrictionBottom = centromereY + constrictionHeight / 2.0;

    cairo_new_path(cr);

    // ---------------------------------------
    // Top Cap
    // ---------------------------------------
    cairo_arc(
        cr,
        x + capRadius,
        y + capRadius,
        capRadius,
        M_PI,
        2 * M_PI
    );


    // -----------------------------------------
    // Right side approaching centromere
    // -----------------------------------------

    cairo_line_to(
        cr,
        x + width,
        constrictionTop
    );

    // Gradually move inward
    cairo_line_to(
        cr,
        x + width - constrictionAmount,
        centromereY
    );

    // Gradually move outward
    cairo_line_to(
        cr,
        x + width,
        constrictionBottom
    );


    // -----------------------------------------
    // Right side below centromere
    // -----------------------------------------

    cairo_line_to(
        cr,
        x + width,
        y + height - capRadius
    );

    // -----------------------------------------
    // Bottom cap
    // -----------------------------------------

    cairo_arc(
        cr,
        x + capRadius,
        y + height - capRadius,
        capRadius,
        0,
        M_PI
    );

    // -----------------------------------------
    // Left side approaching centromere
    // -----------------------------------------

    cairo_line_to(
        cr,
        x,
        constrictionBottom
    );

    // Gradually move inward
    cairo_line_to(
        cr,
        x + constrictionAmount,
        centromereY
    );

    // Gradually move outward
    cairo_line_to(
        cr,
        x,
        constrictionTop
    );

    // Close chromosome
    cairo_close_path(cr);

    // Fill
    cairo_set_source_rgb(
        cr,
        color.r,
        color.g,
        color.b
    );

    cairo_fill_preserve(cr);

    // Outline
    cairo_set_source_rgb(
        cr,
        0.1,
        0.1,
        0.1
    );

    cairo_set_line_width(cr, 1.5);

    cairo_stroke(cr);


    // -----------------------------------------
    // Drawing centromere ellipse
    // -----------------------------------------

    const double centromereStartY = y + height * centromereStart;

    const double centromereEndY = y + height * centromereEnd;

    const double centromereEllipseCenterY = (centromereStartY + centromereEndY) / 2.0;

    const double ellipseWidth = width + 2.0;

    const double ellipseHeight = std::max(5.0, centromereEndY - centromereStartY);

    cairo_save(cr);

    cairo_translate(cr, x + width / 2.0, centromereEllipseCenterY);

    cairo_scale(cr, ellipseWidth / 2.0, ellipseHeight / 2.0);

    cairo_arc(cr, 0.0, 0.0, 1.0, 0.0, 2.0 * M_PI);

    cairo_restore(cr);

    cairo_set_source_rgb(cr, 0.6, 0.6, 0.6); // Gray centromeres

    cairo_fill(cr);

}




void Karyogram::drawDoubleStrandBreakMarker(					// Function to draw black lines as DSB markers on th chromosomes
    cairo_t* cr,
    double x,									// DSB marker X location
    double y,									// DSB marker Y location
    double chromosomeWidth)							// DSB marker will span slightly more than the width of the chromosome for visual clarity.
{

    const double markerWidth = chromosomeWidth + 8.0; 				// Extend DSB markers past the chromosome width for clarity.
    const double markerLineWidth = 1.0;						// Thickness of DSB line marker

    // Black DSB line marker
    cairo_set_source_rgb(
        cr,
        0.0,
        0.0,
        0.0
    );

    cairo_set_line_width(
        cr,
        markerLineWidth
    );

    cairo_move_to(
        cr,
        x - markerWidth / 2.0,
        y
    );

    cairo_line_to(
        cr,
        x + markerWidth / 2.0,
        y
    );

    cairo_stroke(cr);

}










void Karyogram::drawSingleStrandBreakMarker(						// Function to draw white line SSB markers with marker length depending on the number of SSBs in the damage site.
    cairo_t* cr,
    double x,										// SSB marker x pixel position
    double y,										// SSB marker y pixel position
    double markerLength)								// The scaled SSB white line marker lengths
{
    // White SSB marker
    cairo_set_source_rgb(
        cr,
        1.0,
        1.0,
        1.0
    );

    cairo_set_line_width(cr, 1.5);

    cairo_move_to(
        cr,
        x - markerLength / 2.0,
        y
    );

    cairo_line_to(
        cr,
        x + markerLength / 2.0,
        y
    );

    cairo_stroke(cr);
}












void Karyogram::drawSDDsummary(								// Box at top of karyogram to summarize the following details in the associated SDD file:
    cairo_t* cr,
    const std::vector<double>& cellCyclePhase,						// Cell cycle phase from SDD file header to be summarized
    const std::vector<Exposure>& exposures,						// Vector of exposures to summarize the total number of SSBs and DSBs
    const std::vector<double>& doseOrFluence,						// Dose or fluence from SDD file to be summarized
    const std::vector<int>& incidentParticles)						// Incident particles from SDD file header to be summarized
{

    // --------------------------------------------------
    // Summary box position and dimensions
    // --------------------------------------------------
    const double summaryX = 50.0;
    const double summaryY = 25.0;

    const double summaryWidth = 900.0;
    const double summaryHeight = 128.0;


    // --------------------------------------------------
    // Count DSBs and SSBs
    // --------------------------------------------------
    int numberOfDSBs = 0;
    int numberOfSSBs = 0;

    for (const auto& exposure : exposures)
    {
        for (const auto& damage : exposure.damages)
        {
            if (damage.damageType.presenceOfDoubleStrandBreaks == 1)
            {
                ++numberOfDSBs;
            }

            numberOfSSBs += damage.damageType.numSingleBackboneBreaks;
        }
    }


    // -------------------------------------
    // Determine cell cycle phase
    // -------------------------------------
    std::string phaseLabel = cellCyclePhaseMeaning(cellCyclePhase);        	// Helper function in SDDutilities to interpret cell cycle phase code in SDD header.


    // -------------------------------------
    // Determine Dose or Fluence
    // -------------------------------------
    std::string doseOrFluenceSummary = doseOrFluenceMeaning(doseOrFluence);	// Helper function in SDDutilities to interpret dose or fluence in SDD header


    // -------------------------------------
    // Determine incident particle(s)
    // -------------------------------------
    std::string incidentParticlesLabel;
    if (incidentParticles.empty())
    {
        incidentParticlesLabel = "N/A";
    }
    else
    {
        for (std::size_t i = 0; i < incidentParticles.size(); i++)
        {
            incidentParticlesLabel += incidentParticlesMeaning(incidentParticles[i]);
            if (i + 1 < incidentParticles.size())
            {
                incidentParticlesLabel += ", ";
            }
        }
    }


    // --------------------------------------------------
    // Draw summary box
    // --------------------------------------------------

    cairo_set_source_rgb(
        cr,
        0.0,
        0.0,
        0.0
    );

    cairo_set_line_width(cr, 1.5);

    cairo_rectangle(
        cr,
        summaryX,
        summaryY,
        summaryWidth,
        summaryHeight
    );

    cairo_stroke(cr);


    // --------------------------------------------------
    // Draw summarized SDD header field summaries in box
    // --------------------------------------------------

    cairo_set_source_rgb(
        cr,
        0.0,
        0.0,
        0.0
    );

    cairo_select_font_face(
        cr,
        "Sans",
        CAIRO_FONT_SLANT_NORMAL,
        CAIRO_FONT_WEIGHT_NORMAL
    );

    cairo_set_font_size(cr, 20.0);

    // Three rows of text
    const double firstRowY = summaryY + 30.0;
    const double secondRowY = summaryY + 68.0;
    const double thirdRowY = summaryY + 106.0;


    // ------------------------------------------
    // TOP ROW
    // ------------------------------------------

    // Incident Particle(s)
    cairo_move_to(
        cr,
        summaryX + 15.0,
        firstRowY
    );

    cairo_show_text(
        cr,
        ("Incident Particle(s): " + incidentParticlesLabel).c_str()
    );


    // -----------------------------------------
    // SECOND ROW
    // -----------------------------------------

    // Dose / fluence
    cairo_move_to(
        cr,
        summaryX + 15.0,
        secondRowY
    );

    cairo_show_text(
        cr,
        ("Dose/Fluence: " + doseOrFluenceSummary).c_str()
    );

    // Cell cycle phase
    cairo_move_to(
        cr,
        summaryX + 470.0,
        secondRowY
    );

    cairo_show_text(
        cr,
        ("Cell Cycle Phase: " +
         phaseLabel).c_str()
    );


    // --------------------------------------------------
    // THIRD ROW
    // --------------------------------------------------

    // SSBs
    cairo_move_to(
        cr,
        summaryX + 15.0,
        thirdRowY
    );

    cairo_show_text(
        cr,
        ("SSBs: " +
         std::to_string(numberOfSSBs)).c_str()
    );

    // DSBs
    cairo_move_to(
        cr,
        summaryX + 470.0,
        thirdRowY
    );

    cairo_show_text(
        cr,
        ("DSBs: " +
         std::to_string(numberOfDSBs)).c_str()
    );
}








void Karyogram::drawLegend(					// Function to draw a legend of symbols at the bottom of the karyogram
    cairo_t* cr,
    double legendY)						// Height of the legend in pixels, adjusted depending on the height of the image, not hardcoded
{
    // --------------------------------------------------
    // Legend position
    // --------------------------------------------------
    const double legendX = 50.0;				// Hardcoded for now, may decide to adjust based on number of symbols like I did for legendY.


    // --------------------------------------------------
    // Legend border
    // --------------------------------------------------
    const double legendWidth = 900.0;
    const double legendHeight = 70.0;

    const double borderPaddingX = 15.0;
    const double borderPaddingY = 10.0;

    // Center chromosome vertically inside legend
    const double legendTop = legendY - 25.0 - borderPaddingY;
    const double legendCenterY = legendTop + legendHeight / 2.0;

    cairo_set_source_rgb(
        cr,
        0.0,
        0.0,
        0.0
    );

    cairo_set_line_width(cr, 1.5);

    cairo_rectangle(
        cr,
        legendX - borderPaddingX,
        legendY - 25.0 - borderPaddingY,
        legendWidth,
        legendHeight
    );

    cairo_stroke(cr);


    cairo_set_source_rgb(
        cr,
        0.1,
        0.1,
        0.1
    );

    cairo_select_font_face(
        cr,
        "Sans",
        CAIRO_FONT_SLANT_NORMAL,
        CAIRO_FONT_WEIGHT_NORMAL
    );

    cairo_set_font_size(cr, 16.0);


    // Centering legend labels
    cairo_font_extents_t fontExtents;
    cairo_font_extents(cr, &fontExtents);

    const double textBaseline = legendCenterY +
    (fontExtents.ascent - fontExtents.descent) / 2.0;


    // --------------------------------------------------
    // Legend title
    // --------------------------------------------------
    cairo_move_to(
        cr,
        legendX,
        textBaseline
    );

    cairo_show_text(
        cr,
        "Symbols:"
    );


    // -------------------------------------------------
    // Draw generic chromosome symbol
    // -------------------------------------------------
    const double chromosomeWidth = 14.0;
    const double chromosomeHeight = 40.0;

    // Position chromosome by its center
    const double chromosomeCenterX = 150.0;
    const double chromosomeX = chromosomeCenterX - chromosomeWidth / 2.0;
    const double chromosomeY = legendCenterY - chromosomeHeight / 2.0;

    // Generic centromere location
    const double centromereCenter = 0.5;
    const double constrictionAmount = 1.5;
    const double capRadius = chromosomeWidth / 2.0;
    const double centerY = chromosomeY + chromosomeHeight * centromereCenter;
    const double constrictionHeight = 5.0;
    const double constrictionTop = centerY - constrictionHeight / 2.0;
    const double constrictionBottom = centerY + constrictionHeight / 2.0;

    cairo_new_path(cr);

    // Top cap
    cairo_arc(
        cr,
        chromosomeX + capRadius,
        chromosomeY + capRadius,
        capRadius,
        M_PI,
        2.0 * M_PI
    );

    // Right side
    cairo_line_to(
        cr,
        chromosomeX + chromosomeWidth,
        constrictionTop
    );

    cairo_line_to(
        cr,
        chromosomeX + chromosomeWidth - constrictionAmount,
        centerY
    );

    cairo_line_to(
        cr,
        chromosomeX + chromosomeWidth,
        constrictionBottom
    );

    cairo_line_to(
        cr,
        chromosomeX + chromosomeWidth,
        chromosomeY + chromosomeHeight - capRadius
    );

    // Bottom cap
    cairo_arc(
        cr,
        chromosomeX + capRadius,
        chromosomeY + chromosomeHeight - capRadius,
        capRadius,
        0.0,
        M_PI
    );

    // Left side
    cairo_line_to(
        cr,
        chromosomeX,
        constrictionBottom
    );

    cairo_line_to(
        cr,
        chromosomeX + constrictionAmount,
        centerY
    );

    cairo_line_to(
        cr,
        chromosomeX,
        constrictionTop
    );

    cairo_close_path(cr);

    // Generic chromosome colour
    RGB chromosomeColor = generateChromosomeColor(0, 24);       // total number of chromosomes does not matter for generic chromosome drawing

    cairo_set_source_rgb(
        cr,
        chromosomeColor.r,
        chromosomeColor.g,
        chromosomeColor.b
    );

    cairo_fill_preserve(cr);

    // Outline
    cairo_set_source_rgb(
        cr,
        0.1,
        0.1,
        0.1
    );

    cairo_set_line_width(cr, 1.0);

    cairo_stroke(cr);

    // Chromosome label
    cairo_set_source_rgb(
        cr,
        0.0,
        0.0,
        0.0
    );

    cairo_move_to(
        cr,
        chromosomeX + chromosomeWidth + 10.0,
        textBaseline
    );

    cairo_show_text(
        cr,
        "Chromosome"
    );



    // --------------------------------------------------
    // DSB marker symbol
    // --------------------------------------------------
    const double dsbX = 330.0;

    const double dsbY = legendCenterY;

    cairo_set_source_rgb(
        cr,
        0.0,
        0.0,
        0.0
    );

    cairo_set_line_width(cr, 2.0);

    cairo_move_to(
        cr,
        dsbX - 10.0,
        dsbY
    );

    cairo_line_to(
        cr,
        dsbX + 10.0,
        dsbY
    );

    cairo_stroke(cr);

    // DSB label
    cairo_set_source_rgb(
        cr,
        0.0,
        0.0,
        0.0
    );

    cairo_move_to(
        cr,
        dsbX + 14.0,
        textBaseline
    );

    cairo_show_text(
        cr,
        "Double-Strand Break"
    );

    // --------------------------------------------------
    // Single-Strand Break symbol
    // --------------------------------------------------

    const double ssbX = 550.0;
    const double ssbY = legendCenterY;

    const double ssbMarkerWidth = 14.0;
    const double ssbMarkerHeight = 3.0;

    // Draw white rectangle
    cairo_rectangle(
        cr,
        ssbX - ssbMarkerWidth / 2.0,
        ssbY - ssbMarkerHeight / 2.0,
        ssbMarkerWidth,
        ssbMarkerHeight
    );

    cairo_set_source_rgb(
        cr,
        1.0,
        1.0,
        1.0
    );

    cairo_fill_preserve(cr);

    // Black outline
    cairo_set_source_rgb(
        cr,
        0.0,
        0.0,
        0.0
    );

    cairo_set_line_width(cr, 1.0);

    cairo_stroke(cr);

    // SSB label
    cairo_set_source_rgb(
        cr,
        0.0,
        0.0,
        0.0
    );

    cairo_move_to(
        cr,
        ssbX + 10.0,
        textBaseline
    );

    cairo_show_text(
        cr,
        "Single-Strand Break"
    );



    // --------------------------------------------------
    // Centromere marker symbol
    // --------------------------------------------------
    const double centromereX = 770.0;

    const double centromereY =
        legendCenterY;

    const double ellipseWidth = 14.0;
    const double ellipseHeight = 8.0;

    cairo_save(cr);

    cairo_translate(
        cr,
        centromereX,
        centromereY
    );

    cairo_scale(
        cr,
        ellipseWidth / 2.0,
        ellipseHeight / 2.0
    );

    cairo_arc(
        cr,
        0.0,
        0.0,
        1.0,
        0.0,
        2.0 * M_PI
    );

    cairo_restore(cr);

    // Gray fill
    cairo_set_source_rgb(
        cr,
        0.6,
        0.6,
        0.6
    );

    cairo_fill_preserve(cr);

    // Centromere label
    cairo_set_source_rgb(
        cr,
        0.0,
        0.0,
        0.0
    );

    cairo_move_to(
        cr,
        centromereX + 10.0,
        textBaseline
    );

    cairo_show_text(
        cr,
        "Centromere"
    );
}
