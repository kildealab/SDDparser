#include <algorithm>
#include <cmath>
#include <iostream>

#include "Karyogram.h"

bool Karyogram::generateKaryogram(
    const std::vector<double>& chromosomeSizes,
    const std::vector<double>& cellCyclePhase,
    const std::vector<Exposure>& exposures,
    const std::string& outputFilename)
{
    if (chromosomeSizes.empty())
    {
        std::cerr << "ERROR: No chromosome sizes provided.\n";
        return false;
    }

    // Use this vector for plotting DSBs onto karyogram
    std::vector<DSBlocation> dsbs = getDoubleStrandBreaks(exposures);

    // Store the geometry of the chromosomes to draw the DSBs.
    std::vector<ChromosomeGeometry> chromosomeGeometry;

    // The first entry specifies the number of chromosomes.
    const int chromosomeCount =
        static_cast<int>(chromosomeSizes[0]);

    // Store number of homologs for DSB drawing 
    const int homologousSetSize = (chromosomeCount - 2) / 2;

    // Check that the vector contains the expected
    // number of chromosome sizes.
    if (chromosomeSizes.size() !=
        static_cast<size_t>(chromosomeCount + 1))
    {
        std::cerr << "ERROR: Chromosome count does not match "
                  << "the number of chromosome sizes provided.\n";
        return false;
    }

    // Check cell cycle phase to determine number of chromatids per chromosomes to draw.
    bool doubleChromatid = false;

    if (!cellCyclePhase.empty())
    {
        const int phase = static_cast<int>(cellCyclePhase[0]);

        if (phase == 3 || phase == 4 || phase == 5)
        {
            doubleChromatid = true;
    	}
    }

    // We assume the last two chromosomes in the headerr are X and Y.
    if (chromosomeCount < 4 || (chromosomeCount - 2) % 2 != 0)
    {
        std::cerr << "ERROR: Chromosome count is not compatible "
                  << "with the expected homologous-pair + X/Y format.\n";
        return false;
    }

    // Number of homologous chromosome pairs.
    const int homologousPairs = (chromosomeCount - 2) / 2;


    // Find largest chromosome
    double maxBP = 0.0;

    for (size_t i = 1 ; i < chromosomeSizes.size(); i++)
    {
        if (chromosomeSizes[i] > maxBP)
        {
            maxBP = chromosomeSizes[i];
        }
    }

    const int imgWidth = 1000;
//    const int imgHeight = 3100; // was 700, but only showed 12 chroms

    const int columns = 4;
    const double colWidth = static_cast<double>(imgWidth) / columns;

    const double rowHeight = 250.0;
    const double startY = 40.0;
//    const double pairGap = 8.0;  NOT DRAWING HOMOLOGOUS KARYOGRAM YET

    // Draw homologous pairs first
    const int rows = (homologousPairs + columns - 1) / columns;
    const int imgHeight = static_cast<int>(startY + rows * rowHeight);

    const double maxRenderHeight = 180.0;

    const double chromosomeWidth = 14.0;
    const double chromatidGap = 1.5;
    
    cairo_surface_t* surface =
        cairo_image_surface_create(
            CAIRO_FORMAT_ARGB32,
            imgWidth,
            imgHeight
        );

    cairo_t* cr = cairo_create(surface);

    cairo_set_antialias(cr, CAIRO_ANTIALIAS_BEST);

    // White background
    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    cairo_paint(cr);

    // -------------------------------------------
    // Draw homologous chromosome pairs first
    // -------------------------------------------

    for (int i = 0; i < homologousPairs; ++i)
    {
	// First set
	const size_t firstIndex = static_cast<size_t>(i + 1);

	// Second set begins after the first set
	const size_t secondIndex = static_cast<size_t>(i + 1 + homologousPairs);

        int col = i % columns;
        int row = i / columns;

        double groupCenterX =
            (col * colWidth) + (colWidth / 2.0);

        double posY =
            startY + (row * rowHeight);


	// --------------------------------------
	// First Homolog
	// --------------------------------------

	double firstHeight =
            (chromosomeSizes[firstIndex] / maxBP) *
            maxRenderHeight;

        if (firstHeight < 30.0)
        {
            firstHeight = 30.0;
        }


        // ------------------------------------------
        // Second homolog
        // ------------------------------------------

        double secondHeight =
            (chromosomeSizes[secondIndex] / maxBP) *
            maxRenderHeight;

        if (secondHeight < 30.0)
        {
            secondHeight = 30.0;
        }


	// ------------------------------------------
        // Position homologues beside each other
        // ------------------------------------------

        // Generate chromosome colours according to homologous pairs.
    	RGB chromosomeColor = generateChromosomeColor(i);

	const double homologGap = 18.0;

	if (!doubleChromatid)
	{
    	// --------------------------------------------------
    	// Single-chromatid chromosomes
    	// --------------------------------------------------

    	    double leftChromosomeX = groupCenterX - chromosomeWidth - (homologGap / 2.0);

    	    double rightChromosomeX = groupCenterX + (homologGap / 2.0);

    	    drawChromosome(cr, leftChromosomeX, posY, firstHeight, chromosomeColor);

    	    drawChromosome(cr, rightChromosomeX, posY, secondHeight, chromosomeColor);

	    chromosomeGeometry.push_back({i + 1, 1, leftChromosomeX, posY, firstHeight});

	    chromosomeGeometry.push_back({i + homologousSetSize + 1, 1, rightChromosomeX, posY, secondHeight});
	}

	else
	{
    	// --------------------------------------------------
    	// Double-chromatid chromosomes
    	// --------------------------------------------------

            const double homologGap = 40.0;

	    // Homologue 1
    	    double leftChromatid1X = groupCenterX - chromosomeWidth - homologGap / 2.0 - chromatidGap / 2.0;

    	    double leftChromatid2X = leftChromatid1X + chromosomeWidth + chromatidGap;

    	    // Homologue 2
    	    double rightChromatid1X = groupCenterX + homologGap / 2.0;

    	    double rightChromatid2X = rightChromatid1X + chromosomeWidth + chromatidGap;

	    // Draw homologue 1
    	    drawChromosome(cr, leftChromatid1X, posY, firstHeight, chromosomeColor);

    	    drawChromosome(cr, leftChromatid2X, posY, firstHeight, chromosomeColor);

    	    // Draw homologue 2
    	    drawChromosome(cr, rightChromatid1X, posY, secondHeight, chromosomeColor);

    	    drawChromosome(cr, rightChromatid2X, posY, secondHeight, chromosomeColor);

	    // Store drawn chromosome geometries for DSB placement
	    chromosomeGeometry.push_back({i + 1, 1, leftChromatid1X, posY, firstHeight});

	    chromosomeGeometry.push_back({i + 1, 2, leftChromatid2X, posY, firstHeight});

	    chromosomeGeometry.push_back({i + homologousSetSize + 1, 1, rightChromatid1X, posY, secondHeight});

	    chromosomeGeometry.push_back({i + homologousSetSize + 1, 2, rightChromatid2X, posY, secondHeight});

	}

	// ------------------------------------------
        // Label
	// ------------------------------------------

	double labelHeight = std::max(firstHeight, secondHeight);

        cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);

        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);

        cairo_set_font_size(cr, 16.0);

        std::string label =
            "Chr " + std::to_string(i + 1);

        cairo_move_to(cr, groupCenterX - 20.0, posY + labelHeight + 25.0);

        cairo_show_text(cr, label.c_str());
    }

    // --------------------------------------------------
    // X and Y chromosome drawing
    // --------------------------------------------------

    const size_t xIndex = static_cast<size_t>(chromosomeCount);

    const size_t yIndex = static_cast<size_t>(chromosomeCount - 1);

    double xHeight = (chromosomeSizes[xIndex] / maxBP) * maxRenderHeight;

    double yHeight = (chromosomeSizes[yIndex] / maxBP) * maxRenderHeight;

    if (xHeight < 30.0)
    {
        xHeight = 30.0;
    }

    if (yHeight < 30.0)
    {
        yHeight = 30.0;
    }


    // X and Y occupy the remaining positions in the final row.
    const int sexChromosomeRow = (homologousPairs - 1) / columns;

    const int xColumn = 2;
    const int yColumn = 3;

    const double xGroupCenterX = (xColumn * colWidth) + (colWidth / 2.0);

    const double yGroupCenterX = (yColumn * colWidth) + (colWidth / 2.0);

    const double sexChromosomeY = startY + (sexChromosomeRow * rowHeight);

    // Center each chromosome in its column.
    const double xChromosomeX = xGroupCenterX - 7.0;

    const double yChromosomeX = yGroupCenterX - 7.0;

    RGB yColor = generateChromosomeColor(homologousPairs);

    RGB xColor = generateChromosomeColor(homologousPairs + 1);

    if (!doubleChromatid)
    {
    	// Draw X
    	drawChromosome(cr, xChromosomeX, sexChromosomeY, xHeight, xColor);
	// Store geometry of X chromosome for DSB drawing
	chromosomeGeometry.push_back({chromosomeCount, 1, xChromosomeX, sexChromosomeY, xHeight});

        // Draw Y
        drawChromosome(cr, yChromosomeX, sexChromosomeY, yHeight, yColor);
	// Store geometry of Y chromosome for DSB drawing
	chromosomeGeometry.push_back({chromosomeCount - 1, 1, yChromosomeX, sexChromosomeY, yHeight});

    }
    else
    {
	// Draw X chromatids 1 and 2
	drawChromosome(cr, xChromosomeX, sexChromosomeY, xHeight, xColor);
	drawChromosome(cr, xChromosomeX + chromosomeWidth + chromatidGap,
		       sexChromosomeY, xHeight, xColor);
	// Store X chromatids 1 and 2 geometries for DSB drawing
	chromosomeGeometry.push_back({chromosomeCount, 1, xChromosomeX, sexChromosomeY, xHeight});
 	chromosomeGeometry.push_back({chromosomeCount, 2, xChromosomeX + chromosomeWidth + chromatidGap, sexChromosomeY, xHeight});

        // Draw Y chromatids 1 and 2
        drawChromosome(cr, yChromosomeX, sexChromosomeY, yHeight, yColor);
        drawChromosome(cr, yChromosomeX + chromosomeWidth + chromatidGap, 
		       sexChromosomeY, yHeight, yColor);

	// Store Y chromatids 1 and 2 geometries for DSB drawing
	chromosomeGeometry.push_back({chromosomeCount - 1, 1, yChromosomeX, sexChromosomeY, yHeight});
	chromosomeGeometry.push_back({chromosomeCount - 1, 2, yChromosomeX + chromosomeWidth + chromatidGap, sexChromosomeY, yHeight});

    }


    // --------------------------------------------------
    // Draw double-strand break markers
    // --------------------------------------------------

    for (const auto& dsb : dsbs)
    {
        const int chromosomeID = dsb.chromosomeID.chromosomeNumber;

        const int chromatidID = dsb.chromosomeID.chromatidNumber;

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
            // Convert SDD position to fraction
            // ------------------------------------------

            double damageFraction = getDamageFraction(dsb, chromosomeSize);

            // Keep fraction safely within chromosome
            damageFraction = std::clamp(damageFraction, 0.0, 1.0);

            // ------------------------------------------
            // Convert fraction to image coordinates
            // ------------------------------------------

            const double markerX = geometry.x + chromosomeWidth / 2.0;

            const double markerY = geometry.y + damageFraction * geometry.height;

            // ------------------------------------------
            // Draw DSB marker
            // ------------------------------------------

            drawDamageMarker(cr, markerX, markerY);

            // We found the correct geometry, so stop searching.
            break;
        }
    }

    // Y label
    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);

    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);

    cairo_set_font_size(cr, 16.0);

    cairo_move_to(cr, xGroupCenterX - 4.0, sexChromosomeY + xHeight + 25.0);

    cairo_show_text(cr, "X");

    // X label
    cairo_move_to(cr, yGroupCenterX - 4.0, sexChromosomeY + yHeight + 25.0);

    cairo_show_text(cr, "Y");


    cairo_status_t status = cairo_surface_write_to_png(surface, outputFilename.c_str());

    cairo_destroy(cr);
    cairo_surface_destroy(surface);

    if (status != CAIRO_STATUS_SUCCESS)
    {
        std::cerr << "ERROR: Failed to write karyogram image.\n";
        return false;
    }

    return true;
}


std::vector<DSBlocation>
Karyogram::getDoubleStrandBreaks(const std::vector<Exposure>& exposures)
{
    std::vector<DSBlocation> dsbs;

    for (const auto& exposure : exposures)
    {
	for (const auto& damage : exposure.damages)
	{
            if (damage.damageType.presenceOfDoubleStrandBreaks == 1)
            {
		DSBlocation dsb;
		dsb.chromosomeID = damage.chromosomeID;
		dsb.chromosomePosition = damage.chromosomePosition;
                dsbs.push_back(dsb);
            }
	}
    }
    return dsbs;
}


double Karyogram::getDamageFraction(const DSBlocation& dsb, double chromosomeSize)
{
    double position = dsb.chromosomePosition.position;

    // Field 4 is already a fraction from the beginning of the p arm to the end of the q arm.
    if (dsb.chromosomePosition.isFractional)
    {
        return position;
    }

    if (chromosomeSize <= 0.0)
    {
        return 0.0;
    }

    // p arm = top, q arm = bottom, chromosome sizes stored in Mbp
    double chromosomeSizeBP = chromosomeSize * 1000000.0;

    return position / chromosomeSizeBP;
}


void Karyogram::drawChromosome(
    cairo_t* cr,
    double x,
    double y,
    double height,
    RGB color)
{
    const double width = 14.0;

    const double centromereRatio = 0.35;

    const double centromereY =
        y + (height * centromereRatio);

    const double capRadius = width / 2.0;

    cairo_new_path(cr);

    // Top cap
    cairo_arc(
        cr,
        x + capRadius,
        y + capRadius,
        capRadius,
        M_PI,
        2 * M_PI
    );

    cairo_line_to(
        cr,
        x + width - 1.5,
        centromereY - 2.0
    );

    cairo_line_to(
        cr,
        x + width - 1.5,
        centromereY + 2.0
    );

    cairo_line_to(
        cr,
        x + width,
        y + height - capRadius
    );

    // Bottom cap
    cairo_arc(
        cr,
        x + capRadius,
        y + height - capRadius,
        capRadius,
        0,
        M_PI
    );

    cairo_line_to(
        cr,
        x + 1.5,
        centromereY + 2.0
    );

    cairo_line_to(
        cr,
        x + 1.5,
        centromereY - 2.0
    );

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

    // Centromere
    cairo_set_source_rgb(
        cr,
        0.1,
        0.1,
        0.1
    );

    cairo_move_to(
        cr,
        x,
        centromereY
    );

    cairo_line_to(
        cr,
        x + width,
        centromereY
    );

    cairo_stroke(cr);
}


RGB Karyogram::generateChromosomeColor(int chromosomeNumber)
{
    // Evenly distribute hues around the color wheel.
//    double hue =
//        static_cast<double>(chromosomeNumber) /
//        static_cast<double>(totalChromosomes);

//    hue *= 360.0;


    // Golden-ratio spacing gives better visual separation
    // between consecutive chromosome colors.
    const double goldenRatio = 0.618033988749895;

    double hue = std::fmod(chromosomeNumber * goldenRatio, 1.0);

    const double saturation = 0.70;
    const double value = 0.90;

    double c = value * saturation;

    double h = hue * 6.0;

    double x = c * (1.0 - std::abs(std::fmod(h, 2.0) - 1.0));

    double r = 0.0;
    double g = 0.0;
    double b = 0.0;

    if (h < 1.0)
    {
        r = c;
        g = x;
    }
    else if (h < 2.0)
    {
        r = x;
        g = c;
    }
    else if (h < 3.0)
    {
        g = c;
        b = x;
    }
    else if (h < 4.0)
    {
        g = x;
        b = c;
    }
    else if (h < 5.0)
    {
        r = x;
        b = c;
    }
    else
    {
        r = c;
        b = x;
    }

    double m = value - c;

    r += m;
    g += m;
    b += m;

    return {r, g, b};
}

void Karyogram::drawDamageMarker(
    cairo_t* cr,
    double x,
    double y)
{
    const double markerRadius = 2.0;

    // Black marker
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);

    cairo_arc(
        cr,
        x,
        y,
        markerRadius,
        0.0,
        2.0 * M_PI
    );

    cairo_fill(cr);
}
