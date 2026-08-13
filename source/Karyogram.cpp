#include <algorithm>
#include <cmath>
#include <iostream>

#include "Karyogram.h"

bool Karyogram::generateKaryogram(
    const std::vector<double>& chromosomeSizes,
    const std::string& outputFilename)
{
    if (chromosomeSizes.empty())
    {
        std::cerr << "ERROR: No chromosome sizes provided.\n";
        return false;
    }

    // The first entry specifies the number of chromosomes.
    const int chromosomeCount =
        static_cast<int>(chromosomeSizes[0]);

    // Check that the vector contains the expected
    // number of chromosome sizes.
    if (chromosomeSizes.size() !=
        static_cast<size_t>(chromosomeCount + 1))
    {
        std::cerr << "ERROR: Chromosome count does not match "
                  << "the number of chromosome sizes provided.\n";
        return false;
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

        const double homologGap = 8.0;


	double leftChromosomeX = groupCenterX - 14.0 - (homologGap / 2.0);

        double rightChromosomeX = groupCenterX + (homologGap / 2.0);

	// Draw first homolog
        drawChromosome(cr, leftChromosomeX, posY, firstHeight);

	// Draw second homolog
        drawChromosome(cr, rightChromosomeX, posY, secondHeight);


	// ------------------------------------------
        // Label
	// ------------------------------------------

	double labelHeight = std::max(firstHeight, secondHeight);

        cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);

        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);

        cairo_set_font_size(cr, 13.0);

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

    // Draw Y
    drawChromosome(cr, xChromosomeX, sexChromosomeY, xHeight);

    // Draw X
    drawChromosome(cr, yChromosomeX, sexChromosomeY, yHeight);

    // Y label
    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);

    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);

    cairo_set_font_size(cr, 13.0);

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



void Karyogram::drawChromosome(
    cairo_t* cr,
    double x,
    double y,
    double height)
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
        0.8,
        0.8,
        0.8
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
