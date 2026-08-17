#include <algorithm>
#include <cmath>
#include <iostream>

#include "Karyogram.h"

const std::vector<CentromerePosition> humanCentromeres =			// Vector storing the associated chromosome ID centromere positions to draw the centromere location if the human option of '--karyotype' is passed by the user in the command-line
{										// Vector stores the CentromerePosition struct
    {"1", 121535434, 124535434},						// Each vector element contains {'Chromosome Number', Centromere Start Position (bp), Centromere End Position (bp)}
    {"2",  92326171,  95326171},
    {"3",  90504854,  93504854},
    {"4",  49660117,  52660117},
    {"5",  46405641,  49405641},
    {"6",  58830166,  61830166},
    {"7",  58054331,  61054331},
    {"8",  43838887,  46838887},
    {"9",  47367679,  50367679},
    {"10", 39254935,  42245935},
    {"11", 51644205,  54644205},
    {"12", 34856694,  37856694},
    {"13", 16000000,  19000000},
    {"14", 16000000,  19000000},
    {"15", 17000000,  20000000},
    {"16", 35335801,  38335801},
    {"17", 22263006,  25263006},
    {"18", 15460898,  18460898},
    {"19", 24681782,  27681782},
    {"20", 26369569,  29369569},
    {"21", 11288129,  14288129},
    {"22", 13000000,  16000000},
    {"Y", 10100001, 10327000},
    {"X", 58632012, 61632012}
};


const CentromerePosition* Karyogram::getHumanCentromere(			// Function to return the centromere start and end positions for each chromosome number.
    int chromosomeID)
{
    std::string chromosome;

    if (chromosomeID >= 1 && chromosomeID <= 22)
    {
        chromosome = std::to_string(chromosomeID);
    }
    else if (chromosomeID >= 23 && chromosomeID <= 44)				// Chromosome number 23-44 correspond to the homologs of chromosomes 1-22. Y is chromosome 45 and X is chromosome 46 --> MUST MAINTAIN THIS ORDERING
    {
        // Second homologous set.
        chromosome = std::to_string(chromosomeID - 22);
    }
    else if (chromosomeID == 45)						// Chromosome Number 45 = Y chromosome
    {
        chromosome = "Y";
    }
    else if (chromosomeID == 46)						// Chromosome Number 46 = X chromosome
    {
        chromosome = "X";
    }
    else
    {
        return nullptr;								// If other Chromosome Numbers encountered, return nullptr.
    }

    for (const auto& centromere : humanCentromeres)				// Return the centromere start and end locations for each chromosome.
    {
        if (centromere.chromosome == chromosome)
        {
            return &centromere;
        }
    }

    return nullptr;
}




bool Karyogram::generateKaryogram(						// Function that generates the overall Karyogram structure
    const std::vector<double>& chromosomeSizes,					// Requires chromosome sizes in the SDD file header for scaling on the image
    const std::vector<double>& cellCyclePhase,					// Requires cell cycle phase in the SDD file header for determining presence of duplicated chromosomes.
    const std::vector<Exposure>& exposures,					// Accesses the damage locations associated with each chromosome in a given exposure.
    bool humanGenome,								// Boolean to determine if karyogram should be drawn with human centromere ranges or generic centromere locations for non-human genomes.
    const std::string& outputFilename)						// Outputs the image to the outputFilename.
{

    if (chromosomeSizes.empty())						// Make sure SDD header 'Chromosome sizes' Is non-empty.
    {
        std::cerr << "ERROR: No chromosome sizes provided.\n";
        return false;
    }

    std::vector<DSBlocation> dsbVec = getDoubleStrandBreaks(exposures);		// Use this vector for plotting DSBs onto karyogram

    std::vector<ChromosomeGeometry> chromosomeGeometry;				// Store the geometry of the chromosomes in x and y pixel coordinates to draw the DSBs.

    const int chromosomeCount = static_cast<int>(chromosomeSizes[0]);		// The first entry specifies the number of chromosomes.
        
    const int homologousSetSize = (chromosomeCount - 2) / 2;			// Store number of homologs for DSB drawing 

    const double genericCentromereStart = 0.34;					// In case non-human chromosomes specified, draw generic centromere regions
    const double genericCentromereEnd = 0.36;					// Will be replaced later if human genome is specified with actual centromere positions.

    if (chromosomeSizes.size() !=						// Check that the vector contains the expected number of chromosome sizes.
        static_cast<size_t>(chromosomeCount + 1))
    {
        std::cerr << "ERROR: Chromosome count does not match "
                  << "the number of chromosome sizes provided.\n";
        return false;
    }

    bool doubleChromatid = false;						// Check cell cycle phase to determine number of chromatids per chromosomes to draw.

    if (!cellCyclePhase.empty())						// Cell cycle phase of 3, 4, or 5 indicates post-replicated DNA, two chromatids per chromosome.
    {
        const int phase = static_cast<int>(cellCyclePhase[0]);

        if (phase == 3 || phase == 4 || phase == 5)
        {
            doubleChromatid = true;
    	}
    }

    if (chromosomeCount < 4 || (chromosomeCount - 2) % 2 != 0)			// For our purposes the last two chromosome IDs must be Y, then X.
    {
        std::cerr << "ERROR: Chromosome count is not compatible "
                  << "with the expected homologous-pair + X/Y format.\n";
        return false;
    }

    const int homologousPairs = (chromosomeCount - 2) / 2;			// Store number of homologous Pairs for karyogram generation.

    double maxBP = 0.0;								// Find largest chromosome, scale all chromosomes with respect to the largest
    for (size_t i = 1 ; i < chromosomeSizes.size(); i++)
    {
        if (chromosomeSizes[i] > maxBP)
        {
            maxBP = chromosomeSizes[i];
        }
    }

    // ----------------------------------------------------
    // Karyogram dimensions
    // ----------------------------------------------------

    const int imgWidth = 1000;							// Image width in number of pixels
    const int columns = 4;							// Desired number of columns in the karyogram
    const double colWidth = static_cast<double>(imgWidth) / columns;		// Column width also in number of pixels
    const double rowHeight = 250.0;						// Row height in number of pixels
    const double startY = 40.0;							// Choosing the starting Y position for the first row of chromosomes
    const int rows = (homologousPairs + columns - 1) / columns;			// Determine number of rows based on number of chromosomes passed.
    const int imgHeight = static_cast<int>(startY + rows * rowHeight);		// Adjust image height based on the number of rows, also in pixels.
    const double maxRenderHeight = 180.0;					// Limit the maximum height of the image, also in pixels
    const double chromosomeWidth = 14.0;					// Individual chromosomes are 14 pixels wide.
    const double chromatidGap = 1.5;						// The gaps between homologs in a replicated chromosome are 1.5 pixels.
    

    // ----------------------------------------------------
    // Generate the background Karyogram template
    // ----------------------------------------------------
    cairo_surface_t* surface =							// Generate the Karyogram template for the chromosomes to be drawn on.
        cairo_image_surface_create(						// Use the image dimensions previously defined.
            CAIRO_FORMAT_ARGB32,
            imgWidth,
            imgHeight
        );

    cairo_t* cr = cairo_create(surface);

    cairo_set_antialias(cr, CAIRO_ANTIALIAS_BEST);				// Best image quality and curve smoothness.

    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);					// White background
    cairo_paint(cr);

    // -------------------------------------------
    // Draw homologous chromosome pairs first
    // -------------------------------------------

    for (int i = 0; i < homologousPairs; ++i)
    {
	const size_t firstIndex = static_cast<size_t>(i + 1); 			// First homolog

	const size_t secondIndex = static_cast<size_t>(i + 1 + homologousPairs);// Second homolog begins after the first 

	// Initialize centromere start and end locations with the generic locations in case user specified option '--karyogram other'
	double firstCentromereStart = genericCentromereStart;			// First homolog centromere start position (from p arm to q arm)
	double firstCentromereEnd = genericCentromereEnd;			// First homolog centromere end position (from p arm to q arm)
	double secondCentromereStart = genericCentromereStart;			// Second homolog centromere start position (from p arm to q arm)
	double secondCentromereEnd = genericCentromereEnd;			// Second homolog centromere end position (from p arm to q arm)
	
	if (humanGenome)							// Check if humanGenome option passed by the user is true.
	{
    	    const CentromerePosition* centromere = getHumanCentromere(i + 1);

    	    if (centromere != nullptr)
    	    {
		// Scale centromere positions to convert to base pairs from mega base pairs. All centromere locations are stored as fractions of the chromosome length afterwards.
	        firstCentromereStart = static_cast<double>(centromere->start) / (chromosomeSizes[firstIndex] * 1000000.0);

                firstCentromereEnd = static_cast<double>(centromere->end) / (chromosomeSizes[firstIndex] * 1000000.0);

        	secondCentromereStart = static_cast<double>(centromere->start) / (chromosomeSizes[secondIndex] * 1000000.0);

        	secondCentromereEnd = static_cast<double>(centromere->end) / (chromosomeSizes[secondIndex] * 1000000.0);
	    }
	}

        int col = i % columns;
        int row = i / columns;

        double groupCenterX =
            (col * colWidth) + (colWidth / 2.0);

        double posY =
            startY + (row * rowHeight);


	// --------------------------------------
	// First Homolog
	// --------------------------------------

	double firstHeight = (chromosomeSizes[firstIndex] / maxBP) * maxRenderHeight; 	// Scale chromosome heights by the largest chromosome.

        if (firstHeight < 30.0)							      	// Prevent chromosomes from being smaller than 30 pixels
        {
            firstHeight = 30.0;
        }

        // ------------------------------------------
        // Second homolog
        // ------------------------------------------

        double secondHeight = (chromosomeSizes[secondIndex] / maxBP) * maxRenderHeight;	// Scale chromosome heights by the largest chromosome

        if (secondHeight < 30.0)						       	// Prevent chromosomes from being smaller than 30 pixels.
        {
            secondHeight = 30.0;
        }

	// ------------------------------------------
        // Position homologues beside each other
        // ------------------------------------------

    	RGB chromosomeColor = generateChromosomeColor(i);				// Generate chromosome colours according to homologous pairs.

	const double homologGap = 18.0;							// Set gap between homologous pairs to be 18 pixels

	if (!doubleChromatid)								// Depending on SDD header 'Cell cycle phase' entry, doubleChromatid is determined to be true or false.
	{
    	// --------------------------------------------------
    	// Single-chromatid chromosomes
    	// --------------------------------------------------

    	    double leftChromosomeX = groupCenterX - chromosomeWidth - (homologGap / 2.0);// Left homolog position

    	    double rightChromosomeX = groupCenterX + (homologGap / 2.0);		// Right homolog position

	    // Draw homologous pair at the desired coordinates with the desired heights, color, and centromere positions.
	    drawChromosome(cr, leftChromosomeX, posY, firstHeight, chromosomeColor, firstCentromereStart, firstCentromereEnd); 
	    drawChromosome(cr, rightChromosomeX, posY, secondHeight, chromosomeColor, secondCentromereStart, secondCentromereEnd);

	    // Store chromosome geometry (pixel coordinates) to draw corresponding chromosome damages from exposure data.
	    chromosomeGeometry.push_back({i + 1, 1, leftChromosomeX, posY, firstHeight});
	    chromosomeGeometry.push_back({i + homologousSetSize + 1, 1, rightChromosomeX, posY, secondHeight});
	}

	else
	{
    	// --------------------------------------------------
    	// Double-chromatid chromosomes
    	// --------------------------------------------------

            const double homologGap = 40.0;							// Set gap of 40 pixels between successive homologous pairs for post-replicated chromosomes.


	    // Homologue 1 left and right chromatid positions
    	    double leftChromatid1X = groupCenterX - chromosomeWidth - homologGap / 2.0 - chromatidGap / 2.0;
    	    double leftChromatid2X = leftChromatid1X + chromosomeWidth + chromatidGap;

    	    // Homologue 2 left and right chromatid positions
    	    double rightChromatid1X = groupCenterX + homologGap / 2.0;
    	    double rightChromatid2X = rightChromatid1X + chromosomeWidth + chromatidGap;

	    // Draw homologue 1
	    drawChromosome(cr, leftChromatid1X, posY, firstHeight, chromosomeColor, firstCentromereStart, firstCentromereEnd);
	    drawChromosome(cr, leftChromatid2X, posY, firstHeight, chromosomeColor, firstCentromereStart, firstCentromereEnd);

    	    // Draw homologue 2
	    drawChromosome(cr, rightChromatid1X, posY, secondHeight, chromosomeColor, secondCentromereStart, secondCentromereEnd);
	    drawChromosome(cr, rightChromatid2X, posY, secondHeight, chromosomeColor, secondCentromereStart, secondCentromereEnd);

	    // Store drawn chromosome geometries (pixel coordinates) to draw corresponding damages from exposure data.
	    chromosomeGeometry.push_back({i + 1, 1, leftChromatid1X, posY, firstHeight});		// Homolog 1 chromatid 1
	    chromosomeGeometry.push_back({i + 1, 2, leftChromatid2X, posY, firstHeight});		// Homolog 1 chromatid 2
	    chromosomeGeometry.push_back({i + homologousSetSize + 1, 1, rightChromatid1X, posY, secondHeight}); // Homolog 2 chromatid 1
	    chromosomeGeometry.push_back({i + homologousSetSize + 1, 2, rightChromatid2X, posY, secondHeight}); // Homolog 2 chromatid 2

	}

	// ------------------------------------------
        // Label
	// ------------------------------------------
	// Labeling each chromosome 1-22, Y and X. Homologous pairs are given a single chromosome ID.
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

    const size_t xIndex = static_cast<size_t>(chromosomeCount);				// Determining the index of the X chromosome
	
    const size_t yIndex = static_cast<size_t>(chromosomeCount - 1);			// Determining the index of the Y chromosome

    // Set the centromere locations to the generic positions, replaced if user passed '--karyogram human' option.
    double xCentromereStart = genericCentromereStart;
    double xCentromereEnd = genericCentromereEnd;

    double yCentromereStart = genericCentromereStart;
    double yCentromereEnd = genericCentromereEnd;

    if (humanGenome)									// If user specified a human genome '--karyogram human'
    {
	// Get the base pair locations of the centromeres for the X and Y chromosomes
        const CentromerePosition* xCentromere = getHumanCentromere(chromosomeCount);    
        const CentromerePosition* yCentromere = getHumanCentromere(chromosomeCount - 1);

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
    const int sexChromosomeRow = (homologousPairs - 1) / columns;
    // Plot X first then Y next.
    const int xColumn = 2;
    const int yColumn = 3;
    const double xGroupCenterX = (xColumn * colWidth) + (colWidth / 2.0);
    const double yGroupCenterX = (yColumn * colWidth) + (colWidth / 2.0);
    const double sexChromosomeY = startY + (sexChromosomeRow * rowHeight);
    const double xChromosomeX = xGroupCenterX - 7.0;
    const double yChromosomeX = yGroupCenterX - 7.0;


    // Set X and Y chromosome Colors
    RGB yColor = generateChromosomeColor(homologousPairs);
    RGB xColor = generateChromosomeColor(homologousPairs + 1);

    if (!doubleChromatid)							// Depending on if SDD header 'Cell cycle phase' specified a post-replicated chromosome.
    {
    	// Draw X
    	drawChromosome(cr, xChromosomeX, sexChromosomeY, xHeight, xColor, xCentromereStart, xCentromereEnd);
	// Store geometry of X chromosome for DSB drawing
	chromosomeGeometry.push_back({chromosomeCount, 1, xChromosomeX, sexChromosomeY, xHeight});

        // Draw Y
        drawChromosome(cr, yChromosomeX, sexChromosomeY, yHeight, yColor, yCentromereStart, yCentromereEnd);
	// Store geometry of Y chromosome for DSB drawing
	chromosomeGeometry.push_back({chromosomeCount - 1, 1, yChromosomeX, sexChromosomeY, yHeight});

    }
    else
    {
	// Draw X chromatids 1 and 2
	drawChromosome(cr, xChromosomeX, sexChromosomeY, xHeight, xColor, xCentromereStart, xCentromereEnd);
	drawChromosome(cr, xChromosomeX + chromosomeWidth + chromatidGap,
		       sexChromosomeY, xHeight, xColor, xCentromereStart, xCentromereEnd);
	// Store X chromatids 1 and 2 geometries for DSB drawing
	chromosomeGeometry.push_back({chromosomeCount, 1, xChromosomeX, sexChromosomeY, xHeight});
 	chromosomeGeometry.push_back({chromosomeCount, 2, xChromosomeX + chromosomeWidth + chromatidGap, sexChromosomeY, xHeight});

        // Draw Y chromatids 1 and 2
        drawChromosome(cr, yChromosomeX, sexChromosomeY, yHeight, yColor, yCentromereStart, yCentromereEnd);
        drawChromosome(cr, yChromosomeX + chromosomeWidth + chromatidGap, 
		       sexChromosomeY, yHeight, yColor, yCentromereStart, yCentromereEnd);

	// Store Y chromatids 1 and 2 geometries for DSB drawing
	chromosomeGeometry.push_back({chromosomeCount - 1, 1, yChromosomeX, sexChromosomeY, yHeight});
	chromosomeGeometry.push_back({chromosomeCount - 1, 2, yChromosomeX + chromosomeWidth + chromatidGap, sexChromosomeY, yHeight});

    }


    // --------------------------------------------------
    // Draw double-strand break markers
    // --------------------------------------------------

    for (const auto& dsb : dsbVec)					// Loop through the double-strand break vector containing exposure data (damages per chromosome per exposure)
    {
        const int chromosomeID = dsb.chromosomeID.chromosomeNumber;	// Chromosome ID 1-46 for humans

        const int chromatidID = dsb.chromosomeID.chromatidNumber;	// Can be either 1 for unduplicated chromosomes, or 1 or 2 for duplicated chromosomes.

        for (const auto& geometry : chromosomeGeometry)			// Get the associated pixel coordinates for the given chromosome ID and chromatid to draw the double-strand breaks.
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

            drawDamageMarker(cr, markerX, markerY);

            break;
        }
    }

    // X label
    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);

    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);

    cairo_set_font_size(cr, 16.0);

    cairo_move_to(cr, xGroupCenterX - 4.0, sexChromosomeY + xHeight + 25.0);

    cairo_show_text(cr, "X");

    // Y label
    cairo_move_to(cr, yGroupCenterX - 4.0, sexChromosomeY + yHeight + 25.0);

    cairo_show_text(cr, "Y");

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


// ----------------------------------------------------------------------------------- //
// --- Functions to draw individual chromosome shapes and choose chromosome colors --- //
// ----------------------------------------------------------------------------------- //
void Karyogram::drawChromosome(
    cairo_t* cr,
    double x,						// Pixel X coordinate
    double y,						// Pixel Y coordinate
    double height,					// Scaled chromosome pixel height
    RGB color,						// Associated chromosome color
    double centromereStart,				// Centromere start location as fraction of the chromosome length from p arm to q arm
    double centromereEnd)				// Centromere end location as fraction of the chromosome length from p arm to q arm.
{

    // --------------------------------------
    // Centromere position and dimensions
    // --------------------------------------
    const double width = 14.0;
    const double capRadius = width / 2.0;
    const double centromereCenter = (centromereStart + centromereEnd) / 2.0;
    const double centromereY = y + (height * centromereCenter);
    const double constrictionHeight = 8.0;
    const double constrictionAmount = 1.5;
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


RGB Karyogram::generateChromosomeColor(int chromosomeNumber)
{

    // Golden-ratio color hue spacing gives better visual separation between consecutive chromosome colors.
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

    // Black double-strand break marker
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

// ---------------------------------------------------------------------------------------- //
// --------------- Helper functions to convert exposure data to karyogram ----------------- //
// ---------------------------------------------------------------------------------------- //

// Function to return the location in base pairs for a double strand break associated with a given chromosome Number and exposure.
std::vector<DSBlocation> Karyogram::getDoubleStrandBreaks(const std::vector<Exposure>& exposures)
{
    std::vector<DSBlocation> dsbVec;					// Initialize the double-strand break vector to store base pair locations of DSBs per chromosome per exposure.

    for (const auto& exposure : exposures)				// Loop through all exposures
    {
	for (const auto& damage : exposure.damages)			// Loop through the damages in each exposure
	{
            if (damage.damageType.presenceOfDoubleStrandBreaks == 1)	// If double strand break is present, store chromosome ID, base pair location to dsbVec
            {
		DSBlocation dsb;					// Instantiate the DSBlocation object dsb for easy access to DSB location and chromosome ID. 
		dsb.chromosomeID = damage.chromosomeID;
		dsb.chromosomePosition = damage.chromosomePosition;
                dsbVec.push_back(dsb);					// Store in dsbVec the following; {'chromosomeID', chromosomePosition}
            }
	}
    }
    return dsbVec;
}

// Function to convert the damage location stored in the SDD damage block in base pairs to the fractional length of the chromosome. 
double Karyogram::getDamageFraction(const DSBlocation& dsb, double chromosomeSize)
{
    double position = dsb.chromosomePosition.position;

    if (dsb.chromosomePosition.isFractional)			// Check if SDD data field 4 is already in a fractional field format
    {
        return position;					// If so, return the fractional position, no conversion needed.
    }

    if (chromosomeSize <= 0.0)
    {
        return 0.0;
    }

    double chromosomeSizeBP = chromosomeSize * 1000000.0;	// p arm = top, q arm = bottom, chromosome sizes stored in Mbp

    return position / chromosomeSizeBP;				// If SDD data field 4 is not fractional, return scaled position in fractional format. 
}
