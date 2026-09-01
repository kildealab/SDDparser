#include <cairo/cairo.h>
#include <cmath>
#include <algorithm>
#include <map>
#include <string>
#include <iostream>

#include "Karyogram.h"



bool Karyogram::generateSDRkaryogram(					// Function to draw a karyogram of the rearranged chromosome segments on the karyogram
    const SDRmasterHeader& masterHeader,				// Access SDR master header for intact chromosome sizes
    const SDRsubHeader& subHeader, 					// Access SDR subheader for cell-by-cell damage/misrepair data
    bool humanGenome,							// Check if human genome specified
    const std::string& outputFilename)					// Output file name of the karyogram, concatenated with the start of the SDR filename.
{

    // Check for absent data in SDR headers/subHeaders
    if (subHeader.dataRecords.empty())
    {
        std::cerr << "ERROR: No SDR data records to draw for cell " << subHeader.cellID << ".\n";
        return false;
    }

    const std::vector<double>& sizes = masterHeader.intactChromosomeSizes;
    if (sizes.empty())
    {
        std::cerr << "ERROR: SDR master header has no chromosome sizes.\n";
        return false;
    }


    // Group all data records based on their original strand ID (of their first listed fragment) and thus their chromosome region on the karyogram
    std::map<int, std::vector<const SDRdataRecord*>> recordsByOriginalStrand;
    for (const SDRdataRecord& record : subHeader.dataRecords)
    {
        if (record.fragments.empty())
        {
            continue;
        }
        else
        {
            recordsByOriginalStrand[record.fragments[0].oldStrandID].push_back(&record);
        }
    }



    // Determining the chromosome layout, based on how chromosome sizes were passed
    const int chromosomeCount = static_cast<int>(sizes[0]);				// First entry in SDR intact chromosome sizes header is the number of chromosomes listed

    const ChromosomeLayout chromosomeLayout = determineChromosomeLayout(sizes);		// The order of how the chromosome sizes are passed alters the karyogram drawing code.
    const bool hasHomologs = chromosomeLayout != ChromosomeLayout::NON_HOMOLOGOUS;	// Default chromosome layout is non-homologous: [count, 1,2,3,...,N]
    const int homologousPairs = hasHomologs ? (chromosomeCount - 2) / 2 : 0;		// Number of homologous pairs excludes the X and Y sex chromosomes and then divided by two for pairs.
    const int drawableGroups = hasHomologs ? homologousPairs : chromosomeCount - 2;	// Determines the number of chromosome groups to draw in the karyogram, the X and Y chromosomes are excluded and will be drawn at the end.


    // Determine max strand length to scale chromosome sizes on karyogram
    double maxLengthMbp = 0.0;

    for (const SDRdataRecord& record : subHeader.dataRecords)
    {
        double lengthMbp = 0.0;

        for (const SDRfragment& fragment : record.fragments)
        {
            lengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);
        }

        if (lengthMbp > maxLengthMbp)
        {
            maxLengthMbp = lengthMbp;
        }
    }

    if (maxLengthMbp <= 0.0)
    {
        std::cerr << "ERROR: SDR data records have no positive-length fragments.\n";
        return false;
    }


    // --------------------------------------------------
    // Karyogram dimensions
    // --------------------------------------------------

    const int imgWidth = 1000;
    const int columns = 4;
    const double colWidth = static_cast<double>(imgWidth) / columns;
    const double rowHeight = 250.0;
    const double startY = 150.0;
    const int rows = (drawableGroups + columns - 1) / columns;
    int imgHeight;

    // Check if X and Y can fit in last row, if not, put in their own row
    const int sexChromosomeRow = (drawableGroups > 0) ? (drawableGroups - 1) / columns : 0;
    const int lastRowUsedCols = drawableGroups - sexChromosomeRow * columns;
    const bool sexChromosomesFitInLastRow = (columns - lastRowUsedCols) >= 2;
    if (sexChromosomesFitInLastRow)
    {
        imgHeight = static_cast<int>(startY + rows * rowHeight + 60.0);
    }
    else
    {
        imgHeight = static_cast<int>(startY + (rows + 1) * rowHeight + 60.0);
    }

    const double maxRenderHeight = 180.0;
    const double chromosomeWidth = 20.0;
    const double homologGap = 18.0;

    // Draw karyogram layout
    cairo_surface_t* surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, imgWidth, imgHeight);
    cairo_t* cr = cairo_create(surface);

    cairo_set_antialias(cr, CAIRO_ANTIALIAS_BEST);
    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    cairo_paint(cr);

    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12.0);


    // --------------------------------------------------
    // Title
    // --------------------------------------------------
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_move_to(cr, 50.0, 40.0);
    cairo_show_text(cr, ("SDR Chromosome Painting - Cell " + std::to_string(subHeader.cellID)).c_str());


    // --------------------------------------------------
    // Draw each numbered chromosome slot depending on 
    // chromosome sizes layout, X and Y will be drawn later
    // --------------------------------------------------
    for (int i = 0; i < drawableGroups; ++i)
    {
        int firstOldStrandID = 0;
        int secondOldStrandID = -1; // -1 = no homolog

        if (chromosomeLayout == ChromosomeLayout::SPLIT_HOMOLOGS)		// [count, 1,2,3...,1,2,3...]
        {
            firstOldStrandID = i + 1;
            secondOldStrandID = i + 1 + homologousPairs;
        }
        else if (chromosomeLayout == ChromosomeLayout::ADJACENT_HOMOLOGS)	// [count,1,1,2,2,3,3...]
        {
            firstOldStrandID = 2 * i + 1;
            secondOldStrandID = 2 * i + 2;
        }
        else									// [count,1,2,3,...]
        {
            firstOldStrandID = i + 1;
        }


        const int col = i % columns;
        const int row = i / columns;
        const double groupCenterX = (col * colWidth) + (colWidth / 2.0);
        const double posY = startY + (row * rowHeight);

        double labelHeight = maxRenderHeight;                   		// Fallback, if slot has nothing to draw.

	// Homologous chromosome sizes layout specified in SDR header
        if (hasHomologs)
        {

	    // Check if the data record we are looking at is an intact strand record, if so ignore. Only care about misrepaired records.
	    // Check by comparing if every fragments old strand ID is the same as that corresponding data record's new strand ID --> If so, skip.
            auto firstIt = recordsByOriginalStrand.find(firstOldStrandID);
            auto secondIt = recordsByOriginalStrand.find(secondOldStrandID);

            std::vector<const SDRdataRecord*> leftRecords;			// Left homolog
            std::vector<const SDRdataRecord*> rightRecords;			// Right homolog

            if (firstIt != recordsByOriginalStrand.end())
            {
                leftRecords = filterBaselineIfMutated(firstIt->second, chromosomeCount);	// If strand has been mutated, ignore that original strands intact strand data record, draw only the mutated strand, not the original and mutated strands together
            }

            if (secondIt != recordsByOriginalStrand.end())
            {
                rightRecords = filterBaselineIfMutated(secondIt->second, chromosomeCount);	// If strand has been mutated, ignore that original strands intact strand data record, draw only the mutated strand, not the original and mutated strands together
            }

	    // For long deletion mutation which are stacked beside the original chromosome position, add a small gap so they do not overlap.
            const double stackGap = 6.0;
	    // Adjust the gap of the left and right homologs and their corresponding deletions depending on the number of deletions and if there are deletions present.
            const double leftStackWidth = !leftRecords.empty() ? (leftRecords.size() * chromosomeWidth + (leftRecords.size() - 1) * stackGap) : 0.0;
            const double rightStackWidth = !rightRecords.empty() ? (rightRecords.size() * chromosomeWidth + (rightRecords.size() - 1) * stackGap) : 0.0;

            const double leftSlotCenterX = groupCenterX - leftStackWidth / 2.0 - homologGap / 2.0;
            const double rightSlotCenterX = groupCenterX + rightStackWidth / 2.0 + homologGap / 2.0;

	    // After checking for the presence of long deletions drawn beside the original chromosomes, draw them for the left and right homologs.
            if (!leftRecords.empty())
            {
                drawStackedMutations(cr, leftRecords, leftSlotCenterX, posY, chromosomeWidth, maxLengthMbp, maxRenderHeight, humanGenome, masterHeader);
            }

            if (!rightRecords.empty())
            {
                drawStackedMutations(cr, rightRecords, rightSlotCenterX, posY, chromosomeWidth, maxLengthMbp, maxRenderHeight, humanGenome, masterHeader);
            }

	    // Add a label
            if (!leftRecords.empty() || !rightRecords.empty())
            {
                labelHeight = std::max(
                    computeMaxBarHeight(leftRecords, maxLengthMbp, maxRenderHeight),
                    computeMaxBarHeight(rightRecords, maxLengthMbp, maxRenderHeight)
                );
            }

        }

	// Non-homologous chromosome sizes layout specified by the user in the SDR header (simpler, do not need to adjust gap between homologs for possible long deletions)
        else
        {
            auto it = recordsByOriginalStrand.find(firstOldStrandID);

            if (it != recordsByOriginalStrand.end())
            {
		// If strand has been mutated, ignore that original strands intact strand data record, draw only the mutated strand, not the original and mutated strands together
                const std::vector<const SDRdataRecord*> filtered = filterBaselineIfMutated(it->second, chromosomeCount); 

                if (!filtered.empty())		// If data present, draw the mutation and label it
                {
                    drawStackedMutations(cr, filtered, groupCenterX, posY, chromosomeWidth, maxLengthMbp, maxRenderHeight, humanGenome, masterHeader);
                    labelHeight = computeMaxBarHeight(filtered, maxLengthMbp, maxRenderHeight);

                }
            }
        }

        // Group label (chromosome number).
        cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        cairo_set_font_size(cr, 16.0);
        std::string label = "Chr " + std::to_string(i + 1);
        cairo_move_to(cr, groupCenterX - 20.0, posY + labelHeight + 25.0);
        cairo_show_text(cr, label.c_str());

    }

    // ---------------------------------------------------
    // Drawing X and Y groups, placed at the end of last row 
    // if there's room, otherwise they're given their own row.
    // ---------------------------------------------------
    if (chromosomeCount >= 2)
    {
	// Determine positions of X and Y chromosomes depeding on the number of chromosomes specified.
        double sexChromPosY = startY + rows * rowHeight;
        int yChromCol;
        int xChromCol;

        if (sexChromosomesFitInLastRow)					// Put sex chromosomes in last two positions of last row if they fit
        {
            sexChromPosY = startY + sexChromosomeRow * rowHeight;
            yChromCol = lastRowUsedCols;
            xChromCol = lastRowUsedCols + 1;
        }
        else								// Put sex chromosomes in first two positions of a new last row if they don't fit
        {
            sexChromPosY = startY + rows * rowHeight;
            yChromCol = 0;
            xChromCol = 1;
        }

        int yOldStrandID = chromosomeCount - 1;
        int xOldStrandID = chromosomeCount;

        double yCenterPosX = (yChromCol * colWidth) + (colWidth / 2.0);
        double xCenterPosX = (xChromCol * colWidth) + (colWidth / 2.0);

	// Ignore intact strand data records for mutated strand drawing
        auto yIt = recordsByOriginalStrand.find(yOldStrandID);
        auto xIt = recordsByOriginalStrand.find(xOldStrandID);
        if (yIt != recordsByOriginalStrand.end())
        {
	    // If strand has been mutated, ignore that original strands intact strand data record, draw only the mutated strand, not the original and mutated strands together
            const std::vector<const SDRdataRecord*> yRecords = filterBaselineIfMutated(yIt->second, chromosomeCount);

            double yLabelHeight = maxRenderHeight;
            if (!yRecords.empty())
            {
                drawStackedMutations(cr, yRecords, yCenterPosX, sexChromPosY, chromosomeWidth, maxLengthMbp, maxRenderHeight, humanGenome, masterHeader);
                yLabelHeight = computeMaxBarHeight(yRecords, maxLengthMbp, maxRenderHeight);
            }

            cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
            cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
            cairo_set_font_size(cr, 16.0);
            cairo_move_to(cr, yCenterPosX - 10.0, sexChromPosY + yLabelHeight + 25.0);
            cairo_show_text(cr, "Y");
        }

        if (xIt != recordsByOriginalStrand.end())
        {
	    // If strand has been mutated, ignore that original strands intact strand data record, draw only the mutated strand, not the original and mutated strands together
            const std::vector<const SDRdataRecord*> xRecords = filterBaselineIfMutated(xIt->second, chromosomeCount);

            double xLabelHeight = maxRenderHeight;

            if (!xRecords.empty())
            {
                drawStackedMutations(cr, xRecords, xCenterPosX, sexChromPosY, chromosomeWidth, maxLengthMbp, maxRenderHeight, humanGenome, masterHeader);
                xLabelHeight = computeMaxBarHeight(xRecords, maxLengthMbp, maxRenderHeight);
            }

            cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
            cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
            cairo_set_font_size(cr, 16.0);
            cairo_move_to(cr, xCenterPosX - 10.0, sexChromPosY + xLabelHeight + 25.0);
            cairo_show_text(cr, "X");
        }
    }

    cairo_surface_write_to_png(surface, outputFilename.c_str());
    cairo_destroy(cr);
    cairo_surface_destroy(surface);

    return true;
}









void Karyogram::drawPaintedChromosome(			// Function to draw the chromosomes on the karyogram , including the new fragments the chromosome is composed of
    cairo_t* cr,
    double x, 						// X position of the painted segment
    double y, 						// Y position of the painted segment
    double height, 					// Height of the chromosome
    double width, 					// Width of the chromosome
    const std::vector<PaintedSegment>& segments)	// Painted segment vector to add the aberrant strands on the original chromosomes to depict strcutural variations
{
    if (segments.empty())
    {
        return;
    }

    const double capRadius = width / 2.0;



    // --------------------------------------------------
    // Determine centromere locations and draw them
    // --------------------------------------------------
    struct CentromereSpan				// Local struct only available for the SDR chromosome drawings
    {
        double centromereTopY;				// Centromere top height in pixels
        double centromereBottomY;			// Centromere bottom height in pixels
    };

    std::vector<CentromereSpan> centromereSpans;
    for (const PaintedSegment& segment : segments)
    {
        if (segment.hasCentromere)			// Check if fragments/strands have a centromere
        {
            CentromereSpan span{};			// Initialize the CentromereSpan vector span
            span.centromereTopY = y + height * segment.centromereStartFraction;			// Top of centromere in pixels
            span.centromereBottomY = y + height * segment.centromereEndFraction;		// Bottom of centromere in pixels
            centromereSpans.push_back(span);		// Store this centromere's geometry
        }
    }

    // Sort centromere spans in order of height (top before bottom)
    std::sort(centromereSpans.begin(), centromereSpans.end(), [](const CentromereSpan& a, const CentromereSpan& b)
    {
        return a.centromereTopY < b.centromereTopY;
    });

    // Draw constriction of the chromosome near the centromere center.
    const double constrictionHeight = 7.0;
    const double constrictionAmount = 2.0;


    // -----------------------------------------
    // Trace the outer capsule outline - with a 
    // constriction notch at each centromere position
    // along the multiple fragments making up a mutated
    // chromosome
    // -----------------------------------------
    const auto traceCapsulePath = [&]()
    {
        cairo_new_path(cr);

        cairo_arc(cr, x + capRadius, y + capRadius, capRadius, M_PI, 2 * M_PI);

        for (const CentromereSpan& span : centromereSpans)
        {
            const double centromereY = (span.centromereTopY + span.centromereBottomY) / 2.0;
            const double constrictionTop = centromereY - constrictionHeight / 2.0;
            const double constrictionBottom = centromereY + constrictionHeight / 2.0;

            cairo_line_to(cr, x + width, constrictionTop);
            cairo_line_to(cr, x + width - constrictionAmount, centromereY);
            cairo_line_to(cr, x + width, constrictionBottom);
        }

        cairo_line_to(cr, x + width, y + height - capRadius);

        cairo_arc(cr, x + capRadius, y + height - capRadius, capRadius, 0, M_PI);

        for (auto it = centromereSpans.rbegin(); it != centromereSpans.rend(); ++it)
        {
            const double centromereY = (it->centromereTopY + it->centromereBottomY) / 2.0;
            const double constrictionTop = centromereY - constrictionHeight / 2.0;
            const double constrictionBottom = centromereY + constrictionHeight / 2.0;

            cairo_line_to(cr, x, constrictionBottom);
            cairo_line_to(cr, x + constrictionAmount, centromereY);
            cairo_line_to(cr, x, constrictionTop);
        }

        cairo_close_path(cr);
    };

    traceCapsulePath(); // build it once, for the clip


    // --------------------------------------
    // Fill each painted segment, clipped to the capsule shape
    // --------------------------------------
    cairo_save(cr);
    cairo_clip_preserve(cr);

    for (const PaintedSegment& segment : segments)
    {
        const double segmentTop = y + height * segment.startFraction;
        const double segmentBottom = y + height * segment.endFraction;
        const double minSegmentHeightForChevron = 8.0;

        cairo_set_source_rgb(cr, segment.color.r, segment.color.g, segment.color.b);
        cairo_rectangle(cr, x, segmentTop, width, segmentBottom - segmentTop);
        cairo_fill(cr);

        if (segment.isReversed && (segmentBottom - segmentTop) >= minSegmentHeightForChevron)		// Check for balanced inversion
        {
            const double chevronSize = width * 0.5;
            const double chevronCenterY = (segmentTop + segmentBottom) / 2.0;

            drawInversionChevron(cr, x + width / 2.0, chevronCenterY, chevronSize);			// Inversion marker to depict balanced inversions
        }

    }

    cairo_restore(cr);


    // --------------------------------------
    // Segment boundaries and reversed segment markers
    // --------------------------------------
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_set_line_width(cr, 1.0);

    for (std::size_t i = 0; i + 1 < segments.size(); i++)
    {
        const double boundaryY = y + height * segments[i].endFraction;

        cairo_move_to(cr, x, boundaryY);
        cairo_line_to(cr, x + width, boundaryY);
        cairo_stroke(cr);
    }

    cairo_set_line_width(cr, 2.0);

    for (const PaintedSegment& segment : segments)
    {
        if (!segment.isReversed)				// For non-balanced inversion mutations
        {
            continue;
        }

        const double segmentTop = y + height * segment.startFraction;		// Determine the extent of the new fragment to be drawn
        const double segmentBottom = y + height * segment.endFraction;

        cairo_move_to(cr, x, segmentTop);
        cairo_line_to(cr, x + width, segmentTop);
        cairo_stroke(cr);

        cairo_move_to(cr, x, segmentBottom);
        cairo_line_to(cr, x + width, segmentBottom);
        cairo_stroke(cr);

    }


    // --------------------------------------
    // Chromosome outline
    // --------------------------------------

    traceCapsulePath(); // rebuild - the fills above consumed the original path

    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
    cairo_set_line_width(cr, 1.5);
    cairo_stroke(cr);


    // --------------------------------------
    // Draw centromere(s) per strand
    // --------------------------------------
    const double ellipseWidth = width + 2.0;

    for (const CentromereSpan& span : centromereSpans)
    {
        const double centromereEllipseCenterY = (span.centromereTopY + span.centromereBottomY) / 2.0;
        const double ellipseHeight = std::max(5.0, span.centromereBottomY - span.centromereTopY);

        cairo_save(cr);
        cairo_translate(cr, x + width / 2.0, centromereEllipseCenterY);
        cairo_scale(cr, ellipseWidth / 2.0, ellipseHeight / 2.0);
        cairo_arc(cr, 0.0, 0.0, 1.0, 0.0, 2.0 * M_PI);
        cairo_restore(cr);

        cairo_set_source_rgb(cr, 0.6, 0.6, 0.6);
        cairo_fill(cr);
    }
}









// Helper to draw rearranged fragments on their original strands, stacked on top or beside in the same chromosome slot.
void Karyogram::drawStackedMutations(
    cairo_t* cr,
    const std::vector<const SDRdataRecord*>& records, 					// Access SDR data records to determine fragment origins and sizes
    double slotCenterX, 								// Determine the slot center to group the intact chromosomes and their associated mutations
    double posY, 									// Generic Y position variable to determine Y height of the chromosome group
    double chromosomeWidth, 								// Chromosome width shout match fragment width
    double maxLengthMbp, 								// Max length to normalize the heights of all the chromosomes drawn
    double maxRenderHeight, 								// Put a cap on the height of the chromosome heights
    bool humanGenome, 									// Check if user specified a human genome.
    const SDRmasterHeader& masterHeader)						// Check SDR master header for intact chromosome sizes layout
{
    if (records.empty())
    {
        return;
    }

    // Determine width of karyogram
    const double stackGap = 6.0;
    const std::size_t count = records.size();
    double totalWidth = count * chromosomeWidth;

    if (count > 0)
    {
        totalWidth += (count - 1) * stackGap;
    }

    // Determine where to draw the abberrant strands and how large to draw them
    double barX = slotCenterX - totalWidth / 2.0;					// Aberrant strand X position
    for (const SDRdataRecord* record : records)
    {
        double lengthMbp = 0.0;								// Convert Mbp lengths to pixel sizes
        for (const SDRfragment& fragment : record->fragments)
        {
            lengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);// Determine new rearranged strand length
        }

        double barHeight = computeSDRbarHeight(lengthMbp, maxLengthMbp, maxRenderHeight); // Determine the height of the aberrant strand
        const std::vector<PaintedSegment> segments = buildPaintedSegments(*record, humanGenome, masterHeader); // Store all the segments associated with a given data record to be built and drawn

        drawPaintedChromosome(cr, barX, posY, barHeight, chromosomeWidth, segments);	// Draw the chromosome with the rearranged fragments

        // Label new segment with new Strand ID
        cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 10.0);
        cairo_move_to(cr, barX, posY - 6.0);
        cairo_show_text(cr, std::to_string(record->newStrandID).c_str());

        barX += chromosomeWidth + stackGap;						// Adjust the position of the aberrant strand fragment x position
    }
}









// Function to build the painted segment for one new strand after it has been rearranged in the SDR file
std::vector<PaintedSegment> Karyogram::buildPaintedSegments(
    const SDRdataRecord& record, 							// Access the data record for the fragment start and end locations
    bool humanGenome, 									// Check if human genome is passed to adjust centromeres on fragments
    const SDRmasterHeader& masterHeader)						// Access SDR master header for chromosome sizes
{

    std::vector<PaintedSegment> segments;						// Initializing segments vector to hold PaintedSegment objects

    // Determining new strand length by summing the lengths of the fragments that make it up
    double totalLengthMbp = 0.0;
    for (const SDRfragment& fragment : record.fragments)
    {
        totalLengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);
    }
    if (totalLengthMbp <= 0)
    {
        return segments;
    }

    // Temporary variable to hold locations of the segment breaks
    double runningPositionMbp = 0.0;

    for (const SDRfragment& fragment : record.fragments)
    {
        const double fragmentLengthMbp = std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);	// Data field 3 fragment length is the difference in fragment start and end locations in Mbp

        PaintedSegment segment{};							// Instantiate the PaintedSegment object segment{}
        segment.startFraction = runningPositionMbp / totalLengthMbp;			// Find start location of fragment as a fraction of the total chromosome length
        segment.endFraction = (runningPositionMbp + fragmentLengthMbp) / totalLengthMbp;// Find end location of fragment as a fraction of the total chromosome length
        segment.color = getColorForOriginalStrand(fragment.oldStrandID, masterHeader);	// Obtain color of the original chromosome from which this fragment originated
        segment.isReversed = fragment.oldStartPosition > fragment.oldEndPosition;	// If the start location is > the end location of the fragment, an inversion occurred.
        segment.hasCentromere = false;							// By default, segments/fragments do not contain centromeres, unless hasCentromere = true.

        if (fragment.hasCentromere && fragmentLengthMbp > 0.0)				// If fragment does have a centromere
        {
            double centromereStartBP = 0.0;						// Determine centromere range
            double centromereEndBP = 0.0;

            if (getCentromereForOriginalStrand(fragment.oldStrandID, humanGenome, masterHeader, centromereStartBP, centromereEndBP)) // If original strand has centromere
            {
                const double centromereStartMbp = centromereStartBP / 1000000.0;
                const double centromereEndMbp = centromereEndBP / 1000000.0;

		// Ensure the centromere range stays within the span of a given fragment
                const double fragMinMbp = std::min(fragment.oldStartPosition, fragment.oldEndPosition);
                const double fragMaxMbp = std::max(fragment.oldStartPosition, fragment.oldEndPosition);
                // Clamp the centromere span to the portion actually contained within this fragment.
                const double overlapStartMbp = std::max(centromereStartMbp, fragMinMbp);
                const double overlapEndMbp = std::min(centromereEndMbp, fragMaxMbp);

                double localFractionStart;
                double localFractionEnd;

                if (overlapEndMbp > overlapStartMbp)
                {
                    if (!segment.isReversed)						// No balanced inversion
                    {
			// Centromere range on the local fragment as a fraction of the fragment length
                        localFractionStart = (overlapStartMbp - fragMinMbp) / fragmentLengthMbp;
                        localFractionEnd = (overlapEndMbp - fragMinMbp) / fragmentLengthMbp;
                    }
                    else
                    {
                        // Reversed fragment - genomic direction is flipped relative to
                        // the visual (top-to-bottom) direction of the drawn segment.
                        localFractionStart = 1.0 - (overlapEndMbp - fragMinMbp) / fragmentLengthMbp;
                        localFractionEnd = 1.0 - (overlapStartMbp - fragMinMbp) / fragmentLengthMbp;
                    }
                }
                else
                {
                    // hasCentromere is set but our lookup doesn't actually overlap this
                    // fragment's range (e.g. table/data mismatch) - fall back to a
                    // point at the fragment's midpoint rather than drawing nothing.
                    localFractionStart = 0.5;
                    localFractionEnd = 0.5;
                }

		// Adjust the centromere location relative to the fragment length and associate it to the entire mutated chromosome length
		// To be converted to pixels and drawn on the karyogram
                segment.hasCentromere = true;
                segment.centromereStartFraction = segment.startFraction + localFractionStart * (segment.endFraction - segment.startFraction);
                segment.centromereEndFraction = segment.startFraction + localFractionEnd * (segment.endFraction - segment.startFraction);
            }
        }

        segments.push_back(segment);					// Store this segment in the segments vector to be drawn in the karyogram
        runningPositionMbp += fragmentLengthMbp;			// Move on to the next fragment in a given data record
    }

    return segments;							// Return geometrical information of the mutated segments/fragments to be drawn on the karyogram.
}










void Karyogram::drawInversionChevron(					// Function to draw a marker to denote a balanced inversion occurred, at the center of the inverted fragment
    cairo_t* cr, 
    double centerX, 							// X position of the marker in pixels
    double centerY, 							// Y position of the marker in pixels
    double size)							// Size of the chevron adjusted for the size of the balanced inversion
{
    cairo_save(cr);

    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0); 				// White, for contrast against any segment color
    cairo_set_line_width(cr, 1.5);
    cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
    cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);

    cairo_move_to(cr, centerX - size / 2.0, centerY + size / 2.0);
    cairo_line_to(cr, centerX, centerY - size / 2.0);
    cairo_line_to(cr, centerX + size / 2.0, centerY + size / 2.0);

    cairo_stroke(cr);

    cairo_restore(cr);

}









// The same height-scaling formula drawStackedMutations uses per bar -
// pulled out so label positioning can reuse it without duplicating
// the clamp logic.
double Karyogram::computeSDRbarHeight(
    double lengthMbp, 								// Length of the bar fragment
    double maxLengthMbp, 							// Max length of all chromosomes for normalizing chromosome heights
    double maxRenderHeight)							// Max render height of the entire karyogram image
{
    double barHeight = maxRenderHeight * (lengthMbp / maxLengthMbp);		// Calculate segment height and on karyogram in pixels

    if (barHeight < 30.0)							// Prevent segment from being too small
    {
        barHeight = 30.0;
    }

    return barHeight;
}









// Tallest bar among a slot's (possibly stacked) records - used to
// position that slot's label just below whatever actually got drawn.
double Karyogram::computeMaxBarHeight(
    const std::vector<const SDRdataRecord*>& records, 				// Access the fragment sizes
    double maxLengthMbp, 							// Max length in Mbp of the chromosomes for normalization
    double maxRenderHeight)							// Max height of the entire karyogram for scaling
{

    double tallest = 0.0;							// Find tallest portion of a drawable group, including intact and aberrant strands

    for (const SDRdataRecord* record : records)					// Search all data records
    {
        double lengthMbp = 0.0;

        for (const SDRfragment& fragment : record->fragments)
        {
            lengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);
        }

        const double height = computeSDRbarHeight(lengthMbp, maxLengthMbp, maxRenderHeight);

        if (height > tallest)
        {
            tallest = height;
        }
    }

    return tallest;
}









// If a slot contains any genuine mutation-derived record (newStrandID
// >= numOriginalStrands), drop any leftover baseline/unmutated record
// for that same original strand - only the resulting mutated strand(s)
// should be drawn, not the untouched original alongside them.
std::vector<const SDRdataRecord*> Karyogram::filterBaselineIfMutated(
    const std::vector<const SDRdataRecord*>& records,					// Access SDR data records to assess which records contain mutated strands and which are intact
    int numOriginalStrands)								// Use to assess when the mutated SDR data records begin
{

    bool hasMutation = false;

    for (const SDRdataRecord* record : records)
    {
        if (record->newStrandID > numOriginalStrands)					// New strand IDs > 46 for human chromosomes are aberrant/mutated strands
        {
            hasMutation = true;
            break;
        }
    }

    if (!hasMutation)									// If no mutated strands found, return the SDR records of the intact strands
    {
        return records;
    }

    std::vector<const SDRdataRecord*> filtered;						// Vector 'filtered' to store unique data records as being either completely intact, or partially mutated, not both intact and mutated.

    for (const SDRdataRecord* record : records)
    {
        if (record->newStrandID <= numOriginalStrands)					// If new strand ID > 46 for human chromosomes, this is a mutated/aberrant fragment strand
        {
            continue; 									// drop the baseline (intact) data entry now that a real mutation exists for this strand
        }

        filtered.push_back(record);
    }

    return filtered;
}
