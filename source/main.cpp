#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <utility>
#include <map>

#include "SDDparser.h"
#include "SDRparser.h"
#include "Karyogram.h"
#include "KaryogramViewer.h"
#include "SDRutilities.h"

using namespace sddparser;

namespace
{
    void printUsage()
    {
        std::cout
            << "SDDparser\n\n"
            << "Usage:\n"
            << "  ./SDDparser <SDD file>\n"
            << "  ./SDDparser <SDD file> --karyogram human|other [--save-png] [--save-svg]\n"
            << "  ./SDDparser -sdr <SDR file>\n"
            << "  ./SDDparser -sdr <SDR file> --karyogram human|other [--save-png] [--save-svg]\n\n"
            << "Karyograms open in a web viewer in your browser. Add --save-png and/or\n"
            << "--save-svg to write them as PNG and/or SVG files next to the input file instead.\n\n"
            << "Examples:\n"
            << "  ./SDDparser exampleSDD.txt\n"
            << "  ./SDDparser exampleSDD.txt --karyogram human\n"
            << "  ./SDDparser -sdr exampleSDR.txt\n"
            << "  ./SDDparser -sdr exampleSDR.txt --karyogram human\n"
            << "  ./SDDparser -sdr exampleSDR.txt --karyogram human --save-png\n"
            << "  ./SDDparser -sdr exampleSDR.txt --karyogram human --save-png --save-svg\n";
    }

    // Parses an optional "--karyogram human|other [--save-png] [--save-svg]" option starting at
    // argv[startIndex]. If nothing follows the preceding arguments,
    // drawKaryogram is left false and this returns true (nothing to do).
    // Prints an error and returns false on any malformed input.
    bool parseKaryogramOption(
	int argc, 						// Number of command line arguments
	char* argv[], 						// The command line argument at a particular index
	int startIndex, 					// Variable to track the different possible commands
	bool& drawKaryogram, 					// Check if draw karyogram was desired by the user
	std::string& genomeType,				// Check if genome specified was human or other
	bool& savePNG,						// Check if the user wants PNG files instead of the web viewer
	bool& saveSVG)						// Check if the user wants SVG files instead of the web viewer
    {
	drawKaryogram = false;
        genomeType.clear();
        savePNG = false;
        saveSVG = false;

        if (argc <= startIndex)
        {
            return true;
        }

        if (std::string(argv[startIndex]) != "--karyogram")
        {
            std::cerr << "ERROR: Unknown command-line option: "
                      << argv[startIndex] << '\n';
            return false;
        }

        drawKaryogram = true;

	if (argc <= startIndex + 1)
        {
            std::cerr << "ERROR: --karyogram requires "
                      << "either 'human' or 'other'.\n";
            return false;
        }

        genomeType = argv[startIndex + 1];

        if (genomeType != "human" && genomeType != "other")
        {
            std::cerr << "ERROR: Genome type must be "
                      << "'human' or 'other'.\n";
            return false;
        }

        for (int i = startIndex + 2; i < argc; i++)			// --save-png and --save-svg may be given in either order, or together
        {
            const std::string option = argv[i];
            bool& saveFlag = (option == "--save-png") ? savePNG : saveSVG;

            if (option != "--save-png" && option != "--save-svg")
            {
                std::cerr << "ERROR: Unknown command-line option: "
                          << option << '\n';
                return false;
            }

            if (saveFlag)
            {
                std::cerr << "ERROR: " << option << " was given more than once.\n";
                return false;
            }

            saveFlag = true;
        }

	return true;
    }

    // Turns the chromosome labels drawn on a karyogram into clickable areas for the web viewer, each with a summary
    // per homolog (SDR mutations or SDD damages). Only the types that are present are listed.
    std::vector<KaryogramHotspot> buildLabelHotspots(
	const std::vector<KaryogramLabelRegion>& labelRegions,	// From Karyogram::getLabelRegions(), after drawing the cell
	const std::map<int, std::vector<std::pair<std::string, int>>>& countsPerStrand)	// Strand/chromosome ID -> (type, count), nonzero only
    {
        std::vector<KaryogramHotspot> hotspots;

        for (const KaryogramLabelRegion& region : labelRegions)
        {
            KaryogramHotspot hotspot;
            const std::string name = region.label.rfind("Chr ", 0) == 0 ? region.label.substr(4) : region.label;	// "Chr 1" -> "1", "X" stays "X"
            hotspot.title = "Chromosome " + name;
            hotspot.x = region.x;
            hotspot.y = region.y;
            hotspot.width = region.width;
            hotspot.height = region.height;

            // Homologous layouts list both homologs of the pair; non-homologous layouts have one strand per label,
            // which needs no heading.
            for (size_t h = 0; h < region.oldStrandIDs.size(); h++)
            {
                KaryogramPopupSection section;

                if (region.oldStrandIDs.size() > 1)
                {
                    section.heading = "Homolog " + std::to_string(h + 1);
                }

                const auto found = countsPerStrand.find(region.oldStrandIDs[h]);

                if (found != countsPerStrand.end())
                {
                    section.counts = found->second;
                }

                hotspot.sections.push_back(std::move(section));
            }

            hotspots.push_back(std::move(hotspot));
        }

        return hotspots;
    }

    // Base damages, single-strand breaks and double-strand breaks per SDD chromosome number for one exposure, using the
    // same per-chromosome totals as the summary file. Sister chromatids of a homolog share one chromosome number, so
    // their damages are added together. Types with no damage are left out.
    std::map<int, std::vector<std::pair<std::string, int>>> countDamagesPerChromosome(
	const SDDparser& parser,
	const Exposure& exposure)
    {
        std::map<int, std::vector<std::pair<std::string, int>>> result;

        for (const auto& [chromosomeNumber, summary] : parser.summarizeChromosomeDamage(exposure))
        {
            std::vector<std::pair<std::string, int>>& counts = result[chromosomeNumber];

            if (summary.numBaseDamages > 0)        counts.emplace_back("Base damages", summary.numBaseDamages);
            if (summary.numSingleStrandBreaks > 0) counts.emplace_back("Single-strand breaks", summary.numSingleStrandBreaks);
            if (summary.numDoubleStrandBreaks > 0) counts.emplace_back("Double-strand breaks", summary.numDoubleStrandBreaks);
        }

        return result;
    }

    // Writes all karyograms of one input file into a single HTML viewer and opens it in the browser.
    bool showInWebViewer(
	const std::vector<KaryogramImage>& images,		// Karyograms generated in memory, one per cell
	const std::filesystem::path& inputPath,			// Input SDD/SDR file, used to name the viewer
	const std::string& genomeType)				// human or other, kept in the file name like the PNG outputs
    {
        const bool browserAvailable = canOpenBrowser();

        // Where to write the viewer:
        // - With a browser, use the system temp folder so the input folder is not cluttered (the user can still
        //   download PNGs from the viewer).
        // - Without a browser (SSH, HPC cluster, container), the temp folder is on a remote machine the user cannot
        //   easily browse, so write the viewer next to the input file instead, where it is easy to copy back (scp).
        std::filesystem::path viewerFolder = inputPath.parent_path();

        if (browserAvailable)
        {
            std::error_code error;
            const std::filesystem::path tempFolder = std::filesystem::temp_directory_path(error);

            if (!error)						// Keep the input folder if the system has no usable temp folder
            {
                viewerFolder = tempFolder;
            }
        }

        const std::filesystem::path viewerPath = viewerFolder /
            (inputPath.stem().string() + "_karyograms_" + genomeType + ".html");

        if (!writeKaryogramViewer(images, inputPath.filename().string(), viewerPath))
        {
            return false;
        }

        // Always print the path first, so the user can open the viewer by hand whatever happens next.
        std::cout << "Karyogram viewer written to: " << viewerPath << "\n";

        if (!browserAvailable)
        {
            std::cout << "No browser can be opened from this session (for example over SSH).\n"
                      << "Copy the file above to your own computer and open it there, or rerun with --save-png or --save-svg to get image files.\n";
        }
        else if (!openInBrowser(viewerPath))
        {
            std::cout << "Could not open a browser automatically, open the file above in your browser.\n";
        }

        return true;
    }
}










int main(int argc, char* argv[]) 					// Variables in main() brackets allow for the second command line argument to be the SDD file path or -h|--help for usage help. Third argument is an option for plotting the corresponding Karyogram if desired using "--karyogram". The fourth argument specifies either 'human' genome or 'other' genome.
{

    // ---------------------------------------
    // Check Command-line arguments
    // ---------------------------------------

    if(argc < 2)							// Must have at least two command-line arguments, the program and the SDD file path.
    {
	std::cerr << "Usage: ./SDDparser <SDD file> [--karyogram human|other [--save-png] [--save-svg]]\n"
		   << "       ./SDDparser -sdr <SDR file> [--karyogram human|other [--save-png] [--save-svg]]\n";
	return 1;
    }

    std::string firstArg = argv[1];					// If second command-line argument is -h or --help, return usage examples.
    if (firstArg == "-h" || firstArg == "--help")
    {
        printUsage();
        return 0;
    }

    // --------------------------------------
    // SDR MODE
    // --------------------------------------

    // --------------------------------------
    // SDR summary file
    // --------------------------------------
    if (firstArg == "-sdr")
    {
        if (argc < 3)
        {
            std::cerr << "ERROR: -sdr requires an SDR file path.\n";
            return 1;
        }

        std::string filename = argv[2];

        bool drawKaryogram = false;
        std::string genomeType;
        bool savePNG = false;
        bool saveSVG = false;

	if (!parseKaryogramOption(argc, argv, 3, drawKaryogram, genomeType, savePNG, saveSVG))
        {
            return 1;
        }

	bool humanGenome = (genomeType == "human");

        SDRparser parser;

        if (!parser.parseFile(filename))
        {
            std::cout << "Failed to load SDR file " << filename << "\n";
            return 1;
        }

	std::filesystem::path inputPath(filename);
        std::filesystem::path summaryPath = inputPath.parent_path() /
            (inputPath.stem().string() + "_summary.txt");

        if (!parser.writeSummary(summaryPath.string()))
        {
            std::cerr << "Failed to create SDR summary file\n";
            return 1;
        }

        std::cout << "Summary written to: " << summaryPath << "\n";

	// --------------------------------------------- //
	// SDR karyogram drawing
	// --------------------------------------------- //
        if (drawKaryogram)					// If user specifies --karyogram then draw
        {
	    Karyogram karyogram;				// Instantiate the Karyogram object karyogram

	    bool hasFailed = false;

	    const bool saveFiles = savePNG || saveSVG;				// Image files requested, so no web viewer

	    std::vector<KaryogramImage> viewerImages;				// Karyograms kept in memory for the web viewer (unless --save-png/--save-svg)

	    for (const SDRsubHeader& subHeader : parser.getSubHeaders())	// Draw one karyogram per cell subheader present in SDR file
	    {
		std::filesystem::path cellKaryogramPath =			// Modify name of output karyogram
            		inputPath.parent_path() /
            		(inputPath.stem().string() + "_cell" + std::to_string(subHeader.cellID) +
             		"_karyogram_" + genomeType + ".png");

		if (!saveFiles)							// Default: generate in memory and show in the web viewer
		{
		    KaryogramImage image;
		    image.label = "Cell " + std::to_string(subHeader.cellID);
		    image.filename = cellKaryogramPath.filename().string();

		    if (!karyogram.generateSDRkaryogram(parser.getMasterHeader(), subHeader, humanGenome, image.svgBytes, KaryogramFormat::SVG)	// SVG stays sharp when zooming in the viewer
			|| !karyogram.generateSDRkaryogram(parser.getMasterHeader(), subHeader, humanGenome, image.pngBytes))			// PNG for the viewer's PNG download, identical to --save-png
		    {
			std::cerr << "Failed to generate karyogram for cell " << subHeader.cellID << ".\n";
			hasFailed = true;
			continue;
		    }

		    image.hotspots = buildLabelHotspots(karyogram.getLabelRegions(),					// Clickable chromosome labels
			countMutationsPerOriginalStrand(subHeader, parser.getMasterHeader()));
		    image.popupNote = "Events involving several chromosomes are counted on each of them.";
		    viewerImages.push_back(std::move(image));
		    continue;
		}

		// Check if karyogram generated correctly
		if (savePNG)
		{
		    if (!karyogram.generateSDRkaryogram(parser.getMasterHeader(), subHeader, humanGenome, cellKaryogramPath.string()))
		    {
			std::cerr << "Failed to generate karyogram for cell " << subHeader.cellID << ".\n";
			hasFailed = true;
			continue;
		    }

		    std::cout << "Karyogram generated successfully: " << cellKaryogramPath << "\n";
		}

		if (saveSVG)
		{
		    const std::filesystem::path cellKaryogramSVGPath =		// Same name as the PNG, with the .svg extension
			std::filesystem::path(cellKaryogramPath).replace_extension(".svg");

		    if (!karyogram.generateSDRkaryogram(parser.getMasterHeader(), subHeader, humanGenome, cellKaryogramSVGPath.string(), KaryogramFormat::SVG))
		    {
			std::cerr << "Failed to generate SVG karyogram for cell " << subHeader.cellID << ".\n";
			hasFailed = true;
			continue;
		    }

		    std::cout << "Karyogram generated successfully: " << cellKaryogramSVGPath << "\n";
		}

	    }

	    if (!saveFiles && !viewerImages.empty() && !showInWebViewer(viewerImages, inputPath, genomeType))
	    {
		hasFailed = true;
	    }

	    if (hasFailed)
	    {
		return 1;
	    }

        }

        return 0;
    }


    // --------------------------------------
    // SDD MODE
    // --------------------------------------

    std::string filename = firstArg;

    bool drawKaryogram = false;
    std::string genomeType;
    bool savePNG = false;
    bool saveSVG = false;

    if (!parseKaryogramOption(argc, argv, 2, drawKaryogram, genomeType, savePNG, saveSVG))
    {
	return 1;
    }

    bool humanGenome = (genomeType == "human");

    // -------------------------------------
    // Load SDD file
    // -------------------------------------

    SDDparser parser; 							// Instantiate SDDparser object called parser

    if(!parser.load(filename))						// Check if SDD input file loads correctly
    {
        std::cout << "Failed to load SDD file " << filename << "\n";
        return 1;
    }

    // -------------------------------------
    // Create summary file
    // -------------------------------------

    std::filesystem::path inputPath(filename);				// Create Summary Output file to the same path as the input SDD file

    std::filesystem::path summaryPath = inputPath.parent_path() / 	// Create the path to the summary file as the same path as the input SDD file
	(inputPath.stem().string() + "_summary.txt");

    std::ofstream summaryFile(summaryPath);				// Output to desired summary file path

    if(!summaryFile)							// Make sure summary file created successfully.
    {
        std::cerr
            << "Failed to create summary file\n";

        return 1;
    }

    parser.printSummary(summaryFile); 					//Returns Summary of SDD header fields in a summary file at the desired path

    summaryFile.close();

    std::cout << "Summary written to: " << summaryPath << "\n";		// Let user know summary file was created successfully

    // ---------------------------------------------
    // SDD karyogram drawing
    // ---------------------------------------------

    if (drawKaryogram)							// Check if user passed '--karyogram', meaning drawKaryogram == true. 
    {
        const Header& header = parser.getHeader();			// Get header to pass the chromosome sizes and cell cycle phase to the generate Karyogram function.

        Karyogram karyogram;						// Instantiate the Karyogram object called karyogram.

    	bool hasFailed = false;						// default check is false is karyogram drawn correctly.

    	const bool saveFiles = savePNG || saveSVG;			// Image files requested, so no web viewer

    	std::vector<KaryogramImage> viewerImages;			// Karyograms kept in memory for the web viewer (unless --save-png/--save-svg)

    	for (const Exposure& exposure : parser.getExposures())		// Draw a karyogram for each individual cell exposure in a given SDD file
    	{
    	    // Create a vector containing only this exposure.
    	    std::vector<Exposure> singleExposure;
    	    singleExposure.push_back(exposure);

    	    // Create a separate output filename for this cell (0-indexing requires - 1 to exposure ID which is 1-indexed).
    	    std::filesystem::path karyogramPath = inputPath.parent_path() /
            (inputPath.stem().string() + "_cell" + 
	    std::to_string(exposure.exposureID - 1) + "_karyogram_" + 
	    genomeType + ".png");

    	    if (!saveFiles)							// Default: generate in memory and show in the web viewer
    	    {
    	        KaryogramImage image;
    	        image.label = "Cell " + std::to_string(exposure.exposureID - 1);
    	        image.filename = karyogramPath.filename().string();

    	        if (!karyogram.generateSDDkaryogram(
            	    header.chromosome_sizes,
            	    header.cell_cycle_phase,
            	    header.dose_or_fluence,
            	    header.incident_particles,
            	    singleExposure,
            	    humanGenome,
            	    image.svgBytes,
            	    KaryogramFormat::SVG)					// SVG stays sharp when zooming in the viewer
    	            || !karyogram.generateSDDkaryogram(				// PNG for the viewer's PNG download, identical to --save-png
            	    header.chromosome_sizes,
            	    header.cell_cycle_phase,
            	    header.dose_or_fluence,
            	    header.incident_particles,
            	    singleExposure,
            	    humanGenome,
            	    image.pngBytes))
    	        {
            	    std::cerr << "Failed to generate karyogram for cell "
                      << exposure.exposureID << ".\n";

            	    hasFailed = true;
            	    continue;
    	        }

    	        image.hotspots = buildLabelHotspots(karyogram.getLabelRegions(),			// Clickable chromosome labels
    	            countDamagesPerChromosome(parser, exposure));
    	        image.popupEmptyText = "No damage detected";					// SDD files record damage, not mutations
    	        viewerImages.push_back(std::move(image));
    	        continue;
    	    }

    	    if (savePNG)
    	    {
    	        if (!karyogram.generateSDDkaryogram(				// Check if generateSDDkaryogram correctly runs
            	    header.chromosome_sizes,
            	    header.cell_cycle_phase,
            	    header.dose_or_fluence,
            	    header.incident_particles,
            	    singleExposure,
            	    humanGenome,
            	    karyogramPath.string()))
    	        {
            	    std::cerr << "Failed to generate karyogram for cell "
                      << exposure.exposureID << ".\n";

            	    hasFailed = true;
            	    continue;
    	        }

                std::cout << "Karyogram generated successfully: "
                      << karyogramPath << "\n";
    	    }

    	    if (saveSVG)
    	    {
    	        const std::filesystem::path karyogramSVGPath =			// Same name as the PNG, with the .svg extension
    	            std::filesystem::path(karyogramPath).replace_extension(".svg");

    	        if (!karyogram.generateSDDkaryogram(
            	    header.chromosome_sizes,
            	    header.cell_cycle_phase,
            	    header.dose_or_fluence,
            	    header.incident_particles,
            	    singleExposure,
            	    humanGenome,
            	    karyogramSVGPath.string(),
            	    KaryogramFormat::SVG))
    	        {
            	    std::cerr << "Failed to generate SVG karyogram for cell "
                      << exposure.exposureID << ".\n";

            	    hasFailed = true;
            	    continue;
    	        }

                std::cout << "Karyogram generated successfully: "
                      << karyogramSVGPath << "\n";
    	    }
    	}

    	if (!saveFiles && !viewerImages.empty() && !showInWebViewer(viewerImages, inputPath, genomeType))
    	{
    	    hasFailed = true;
    	}

    	if (hasFailed)
    	{
    	    return 1;
    	}

    }
    return 0;
}

