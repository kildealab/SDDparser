#include <iostream>
#include <fstream>
#include <filesystem>
#include <functional>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include <cstdio>

#if defined(_WIN32)
#include <io.h>                                                 // _isatty, _fileno
#else
#include <unistd.h>                                             // isatty, fileno
#endif

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
            << "  ./SDDparser <SDD file> --karyogram human|other [--save-png] [--save-svg] [--save-zip]\n"
            << "  ./SDDparser -sdr <SDR file>\n"
            << "  ./SDDparser -sdr <SDR file> --karyogram human|other [--save-png] [--save-svg] [--save-zip]\n\n"
            << "Karyograms open in a web viewer in your browser. Add --save-png and/or\n"
            << "--save-svg to write them as PNG and/or SVG files next to the input file instead.\n"
            << "Add --save-zip to put all karyograms in one zip file next to the input file instead\n"
            << "of one file per cell (PNG by default; add --save-svg for SVG, or both flags for both).\n\n"
            << "Examples:\n"
            << "  ./SDDparser exampleSDD.txt\n"
            << "  ./SDDparser exampleSDD.txt --karyogram human\n"
            << "  ./SDDparser -sdr exampleSDR.txt\n"
            << "  ./SDDparser -sdr exampleSDR.txt --karyogram human\n"
            << "  ./SDDparser -sdr exampleSDR.txt --karyogram human --save-png\n"
            << "  ./SDDparser -sdr exampleSDR.txt --karyogram human --save-png --save-svg\n"
            << "  ./SDDparser -sdr exampleSDR.txt --karyogram human --save-zip\n";
    }

    // Parses an optional "--karyogram human|other [--save-png] [--save-svg] [--save-zip]" option starting at
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
	bool& saveSVG,						// Check if the user wants SVG files instead of the web viewer
	bool& saveZip)						// Check if the user wants all karyograms in one zip file instead of the web viewer
    {
	drawKaryogram = false;
        genomeType.clear();
        savePNG = false;
        saveSVG = false;
        saveZip = false;

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

        for (int i = startIndex + 2; i < argc; i++)			// --save-png, --save-svg and --save-zip may be given in any order, or together
        {
            const std::string option = argv[i];
            bool& saveFlag = (option == "--save-png") ? savePNG : (option == "--save-svg") ? saveSVG : saveZip;

            if (option != "--save-png" && option != "--save-svg" && option != "--save-zip")
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
                      << "Copy the file above to your own computer and open it there, or rerun with --save-png, --save-svg or --save-zip to get image files.\n";
        }
        else if (!openInBrowser(viewerPath))
        {
            std::cout << "Could not open a browser automatically, open the file above in your browser.\n";
        }

        return true;
    }


    // Shows "Processing cell 18/100 (Cell ID 17)..." while the cells are drawn. In a terminal the line is rewritten in place, so a file with
    // 100 cells does not fill the screen; when the output goes to a file or a pipe (e.g. a cluster job log), each cell gets its own line instead.
    class ProgressLine
    {
    public:
        ProgressLine()
#if defined(_WIN32)
            : inTerminal(_isatty(_fileno(stdout)) != 0)
#else
            : inTerminal(isatty(fileno(stdout)) != 0)
#endif
        {
        }

        void show(std::size_t position, std::size_t total, int cellID)
        {
            const std::string text = "Processing cell " + std::to_string(position) + "/" + std::to_string(total) + " (Cell ID " + std::to_string(cellID) + ")...";

            if (inTerminal)
            {
                std::cout << "\r" << text << std::string(lastLength > text.size() ? lastLength - text.size() : 0, ' ') << std::flush;
                lastLength = text.size();
            }
            else
            {
                std::cout << text << "\n";
            }
        }

        // Ends the line being rewritten, so the next message (or an error) starts on its own line.
        void finish()
        {
            if (inTerminal && lastLength > 0)
            {
                std::cout << "\n";
                lastLength = 0;
            }
        }

    private:
        bool inTerminal;
        std::size_t lastLength = 0;
    };


    // Karyograms of one input file, collected cell by cell, and how they are delivered: in the web viewer (default), as one PNG/SVG file per cell
    // (--save-png/--save-svg), or together in one zip file (--save-zip).
    struct KaryogramOutput
    {
        bool savePNG = false;
        bool saveSVG = false;
        bool saveZip = false;

        std::vector<KaryogramImage> viewerImages;						// Web viewer
        std::vector<std::pair<std::string, std::vector<unsigned char>>> zipEntries;	// --save-zip: {file name, contents}
        int pngFilesWritten = 0;
        int svgFilesWritten = 0;

        bool saveFiles() const { return savePNG || saveSVG || saveZip; }			// Any save option replaces the web viewer
        bool zipPNG() const { return saveZip && (savePNG || !saveSVG); }		// --save-zip holds PNGs unless only --save-svg was asked for
        bool zipSVG() const { return saveZip && saveSVG; }
    };


    bool writeBytesToFile(const std::filesystem::path& path, const std::vector<unsigned char>& bytes)
    {
        std::ofstream file(path, std::ios::binary);
        file.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
        return static_cast<bool>(file);
    }


    // Saves one cell's karyogram as asked: PNG/SVG files next to the input file, or entries for the zip file. generate draws the cell into the
    // given bytes in the given format. Returns false (after printing what failed) if a karyogram could not be drawn or written.
    bool saveCellKaryogram(
	KaryogramOutput& output,
	const std::filesystem::path& pngPath,						// e.g. <input>_cell17_karyogram_other.png; the SVG uses the same name
	const std::function<bool(std::vector<unsigned char>&, KaryogramFormat)>& generate)
    {
        const bool wantPNG = output.saveZip ? output.zipPNG() : output.savePNG;
        const bool wantSVG = output.saveZip ? output.zipSVG() : output.saveSVG;

        for (const KaryogramFormat format : {KaryogramFormat::PNG, KaryogramFormat::SVG})
        {
            if ((format == KaryogramFormat::PNG && !wantPNG) || (format == KaryogramFormat::SVG && !wantSVG))
            {
                continue;
            }

            const std::filesystem::path path = (format == KaryogramFormat::PNG) ? pngPath : std::filesystem::path(pngPath).replace_extension(".svg");
            std::vector<unsigned char> bytes;

            if (!generate(bytes, format))
            {
                return false;
            }

            if (output.saveZip)
            {
                output.zipEntries.emplace_back(path.filename().string(), std::move(bytes));
            }
            else if (!writeBytesToFile(path, bytes))
            {
                std::cerr << "ERROR: Failed to write karyogram file " << path << "\n";
                return false;
            }
            else
            {
                ++(format == KaryogramFormat::PNG ? output.pngFilesWritten : output.svgFilesWritten);
            }
        }

        return true;
    }


    // Delivers the collected karyograms: opens the web viewer, writes the zip file, or reports where the PNG/SVG files went.
    bool deliverKaryograms(KaryogramOutput& output, const std::filesystem::path& inputPath, const std::string& genomeType)
    {
        if (!output.saveFiles())
        {
            return output.viewerImages.empty() || showInWebViewer(output.viewerImages, inputPath, genomeType);
        }

        if (output.saveZip)
        {
            if (output.zipEntries.empty())
            {
                return true;
            }

            const std::filesystem::path zipPath = inputPath.parent_path() / (inputPath.stem().string() + "_karyograms_" + genomeType + ".zip");

            if (!writeZipArchive(output.zipEntries, zipPath))
            {
                return false;
            }

            std::cout << "Karyograms written to: " << zipPath << " (" << output.zipEntries.size() << " files)\n";
            return true;
        }

        const std::filesystem::path folder = std::filesystem::absolute(inputPath).parent_path();

        if (output.pngFilesWritten > 0)
        {
            std::cout << "PNG karyograms written to: " << folder << " (" << output.pngFilesWritten << " files)\n";
        }

        if (output.svgFilesWritten > 0)
        {
            std::cout << "SVG karyograms written to: " << folder << " (" << output.svgFilesWritten << " files)\n";
        }

        return true;
    }


    // "17, 24, 28"
    std::string joinCellIDs(const std::vector<int>& cellIDs)
    {
        std::string text;

        for (std::size_t i = 0; i < cellIDs.size(); ++i)
        {
            text += (i > 0 ? ", " : "") + std::to_string(cellIDs[i]);
        }

        return text;
    }


    // Last lines of every run: the cells without karyograms and why, then either the success message or what went wrong.
    // Returns the program's exit code (0 = success).
    int reportFinalStatus(std::size_t cellCount, const std::vector<int>& absentDataCells, const std::vector<int>& failedCells, bool otherError)
    {
        if (!absentDataCells.empty())
        {
            std::cout << "Karyograms not generated for the following cells due to absent SDR data: " << joinCellIDs(absentDataCells) << "\n";
        }

        if (!failedCells.empty())
        {
            std::cerr << "ERROR: Karyograms could not be generated for the following cells: " << joinCellIDs(failedCells) << "\n";
        }

        if (!failedCells.empty() || otherError)
        {
            std::cerr << "Finished with errors, see the messages above.\n";
            return 1;
        }

        std::cout << "All " << cellCount << " cells processed successfully.\n";
        return 0;
    }
}










int main(int argc, char* argv[]) 					// Variables in main() brackets allow for the second command line argument to be the SDD file path or -h|--help for usage help. Third argument is an option for plotting the corresponding Karyogram if desired using "--karyogram". The fourth argument specifies either 'human' genome or 'other' genome.
{

    // ---------------------------------------
    // Check Command-line arguments
    // ---------------------------------------

    if(argc < 2)							// Must have at least two command-line arguments, the program and the SDD file path.
    {
	std::cerr << "Usage: ./SDDparser <SDD file> [--karyogram human|other [--save-png] [--save-svg] [--save-zip]]\n"
		   << "       ./SDDparser -sdr <SDR file> [--karyogram human|other [--save-png] [--save-svg] [--save-zip]]\n";
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
        KaryogramOutput output;

	if (!parseKaryogramOption(argc, argv, 3, drawKaryogram, genomeType, output.savePNG, output.saveSVG, output.saveZip))
        {
            return 1;
        }

	bool humanGenome = (genomeType == "human");

        std::cout << "SDR file processing: " << filename << "...\n";

        SDRparser parser;

        if (!parser.parseFile(filename))
        {
            std::cerr << "ERROR: Failed to load SDR file " << filename << "\n";
            return 1;
        }

	std::filesystem::path inputPath(filename);
        const std::vector<SDRsubHeader>& subHeaders = parser.getSubHeaders();

        std::vector<int> absentDataCells;				// Cells with no SDR data records: nothing to draw, not an error
        std::vector<int> failedCells;					// Cells whose karyogram could not be drawn or written: an error

	// --------------------------------------------- //
	// SDR karyogram drawing, one cell at a time
	// --------------------------------------------- //
        if (drawKaryogram)						// If user specifies --karyogram then draw
        {
	    Karyogram karyogram;					// Instantiate the Karyogram object karyogram
	    ProgressLine progress;

	    for (std::size_t i = 0; i < subHeaders.size(); ++i)		// Draw one karyogram per cell subheader present in SDR file
	    {
		const SDRsubHeader& subHeader = subHeaders[i];
		progress.show(i + 1, subHeaders.size(), subHeader.cellID);

		// A cell without data records has nothing to draw. It is listed at the end rather than reported as an error,
		// since the file itself is fine (generateSDRkaryogram() would otherwise refuse it with an error).
		if (subHeader.dataRecords.empty())
		{
		    absentDataCells.push_back(subHeader.cellID);
		    continue;
		}

		std::filesystem::path cellKaryogramPath =			// Modify name of output karyogram
            		inputPath.parent_path() /
            		(inputPath.stem().string() + "_cell" + std::to_string(subHeader.cellID) +
             		"_karyogram_" + genomeType + ".png");

		const auto generate = [&](std::vector<unsigned char>& bytes, KaryogramFormat format)
		{
		    return karyogram.generateSDRkaryogram(parser.getMasterHeader(), subHeader, humanGenome, bytes, format);
		};

		if (!output.saveFiles())					// Default: generate in memory and show in the web viewer
		{
		    KaryogramImage image;
		    image.label = "Cell " + std::to_string(subHeader.cellID);
		    image.filename = cellKaryogramPath.filename().string();

		    if (!generate(image.svgBytes, KaryogramFormat::SVG)	// SVG stays sharp when zooming in the viewer
			|| !generate(image.pngBytes, KaryogramFormat::PNG))	// PNG for the viewer's PNG download, identical to --save-png
		    {
			failedCells.push_back(subHeader.cellID);
			continue;
		    }

		    image.hotspots = buildLabelHotspots(karyogram.getLabelRegions(),					// Clickable chromosome labels
			countMutationsPerOriginalStrand(subHeader, parser.getMasterHeader()));
		    image.popupNote = "Events involving several chromosomes are counted on each of them.";
		    output.viewerImages.push_back(std::move(image));
		    continue;
		}

		if (!saveCellKaryogram(output, cellKaryogramPath, generate))	// --save-png / --save-svg / --save-zip
		{
		    failedCells.push_back(subHeader.cellID);
		}
	    }

	    progress.finish();
        }
        else
        {
            std::cout << "Summarizing " << subHeaders.size() << " cells...\n";
        }

	// --------------------------------------
	// SDR summary file
	// --------------------------------------
        std::filesystem::path summaryPath = inputPath.parent_path() /
            (inputPath.stem().string() + "_summary.txt");

        bool otherError = false;

        if (!parser.writeSummary(summaryPath.string()))
        {
            std::cerr << "ERROR: Failed to create SDR summary file " << summaryPath << "\n";
            otherError = true;
        }
        else
        {
            std::cout << "Summary written to: " << summaryPath << "\n";
        }

        if (drawKaryogram && !deliverKaryograms(output, inputPath, genomeType))
        {
            otherError = true;
        }

        return reportFinalStatus(subHeaders.size(), absentDataCells, failedCells, otherError);
    }


    // --------------------------------------
    // SDD MODE
    // --------------------------------------

    std::string filename = firstArg;

    bool drawKaryogram = false;
    std::string genomeType;
    KaryogramOutput output;

    if (!parseKaryogramOption(argc, argv, 2, drawKaryogram, genomeType, output.savePNG, output.saveSVG, output.saveZip))
    {
	return 1;
    }

    bool humanGenome = (genomeType == "human");

    // -------------------------------------
    // Load SDD file
    // -------------------------------------

    std::cout << "SDD file processing: " << filename << "...\n";

    SDDparser parser; 							// Instantiate SDDparser object called parser

    if(!parser.load(filename))						// Check if SDD input file loads correctly
    {
        std::cerr << "ERROR: Failed to load SDD file " << filename << "\n";
        return 1;
    }

    std::filesystem::path inputPath(filename);				// Outputs go next to the input SDD file
    const std::vector<Exposure>& exposures = parser.getExposures();

    std::vector<int> failedCells;					// Cells whose karyogram could not be drawn or written: an error

    // ---------------------------------------------
    // SDD karyogram drawing, one cell at a time
    // ---------------------------------------------

    if (drawKaryogram)							// Check if user passed '--karyogram', meaning drawKaryogram == true.
    {
        const Header& header = parser.getHeader();			// Get header to pass the chromosome sizes and cell cycle phase to the generate Karyogram function.

        Karyogram karyogram;						// Instantiate the Karyogram object called karyogram.
        ProgressLine progress;

    	for (std::size_t i = 0; i < exposures.size(); ++i)		// Draw a karyogram for each individual cell exposure in a given SDD file
    	{
    	    const Exposure& exposure = exposures[i];
    	    const int cellID = exposure.exposureID - 1;		// Cells are numbered from 0, exposure IDs from 1
    	    progress.show(i + 1, exposures.size(), cellID);

    	    // Create a vector containing only this exposure.
    	    std::vector<Exposure> singleExposure;
    	    singleExposure.push_back(exposure);

    	    // Create a separate output filename for this cell.
    	    std::filesystem::path karyogramPath = inputPath.parent_path() /
            (inputPath.stem().string() + "_cell" +
	    std::to_string(cellID) + "_karyogram_" +
	    genomeType + ".png");

    	    const auto generate = [&](std::vector<unsigned char>& bytes, KaryogramFormat format)
    	    {
    	        return karyogram.generateSDDkaryogram(
            	    header.chromosome_sizes,
            	    header.cell_cycle_phase,
            	    header.dose_or_fluence,
            	    header.incident_particles,
            	    singleExposure,
            	    humanGenome,
            	    bytes,
            	    format);
    	    };

    	    if (!output.saveFiles())					// Default: generate in memory and show in the web viewer
    	    {
    	        KaryogramImage image;
    	        image.label = "Cell " + std::to_string(cellID);
    	        image.filename = karyogramPath.filename().string();

    	        if (!generate(image.svgBytes, KaryogramFormat::SVG)		// SVG stays sharp when zooming in the viewer
    	            || !generate(image.pngBytes, KaryogramFormat::PNG))	// PNG for the viewer's PNG download, identical to --save-png
    	        {
            	    failedCells.push_back(cellID);
            	    continue;
    	        }

    	        image.hotspots = buildLabelHotspots(karyogram.getLabelRegions(),			// Clickable chromosome labels
    	            countDamagesPerChromosome(parser, exposure));
    	        image.popupEmptyText = "No damage detected";					// SDD files record damage, not mutations
    	        output.viewerImages.push_back(std::move(image));
    	        continue;
    	    }

    	    if (!saveCellKaryogram(output, karyogramPath, generate))	// --save-png / --save-svg / --save-zip
    	    {
    	        failedCells.push_back(cellID);
    	    }
    	}

    	progress.finish();
    }
    else
    {
        std::cout << "Summarizing " << exposures.size() << " cells...\n";
    }

    // -------------------------------------
    // Create summary file
    // -------------------------------------

    std::filesystem::path summaryPath = inputPath.parent_path() / 	// Create the path to the summary file as the same path as the input SDD file
	(inputPath.stem().string() + "_summary.txt");

    bool otherError = false;
    std::ofstream summaryFile(summaryPath);				// Output to desired summary file path

    if(!summaryFile)							// Make sure summary file created successfully.
    {
        std::cerr << "ERROR: Failed to create summary file " << summaryPath << "\n";
        otherError = true;
    }
    else
    {
        parser.printSummary(summaryFile); 				// Returns Summary of SDD header fields in a summary file at the desired path
        summaryFile.close();
        std::cout << "Summary written to: " << summaryPath << "\n";	// Let user know summary file was created successfully
    }

    if (drawKaryogram && !deliverKaryograms(output, inputPath, genomeType))
    {
        otherError = true;
    }

    return reportFinalStatus(exposures.size(), {}, failedCells, otherError);
}
