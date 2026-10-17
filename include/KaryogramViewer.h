#ifndef KARYOGRAM_VIEWER_H
#define KARYOGRAM_VIEWER_H

#include <string>
#include <vector>
#include <filesystem>

namespace sddparser
{

struct KaryogramImage                           // One karyogram to show in the web viewer.
{
    std::string label;                          // Shown in the cell list, e.g. "Cell 0".
    std::string filename;                       // Suggested name when the user downloads the PNG from the browser.
    std::vector<unsigned char> pngBytes;        // Encoded PNG, as returned by the in-memory generateSDD/SDRkaryogram overloads.
};


// Writes a single self-contained HTML file with every karyogram embedded in it (no server or extra files needed).
bool writeKaryogramViewer(
    const std::vector<KaryogramImage>& images,  // Karyograms to show, in order.
    const std::string& title,                   // Page title, e.g. the input file name.
    const std::filesystem::path& htmlPath       // Where to write the HTML file.
);


// Returns false when no browser can be shown on this machine, e.g. over SSH, on an HPC cluster or in a container,
// or in WSL with Windows interop switched off. Used to decide where to write the viewer and whether to open it.
bool canOpenBrowser();



// Opens a file in the system's default web browser (xdg-open on Linux, open on macOS, start on Windows).
bool openInBrowser(const std::filesystem::path& htmlPath);

} // namespace sddparser

#endif
