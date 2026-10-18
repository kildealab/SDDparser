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
    std::string filename;                       // Name used for downloads from the viewer, e.g. "cell0_karyogram.png" (the extension is swapped to .svg for SVG downloads).
    std::vector<unsigned char> svgBytes;        // Vector image (KaryogramFormat::SVG): shown in the viewer, stays sharp at any zoom level, offered as SVG download.
    std::vector<unsigned char> pngBytes;        // Raster image (KaryogramFormat::PNG): offered as PNG download, and shown only when svgBytes is empty.
                                                // Either one may be left empty; the viewer then only offers the format it has.
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

// Opens a file in the system's default web browser (xdg-open on Linux, open on macOS, start on Windows, the Windows
// browser from WSL). Returns false if it could not, in which case the caller should tell the user where the file is.
bool openInBrowser(const std::filesystem::path& htmlPath);

} // namespace sddparser

#endif
