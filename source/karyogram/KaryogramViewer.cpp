#include "KaryogramViewer.h"
#include "generalUtilities.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>




namespace sddparser
{




namespace
{
    // Encode PNG bytes as base64 so they can be embedded in the HTML as data: URLs.
    std::string base64Encode(const std::vector<unsigned char>& bytes)
    {
        static const char table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

        std::string encoded;
        encoded.reserve(((bytes.size() + 2) / 3) * 4);

        size_t i = 0;
        for (; i + 2 < bytes.size(); i += 3)                                    // Every 3 input bytes become 4 output characters.
        {
            const unsigned int chunk = (bytes[i] << 16) | (bytes[i + 1] << 8) | bytes[i + 2];
            encoded += table[(chunk >> 18) & 0x3F];
            encoded += table[(chunk >> 12) & 0x3F];
            encoded += table[(chunk >> 6) & 0x3F];
            encoded += table[chunk & 0x3F];
        }

        const size_t remaining = bytes.size() - i;                             // Pad the last 1 or 2 bytes with '='.
        if (remaining == 1)
        {
            const unsigned int chunk = bytes[i] << 16;
            encoded += table[(chunk >> 18) & 0x3F];
            encoded += table[(chunk >> 12) & 0x3F];
            encoded += "==";
        }
        else if (remaining == 2)
        {
            const unsigned int chunk = (bytes[i] << 16) | (bytes[i + 1] << 8);
            encoded += table[(chunk >> 18) & 0x3F];
            encoded += table[(chunk >> 12) & 0x3F];
            encoded += table[(chunk >> 6) & 0x3F];
            encoded += '=';
        }

        return encoded;
    }

    // Escape text placed inside the HTML (file names can contain characters such as & or <).
    std::string htmlEscape(const std::string& text)
    {
        std::string escaped;
        for (char c : text)
        {
            switch (c)
            {
                case '&':  escaped += "&amp;";  break;
                case '<':  escaped += "&lt;";   break;
                case '>':  escaped += "&gt;";   break;
                case '"':  escaped += "&quot;"; break;
                case '\'': escaped += "&#39;";  break;
                default:   escaped += c;
           }
        }
        return escaped;
    }













#if !defined(_WIN32) && !defined(__APPLE__)
    // Detect Windows Subsystem for Linux, where xdg-open usually has no browser and the Windows browser must be used instead.
    bool isRunningInWSL()
    {
        if (std::getenv("WSL_DISTRO_NAME") != nullptr)
        {
            return true;
        }

        std::ifstream version("/proc/version");
        std::string line;
        std::getline(version, line);
        return toLower(line).find("microsoft") != std::string::npos;
    }



    // Check that WSL can launch Windows programs.
    // Why: "Windows interop" can be switched off in WSL (common on locked-down work laptops). Then wslpath or
    // explorer.exe are missing, and we must fall back to telling the user where the file is.
    bool wslCanLaunchWindowsPrograms()
    {
        return std::system("command -v wslpath > /dev/null 2>&1 && command -v explorer.exe > /dev/null 2>&1") == 0;
    }
#endif




#if defined(__APPLE__)
    // Check whether this session is connected over SSH (only needed on macOS, Linux checks for a display instead).
    // Why: on a remote server, HPC cluster or container reached over SSH, the browser would have to open on the
    // remote machine, which has no screen. SSH sets these variables in every remote session.
    bool isRemoteSSHsession()
    {
        return std::getenv("SSH_CONNECTION") != nullptr || std::getenv("SSH_TTY") != nullptr;
    }
#endif





#if !defined(_WIN32)
    // Quote a path for a POSIX shell (Linux, macOS, WSL).
    // Why: the path goes into a std::system() command line. Input file names can contain spaces, $, quotes or other
    // characters the shell would otherwise interpret, which could break the command or run something unintended.
    // Inside single quotes the shell treats everything literally; an embedded ' is written as '\'' (close, escaped
    // quote, reopen).
    std::string shellQuote(const std::string& path)
    {
        std::string quoted = "'";
        for (char c : path)
        {
            quoted += (c == '\'') ? std::string("'\\''") : std::string(1, c);
        }
        quoted += "'";
        return quoted;
    }
#endif

}








bool writeKaryogramViewer(
    const std::vector<KaryogramImage>& images,
    const std::string& title,
    const std::filesystem::path& htmlPath)
{
    if (images.empty())
    {
        std::cerr << "ERROR: No karyograms to show in the web viewer.\n";
        return false;
    }

    std::ofstream html(htmlPath, std::ios::binary);

    if (!html)
    {
        std::cerr << "ERROR: Failed to create web viewer file " << htmlPath << "\n";
        return false;
    }

    const std::string safeTitle = htmlEscape(title);

    // ------------------------------------------
    // Page header and styling
    // ------------------------------------------
    html << "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
         << "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
         << "<title>" << safeTitle << " karyograms</title>\n"
         << "<style>\n"
         << "body { margin: 0; font-family: sans-serif; background: #f4f4f4; color: #222; }\n"
         << "header { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; padding: 10px 16px; background: #fff; border-bottom: 1px solid #ddd; }\n"
         << "header h1 { font-size: 16px; margin: 0 12px 0 0; }\n"
         << "button, select { font-size: 14px; padding: 4px 10px; }\n"
         << "#counter { color: #666; font-size: 14px; }\n"
         << "main { padding: 16px; text-align: center; }\n"
         << "main img { max-width: 100%; background: #fff; box-shadow: 0 1px 4px rgba(0,0,0,0.15); }\n"
         << "</style>\n</head>\n<body>\n";

    // ------------------------------------------
    // Toolbar to switch between karyograms
    // ------------------------------------------
    html << "<header>\n<h1>" << safeTitle << "</h1>\n"
         << "<button id=\"prev\">&larr; Previous</button>\n"
         << "<select id=\"cell\">\n";

    for (size_t i = 0; i < images.size(); i++)
    {
        html << "<option value=\"" << i << "\">" << htmlEscape(images[i].label) << "</option>\n";
    }

    html << "</select>\n"
         << "<button id=\"next\">Next &rarr;</button>\n"
         << "<span id=\"counter\"></span>\n"
         << "</header>\n"
         << "<main><img id=\"karyogram\" alt=\"Karyogram\"></main>\n";

    // ------------------------------------------
    // Embedded karyogram images and switching logic
    // ------------------------------------------
    html << "<script>\nconst karyograms = [\n";

    for (const KaryogramImage& image : images)
    {
        html << "  { filename: \"" << htmlEscape(image.filename) << "\", src: \"data:image/png;base64,"
             << base64Encode(image.pngBytes) << "\" },\n";
    }

    html << "];\n"
         << "const img = document.getElementById('karyogram');\n"
         << "const select = document.getElementById('cell');\n"
         << "const counter = document.getElementById('counter');\n"
         << "let current = 0;\n"
         << "function show(index) {\n"
         << "  current = (index + karyograms.length) % karyograms.length;\n"
         << "  img.src = karyograms[current].src;\n"
         << "  select.value = current;\n"
         << "  counter.textContent = (current + 1) + ' of ' + karyograms.length;\n"
         << "}\n"
         << "document.getElementById('prev').onclick = () => show(current - 1);\n"
         << "document.getElementById('next').onclick = () => show(current + 1);\n"
         << "select.onchange = () => show(Number(select.value));\n"
         << "show(0);\n"
         << "</script>\n</body>\n</html>\n";

    html.close();

    if (!html)
    {
        std::cerr << "ERROR: Failed to write web viewer file " << htmlPath << "\n";
        return false;
    }

    return true;
}









bool canOpenBrowser()
{
#if defined(_WIN32)
    return true;                                                                // Native Windows always has a desktop and a default browser.
#elif defined(__APPLE__)
    return !isRemoteSSHsession();                                               // A Mac reached over SSH would open the browser on its own (unseen) screen.
#else
    if (isRunningInWSL())
    {
        return wslCanLaunchWindowsPrograms();                                   // WSL opens the Windows browser, which needs Windows interop.
    }

    // Plain Linux: a browser can only appear if there is a graphical display. Desktops set DISPLAY (X11) or
    // WAYLAND_DISPLAY (Wayland); SSH with X forwarding (ssh -X) also sets DISPLAY, so that case still works.
    // Why: without a display, xdg-open either fails or, worse, starts a text-mode browser (such as lynx or w3m)
    // inside the terminal, which takes over the terminal and blocks the program until the user quits it.
    return std::getenv("DISPLAY") != nullptr || std::getenv("WAYLAND_DISPLAY") != nullptr;
#endif
}












bool openInBrowser(const std::filesystem::path& htmlPath)
{
    if (!canOpenBrowser())                                                      // Safeguard for headless machines, see canOpenBrowser().
    {
        return false;
    }

    const std::string path = std::filesystem::absolute(htmlPath).string();

#if defined(_WIN32)
    if (path.find('"') != std::string::npos)                                    // Windows paths cannot contain quotes, refuse anything odd.
    {
        return false;
    }
    const std::string command = "start \"\" \"" + path + "\"";			// The empty "" is the window title, start would otherwise treat the quoted path as the title.
#elif defined(__APPLE__)
    const std::string command = "open " + shellQuote(path)
#else
    if (isRunningInWSL())                                                       // wslpath -w turns /tmp/x.html into \\wsl.localhost\<distro>\tmp\x.html, a path Windows programs can read,
    {										 // and explorer.exe opens it with the user's default Windows browser.
	const std::string wslCommand = "explorer.exe \"$(wslpath -w " + shellQuote(path) + ")\" > /dev/null 2>&1";

	// explorer.exe returns 1 even when it succeeds, so its exit code cannot tell us anything. It is stored in a
        // variable only because some compilers warn when the result of std::system() is ignored.

	const int explorerStatus = std::system(wslCommand.c_str());
        (void)explorerStatus;
        return true;
    }

    const std::string command = "xdg-open " + shellQuote(path) + " > /dev/null 2>&1";	// Output is hidden because xdg-open can print browser noise into the CLI output.
#endif

    return std::system(command.c_str()) == 0;
}




} // namespace sddparser
