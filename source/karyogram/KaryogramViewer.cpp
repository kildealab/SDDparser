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

    // Write text as a JavaScript string literal (with quotes) for use inside the viewer's <script>.
    // Why: file names and titles come from the input file name, and could contain quotes, backslashes or "</script>",
    // which would otherwise break the script. HTML escaping does not apply inside <script>, so JavaScript escaping is used.
    std::string jsString(const std::string& text)
    {
        std::string quoted = "\"";
        for (char c : text)
        {
            switch (c)
            {
                case '\\': quoted += "\\\\"; break;
                case '"':  quoted += "\\\""; break;
                case '<':  quoted += "\\u003c"; break;
                case '\n': quoted += "\\n"; break;
                case '\r': quoted += "\\r"; break;
                default:   quoted += c;
            }
        }
        return quoted + "\"";
    }

#if !defined(_WIN32) && !defined(__APPLE__)
    // Detect Windows Subsystem for Linux (WSL).
    // Why: a program built inside WSL is a Linux program, so it would normally use xdg-open. But WSL usually has no
    // Linux browser installed, so xdg-open fails silently and nothing appears. In WSL we instead hand the file to the
    // Windows browser. WSL sets WSL_DISTRO_NAME, and its kernel version string mentions "microsoft" as a backup check.
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
    // The karyogram sits inside a scrollable "viewport". Zooming changes the image's real width and height (not a CSS
    // scale transform), so the browser redraws the SVG at the new size and lines and text stay sharp at every zoom level.
    html << "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
         << "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
         << "<title>" << safeTitle << " karyograms</title>\n"
         << R"HTML(<style>
html, body { height: 100%; margin: 0; }
body { display: flex; flex-direction: column; font-family: sans-serif; background: #e9e9e9; color: #222; }
header { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; padding: 8px 16px; background: #fff; border-bottom: 1px solid #ccc; }
header h1 { font-size: 16px; margin: 0 12px 0 0; }
.group { display: flex; gap: 4px; align-items: center; }
.spacer { flex: 1; }
button, select { font-size: 14px; padding: 4px 10px; cursor: pointer; }
#counter, #zoomLevel { color: #555; font-size: 14px; min-width: 52px; text-align: center; }
#viewport { flex: 1; overflow: auto; position: relative; cursor: grab; user-select: none; touch-action: none; }
#viewport.dragging { cursor: grabbing; }
#stage { position: relative; }
#karyogram { position: absolute; display: block; background: #fff; box-shadow: 0 1px 6px rgba(0,0,0,0.25); }
#hint { padding: 4px 16px; font-size: 12px; color: #666; background: #fff; border-top: 1px solid #ccc; }
</style>
</head>
<body>
)HTML";

    // ------------------------------------------
    // Toolbar: switch between karyograms, zoom controls
    // ------------------------------------------
    html << "<header>\n<h1>" << safeTitle << "</h1>\n"
         << "<div class=\"group\">\n"
         << "<button id=\"prev\" title=\"Previous karyogram\">&larr; Previous</button>\n"
         << "<select id=\"cell\" title=\"Choose a karyogram\">\n";

    for (size_t i = 0; i < images.size(); i++)
    {
        html << "<option value=\"" << i << "\">" << htmlEscape(images[i].label) << "</option>\n";
    }

    html << R"HTML(</select>
<button id="next" title="Next karyogram">Next &rarr;</button>
<span id="counter"></span>
</div>
<div class="spacer"></div>
<div class="group">
<button id="zoomOut" title="Zoom out (-)">&minus;</button>
<span id="zoomLevel"></span>
<button id="zoomIn" title="Zoom in (+)">+</button>
<button id="fit" title="Fit the whole karyogram in the window (0)">Fit</button>
<button id="actual" title="Show at original size (1)">100%</button>
</div>
<div class="group">
<select id="format" title="File format for downloads">
<option value="png">PNG</option>
<option value="svg">SVG</option>
</select>
<button id="download" title="Download the karyogram shown">Download</button>
<button id="downloadAll" title="Download every karyogram in one .zip file">Download all (.zip)</button>
</div>
</header>
<div id="viewport"><div id="stage"><img id="karyogram" alt="Karyogram" draggable="false"></div></div>
<div id="hint">Scroll or pinch to zoom at the cursor &middot; drag to pan &middot; double-click to zoom in &middot; keys: + &minus; 0 (fit) 1 (100%)</div>
)HTML";

    // ------------------------------------------
    // Embedded karyogram images
    // ------------------------------------------
    // Each image is embedded as a base64 data: URL so the HTML file works on its own, with no other files or server.
    // The SVG is shown and both formats are kept for the download buttons.
    html << "<script>\nconst TITLE = " << jsString(title) << ";\nconst karyograms = [\n";

    for (const KaryogramImage& image : images)
    {
        html << "  { filename: " << jsString(image.filename)
             << ",\n    svg: \"" << (image.svgBytes.empty() ? "" : "data:image/svg+xml;base64," + base64Encode(image.svgBytes))
             << "\",\n    png: \"" << (image.pngBytes.empty() ? "" : "data:image/png;base64," + base64Encode(image.pngBytes))
             << "\" },\n";
    }

    // ------------------------------------------
    // Switching, zoom and pan logic
    // ------------------------------------------
    html << R"HTML(];

const viewport = document.getElementById('viewport');
const stage = document.getElementById('stage');
const img = document.getElementById('karyogram');
const select = document.getElementById('cell');
const counter = document.getElementById('counter');
const zoomLevel = document.getElementById('zoomLevel');

const MIN_SCALE = 0.1, MAX_SCALE = 20, STEP = 1.25, MARGIN = 16;
let current = 0;
let scale = 1;
let fitMode = true;            // While true the image refits when the window is resized.

// Scale at which the whole karyogram fits inside the window.
function fitScale() {
  const w = img.naturalWidth || 1, h = img.naturalHeight || 1;
  return Math.min((viewport.clientWidth - 2 * MARGIN) / w, (viewport.clientHeight - 2 * MARGIN) / h);
}

// Offset of the image inside the stage: centred when smaller than the window, otherwise a small margin.
function offsets(s) {
  const w = img.naturalWidth * s, h = img.naturalHeight * s;
  return {
    x: Math.max(MARGIN, (viewport.clientWidth - w) / 2),
    y: Math.max(MARGIN, (viewport.clientHeight - h) / 2),
  };
}

// Apply a new zoom level while keeping the image point under (screenX, screenY) fixed on screen.
function setScale(newScale, screenX, screenY) {
  newScale = Math.min(MAX_SCALE, Math.max(MIN_SCALE, newScale));
  if (screenX === undefined) { screenX = viewport.clientWidth / 2; screenY = viewport.clientHeight / 2; }

  const before = offsets(scale);
  const imageX = (viewport.scrollLeft + screenX - before.x) / scale;    // Point under the cursor, in image pixels.
  const imageY = (viewport.scrollTop + screenY - before.y) / scale;

  scale = newScale;
  const after = offsets(scale);
  const w = img.naturalWidth * scale, h = img.naturalHeight * scale;
  img.style.width = w + 'px';              // Real layout size, so the SVG is redrawn sharply at this size.
  img.style.height = h + 'px';
  img.style.left = after.x + 'px';
  img.style.top = after.y + 'px';
  stage.style.width = Math.max(viewport.clientWidth, w + 2 * after.x) + 'px';
  stage.style.height = Math.max(viewport.clientHeight, h + 2 * after.y) + 'px';

  viewport.scrollLeft = imageX * scale + after.x - screenX;
  viewport.scrollTop = imageY * scale + after.y - screenY;
  zoomLevel.textContent = Math.round(scale * 100) + '%';
}

function fit() { fitMode = true; setScale(fitScale()); }
function zoomBy(factor, x, y) { fitMode = false; setScale(scale * factor, x, y); }

function show(index) {
  current = (index + karyograms.length) % karyograms.length;
  const first = !img.src;
  img.onload = () => { if (first || fitMode) fit(); else setScale(scale); };   // Keep the zoom when switching cells, to compare them.
  img.src = karyograms[current].svg || karyograms[current].png;      // Prefer the sharp vector image.
  select.value = current;
  counter.textContent = (current + 1) + ' of ' + karyograms.length;
}

// Mouse wheel and trackpad pinch (sent as wheel events with ctrlKey) zoom around the cursor.
viewport.addEventListener('wheel', (e) => {
  e.preventDefault();
  const r = viewport.getBoundingClientRect();
  const delta = e.deltaMode === 1 ? e.deltaY * 33 : e.deltaY;           // Line-based scrolling (Firefox) to pixels.
  zoomBy(Math.exp(-delta * (e.ctrlKey ? 0.01 : 0.002)), e.clientX - r.left, e.clientY - r.top);
}, { passive: false });

// Drag with the mouse (or one finger) to pan.
let drag = null;
viewport.addEventListener('pointerdown', (e) => {
  if (e.button !== 0) return;
  drag = { x: e.clientX, y: e.clientY, left: viewport.scrollLeft, top: viewport.scrollTop };
  viewport.setPointerCapture(e.pointerId);
  viewport.classList.add('dragging');
});
viewport.addEventListener('pointermove', (e) => {
  if (!drag) return;
  viewport.scrollLeft = drag.left - (e.clientX - drag.x);
  viewport.scrollTop = drag.top - (e.clientY - drag.y);
});
const endDrag = () => { drag = null; viewport.classList.remove('dragging'); };
viewport.addEventListener('pointerup', endDrag);
viewport.addEventListener('pointercancel', endDrag);

viewport.addEventListener('dblclick', (e) => {
  const r = viewport.getBoundingClientRect();
  zoomBy(2, e.clientX - r.left, e.clientY - r.top);
});

document.getElementById('zoomIn').onclick = () => zoomBy(STEP);
document.getElementById('zoomOut').onclick = () => zoomBy(1 / STEP);
document.getElementById('fit').onclick = fit;
document.getElementById('actual').onclick = () => { fitMode = false; setScale(1); };
document.getElementById('prev').onclick = () => show(current - 1);
document.getElementById('next').onclick = () => show(current + 1);
select.onchange = () => show(Number(select.value));

document.addEventListener('keydown', (e) => {
  if (e.target === select || e.ctrlKey || e.metaKey || e.altKey) return;
  if (e.key === '+' || e.key === '=') zoomBy(STEP);
  else if (e.key === '-' || e.key === '_') zoomBy(1 / STEP);
  else if (e.key === '0') fit();
  else if (e.key === '1') { fitMode = false; setScale(1); }
});

// ------------------------------------------------------------
// Downloads
// ------------------------------------------------------------
const formatSelect = document.getElementById('format');
for (const option of formatSelect.options) {                       // Only offer formats that were embedded.
  option.disabled = !karyograms.some((k) => k[option.value]);
}
if (formatSelect.selectedOptions[0].disabled) formatSelect.value = 'svg';

// "cell0_karyogram_human.png" -> "cell0_karyogram_human.svg" for SVG downloads.
function downloadName(k, format) {
  return k.filename.replace(/\.(png|svg)$/i, '') + '.' + format;
}

// Decode an embedded base64 data: URL back into the original file bytes.
function dataUrlBytes(dataUrl) {
  const binary = atob(dataUrl.slice(dataUrl.indexOf(',') + 1));
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
  return bytes;
}

// Hand a file to the browser's normal download (it asks where to save, or uses the Downloads folder).
function saveFile(bytes, name, type) {
  const url = URL.createObjectURL(new Blob([bytes], { type }));
  const link = document.createElement('a');
  link.href = url;
  link.download = name;
  document.body.appendChild(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

const MIME = { png: 'image/png', svg: 'image/svg+xml' };

document.getElementById('download').onclick = () => {
  const k = karyograms[current], format = formatSelect.value;
  if (!k[format]) { alert('This karyogram has no ' + format.toUpperCase() + ' version.'); return; }
  saveFile(dataUrlBytes(k[format]), downloadName(k, format), MIME[format]);
};

// Build a .zip file in the browser (stored, not compressed: PNGs are already compressed and this keeps the code short).
// Why a zip: browsers block or prompt for many separate downloads at once, a single file avoids that.
const CRC_TABLE = (() => {
  const table = new Uint32Array(256);
  for (let n = 0; n < 256; n++) {
    let c = n;
    for (let k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320 ^ (c >>> 1) : c >>> 1;
    table[n] = c >>> 0;
  }
  return table;
})();

function crc32(bytes) {
  let c = 0xFFFFFFFF;
  for (let i = 0; i < bytes.length; i++) c = CRC_TABLE[(c ^ bytes[i]) & 0xFF] ^ (c >>> 8);
  return (c ^ 0xFFFFFFFF) >>> 0;
}

function makeZip(files) {
  const encoder = new TextEncoder();
  const chunks = [], central = [];
  let offset = 0;
  const now = new Date();
  const dosTime = (now.getHours() << 11) | (now.getMinutes() << 5) | (now.getSeconds() >> 1);
  const dosDate = ((now.getFullYear() - 1980) << 9) | ((now.getMonth() + 1) << 5) | now.getDate();

  for (const file of files) {
    const name = encoder.encode(file.name), crc = crc32(file.bytes), size = file.bytes.length;

    const local = new DataView(new ArrayBuffer(30));
    local.setUint32(0, 0x04034B50, true);   // Local file header signature
    local.setUint16(4, 20, true);           // Version needed to extract
    local.setUint16(6, 0x0800, true);       // Flags: file name is UTF-8
    local.setUint16(8, 0, true);            // Method: stored
    local.setUint16(10, dosTime, true);
    local.setUint16(12, dosDate, true);
    local.setUint32(14, crc, true);
    local.setUint32(18, size, true);        // Compressed size
    local.setUint32(22, size, true);        // Uncompressed size
    local.setUint16(26, name.length, true);
    local.setUint16(28, 0, true);           // Extra field length
    chunks.push(new Uint8Array(local.buffer), name, file.bytes);

    const entry = new DataView(new ArrayBuffer(46));
    entry.setUint32(0, 0x02014B50, true);   // Central directory signature
    entry.setUint16(4, 20, true);           // Version made by
    entry.setUint16(6, 20, true);           // Version needed
    entry.setUint16(8, 0x0800, true);
    entry.setUint16(10, 0, true);
    entry.setUint16(12, dosTime, true);
    entry.setUint16(14, dosDate, true);
    entry.setUint32(16, crc, true);
    entry.setUint32(20, size, true);
    entry.setUint32(24, size, true);
    entry.setUint16(28, name.length, true);
    entry.setUint32(42, offset, true);      // Where this file's local header starts
    central.push(new Uint8Array(entry.buffer), name);

    offset += 30 + name.length + size;
  }

  const centralSize = central.reduce((total, part) => total + part.length, 0);
  const end = new DataView(new ArrayBuffer(22));
  end.setUint32(0, 0x06054B50, true);       // End of central directory signature
  end.setUint16(8, files.length, true);
  end.setUint16(10, files.length, true);
  end.setUint32(12, centralSize, true);
  end.setUint32(16, offset, true);
  return new Blob([...chunks, ...central, new Uint8Array(end.buffer)], { type: 'application/zip' });
}

document.getElementById('downloadAll').onclick = () => {
  const format = formatSelect.value;
  const files = karyograms.filter((k) => k[format])
    .map((k) => ({ name: downloadName(k, format), bytes: dataUrlBytes(k[format]) }));
  const zipName = TITLE.replace(/\.[^.]*$/, '') + '_karyograms_' + format + '.zip';
  saveFile(makeZip(files), zipName, 'application/zip');
};

window.addEventListener('resize', () => { if (fitMode) fit(); else setScale(scale); });

show(0);
</script>
</body>
</html>
)HTML";

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
    // Windows file names cannot contain double quotes, so a quote here means something is wrong with the path.
    // Refusing it keeps the quoting of the "start" command below intact.
    if (path.find('"') != std::string::npos)
    {
        return false;
    }
    const std::string command = "start \"\" \"" + path + "\"";                    // The empty "" is the window title, start would otherwise treat the quoted path as the title.
#elif defined(__APPLE__)
    const std::string command = "open " + shellQuote(path);
#else
    if (isRunningInWSL())
    {
        // wslpath -w turns /tmp/x.html into \\wsl.localhost\<distro>\tmp\x.html, a path Windows programs can read,
        // and explorer.exe opens it with the user's default Windows browser.
        const std::string wslCommand =
            "explorer.exe \"$(wslpath -w " + shellQuote(path) + ")\" > /dev/null 2>&1";

        // explorer.exe returns 1 even when it succeeds, so its exit code cannot tell us anything. It is stored in a
        // variable only because some compilers warn when the result of std::system() is ignored.
        const int explorerStatus = std::system(wslCommand.c_str());
        (void)explorerStatus;
        return true;
    }

    const std::string command = "xdg-open " + shellQuote(path) + " > /dev/null 2>&1";   // Output is hidden because xdg-open can print browser noise into the CLI output.
#endif

    return std::system(command.c_str()) == 0;
}

} // namespace sddparser
