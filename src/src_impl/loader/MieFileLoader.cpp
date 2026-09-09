#include <loader/MieFileLoader.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace exaero {

namespace {

[[noreturn]] void loader_fatal(const std::string& msg) {
    throw std::runtime_error("FATAL ERROR: MieFileLoader: " + msg);
}

// Normalize a Fortran-style exponent token ("1.0d-05", "2D+3") to C form ("1.0e-05").
// Only a d/D that follows digit-or-dot-and-precedes-sign-or-digit is rewritten.
std::string normalize_exponent(const std::string& tok) {
    std::string s = tok;
    for (std::size_t i = 1; i + 1 < s.size(); ++i) {
        if ((s[i] == 'd' || s[i] == 'D') &&
            (std::isdigit(static_cast<unsigned char>(s[i - 1])) || s[i - 1] == '.') &&
            (std::isdigit(static_cast<unsigned char>(s[i + 1])) || s[i + 1] == '+' || s[i + 1] == '-')) {
            s[i] = 'e';
        }
    }
    return s;
}

double parse_double(const std::string& tok, const std::string& where) {
    std::string s = normalize_exponent(tok);
    try {
        std::size_t pos = 0;
        double v = std::stod(s, &pos);
        if (pos != s.size() || !std::isfinite(v)) loader_fatal("malformed number '" + tok + "' in " + where);
        return v;
    } catch (const std::exception&) {
        loader_fatal("malformed number '" + tok + "' in " + where);
    }
}

// Split a line into whitespace-separated tokens.
std::vector<std::string> tokens(const std::string& line) {
    std::vector<std::string> out;
    std::istringstream ss(line);
    std::string t;
    while (ss >> t) out.push_back(t);
    return out;
}

// Consume the rest of a directive line as a free-form string (citation / version).
std::string rest_of(const std::string& line, std::size_t from) {
    std::string s = line.substr(std::min(from, line.size()));
    auto a = s.find_first_not_of(" \t");
    if (a == std::string::npos) return "";
    auto b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

std::size_t axis_length(const SpeciesCurve& c, const std::string& name) {
    if (name == "radius") return c.radius.size();
    if (name == "rh") return c.rh.size();
    if (name == "lambda") return c.lambda.size();
    if (name == "pol") return c.pol.size();
    if (name == "moment") return c.moment.size();
    return 0;
}

} // namespace

bool MieFileLoader::path_is_safe(const std::string& path) {
    if (path.empty() || path.size() > 4096) return false;
    if (path.find('\0') != std::string::npos) return false;
    // Reject directory traversal (Security baseline: path sanitization). Absolute and
    // relative paths are both acceptable — the file is opened directly, never via a shell.
    if (path.find("..") != std::string::npos) return false;
    for (char ch : path) {
        unsigned char c = static_cast<unsigned char>(ch);
        if (std::isalnum(c) || ch == '_' || ch == '-' || ch == '.' || ch == '/' || ch == '+') continue;
        return false; // spaces, shell metacharacters, control bytes
    }
    return true;
}

std::vector<SpeciesCurve> MieFileLoader::load(const std::string& path) {
    if (!path_is_safe(path)) loader_fatal("rejected unsafe path: " + path);
    std::ifstream file(path);
    if (!file.is_open()) loader_fatal("failed to open file: " + path);

    std::vector<SpeciesCurve> curves;
    std::string format;
    SpeciesCurve* cur = nullptr;
    std::string version_global, citation_global;

    std::string line;
    while (std::getline(file, line)) {
        // Strip trailing CR (CRLF tolerance) and leading whitespace.
        while (!line.empty() && (line.back() == '\r' || line.back() == ' ' || line.back() == '\t')) line.pop_back();
        auto first = line.find_first_not_of(" \t");
        if (first == std::string::npos) continue; // blank
        // Pure comment lines (SPEC-CRTM-004): '#' or '!' not introducing a directive.
        if ((line[first] == '#' || line[first] == '!')) {
            if (line[first] != '#' || line.find("@") == std::string::npos) continue;
        }

        std::string body = (line[first] == '#') ? line.substr(first + 1) : line.substr(first);
        auto tk = tokens(body);
        if (tk.empty()) continue;
        const std::string& kw = tk[0];

        if (kw == "@format") {
            if (tk.size() < 2) loader_fatal("@format missing value");
            format = tk[1];
            if (format != "exaero_mie_v1") loader_fatal("unsupported format '" + format + "'");
        } else if (kw == "@source_version") {
            version_global = rest_of(body, 15);
            if (version_global.empty()) loader_fatal("@source_version must be non-empty (FR-016)");
            if (cur) cur->source_version = version_global;
        } else if (kw == "@citation") {
            citation_global = rest_of(body, 9);
            if (cur) cur->citation = citation_global;
        } else if (kw == "@species") {
            if (tk.size() < 2) loader_fatal("@species missing name");
            curves.emplace_back();
            cur = &curves.back();
            cur->species_name = tk[1];
            cur->source_label = tk[1];
            cur->source_version = version_global;
            cur->citation = citation_global;
        } else if (kw == "@axis") {
            if (!cur) loader_fatal("@axis outside @species block");
            if (tk.size() < 4) loader_fatal("@axis needs <name> <units> <n>");
            const std::string name = tk[1];
            const std::string units = tk[2];
            const int n = static_cast<int>(parse_double(tk[3], "@axis count"));
            if (n <= 0) loader_fatal("@axis " + name + ": non-positive count");
            std::vector<double> vals;
            for (std::size_t i = 4; i < tk.size(); ++i) vals.push_back(parse_double(tk[i], "@axis values"));
            // Values may continue on following plain lines until n are collected.
            while (static_cast<int>(vals.size()) < n && std::getline(file, line)) {
                for (const auto& t : tokens(line)) vals.push_back(parse_double(t, "@axis values"));
            }
            if (static_cast<int>(vals.size()) != n)
                loader_fatal("@axis " + name + ": truncated (" + std::to_string(vals.size()) + " of " + std::to_string(n) + ")");
            auto& dst = (name == "radius") ? cur->radius : (name == "rh") ? cur->rh
                      : (name == "lambda") ? cur->lambda : (name == "pol") ? cur->pol : cur->moment;
            dst = std::move(vals);
            if (name == "radius") {
                cur->radius_resamplable = cur->radius.size() >= 2 &&
                    std::adjacent_find(cur->radius.begin(), cur->radius.end(),
                                       [](double a, double b) { return !(b > a); }) == cur->radius.end();
            }
            (void)units; // axis units validated by the store's curve validation
        } else if (kw == "@field") {
            if (!cur) loader_fatal("@field outside @species block");
            if (tk.size() < 4) loader_fatal("@field needs <name> <dims> <units>");
            const std::string name = tk[1];
            // dims: comma-joined axis names, e.g. radius,rh,lambda
            std::vector<std::string> dims;
            {
                std::stringstream dss(tk[2]);
                std::string d;
                while (std::getline(dss, d, ',')) dims.push_back(d);
            }
            const std::string units = tk[3];
            std::size_t count = 1;
            int rank = 1;
            for (const auto& d : dims) {
                std::size_t len = axis_length(*cur, d);
                if (len == 0) loader_fatal("@field " + name + ": unknown or empty axis '" + d + "'");
                count *= len;
            }
            rank = static_cast<int>(dims.size());
            CurveField f;
            f.unit = units;
            f.rank = rank;
            f.delivery = DeliverySource::RuntimeFile;
            f.dims = dims;
            f.category = (name == "pmom") ? AttributeCategory::PolarizedMoment
                       : (dims.size() >= 3 || name == "qext" || name == "qsca" || name == "bext" ||
                          name == "bsca" || name == "bbck" || name == "g" || name == "refreal" ||
                          name == "refimag" || name == "qabs" || name == "ssa" || name == "lidar_ratio")
                           ? AttributeCategory::SpectralOptical
                           : AttributeCategory::Microphysical;
            std::size_t have = 0;
            while (have < count && std::getline(file, line)) {
                for (const auto& t : tokens(line)) {
                    if (have >= count) loader_fatal("@field " + name + ": too many values");
                    f.values.push_back(parse_double(t, "@field " + name));
                    ++have;
                }
            }
            if (have < count)
                loader_fatal("@field " + name + ": truncated (" + std::to_string(have) + " of " + std::to_string(count) + ")");
            cur->fields[name] = std::move(f);
        } else {
            loader_fatal("unknown directive '" + kw + "'");
        }
    }

    if (format.empty()) loader_fatal("missing @format directive");
    if (curves.empty()) loader_fatal("file declares no @species blocks");
    if (version_global.empty()) loader_fatal("missing @source_version (FR-016)");
    for (auto& c : curves) {
        if (c.source_version.empty()) c.source_version = version_global;
    }
    return curves;
}

} // namespace exaero
