#include <loader/MieTableStore.hpp>
#include <loader/MieFileLoader.hpp>

#include <data/generated/MieTableData.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace exaero {

namespace {

// Fail-loud diagnostic shared by every validation/query guard (FR-009).
[[noreturn]] void fatal(const std::string& msg) {
    throw std::runtime_error("FATAL ERROR: MieTableStore: " + msg);
}

constexpr double kPi = 3.14159265358979323846;

// Slack for float32-sourced grid values (the pinned source is float32 widened to float64).
constexpr double kF32Slack = 1e-6;

// Map a public (category, index) to the canonical field name resolved by the store.
// Returns "" for an index with no backing field (=> NotInSource, never a silent 0).
std::string field_for(AttributeCategory cat, int idx) {
    using namespace exaero;
    if (cat == AttributeCategory::Microphysical) {
        switch (idx) {
            case microphysical_indices::WET_PARTICLE_DENSITY: return "wet_particle_density";
            case microphysical_indices::GROWTH_FACTOR:        return "growth_factor";
            case microphysical_indices::EFFECTIVE_RADIUS:     return "rEff";
            case microphysical_indices::MASS_MEAN_RADIUS:     return "rMass_mean"; // no source var
            case microphysical_indices::BIN_LOWER_RADIUS:     return "rLow";
            case microphysical_indices::BIN_UPPER_RADIUS:     return "rUp";
            case microphysical_indices::VOLUME_PER_MASS:      return "volume_per_mass";
            case microphysical_indices::AREA_PER_MASS:        return "area_per_mass";
            case microphysical_indices::PARTICLE_MASS:        return "particle_mass";
            default: return "";
        }
    }
    if (cat == AttributeCategory::SpectralOptical) {
        switch (idx) {
            case spectral_optical_indices::EXTINCTION_EFFICIENCY:  return "qext";
            case spectral_optical_indices::SCATTERING_EFFICIENCY:  return "qsca";
            case spectral_optical_indices::ABSORPTION_EFFICIENCY:  return "qabs";
            case spectral_optical_indices::MASS_EXTINCTION:        return "bext";
            case spectral_optical_indices::MASS_SCATTERING:        return "bsca";
            case spectral_optical_indices::MASS_BACKSCATTER:       return "bbck";
            case spectral_optical_indices::LIDAR_RATIO:            return "lidar_ratio";
            case spectral_optical_indices::ASYMMETRY_FACTOR:       return "g";
            case spectral_optical_indices::SINGLE_SCATTERING_ALBEDO: return "ssa";
            case spectral_optical_indices::REFRACTIVE_INDEX_REAL:  return "refreal";
            case spectral_optical_indices::REFRACTIVE_INDEX_IMAG:  return "refimag";
            default: return "";
        }
    }
    if (cat == AttributeCategory::PolarizedMoment) {
        return "pmom"; // supplied only via runtime file (FR-017); absent => NotInSource
    }
    return "";
}

// Canonical unit per field (data-model.md). Used for provenance + unit validation.
// The public static MieTableStore::unit_for_field forwards here (single source of truth).
std::string unit_for_field_impl(const std::string& field) {
    static const std::map<std::string, std::string> u = {
        {"wet_particle_density", "kg m^-3"}, {"growth_factor", "1"},
        {"rEff", "m"}, {"rMass_mean", "m"}, {"rLow", "m"}, {"rUp", "m"},
        {"volume_per_mass", "m^3 kg^-1"}, {"area_per_mass", "m^2 kg^-1"},
        {"particle_mass", "kg"},
        {"qext", "1"}, {"qsca", "1"}, {"qabs", "1"},
        {"bext", "m^2 (kg dry mass)^-1"}, {"bsca", "m^2 (kg dry mass)^-1"},
        {"bbck", "m^2 (kg dry mass)^-1 sr^-1"}, {"lidar_ratio", "sr"},
        {"g", "1"}, {"ssa", "1"}, {"refreal", "1"}, {"refimag", "1"},
        {"pmom", "1"},
    };
    auto it = u.find(field);
    return it == u.end() ? std::string("1") : it->second;
}

// Linear search for the bracketing interval on a strictly-increasing coordinate axis.
// Deterministic tie rule: exact grid hit returns (i, i, 0.0); below-min clamps to first,
// above-max clamps to last (documented, no extrapolation beyond declared edges).
void locate(const std::vector<double>& axis, double v, int& i0, int& i1, double& w) {
    const int n = static_cast<int>(axis.size());
    if (n == 1) { i0 = i1 = 0; w = 0.0; return; }
    if (v <= axis.front()) {
        i0 = 0;
        w = (v == axis.front()) ? 0.0 : 0.0;
        i1 = (v == axis.front()) ? 0 : 1;
        if (v == axis.front()) { i1 = 0; }
        return;
    }
    if (v >= axis.back()) { i0 = n - 1; i1 = n - 1; w = 0.0; return; }
    // binary search upper bound
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (axis[mid] <= v) lo = mid; else hi = mid;
    }
    i0 = lo; i1 = hi;
    double span = axis[i1] - axis[i0];
    w = (span > 0.0) ? (v - axis[i0]) / span : 0.0;
}

} // namespace

const char* MieTableStore::unit_for_field(const std::string& field) {
    // Stable storage: the canonical unit set is process-lifetime static (FR-004).
    static const std::map<std::string, std::string> units = {
        {"wet_particle_density", "kg m^-3"}, {"growth_factor", "1"},
        {"rEff", "m"}, {"rMass_mean", "m"}, {"rLow", "m"}, {"rUp", "m"},
        {"volume_per_mass", "m^3 kg^-1"}, {"area_per_mass", "m^2 kg^-1"},
        {"particle_mass", "kg"},
        {"qext", "1"}, {"qsca", "1"}, {"qabs", "1"},
        {"bext", "m^2 (kg dry mass)^-1"}, {"bsca", "m^2 (kg dry mass)^-1"},
        {"bbck", "m^2 (kg dry mass)^-1 sr^-1"}, {"lidar_ratio", "sr"},
        {"g", "1"}, {"ssa", "1"}, {"refreal", "1"}, {"refimag", "1"},
        {"pmom", "1"},
    };
    auto it = units.find(field);
    return it == units.end() ? "1" : it->second.c_str();
}

// --- Baked-in registry: label -> raw array pointers (dims read from data, not hardcoded) ---
namespace {
struct BakedArrays {
    const double* radius; const double* rh; const double* lambda;
    const double* rLow; const double* rUp; const double* rEff; const double* rMass;
    const double* qext; const double* qsca; const double* bext; const double* bsca;
    const double* bbck; const double* g; const double* refreal; const double* refimag;
};

template <const double* A, const double* B>
BakedArrays unused_guard(); // never called; keeps -Wunused happy pattern out of the way

#define EXAERO_BAKED(lbl) BakedArrays{ \
    exaero::mie_baked::lbl##_radius, exaero::mie_baked::lbl##_rh, exaero::mie_baked::lbl##_lambda, \
    exaero::mie_baked::lbl##_rLow, exaero::mie_baked::lbl##_rUp, exaero::mie_baked::lbl##_rEff, \
    exaero::mie_baked::lbl##_rMass, exaero::mie_baked::lbl##_qext, exaero::mie_baked::lbl##_qsca, \
    exaero::mie_baked::lbl##_bext, exaero::mie_baked::lbl##_bsca, exaero::mie_baked::lbl##_bbck, \
    exaero::mie_baked::lbl##_g, exaero::mie_baked::lbl##_refreal, exaero::mie_baked::lbl##_refimag }

const BakedArrays& baked_for(int species_i) {
    // Order matches mie_baked::SPECIES_LABELS (DU, SS, SU, BC, OC, NI).
    static const BakedArrays table[] = {
        EXAERO_BAKED(DU), EXAERO_BAKED(SS), EXAERO_BAKED(SU),
        EXAERO_BAKED(BC), EXAERO_BAKED(OC), EXAERO_BAKED(NI),
    };
    static_assert(sizeof(table) / sizeof(table[0]) == exaero::mie_baked::NUM_SPECIES,
                  "baked registry must match generated NUM_SPECIES");
    return table[species_i];
}
#undef EXAERO_BAKED

std::vector<double> slice(const double* base, int n) {
    return std::vector<double>(base, base + n);
}

} // namespace

MieTableStore& MieTableStore::instance() {
    static MieTableStore store;
    return store;
}

MieTableStore::~MieTableStore() = default;

void MieTableStore::ensure_baked_in_loaded() {
    if (baked_loaded_) return;
    curves_.clear();
    for (int i = 0; i < mie_baked::NUM_SPECIES; ++i) {
        const std::string label = mie_baked::SPECIES_LABELS[i];
        const auto& dims = mie_baked::SPECIES_DIMS[i];
        const BakedArrays& a = baked_for(i);
        SpeciesCurve c;
        c.source_label = label;
        c.species_name = label; // default config name == label until config binds
        c.source_version = mie_baked::SOURCE_VERSION; // FR-016
        c.citation = mie_baked::CITATION;
        c.radius = slice(a.radius, dims.n_radius);
        c.rh = slice(a.rh, dims.n_rh);
        c.lambda = slice(a.lambda, dims.n_lambda);
        c.radius_resamplable = std::adjacent_find(c.radius.begin(), c.radius.end(),
                                                 [](double a, double b) { return !(b > a); }) == c.radius.end() &&
                               c.radius.size() >= 2;
        auto put = [&](const std::string& name, std::vector<double> v, int rank) {
            CurveField f; f.values = std::move(v); f.rank = rank; f.unit = unit_for_field_impl(name);
            c.fields[name] = std::move(f);
        };
        int nR = dims.n_radius, nH = dims.n_rh, nL = dims.n_lambda;
        put("rLow", slice(a.rLow, nR), 1);
        put("rUp", slice(a.rUp, nR), 1);
        put("rEff", slice(a.rEff, nR * nH), 2);
        put("particle_mass", slice(a.rMass, nR * nH), 2);
        put("qext", slice(a.qext, nR * nH * nL), 3);
        put("qsca", slice(a.qsca, nR * nH * nL), 3);
        put("bext", slice(a.bext, nR * nH * nL), 3);
        put("bsca", slice(a.bsca, nR * nH * nL), 3);
        put("bbck", slice(a.bbck, nR * nH * nL), 3);
        put("g", slice(a.g, nR * nH * nL), 3);
        put("refreal", slice(a.refreal, nR * nH * nL), 3);
        put("refimag", slice(a.refimag, nR * nH * nL), 3);
        derive_fields(c);
        validate_curve(c);
        curves_[label] = std::move(c);
    }
    activation_mask_ = attribute_category_bit(AttributeCategory::Microphysical) |
                       attribute_category_bit(AttributeCategory::SpectralOptical);
    baked_loaded_ = true;
}

// True when a field is present and came from a layer above baked (file/config): then a
// derived recomputation must NOT clobber the supplied values (FR-012 precedence).
static bool supplied_above_baked(const SpeciesCurve& c, const std::string& name) {
    auto it = c.fields.find(name);
    return it != c.fields.end() && it->second.delivery != DeliverySource::BakedIn;
}

void MieTableStore::derive_fields(SpeciesCurve& c) {
    const int nR = c.n_radius(), nH = c.n_rh(), nL = c.n_lambda();
    if (!c.fields.count("rEff") || !c.fields.count("particle_mass")) return; // need raw inputs
    auto& rEff = c.fields["rEff"].values;
    auto& rMass = c.fields["particle_mass"].values;

    // Derived microphysical (rank 2: radius, rh).
    CurveField gf; gf.rank = 2; gf.unit = unit_for_field_impl("growth_factor");
    gf.values.assign(nR * nH, 0.0);
    CurveField wd; wd.rank = 2; wd.unit = unit_for_field_impl("wet_particle_density");
    wd.values.assign(nR * nH, 0.0);
    CurveField vm; vm.rank = 2; vm.unit = unit_for_field_impl("volume_per_mass");
    vm.values.assign(nR * nH, 0.0);
    CurveField am; am.rank = 2; am.unit = unit_for_field_impl("area_per_mass");
    am.values.assign(nR * nH, 0.0);
    for (int b = 0; b < nR; ++b) {
        double rEff_dry = rEff[b * nH + 0]; // rh index 0 == dry
        for (int h = 0; h < nH; ++h) {
            int i = b * nH + h;
            double re = rEff[i], m = rMass[i];
            double vol = (4.0 / 3.0) * kPi * re * re * re;
            gf.values[i] = (rEff_dry > 0.0) ? re / rEff_dry : 1.0;
            wd.values[i] = (vol > 0.0) ? m / vol : 0.0;
            vm.values[i] = (m > 0.0) ? vol / m : 0.0;
            am.values[i] = (m > 0.0) ? (kPi * re * re) / m : 0.0;
        }
    }
    if (!supplied_above_baked(c, "growth_factor")) c.fields["growth_factor"] = std::move(gf);
    if (!supplied_above_baked(c, "wet_particle_density")) c.fields["wet_particle_density"] = std::move(wd);
    if (!supplied_above_baked(c, "volume_per_mass")) c.fields["volume_per_mass"] = std::move(vm);
    if (!supplied_above_baked(c, "area_per_mass")) c.fields["area_per_mass"] = std::move(am);

    if (nL > 0) {
        auto& qext = c.fields["qext"].values;
        auto& qsca = c.fields["qsca"].values;
        auto& bext = c.fields["bext"].values;
        auto& bbck = c.fields["bbck"].values;
        const std::size_t N = qext.size();
        CurveField qabs; qabs.rank = 3; qabs.unit = "1"; qabs.values.assign(N, 0.0);
        CurveField ssa; ssa.rank = 3; ssa.unit = "1"; ssa.values.assign(N, 0.0);
        CurveField lr; lr.rank = 3; lr.unit = unit_for_field_impl("lidar_ratio"); lr.values.assign(N, 0.0);
        for (std::size_t k = 0; k < N; ++k) {
            qabs.values[k] = qext[k] - qsca[k];
            ssa.values[k] = (qext[k] > 1e-30) ? qsca[k] / qext[k] : 0.0; // guard div0
            lr.values[k] = (bbck[k] > 1e-30) ? bext[k] / bbck[k] : 0.0;  // guard div0
        }
        if (!supplied_above_baked(c, "qabs")) c.fields["qabs"] = std::move(qabs);
        if (!supplied_above_baked(c, "ssa")) c.fields["ssa"] = std::move(ssa);
        if (!supplied_above_baked(c, "lidar_ratio")) c.fields["lidar_ratio"] = std::move(lr);
    }
}

void MieTableStore::validate_curve(const SpeciesCurve& c) const {
    if (c.source_version.empty())
        fatal("curve " + c.species_name + " has empty source_version (FR-016)");
    if (c.radius.empty() || c.rh.empty()) fatal("curve " + c.species_name + " has empty axis");
    for (double r : c.radius) if (!(r > 0.0) || !std::isfinite(r))
        fatal("curve " + c.species_name + ": radius must be finite and > 0");
    for (double h : c.rh) if (!(h >= -kF32Slack && h <= 0.99 + kF32Slack))
        fatal("curve " + c.species_name + ": rh out of [0,0.99]");
    for (int i = 1; i < static_cast<int>(c.rh.size()); ++i)
        if (!(c.rh[i] > c.rh[i - 1])) fatal("curve " + c.species_name + ": rh not strictly increasing");

    auto check = [&](const std::string& name, auto&& pred, const char* what) {
        auto it = c.fields.find(name);
        if (it == c.fields.end()) return;
        for (double v : it->second.values) if (!pred(v))
            fatal("curve " + c.species_name + ": " + name + " " + what);
    };
    check("rEff", [](double v) { return v > 0.0 && std::isfinite(v); }, "must be > 0");
    check("particle_mass", [](double v) { return v > 0.0 && std::isfinite(v); }, "must be > 0");
    check("growth_factor", [](double v) { return v >= 1.0 - kF32Slack && std::isfinite(v); }, "must be >= 1");
    check("wet_particle_density", [](double v) { return v > 0.0 && std::isfinite(v); }, "must be > 0");
    check("qext", [](double v) { return v >= -kF32Slack && std::isfinite(v); }, "must be >= 0");
    check("qsca", [](double v) { return v >= -kF32Slack && std::isfinite(v); }, "must be >= 0");
    check("ssa", [](double v) { return v >= -kF32Slack && v <= 1.0 + kF32Slack; }, "outside [0,1]");
    check("g", [](double v) { return v >= -1.0 - kF32Slack && v <= 1.0 + kF32Slack; }, "outside [-1,1]");
    // qsca <= qext
    auto qe = c.fields.find("qext"), qs = c.fields.find("qsca");
    if (qe != c.fields.end() && qs != c.fields.end())
        for (std::size_t k = 0; k < qe->second.values.size(); ++k)
            if (qs->second.values[k] > qe->second.values[k] + 1e-5)
                fatal("curve " + c.species_name + ": qsca > qext");
}

void MieTableStore::reset_to_baked_in() {
    baked_loaded_ = false;
    ensure_baked_in_loaded();
}

void MieTableStore::apply_runtime_file(const std::string& path) {
    if (path.empty()) return;
    // MieFileLoader parses + validates the portable file (SPEC-CRTM-004 pattern) and aborts
    // with FATAL ERROR on any failure (FR-009). File values override baked per key (FR-012).
    std::vector<SpeciesCurve> loaded = MieFileLoader::load(path);
    for (auto& c : loaded) {
        for (auto& [name, f] : c.fields) f.delivery = DeliverySource::RuntimeFile;
        auto it = curves_.find(c.species_name);
        if (it == curves_.end()) {
            derive_fields(c);   // file supplies raw curves; derived fields follow (R9)
            validate_curve(c);
            curves_[c.species_name] = std::move(c);
        } else {
            SpeciesCurve& dst = it->second;
            if (!c.lambda.empty()) { dst.lambda = c.lambda; }
            for (auto& [name, f] : c.fields) dst.fields[name] = std::move(f);
            // New axes may change monotonicity; recompute resamplability from data (R9).
            dst.radius_resamplable = dst.radius.size() >= 2 &&
                std::adjacent_find(dst.radius.begin(), dst.radius.end(),
                                   [](double a, double b) { return !(b > a); }) == dst.radius.end();
            derive_fields(dst);
            validate_curve(dst);
        }
    }
}

void MieTableStore::apply_curve_config(const std::vector<SpeciesCurveConfig>& curves) {
    for (const auto& cc : curves) {
        if (cc.interpolate != "linear")
            fatal("unknown interpolate strategy '" + cc.interpolate + "' for species " + cc.species_name);
        if (cc.source_label.empty()) continue; // no curve bound => stays NotInSource
        auto base = curves_.find(cc.source_label);
        if (base == curves_.end())
            fatal("config species " + cc.species_name + " binds unknown source '" + cc.source_label + "'");

        SpeciesCurve c = base->second; // copy the merged (baked+file) curve
        c.species_name = cc.species_name;
        if (cc.solver_radius_node >= 0) c.solver_radius_node = cc.solver_radius_node;

        // A pure alias ({source: X}, no nodes/overrides) keeps the source delivery so
        // grid-point queries stay bit-identical to the golden references (C11). Anything
        // else is a config-modified curve (R9).
        if (!cc.radius_nodes.empty() || !cc.overrides.empty()) c.config_modified = true;

        // Radius-only resample onto config nodes (R9). Disabled when the source radius axis
        // is not strictly increasing (multi-mode BC/OC): duplicate coordinates cannot define
        // a curve, so we keep the source index axis and flag the species non-resamplable.
        if (!cc.radius_nodes.empty()) {
            const auto& nodes = cc.radius_nodes;
            for (int i = 1; i < static_cast<int>(nodes.size()); ++i)
                if (!(nodes[i] > nodes[i - 1]))
                    fatal("config species " + cc.species_name + ": radius_nodes not strictly increasing");
            if (!c.radius_resamplable) {
                // Explicit, non-silent: we do NOT interpolate over duplicate radii.
                c.radius = nodes;
                c.radius_resamplable = false;
                // Resample is a no-op for non-curve axes; fields retain source index order.
            } else {
                resample_radius(c, nodes);
                c.radius = nodes;
            }
        }

        // Apply per-attribute overrides at the config radius nodes (delivery=config, R9).
        for (const auto& ov : cc.overrides) {
            std::string fname = field_for(ov.category, ov.attribute_index);
            if (fname.empty())
                fatal("config species " + cc.species_name + ": override index has no field");
            int nR = c.n_radius();
            if (static_cast<int>(ov.values.size()) != nR)
                fatal("config species " + cc.species_name + ": override '" + fname +
                      "' length " + std::to_string(ov.values.size()) + " != radius nodes " + std::to_string(nR));
            CurveField f;
            f.rank = (ov.category == AttributeCategory::Microphysical) ? 2 : 3;
            f.unit = unit_for_field_impl(fname);
            f.delivery = DeliverySource::Config;
            int nH = c.n_rh(), nL = std::max(1, c.n_lambda());
            f.values.assign(static_cast<std::size_t>(nR) * nH * nL, 0.0);
            for (int b = 0; b < nR; ++b) {
                for (int h = 0; h < nH; ++h) {
                    for (int l = 0; l < nL; ++l) {
                        std::size_t idx = (static_cast<std::size_t>(b) * nH + h) * nL + l;
                        f.values[idx] = ov.values[b]; // constant across rh/lambda at that node
                    }
                }
            }
            c.fields[fname] = std::move(f);
        }

        // Recompute derived fields from the (possibly resampled/overridden) raw inputs so
        // growth/density/volume/area/qabs/ssa/lidar stay consistent with the curve. Fields
        // supplied by file or config (delivery != baked) are never clobbered (FR-012, R9).
        derive_fields(c);
        validate_curve(c);
        curves_[cc.species_name] = std::move(c);
    }
}

// Resample every rank-2/rank-3 field in a curve along radius onto new nodes.
void MieTableStore::resample_radius(SpeciesCurve& c, const std::vector<double>& nodes) {
    const int nOld = c.n_radius();
    const int nH = c.n_rh(), nL = std::max(1, c.n_lambda());
    std::vector<double>& oldR = c.radius;
    for (auto& [name, f] : c.fields) {
        if (f.rank < 2) continue;
        std::vector<double> src = f.values;
        int nCol = (f.rank == 2) ? 1 : nL;
        std::size_t newN = static_cast<std::size_t>(nodes.size()) * nH * nCol;
        f.values.assign(newN, 0.0);
        for (int nb = 0; nb < static_cast<int>(nodes.size()); ++nb) {
            int b0, b1; double wr; locate(oldR, nodes[nb], b0, b1, wr);
            for (int h = 0; h < nH; ++h) {
                for (int l = 0; l < nCol; ++l) {
                    std::size_t i0 = (static_cast<std::size_t>(b0) * nH + h) * nCol + l;
                    std::size_t i1 = (static_cast<std::size_t>(b1) * nH + h) * nCol + l;
                    std::size_t o = (static_cast<std::size_t>(nb) * nH + h) * nCol + l;
                    f.values[o] = src[i0] + wr * (src[i1] - src[i0]);
                }
            }
        }
        (void)nOld;
    }
}

bool MieTableStore::is_activated(const std::string& species, AttributeCategory category) const {
    if (!(activation_mask_ & attribute_category_bit(category))) return false;
    if (!activated_species_.empty() &&
        std::find(activated_species_.begin(), activated_species_.end(), species) == activated_species_.end())
        return false;
    return true;
}

void MieTableStore::set_activation(const std::vector<std::string>& species, int categories_mask) {
    activation_mask_ = categories_mask;
    activated_species_ = species; // empty => all bound species
}

AttributeStatus MieTableStore::query(const std::string& species, AttributeCategory category,
                                     int attribute_index, double rh, double wavelength_m,
                                     int radius_index, double* value_out, Provenance* prov_out) const {
    auto it = curves_.find(species);
    if (it == curves_.end()) return AttributeStatus::NotInSource;
    const SpeciesCurve& c = it->second;
    if (!is_activated(species, category)) return AttributeStatus::NotActivated;

    std::string fname = field_for(category, attribute_index);
    if (fname.empty()) return AttributeStatus::NotInSource;
    auto fit = c.fields.find(fname);
    if (fit == c.fields.end()) return AttributeStatus::NotInSource; // e.g. mass_mean_radius, pmom
    const CurveField& f = fit->second;

    if (!std::isfinite(rh)) return AttributeStatus::NotInSource; // reject NaN/Inf (FR-007, never silent 0)
    double rh_clamped = std::clamp(rh, 0.0, 0.99);

    int nR = c.n_radius();
    int b = (radius_index < 0) ? 0 : radius_index;
    if (b >= nR) return AttributeStatus::NotInSource;

    const int nH = c.n_rh(), nL = std::max(1, c.n_lambda());
    int h0, h1; double wh; locate(c.rh, rh_clamped, h0, h1, wh);

    double value = 0.0;
    bool interp = false;
    if (f.rank == 1) {
        value = f.values[b];
        interp = false;
    } else if (f.rank == 2) {
        double v0 = f.values[static_cast<std::size_t>(b) * nH + h0];
        double v1 = f.values[static_cast<std::size_t>(b) * nH + h1];
        value = v0 + wh * (v1 - v0);
        interp = wh > 0.0 && h0 != h1;
    } else if (f.rank == 5) {
        // Polarized moment (rank 5). Minimal correct grid-point read; (element, moment)
        // decomposition + interpolation refined in T032. radius_index addresses a flat slot.
        if (b < 0 || static_cast<std::size_t>(b) >= f.values.size())
            return AttributeStatus::NotInSource;
        value = f.values[static_cast<std::size_t>(b)];
        interp = false;
    } else { // rank 3: interpolate rh and lambda
        if (!std::isfinite(wavelength_m)) return AttributeStatus::NotInSource; // spectral needs a real band (never silent 0)
        int l0, l1; double wl; locate(c.lambda, wavelength_m, l0, l1, wl);
        auto at = [&](int bb, int hh, int ll) {
            return f.values[((static_cast<std::size_t>(bb) * nH) + hh) * nL + ll];
        };
        double a = at(b, h0, l0), bb = at(b, h1, l0);
        double c2 = at(b, h0, l1), d = at(b, h1, l1);
        double e = a + wh * (bb - a);
        double g2 = c2 + wh * (d - c2);
        value = e + wl * (g2 - e);
        interp = (wh > 0.0 && h0 != h1) || (wl > 0.0 && l0 != l1);
    }

    if (value_out) *value_out = value;
    if (prov_out) {
        prov_out->species = species;
        prov_out->attribute = fname;
        prov_out->unit = f.unit;
        prov_out->source_version = c.source_version;
        prov_out->citation = c.citation;
        prov_out->delivery = f.delivery;
        prov_out->interpolated = interp;
        prov_out->status = interp ? AttributeStatus::Interpolated
            : (f.delivery == DeliverySource::BakedIn ? AttributeStatus::Available
               : f.delivery == DeliverySource::RuntimeFile ? AttributeStatus::AvailableFile
               : AttributeStatus::AvailableConfig);
    }
    if (interp) return AttributeStatus::Interpolated;
    switch (f.delivery) {
        case DeliverySource::BakedIn: return AttributeStatus::Available;
        case DeliverySource::RuntimeFile: return AttributeStatus::AvailableFile;
        case DeliverySource::Config: return AttributeStatus::AvailableConfig;
    }
    return AttributeStatus::Available;
}

const std::vector<std::string>& MieTableStore::serialized_field_order() {
    // Fixed, documented order shared by the solver reader (R10). Axes serialized separately.
    static const std::vector<std::string> order = {
        "rEff", "particle_mass", "rLow", "rUp", "growth_factor",
        "wet_particle_density", "volume_per_mass", "area_per_mass",
        "qext", "qsca", "qabs", "bext", "bsca", "bbck", "lidar_ratio", "g",
        "ssa", "refreal", "refimag",
    };
    return order;
}

void MieTableStore::serialize_activated(const std::vector<std::string>& species_order,
                                         std::vector<double>& pool_out,
                                         std::vector<CurvePoolDescriptor>& desc_out) const {
    pool_out.clear();
    desc_out.clear();
    const auto& order = serialized_field_order();
    for (int si = 0; si < static_cast<int>(species_order.size()); ++si) {
        const std::string& name = species_order[si];
        auto it = curves_.find(name);
        CurvePoolDescriptor d;
        d.species_index = si;
        d.source_label = (it == curves_.end()) ? "" : it->second.source_label;
        d.field_names = order;
        if (it == curves_.end()) {
            d.n_radius = d.n_rh = d.n_lambda = 0;
            d.field_offsets.assign(order.size(), -1);
            d.curve_offset = static_cast<int>(pool_out.size());
            desc_out.push_back(std::move(d));
            continue;
        }
        const SpeciesCurve& c = it->second;
        d.n_radius = c.n_radius(); d.n_rh = c.n_rh(); d.n_lambda = c.n_lambda();
        d.curve_offset = static_cast<int>(pool_out.size());
        d.field_offsets.assign(order.size(), -1);
        for (std::size_t fi = 0; fi < order.size(); ++fi) {
            auto fit = c.fields.find(order[fi]);
            if (fit == c.fields.end()) continue;
            d.field_offsets[fi] = static_cast<std::int64_t>(pool_out.size()) - d.curve_offset;
            pool_out.insert(pool_out.end(), fit->second.values.begin(), fit->second.values.end());
        }
        desc_out.push_back(std::move(d));
    }
}

const SpeciesCurve* MieTableStore::find_curve(const std::string& species) const {
    auto it = curves_.find(species);
    return it == curves_.end() ? nullptr : &it->second;
}

} // namespace exaero
