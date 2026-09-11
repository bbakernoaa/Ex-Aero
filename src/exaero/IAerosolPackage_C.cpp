#include <cstring>
#include <exaero/AerosolIndices.hpp>
#include <exaero/Environment.hpp>
#include <exaero/IAerosolPackage_C.h>
#include <exception>
#include <gocart/GocartPackage.hpp>

extern "C" {

exaero_package_t exaero_create_gocart_package() {
  return static_cast<exaero_package_t>(new exaero::GocartPackage());
}

void exaero_free_package(exaero_package_t pkg) {
  if (pkg) {
    delete static_cast<exaero::GocartPackage *>(pkg);
  }
}

void exaero_initialize_package(exaero_package_t pkg, const char *config_yaml,
                               char *errmsg, int *errflg) {
  if (!pkg) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg, "EX-aero Error: Null package handle");
    return;
  }
  try {
    static_cast<exaero::GocartPackage *>(pkg)->initialize(
        std::string(config_yaml));
    if (errflg)
      *errflg = 0;
    if (errmsg)
      errmsg[0] = '\0';
  } catch (const std::exception &e) {
    if (errflg)
      *errflg = 1;
    if (errmsg) {
      std::strncpy(errmsg, e.what(), 255);
      errmsg[255] = '\0';
    }
  } catch (...) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(
          errmsg,
          "EX-aero Error: Unknown exception occurred during initialization");
  }
}

void exaero_compute_diagnostics(exaero_package_t pkg, int num_cells,
                                int num_levels, int num_species,
                                const double *temp_ptr, const double *pres_ptr,
                                const double *dens_ptr, const double *rh_ptr,
                                const double *thick_ptr,
                                const double *state_ptr, double *diags_ptr,
                                char *errmsg, int *errflg) {

  if (!pkg) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg, "EX-aero Error: Null package handle");
    return;
  }
  try {
    auto *package = static_cast<exaero::GocartPackage *>(pkg);

    // Map incoming flat Fortran/C raw pointers to standard layout C++20 mdspan
    // Views on-the-fly (LayoutLeft)
    exaero::View2D<const double> temperature(temp_ptr, num_cells, num_levels);
    exaero::View2D<const double> pressure(pres_ptr, num_cells, num_levels);
    exaero::View2D<const double> air_density(dens_ptr, num_cells, num_levels);
    exaero::View2D<const double> relative_humidity(rh_ptr, num_cells,
                                                   num_levels);
    exaero::View2D<const double> layer_thickness(thick_ptr, num_cells,
                                                 num_levels);

    exaero::EnvironmentalStateView env{temperature, pressure, air_density,
                                       relative_humidity, layer_thickness};

    exaero::View3D<const double> state(state_ptr, num_cells, num_levels,
                                       num_species);
    exaero::View3D<double> diagnostics_out(
        diags_ptr, num_cells, num_levels,
        exaero::diagnostic_indices::NUM_DIAGNOSTICS);

    // Execute dynamic diagnostics solver
    package->computeDerivedDiagnostics(env, state, diagnostics_out);

    if (errflg)
      *errflg = 0;
    if (errmsg)
      errmsg[0] = '\0';
  } catch (const std::exception &e) {
    if (errflg)
      *errflg = 1;
    if (errmsg) {
      std::strncpy(errmsg, e.what(), 255);
      errmsg[255] = '\0';
    }
  } catch (...) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg, "EX-aero Error: Unknown exception occurred during "
                          "diagnostics calculation");
  }
}

void exaero_compute_emissions(exaero_package_t pkg, int num_cells,
                              int num_levels, int num_raw_species,
                              int num_target_species, int flux_type_code,
                              const double *thick_ptr,
                              const double *raw_emissions_ptr,
                              double *target_emissions_out_ptr, char *errmsg,
                              int *errflg) {

  if (!pkg) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg, "EX-aero Error: Null package handle");
    return;
  }
  try {
    auto *package = static_cast<exaero::IAerosolPackage *>(pkg);

    // Stub environmental view (only layer_thickness is required for emissions
    // scaling)
    exaero::View2D<const double> temperature(nullptr, 0, 0);
    exaero::View2D<const double> pressure(nullptr, 0, 0);
    exaero::View2D<const double> air_density(nullptr, 0, 0);
    exaero::View2D<const double> relative_humidity(nullptr, 0, 0);
    exaero::View2D<const double> layer_thickness(thick_ptr, num_cells,
                                                 num_levels);

    exaero::EnvironmentalStateView env{temperature, pressure, air_density,
                                       relative_humidity, layer_thickness};

    exaero::View3D<const double> flux(raw_emissions_ptr, num_cells, num_levels,
                                      num_raw_species);
    exaero::FluxType type = (flux_type_code == 1)
                                ? exaero::FluxType::AREA_FLUX
                                : exaero::FluxType::MASS_CONCENTRATION_RATE;
    exaero::EmissionsInputView emissions_in{flux, type};

    exaero::View3D<double> emissions_out(target_emissions_out_ptr, num_cells,
                                         num_levels, num_target_species);

    package->computeEmissions(env, emissions_in, emissions_out);

    if (errflg)
      *errflg = 0;
    if (errmsg)
      errmsg[0] = '\0';
  } catch (const std::exception &e) {
    if (errflg)
      *errflg = 1;
    if (errmsg) {
      std::strncpy(errmsg, e.what(), 255);
      errmsg[255] = '\0';
    }
  } catch (...) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg, "EX-aero Error: Unknown exception occurred during "
                          "emissions calculation");
  }
}

void exaero_compute_optics(exaero_package_t pkg, int num_cells, int num_levels,
                           int num_bands, int num_species,
                           const double *wavelengths_ptr,
                           const double *temp_ptr, const double *pres_ptr,
                           const double *dens_ptr, const double *rh_ptr,
                           const double *thick_ptr, const double *state_ptr,
                           double *optics_ptr, char *errmsg, int *errflg) {

  if (!pkg) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg, "EX-aero Error: Null package handle");
    return;
  }
  try {
    auto *package = static_cast<exaero::GocartPackage *>(pkg);

    exaero::View2D<const double> temperature(temp_ptr, num_cells, num_levels);
    exaero::View2D<const double> pressure(pres_ptr, num_cells, num_levels);
    exaero::View2D<const double> air_density(dens_ptr, num_cells, num_levels);
    exaero::View2D<const double> relative_humidity(rh_ptr, num_cells,
                                                   num_levels);
    exaero::View2D<const double> layer_thickness(thick_ptr, num_cells,
                                                 num_levels);

    exaero::EnvironmentalStateView env{temperature, pressure, air_density,
                                       relative_humidity, layer_thickness};

    exaero::View1D<const double> wavelengths(wavelengths_ptr, num_bands);
    exaero::View3D<const double> state(state_ptr, num_cells, num_levels,
                                       num_species);
    exaero::View4D<double> optics_out(optics_ptr, num_cells, num_levels,
                                      num_bands,
                                      exaero::optical_indices::NUM_OPTICS);

    // Execute dynamic optics solver
    package->computeOptics(env, state, wavelengths, optics_out);

    if (errflg)
      *errflg = 0;
    if (errmsg)
      errmsg[0] = '\0';
  } catch (const std::exception &e) {
    if (errflg)
      *errflg = 1;
    if (errmsg) {
      std::strncpy(errmsg, e.what(), 255);
      errmsg[255] = '\0';
    }
  } catch (...) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg, "EX-aero Error: Unknown exception occurred during "
                          "optics calculation");
  }
}

void exaero_compute_ccn(exaero_package_t pkg, int num_cells, int num_levels,
                        int num_ss, int num_species, const double *ss_ptr,
                        const double *temp_ptr, const double *rh_ptr,
                        const double *state_ptr, double *ccn_ptr, char *errmsg,
                        int *errflg) {

  if (!pkg) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg, "EX-aero Error: Null package handle");
    return;
  }
  try {
    auto *package = static_cast<exaero::GocartPackage *>(pkg);

    exaero::View2D<const double> temperature(temp_ptr, num_cells, num_levels);
    exaero::View2D<const double> relative_humidity(rh_ptr, num_cells,
                                                   num_levels);

    // CCCN does not need air_density, pressure or thickness, so we use
    // lightweight dummy views
    double dummy_raw[1] = {1.2};
    exaero::View2D<const double> air_density_view(dummy_raw, num_cells,
                                                  num_levels);
    exaero::View2D<const double> pressure_view(dummy_raw, num_cells,
                                               num_levels);
    exaero::View2D<const double> layer_thickness_view(dummy_raw, num_cells,
                                                      num_levels);

    exaero::EnvironmentalStateView env{temperature, pressure_view,
                                       air_density_view, relative_humidity,
                                       layer_thickness_view};

    exaero::View1D<const double> supersaturations(ss_ptr, num_ss);
    exaero::View3D<const double> state(state_ptr, num_cells, num_levels,
                                       num_species);
    exaero::View4D<double> ccn_out(ccn_ptr, num_cells, num_levels, num_ss, 1);

    // Execute dynamic cloud activation solver
    package->computeCCN(env, state, supersaturations, ccn_out);

    if (errflg)
      *errflg = 0;
    if (errmsg)
      errmsg[0] = '\0';
  } catch (const std::exception &e) {
    if (errflg)
      *errflg = 1;
    if (errmsg) {
      std::strncpy(errmsg, e.what(), 255);
      errmsg[255] = '\0';
    }
  } catch (...) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg, "EX-aero Error: Unknown exception occurred during "
                          "cloud CCN calculation");
  }
}

// --- Dynamic Species-to-Index Queries (Hole 4) ---
int exaero_get_species_index(exaero_package_t pkg, const char *name) {
  if (!pkg || !name)
    return -1;
  return static_cast<exaero::GocartPackage *>(pkg)->getSpeciesIndex(
      std::string(name));
}

void exaero_get_species_name(exaero_package_t pkg, int index, char *name_out,
                             int max_len) {
  if (!pkg || !name_out || max_len <= 0)
    return;
  try {
    std::string name =
        static_cast<exaero::GocartPackage *>(pkg)->getSpeciesName(index);
    int len = static_cast<int>(name.length());
    if (len >= max_len) {
      len = max_len - 1;
    }
    for (int i = 0; i < len; ++i) {
      name_out[i] = name[i];
    }
    name_out[len] = '\0'; // ensure null-termination
  } catch (...) {
    name_out[0] = '\0';
  }
}

// --- GEOSmie MIE attribute surface (contract §3, CCPP errmsg/errflg) ---
void exaero_set_attribute_activation(exaero_package_t pkg,
                                     const char *const *species_names,
                                     int num_species, int categories_mask,
                                     const char *runtime_file_path,
                                     char *errmsg, int *errflg) {

  if (!pkg) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg, "EX-aero Error: Null package handle");
    return;
  }
  try {
    std::vector<std::string> names;
    if (species_names && num_species > 0) {
      names.reserve(num_species);
      for (int i = 0; i < num_species; ++i)
        names.emplace_back(species_names[i] ? species_names[i] : "");
    }
    std::string file =
        runtime_file_path ? std::string(runtime_file_path) : std::string();
    static_cast<exaero::IAerosolPackage *>(pkg)->setAttributeActivation(
        names, categories_mask, file);
    if (errflg)
      *errflg = 0;
    if (errmsg)
      errmsg[0] = '\0';
  } catch (const std::exception &e) {
    if (errflg)
      *errflg = 1; // FATAL on bad file, no fallback (FR-009)
    if (errmsg) {
      std::strncpy(errmsg, e.what(), 255);
      errmsg[255] = '\0';
    }
  } catch (...) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg, "EX-aero Error: Unknown exception during activation");
  }
}

void exaero_set_species_curve_config(
    exaero_package_t pkg, const char *const *species_names,
    const char *const *source_labels, const double *const *radius_nodes,
    const int *num_radius_nodes, const double *const *override_values,
    const int *num_override_attributes,
    const int *const *override_attribute_ids, const char *const *interpolate,
    int num_curves, char *errmsg, int *errflg) {

  if (!pkg) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg, "EX-aero Error: Null package handle");
    return;
  }
  try {
    std::vector<exaero::SpeciesCurveConfig> curves;
    for (int i = 0; i < num_curves; ++i) {
      exaero::SpeciesCurveConfig cc;
      if (species_names && species_names[i])
        cc.species_name = species_names[i];
      if (source_labels && source_labels[i])
        cc.source_label = source_labels[i];
      if (interpolate && interpolate[i] && interpolate[i][0])
        cc.interpolate = interpolate[i];
      const int nR = (num_radius_nodes && num_radius_nodes[i] > 0)
                         ? num_radius_nodes[i]
                         : 0;
      if (nR > 0 && radius_nodes && radius_nodes[i]) {
        cc.radius_nodes.assign(radius_nodes[i], radius_nodes[i] + nR);
      }
      const int nAttr =
          (num_override_attributes && num_override_attributes[i] > 0)
              ? num_override_attributes[i]
              : 0;
      if (nAttr > 0 && override_values && override_values[i] &&
          override_attribute_ids && override_attribute_ids[i]) {
        const double *row = override_values[i];
        for (int a = 0; a < nAttr; ++a) {
          const int code = override_attribute_ids[i][a];
          exaero::SpeciesCurveConfig::Override ov;
          ov.category = static_cast<exaero::AttributeCategory>(code >> 8);
          ov.attribute_index = code & 0xff;
          if (nR > 0)
            ov.values.assign(row + static_cast<std::size_t>(a) * nR,
                             row + static_cast<std::size_t>(a + 1) * nR);
          cc.overrides.push_back(std::move(ov));
        }
      }
      curves.push_back(std::move(cc));
    }
    static_cast<exaero::IAerosolPackage *>(pkg)->setSpeciesCurveConfig(curves);
    if (errflg)
      *errflg = 0;
    if (errmsg)
      errmsg[0] = '\0';
  } catch (const std::exception &e) {
    if (errflg)
      *errflg = 1; // FATAL on bad config, no fallback (FR-009)
    if (errmsg) {
      std::strncpy(errmsg, e.what(), 255);
      errmsg[255] = '\0';
    }
  } catch (...) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg,
                  "EX-aero Error: Unknown exception during curve config");
  }
}

void exaero_compute_attributes(
    exaero_package_t pkg, int num_cells, int num_levels, const double *temp_ptr,
    const double *pres_ptr, const double *dens_ptr, const double *rh_ptr,
    const double *thick_ptr, const double *state_ptr, int species_index,
    int category, int num_bands, const double *wavelengths_ptr,
    double *attributes_out, int *status_out, char *errmsg, int *errflg) {

  if (!pkg) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg, "EX-aero Error: Null package handle");
    return;
  }
  try {
    auto *package = static_cast<exaero::IAerosolPackage *>(pkg);
    exaero::View2D<const double> temperature(temp_ptr, num_cells, num_levels);
    exaero::View2D<const double> pressure(pres_ptr, num_cells, num_levels);
    exaero::View2D<const double> air_density(dens_ptr, num_cells, num_levels);
    exaero::View2D<const double> relative_humidity(rh_ptr, num_cells,
                                                   num_levels);
    exaero::View2D<const double> layer_thickness(thick_ptr, num_cells,
                                                 num_levels);
    exaero::EnvironmentalStateView env{temperature, pressure, air_density,
                                       relative_humidity, layer_thickness};

    exaero::View3D<const double> state(state_ptr, num_cells, num_levels, 1);
    exaero::View1D<const double> wavelengths(wavelengths_ptr, num_bands);

    // num_attributes derived from category (never a fixed literal, R10).
    // Spectral multi-band output is band-major: num_bands * num_attributes
    // slots (T027). PolarizedMoment slots enumerate the species' (element,
    // moment) pairs from its curve (data, FR-003); a species without moments
    // yields 0 slots (FR-008).
    int num_attr = 0;
    switch (static_cast<exaero::AttributeCategory>(category)) {
    case exaero::AttributeCategory::Microphysical:
      num_attr = exaero::microphysical_indices::NUM_ATTRIBUTES;
      break;
    case exaero::AttributeCategory::SpectralOptical:
      num_attr = exaero::spectral_optical_indices::NUM_ATTRIBUTES;
      break;
    case exaero::AttributeCategory::PolarizedMoment: {
      int n_pol = 0, n_moment = 0;
      package->momentCounts(species_index, &n_pol, &n_moment);
      num_attr = n_pol * n_moment;
      break;
    }
    default:
      num_attr = 0;
    }
    const int bands = (static_cast<exaero::AttributeCategory>(category) ==
                           exaero::AttributeCategory::SpectralOptical &&
                       num_bands > 0)
                          ? num_bands
                          : 1;
    const int num_slots = bands * num_attr;
    exaero::View3D<double> attrs_out(attributes_out, num_cells, num_levels,
                                     num_slots);
    exaero::View3D<int> stat_out;
    exaero::View3D<int> *stat_ptr = nullptr;
    if (status_out) {
      stat_out =
          exaero::View3D<int>(status_out, num_cells, num_levels, num_slots);
      stat_ptr = &stat_out;
    }
    package->computeAttributes(env, state, species_index,
                               static_cast<exaero::AttributeCategory>(category),
                               wavelengths, attrs_out, stat_ptr);
    if (errflg)
      *errflg = 0;
    if (errmsg)
      errmsg[0] = '\0';
  } catch (const std::exception &e) {
    if (errflg)
      *errflg = 1;
    if (errmsg) {
      std::strncpy(errmsg, e.what(), 255);
      errmsg[255] = '\0';
    }
  } catch (...) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg,
                  "EX-aero Error: Unknown exception during compute_attributes");
  }
}

void exaero_query_attribute(exaero_package_t pkg, int species_index,
                            int category, int attribute_index, double rh,
                            double wavelength_m, double *value_out,
                            char *unit_out, int unit_max, char *version_out,
                            int version_max, int *status_out, char *errmsg,
                            int *errflg) {

  if (!pkg) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg, "EX-aero Error: Null package handle");
    return;
  }
  try {
    auto *package = static_cast<exaero::IAerosolPackage *>(pkg);
    exaero::ProvenanceInfo prov{};
    double value = 0.0;
    const exaero::AttributeStatus st = package->queryAttribute(
        species_index, static_cast<exaero::AttributeCategory>(category),
        attribute_index, rh, wavelength_m, &value, &prov);
    // FR-008: write the caller's value only on an available/interpolated result
    // -- never a silent 0 for not-activated / not-in-source.
    const bool ok = st == exaero::AttributeStatus::Available ||
                    st == exaero::AttributeStatus::AvailableFile ||
                    st == exaero::AttributeStatus::AvailableConfig ||
                    st == exaero::AttributeStatus::Interpolated;
    if (value_out && ok)
      *value_out = value;
    if (status_out)
      *status_out = static_cast<int>(st);
    if (unit_out && unit_max > 0) {
      std::strncpy(unit_out, prov.unit, unit_max - 1);
      unit_out[unit_max - 1] = '\0';
    }
    if (version_out && version_max > 0) {
      std::strncpy(version_out, prov.source_version, version_max - 1);
      version_out[version_max - 1] = '\0';
    }
    if (errflg)
      *errflg = 0;
    if (errmsg)
      errmsg[0] = '\0';
  } catch (const std::exception &e) {
    if (errflg)
      *errflg = 1;
    if (errmsg) {
      std::strncpy(errmsg, e.what(), 255);
      errmsg[255] = '\0';
    }
  } catch (...) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg,
                  "EX-aero Error: Unknown exception during query_attribute");
  }
}

int exaero_get_moment_counts(exaero_package_t pkg, int species_index,
                             int *num_pol_out, int *num_moment_out,
                             char *errmsg, int *errflg) {

  if (!pkg) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg, "EX-aero Error: Null package handle");
    return 1;
  }
  try {
    auto *package = static_cast<exaero::IAerosolPackage *>(pkg);
    package->momentCounts(species_index, num_pol_out, num_moment_out);
    if (errflg)
      *errflg = 0;
    if (errmsg)
      errmsg[0] = '\0';
    return 0;
  } catch (const std::exception &e) {
    if (errflg)
      *errflg = 1;
    if (errmsg) {
      std::strncpy(errmsg, e.what(), 255);
      errmsg[255] = '\0';
    }
    return 1;
  } catch (...) {
    if (errflg)
      *errflg = 1;
    if (errmsg)
      std::strcpy(errmsg,
                  "EX-aero Error: Unknown exception during get_moment_counts");
    return 1;
  }
}

void exaero_init_environment() { exaero::initialize_environment(); }

void exaero_finalize_environment() { exaero::finalize_environment(); }

} // extern "C"
