#include <utils/LutGenerator.hpp>
#include <yaml-cpp/yaml.h>
#include <stdexcept>

namespace exaero {

    LutGenerator::LutGenerator(const std::string& config_yaml) {
        YAML::Node config = YAML::Load(config_yaml);
        if (!config["generation_grid"]) {
            throw std::runtime_error("Config missing 'generation_grid' field");
        }
        
        auto grid = config["generation_grid"];
        auto rh_node = grid["rh_bins"];
        auto wave_node = grid["wavelengths"];

        if (!rh_node || !rh_node.IsSequence()) {
            throw std::runtime_error("Config 'rh_bins' is missing or not a sequence");
        }
        if (!wave_node || !wave_node.IsSequence()) {
            throw std::runtime_error("Config 'wavelengths' is missing or not a sequence");
        }

        num_rh_ = rh_node.size();
        num_bands_ = wave_node.size();
        num_species_ = 1; // single species generator for now
        num_moments_ = grid["legendre_moments"] ? grid["legendre_moments"].as<int>() : 16;

        for (int i = 0; i < num_rh_; ++i) rh_bins_.push_back(rh_node[i].as<double>());
        for (int i = 0; i < num_bands_; ++i) wavelengths_.push_back(wave_node[i].as<double>());

        // Allocate memory buffers
        ext_data_.resize(num_rh_ * num_bands_ * num_species_, 0.0);
        sca_data_.resize(num_rh_ * num_bands_ * num_species_, 0.0);
        moments_data_.resize(num_rh_ * num_bands_ * num_species_ * num_moments_, 0.0);

        // Wrap standard mdspan views zero-copy
        ext_coeff_view_ = View3D<double>(ext_data_.data(), num_rh_, num_bands_, num_species_);
        sca_coeff_view_ = View3D<double>(sca_data_.data(), num_rh_, num_bands_, num_species_);
        moments_view_ = View4D<double>(moments_data_.data(), num_rh_, num_bands_, num_species_, num_moments_);
    }

    void LutGenerator::generate(const std::string& custom_spheroid_json_path) {
        // Stub for compilation verification
    }

    std::string LutGenerator::to_yaml() const {
        return "species:\n  - name: \"Generated_Lut\"\n";
    }

} // namespace exaero
