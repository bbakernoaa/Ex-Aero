#include <loader/Provenance.hpp>

#include <cstring>

namespace exaero {

void copy_provenance_field(char* dst, std::size_t dst_size, const std::string& src) {
    if (dst == nullptr || dst_size == 0) return;
    const std::size_t n = std::min(src.size(), dst_size - 1);
    std::memcpy(dst, src.data(), n);
    dst[n] = '\0';
}

ProvenanceInfo Provenance::to_public(int num_radius, int num_rh, int num_lambda,
                                     int num_pol, int num_moment) const {
    ProvenanceInfo info{};
    copy_provenance_field(info.species, sizeof(info.species), species);
    copy_provenance_field(info.unit, sizeof(info.unit), unit);
    copy_provenance_field(info.source_version, sizeof(info.source_version), source_version);
    copy_provenance_field(info.citation, sizeof(info.citation), citation);
    info.delivery_source = static_cast<int>(delivery);
    info.interpolated = interpolated ? 1 : 0;
    info.status = static_cast<int>(status);
    info.num_radius = num_radius;
    info.num_rh = num_rh;
    info.num_lambda = num_lambda;
    info.num_pol = num_pol;
    info.num_moment = num_moment;
    return info;
}

} // namespace exaero
