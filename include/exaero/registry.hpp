#ifndef EXAERO_REGISTRY_HPP
#define EXAERO_REGISTRY_HPP

#include <string>
#include "exaero/state.hpp"

namespace exaero {
class Registry {
public:
    static bool has_mechanism(const std::string& name);
    static void dispatch(ExaeroContext& ctx, UnmanagedDeviceState& state, double dt);
};
}
#endif
