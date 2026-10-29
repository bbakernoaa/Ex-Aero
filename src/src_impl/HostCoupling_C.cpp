/// @file HostCoupling_C.cpp
/// @brief Atomic metadata conversion into caller-owned C/Fortran buffers.
#include "HostCouplingInternal.hpp"
#include <exaero/HostCoupling_C.h>
using namespace exaero::host_detail;
extern "C" {
void exaero_host_coupling_counts(exaero_package_t h,int *rep,int *nf,int *np,int *ns,char *error,int *flag) {
  guarded(error,flag,[&]{require(rep&&nf&&np&&ns,"null metadata output");const auto v=exaero::hostCouplingInfo(package(h));
    *rep=int(v.representation);*nf=v.field_count;*np=int(v.populations.size());*ns=int(v.constituents.size());});
}
void exaero_host_population_info(exaero_package_t h,int i,int *index,int *kind,int *field,char *name,int capacity,char *error,int *flag) {
  guarded(error,flag,[&]{require(index&&kind&&field,"null population metadata output");const auto info=exaero::hostCouplingInfo(package(h));
    require(i>=0&&std::size_t(i)<info.populations.size(),"invalid population index");const auto &v=info.populations[i];nameBuffer(v.name,name,capacity);
    std::memcpy(name,v.name.c_str(),v.name.size()+1);*index=v.index;*kind=int(v.number_kind);*field=v.number_field;});
}
void exaero_host_constituent_info(exaero_package_t h,int i,int *pop,int *species,int *field,int *permitted,char *name,int capacity,char *error,int *flag) {
  guarded(error,flag,[&]{require(pop&&species&&field&&permitted,"null constituent metadata output");const auto info=exaero::hostCouplingInfo(package(h));
    require(i>=0&&std::size_t(i)<info.constituents.size(),"invalid constituent index");const auto &v=info.constituents[i];nameBuffer(v.name,name,capacity);
    std::memcpy(name,v.name.c_str(),v.name.size()+1);*pop=v.population_index;*species=v.species_index;*field=v.mass_field;*permitted=int(v.permitted);});
}
void exaero_supports_host_operation(exaero_package_t h,int op,int *supported,char *error,int *flag) {
  guarded(error,flag,[&]{require(supported,"null capability output");const bool v=exaero::supportsHostOperation(package(h),static_cast<exaero::HostOperation>(op));*supported=int(v);});
}
}
