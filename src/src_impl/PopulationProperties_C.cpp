/// @file PopulationProperties_C.cpp
/// @brief Exception-safe physical property buffers for C and Fortran hosts.
#include "HostCouplingInternal.hpp"
#include <exaero/PopulationProperties.hpp>
#include <exaero/PopulationProperties_C.h>
using namespace exaero;
using namespace exaero::host_detail;
extern "C" void exaero_population_properties(exaero_package_t h,int nc,int nl,int nf,int np,
    const double *t,const double *rh,const double *state,double *values,int *available,int *status,char *error,int *flag) {
  guarded(error,flag,[&]{auto &p=package(h);constexpr int count=int(PopulationProperty::Count);
    dimensions(p,nc,nl,nf,np,count);require(state&&values&&available&&status,"null property buffer");
    View4D<double> out(values,nc,nl,np,count);View4D<int> codes(available,nc,nl,np,count);View3D<int> st(status,nc,nl,np);
    populationProperties(p,environment(nc,nl,t,rh),View3D<const double>(state,nc,nl,nf),out,codes,st);});
}
