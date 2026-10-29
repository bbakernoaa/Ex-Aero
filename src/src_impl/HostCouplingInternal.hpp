#pragma once
/// @file HostCouplingInternal.hpp
/// @brief Private checked host-array and C error adapters; no I/O or owning raw pointers.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <exaero/HostCoupling.hpp>
#include <exaero/IAerosolPackage_C.h>
#include <limits>
#include <stdexcept>
namespace exaero::host_detail {
inline void require(bool ok, const char *why) {
  if (!ok) throw std::invalid_argument(std::string("FATAL ERROR: ")+why);
}
template<class V> void shape(const V &v, std::initializer_list<std::size_t> dims) {
  require(v.data_handle()!=nullptr,"null host coupling array");
  std::size_t i=0;
  for(auto n:dims) require(v.extent(i++)==n,"host coupling array shape mismatch");
}
inline void separate(const void *a,std::size_t an,const void *b,std::size_t bn) {
  const auto ap=reinterpret_cast<std::uintptr_t>(a),bp=reinterpret_cast<std::uintptr_t>(b);
  const auto max=std::numeric_limits<std::uintptr_t>::max();
  require(ap<=max-an && bp<=max-bn && (ap+an<=bp || bp+bn<=ap),"overlapping host coupling arrays");
}
template<class A,class B> void disjoint(const A &a,const B &b) {
  separate(a.data_handle(),a.size()*sizeof(*a.data_handle()),b.data_handle(),b.size()*sizeof(*b.data_handle()));
}
template<class F> void guarded(char *error,int *flag,F action) noexcept {
  try { action(); if(error) error[0]='\0'; if(flag) *flag=0; }
  catch(const std::exception &e) { if(error){std::strncpy(error,e.what(),255);error[255]='\0';} if(flag)*flag=1; }
  catch(...) { if(error)std::strcpy(error,"FATAL ERROR: host coupling failed"); if(flag)*flag=1; }
}
inline IAerosolPackage &package(exaero_package_t handle) {
  require(handle!=nullptr,"null host coupling package"); return *static_cast<IAerosolPackage*>(handle);
}
inline HostCouplingInfo dimensions(IAerosolPackage &p,int nc,int nl,int nf,int np,int extra=1) {
  require(nc>0&&nl>0&&nf>0&&np>0&&extra>0,"invalid host coupling dimensions");
  const auto extent=std::max({nf,np,extra});
  require(nc<=std::numeric_limits<int>::max()/nl/extent/extra,"host coupling dimension overflow");
  auto info=hostCouplingInfo(p);
  require(nf==info.field_count&&std::size_t(np)==info.populations.size(),"host coupling population/field count mismatch"); return info;
}
inline EnvironmentalStateView environment(int nc,int nl,const double *temp,const double *rh) {
  require(temp&&rh,"null host coupling environment");
  return {View2D<const double>(temp,nc,nl),View2D<const double>(nullptr,nc,nl),View2D<const double>(nullptr,nc,nl),View2D<const double>(rh,nc,nl),View2D<const double>(nullptr,nc,nl)};
}
inline void nameBuffer(const std::string &name,char *out,int capacity) {
  require(out&&capacity>0&&name.size()<std::size_t(capacity),"host metadata name buffer too small");
}
} // namespace exaero::host_detail
