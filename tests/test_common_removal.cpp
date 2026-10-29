#ifdef NDEBUG
#undef NDEBUG
#endif
#include "host_coupling_fixture.hpp"
#include <exaero/HostCoupling.hpp>
#include <exaero/Removal.hpp>
#include <exaero/Removal_C.h>
#include <limits>
int main() {
  host_test::Runtime runtime;
  using namespace exaero;
  for(bool modal:{false,true}) {
    auto p=host_test::make(modal);const int nc=2,nl=3,nf=modal?6:2,np=2;
    std::vector<double>x(nc*nl*nf),rate(x.size()),mass(x.size()),number(nc*nl*np),nr(number.size(),1e-4),vel(number.size(),.01),dz(nc*nl,10);
    View3D<double> xv(x.data(),nc,nl,nf),rv(rate.data(),nc,nl,nf),mv(mass.data(),nc,nl,nf),nv(number.data(),nc,nl,np);
    View3D<const double> input(x.data(),nc,nl,nf),rates(rate.data(),nc,nl,nf),nrate(nr.data(),nc,nl,np),velocity(vel.data(),nc,nl,np);
    for(int k=0;k<nl;++k)for(int c=0;c<nc;++c){
      xv(c,k,0)=2e-9*(c+k+1);xv(c,k,1)=3e-9*(c+k+1);rv(c,k,0)=5e-4;rv(c,k,1)=2e-4;
      if(modal){xv(c,k,2)=5e7; xv(c,k,3)=1e-9; xv(c,k,5)=1e7; rv(c,k,3)=1e-3;}
    }
    const auto original=x;
    removalRateTendencies(*p,input,rates,modal?&nrate:nullptr,mv,nv);
    assert(host_test::near(mv(0,0,0),-1e-12)&&host_test::near(mv(0,0,1),-6e-13));
    if(modal)assert(nv(0,0,0)==-5e3&&mv(0,0,2)==0&&mv(0,0,5)==0);
    auto expected_mass=mass,expected_number=number;int flag;char error[256];
    exaero_removal_rate_tendencies(p.get(),nc,nl,nf,np,x.data(),rate.data(),modal?nr.data():nullptr,mass.data(),number.data(),error,&flag);
    assert(flag==0&&mass==expected_mass&&number==expected_number);
    View2D<const double> depth(dz.data(),nc,nl);
    removalVelocityTendencies(*p,depth,input,velocity,mv,nv);
    assert(host_test::near(mv(1,2,0),-.001*xv(1,2,0)));
    if(modal)assert(host_test::near(nv(1,2,0),-5e4));
    exaero_removal_velocity_tendencies(p.get(),nc,nl,nf,np,x.data(),vel.data(),dz.data(),mass.data(),number.data(),error,&flag);
    assert(flag==0&&host_test::near(mv(1,2,0),-.001*xv(1,2,0))&&x==original);
    std::fill(mass.begin(),mass.end(),77);std::fill(number.begin(),number.end(),88);
    host_test::rejects([&]{removalRateTendencies(*p,input,rates,modal?nullptr:&nrate,mv,nv);});
    rate.back()=std::numeric_limits<double>::quiet_NaN();
    exaero_removal_rate_tendencies(p.get(),nc,nl,nf,np,x.data(),rate.data(),modal?nr.data():nullptr,mass.data(),number.data(),error,&flag);
    assert(flag!=0&&mass.front()==77&&mass.back()==77&&number.back()==88);
    exaero_removal_velocity_tendencies(p.get(),nc,nl,nf,np+1,x.data(),vel.data(),dz.data(),mass.data(),number.data(),error,&flag);
    assert(flag!=0&&mass.front()==77&&number.back()==88);
    assert(supportsHostOperation(*p,HostOperation::RemovalRates)&&supportsHostOperation(*p,HostOperation::RemovalVelocity));
  }
}
