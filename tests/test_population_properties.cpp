#ifdef NDEBUG
#undef NDEBUG
#endif
#include "host_coupling_fixture.hpp"
#include <exaero/HostCoupling.hpp>
#include <exaero/ModalRemovalMath.hpp>
#include <exaero/PopulationProperties.hpp>
#include <exaero/PopulationProperties_C.h>
#include <limits>
int main(){
  host_test::Runtime runtime;using namespace exaero;
  for(bool modal:{false,true}){
    auto p=host_test::make(modal);host_test::Grid g(modal);
    const int count=int(PopulationProperty::Count),n=g.nc*g.nl*g.np;
    std::vector<double>values(n*count);std::vector<int>available(n*count),status(n);
    View4D<double> v(values.data(),g.nc,g.nl,g.np,count);
    View4D<int> a(available.data(),g.nc,g.nl,g.np,count);View3D<int> st(status.data(),g.nc,g.nl,g.np);
    auto get=[&](PopulationProperty prop){return v(0,0,0,int(prop));};
    if(modal){View3D<double>x(g.state.data(),g.nc,g.nl,g.nf);for(int k=0;k<g.nl;++k)for(int c=0;c<g.nc;++c){x(c,k,1)=0;x(c,k,2)=x(c,k,0)/1000/modalNumberVolume(1e-7,1.5);}}
    populationProperties(*p,g.env(),g.input(),v,a,st);
    const double volume=g.state[0]/1000,rg=1e-7,l=std::log(1.5),number=volume/modalNumberVolume(rg,1.5);
    assert(host_test::near(get(PopulationProperty::DryVolume),volume));
    assert(host_test::near(get(PopulationProperty::Number),number));
    assert(host_test::near(get(PopulationProperty::DryNumberMedianRadius),rg));
    assert(host_test::near(get(PopulationProperty::DrySurfaceArea),4*std::acos(-1.)*number*rg*rg*std::exp(2*l*l)));
    assert(host_test::near(get(PopulationProperty::EffectiveDryRadius),rg*std::exp(2.5*l*l)));
    assert(a(0,0,0,int(PopulationProperty::LowerDryRadius))==int(PropertyAvailability::Unavailable));
    if(!modal)assert(host_test::near(get(PopulationProperty::DryScavengingDiameter),8e-7));
    else {
      View3D<double>x(g.state.data(),g.nc,g.nl,g.nf);
      const double total=4e-12;
      x(0,0,0)=.75*total*1000;x(0,0,1)=.25*total*2000;
      populationProperties(*p,g.env(),g.input(),v,a,st);
      assert(host_test::near(get(PopulationProperty::Kappa),.45));const auto dry=get(PopulationProperty::DryNumberMedianRadius),wet=get(PopulationProperty::WetNumberMedianRadius);
      x(0,0,0)=.25*total*1000;x(0,0,1)=.75*total*2000;
      populationProperties(*p,g.env(),g.input(),v,a,st);
      assert(host_test::near(get(PopulationProperty::Kappa),.15));assert(host_test::near(get(PopulationProperty::DryNumberMedianRadius),dry)&&get(PopulationProperty::WetNumberMedianRadius)<wet);
    }
    const auto expected=values;int flag;char error[256];
    exaero_population_properties(p.get(),g.nc,g.nl,g.nf,g.np,g.temp.data(),g.rh.data(),g.state.data(),values.data(),available.data(),status.data(),error,&flag);
    assert(flag==0&&values==expected);
    std::fill(g.state.begin(),g.state.end(),0);populationProperties(*p,g.env(),g.input(),v,a,st);
    assert(get(PopulationProperty::Mass)==0&&get(PopulationProperty::Number)==0&&status[0]==int(PopulationStatus::Empty));
    if(modal)assert(a(0,0,0,int(PopulationProperty::Kappa))==int(PropertyAvailability::Unavailable));
    std::fill(values.begin(),values.end(),77);std::fill(available.begin(),available.end(),88);std::fill(status.begin(),status.end(),99);g.rh.back()=std::numeric_limits<double>::quiet_NaN();
    exaero_population_properties(p.get(),g.nc,g.nl,g.nf,g.np,g.temp.data(),g.rh.data(),g.state.data(),values.data(),available.data(),status.data(),error,&flag);
    assert(flag!=0&&values.front()==77&&values.back()==77&&available.back()==88&&status.back()==99);
    assert(supportsHostOperation(*p,HostOperation::PopulationProperties));
  }
}
