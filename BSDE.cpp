#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <iostream>
#include "BSDE2.h"

using namespace std;

double call(double x,double strike){
  if (x-strike>0){
    return x-strike;}
  else{
    return 0;
  }
};

double norm_pdf(const double& x) {
    return (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x);
}

double norm_cdf(const double& x) {
    double k = 1.0/(1.0 + 0.2316419*x);
    double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));

    if (x >= 0.0) {
        return (1.0 - (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x) * k_sum);
    } else {
        return 1.0 - norm_cdf(-x);
    }
}

double d_j(const int& j, const double& S, const double& K, const double& r, const double& v, const double& T) {
    return (log(S/K) + (r + (pow(-1,j-1))*0.5*v*v)*T)/(v*(pow(T,0.5)));
}

double call_price(const double& S, const double& K, const double& r, const double& v, const double& T) {
    return S * norm_cdf(d_j(1, S, K, r, v, T))-K*exp(-r*T) * norm_cdf(d_j(2, S, K, r, v, T));
}

int main(){


  // std::function<double(const double &,const double &, const double &, const double &)> test_fixe = [=] (double time, double x, double y, double z){ return  (y+2)/(y+1) ;};
  //
  // double solution = point_fixe<30>(1.,1.,0.,1.,0.,test_fixe);
  // std::cout << solution << std::endl;

  random_device rd;
  random_device rd1;
  auto seed = rd();
  mt19937_64 gen(seed);


  double S0 = 100.;
  double sigma = 0.2;
  double T = 1.;
  double r = 0.05;
  double strike = 100.;
  double mu = r - sigma * sigma / 2. ;

  std::function<double(const double &)> payoff = [=] (double x){return call(x,strike);};
  std::function<double(const double &,const double &, const double &, const double &)> generator = [=] (double time, double x, double y, double z){ return - r * y ;};

  BS_traj Bs_test(S0, sigma, T,mu);
  BSDE<BS_traj> bsde_test(Bs_test,T,payoff,generator);

  // 
  // std::cout << " methode explicite 1 = " << std::endl;
  // std::cout << bsde_test.resolution_explicite(gen) << std::endl;
  // std::cout << std::endl;
  // std::cout << " methode explicite 2 = " << std::endl;
  // std::cout << bsde_test.resolution_explicite_2(gen) << std::endl;
  // std::cout << std::endl;
  std::cout << " methode implicite = " << std::endl;
  std::cout << bsde_test.resolution_implicite(gen) << std::endl;

  std::cout << std::endl;
  std::cout << call_price(S0,strike,r,sigma,T) << std::endl;
  return 0;
}
