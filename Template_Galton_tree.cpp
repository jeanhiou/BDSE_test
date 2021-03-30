#include<iostream>
#include "Template_Galton_tree.h"

using namespace std;

double positive_part(const double& x){
  if (x>=0){
    return x;
  }
  else{
    return 0;
  };
};

int main(){

  random_device rd;
  mt19937_64 gen(rd());


  double lambda = 2.;
  double r = 0.05;
  double sigma = 0.2;
  double S0 = 1.;
  state<double> ini_position = { 0.,0.};
  int N = 20;
  Brownien_geo B_geo(r,sigma);
  Brownien B_si;

  exponentiel_distribution Life_time(lambda);

  std::function<double(double const & )> payoff = [=] (double const &x){return positive_part(x) ;};
  double Maturity = 1.;
  double spot = S0  ;

  int N_simulations = 1000;
  Branch_diffusion_simple<Brownien,exponentiel_distribution> B_test(B_si,Life_time,payoff,Maturity,spot);


  timer t;
  t.reset();
  auto result = monte_carlo(B_test, gen, N_simulations);
  cout << t << "\t" << result << endl;



  return 0;
};
