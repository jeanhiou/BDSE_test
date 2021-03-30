#include<iostream>
#include "Template_Galton_tree.h"

using namespace std;

int main(){

  random_device rd;
  mt19937_64 gen(rd());


  double lambda = 2.;
  double r = 0.05;
  double sigma = 0.2;
  double S0 = 1.;
  state<double> ini_position = { 0.,S0};
  int N = 20;
  double t = 1.;
  Brownien_geo B_geo(r,sigma);

  exponentiel_distribution Life_time(lambda);

  std::function<double(double const & )> payoff = [=] (double const &x){return x*x ;};
  double Maturity = 2.;
  double spot = S0  ;

  int N_simulations = 10000;
  std::cout << B_geo(gen,ini_position,1.,N) << std::endl;

  Node* root;
  root=create(ini_position,Maturity,N,B_geo,Life_time,gen);
  ofstream of1("particles.txt");
  SaveNodes(root,of1);
  of1.close();



  // Branch_diffusion_simple<Brownien_geo,exponentiel_distribution> B_test(B_geo,Life_time,payoff,Maturity,spot);
  // std::cout << B_test(N_simulations,gen) << std::endl;



  return 0;
};
