#include<ostream>
#include "N_ary_tree.h"
#include "Monte_carlo_template/random_va.h"
#include "Monte_carlo_template/timer.hpp"
#include "Monte_carlo_template/compose.hpp"
#include "Monte_carlo_template/quantile.hpp"
#include "Monte_carlo_template/mc.hpp"
#include "Monte_carlo_template/monte.hpp"

using namespace std;

double positive_part(double x){
  return (x>0) ? x : 0 ;
};

int main(){

  // random_generator //
  random_device rd;
  mt19937_64 gen(rd());

  double lambda = 0.05;
  double r = 0.;
  double sigma = 0.2;

  Brownien_geo B_geo(r,sigma);
  Brownien B_si ;

  // Life time for particles //
  exponentiel_distribution Life_time(lambda);

  // payoff //
  std::function<double(double const & )> payoff = [=] (double const &x){return positive_part(x - 1 ) ;};

  // maturity and spot value //
  double Maturity = 10.;
  double spot = 1.  ;
  //
  // vector<double> a = {0.0580,0.5,0.8164,0.,0.4043};  //
  // // vector<double> a_test = {0.,-1/3,1./3.,-1./3.};

  // // structure de la solution //
  vector<double> a = { 0, 1.};
  // vector<double> proba = random_p_choice(gen,a);
  vector<double> proba = {0.,1.};
  // vector<double> proba_uni = { 0. , 0.33 , 0.33 , 0.33};
  Children childs(proba);

  Branch_diffusion<Brownien_geo,exponentiel_distribution> B_cva(B_geo,Life_time,payoff,Maturity,spot,a,proba,childs);

  // Monte_carlo_simulation //
  // // std::normal_distribution<double> G(2,3);
  //
  int N_simulations = 100000;
  timer t;
  t.reset();
  auto result = monte_carlo(B_cva, gen, N_simulations);
  cout << t << "\t" << result << endl;
  cout << "intervalle de confiance taille" << result.ic_size() << std::endl;


  return 0;
}
