#include<ostream>
#include "N_ary_tree.h"
#include "random_va.h"

using namespace std;

int main(){

  vector<double> proba =  {0.2 , 0.5, 0.3};

  Children Childs( cumul_proba(proba));

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

  std::function<double(double const & )> payoff = [=] (double const &x){return x ;};
  double Maturity = 1.;
  double spot = S0  ;

  int N_simulations = 1000;



  Node* root = create(ini_position,Maturity,N,B_si,Life_time,Childs,gen);
  ofstream of("ultra_particles.txt");
  SaveNodes(root,of);
  of.close();


  return 0;
}
