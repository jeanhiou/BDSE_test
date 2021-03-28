#include<iostream>
#include "Galton_tree.h"

using namespace std;

int main(){

  random_device rd;
  mt19937_64 gen(rd());

  double lambda = 2.;
  Normalized_Poisson_process P(lambda) ;
  Normalized_Poisson_process P2(3.);
  auto m = Galton(P,P2,gen);
  exponentiel_distribution E(1);

  state<double> ini_position = {0.,0.};

  path<double> path = path_sim(ini_position,10);

  Node* root;
  root=create(ini_position,3.);
  std::cout << getLeafCount(root) << std::endl;

  ofstream of("particles.txt");
  SaveNodes(root,of);
  of.close();

  // vector<State<double>> try_test = read_vect("particles.txt");
  // std::cout << try_test << std::endl;


  return 0;
};
