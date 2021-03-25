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

  State<double> ini_position = {0.,0.};
  Node *root;
  root=create(ini_position,3.);
  std::cout << std::endl;
  preorder(root);
  std::cout << " number of particules emmited = " << getLeafCount(root) << std::endl;
  std::cout << " max_particules_emises par le parent = " << maxDepth(root) << std::endl;
  ofstream of("particles.txt");
  SaveNodes(root,of);
  of.close();
  vector<State<double>> try_test = read_vect("particles.txt");
  std::cout << try_test << std::endl;

  return 0;
};
