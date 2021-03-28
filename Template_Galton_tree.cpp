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

  Node* root;
  root=create(ini_position,1.,30);
  std::cout << endl;
  std::cout << getLeafCount(root) << std::endl;
  std::cout << endl;

  ofstream of("particles.txt");
  SaveNodes(root,of);
  of.close();

  ofstream of1("last_leaves.txt");
  LastLeaves(root,of1);
  of1.close();

  vector<state<double>> read_vec_T = read_vect_particles_T("last_leaves.txt",30,1.);
  std::cout << " le vecteur des temps = " << endl;
  std::cout << read_vec_T << std::endl;

  // vector<path<double>> try_test = read_vect ("particles.txt",20);
  // std::cout << try_test[0] << std::endl;


  return 0;
};
