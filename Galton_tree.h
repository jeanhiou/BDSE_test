#include<iostream>
#include<random>
#include<vector>
#include <bits/stdc++.h>
#include <fstream>
#include "var.h"

#define COUNT 10

using namespace std;


template <typename T>
struct state {
    double time;
    T value;
};

template <typename T>
std::ostream & operator<<(std::ostream & o, state<T> const & s) {
    return o << s.time << "\t" << s.value;
}

template <typename T>
struct path : protected std::vector<state<T>> {
    using vec = std::vector<state<T>>;  // alias de nom
    using vec::vec;             // constructeur de la classe vector utilisable
    using vec::operator[];      // opérateur [] utilisable (public)
    using vec::begin;           // itérateurs utilisables (for-range loop)
    using vec::end;
    using vec::size;            // utile !
};

template <typename T>
std::ostream & operator<<(std::ostream & o, path<T> const & p) {
    for (auto const & st : p)
        o << st << std::endl;
    return o << std::endl;
}

template<typename TDistrib1,typename TDistrib2, typename TGen>
struct Galton_Watson{

  Galton_Watson(const TDistrib1 & X,const TDistrib2 & Z, const TGen &gen1 ): X(X),Z_0(Z),gen(gen1),membres(0) {};

  int Taille_population(){
    return membres;
  };

  void start(){
    membres += Z_0(gen);
  };

  void operator()(){
    int somme = 0;
    for (int i = 0;i< membres ;i++){
      somme += X(gen);
    };
    membres += somme;

  };

private:
  TDistrib1 X;
  TDistrib2 Z_0;
  TGen gen;
  int membres;
};

template <typename TDistrib1,typename TDistrib2,typename TGen>
inline Galton_Watson<TDistrib1,TDistrib2,TGen> Galton(const TDistrib1 & X,const TDistrib2 & Z, const TGen &gen ){ return Galton_Watson<TDistrib1,TDistrib2,TGen>(X,Z,gen);};


path<double> path_sim(state<double> ini_position,int N){
  random_device rd;
  mt19937_64 gen(rd());
  exponentiel_distribution E(1);
  std::normal_distribution<double> G(0,1);
  double left_time = E(gen);
  double delta_t = left_time/N ;
  path<double> Traj(N+1);
  Traj[0].time = ini_position.time;
  Traj[0].value = ini_position.value;
  for (int i = 0;i<N;i++){
    Traj[i+1].value = Traj[i].time + sqrt(delta_t) *G(gen);
    Traj[i+1].time  = Traj[i].time + delta_t;
  };
  return Traj;
};

struct Node
{
    path<double> key;
    Node* left;
    Node* right;
};

struct Node* newNode(path<double> data) {
  struct Node* tree = new (struct Node);
  tree->key = data;
  tree->left = NULL;
  tree->right = NULL;
  return tree;
}

struct Node* create(state<double> init_position,double T)
   {
    Node* p;
   	path<double> ini_recur = path_sim(init_position,4);
    if (ini_recur[4].time > T){
      return NULL;
    }
    else
    {
    p= newNode(ini_recur);
   	p->left=create(ini_recur[4],T);
   	p->right=create(ini_recur[4],T);
  }

  return p;
  };


void preorder(Node* t)		//address of root node is passed in t
{
	if(t!=NULL)
	{
		preorder(t->left);		//preorder traversal on left subtree
		preorder(t->right);		//preorder traversal om right subtree
	}
}

//
unsigned int getLeafCount(struct Node* node)
{
    if(node == NULL)
        return 0;
    else{
        return 1 + getLeafCount(node->left)+ getLeafCount(node->right);
      }
};


int maxDepth(Node* node)
{
    if (node == NULL)
        return 0;
    else
    {
        /* compute the depth of each subtree */
        int lDepth = maxDepth(node->left);
        int rDepth = maxDepth(node->right);

        /* use the larger one */
        if (lDepth > rDepth)
            return(lDepth + 1);
        else return(rDepth + 1);
    }
}

void LastLeaves(Node* root){
  if(root == NULL)
  {
    return;
  }
  if (root-> left != NULL | root -> right != NULL )
  {
      LastLeaves(root->right);
      LastLeaves(root->left);
  }
  else
  {
    std::cout << root-> key << std::endl;
    std::cout << endl;
    return;
  };
};

void SaveNodes(Node* root,ofstream& of){
  if(root == NULL)
  {
    return;
  }
  else
  {   of << root-> key << std::endl;
      of << endl;
      SaveNodes(root->right,of);
      SaveNodes(root->left,of);
  }
};

// //
// vector<path<double>> read_vect (const char *Nomfich)
// {
//   ifstream fileIn(Nomfich); // fileName is not a good name for a file!
//
//   double data;
//   vector<double> myVec ;
//   while (fileIn >> data)
//   {
//     myVec.push_back(data);
//   }
//   int N = myVec.size()/2;
//   vector<path<double>> mes_petites_particules(N);
//   for (int i = 0;i<N;i++){
//     mes_petites_particules[i] = {myVec[2*i],myVec[2*i+1]};
//   };
//   return mes_petites_particules;
// };
//
//
std::ostream & operator<<(std::ostream &o,const Node& n){
  o << n.key << std::endl;
  return o << std::endl;
};

struct Trunk
{
    Trunk *prev;
    string str;

    Trunk(Trunk *prev, string str)
    {
        this->prev = prev;
        this->str = str;
    }
};

// Helper function to print branches of the binary tree
void showTrunks(Trunk *p)
{
    if (p == nullptr) {
        return;
    }

    showTrunks(p->prev);
    cout << p->str;
}

// Recursive function to print a binary tree.
// It uses the inorder traversal.
void printTree(Node* root, Trunk *prev, bool isLeft)
{
    if (root == nullptr) {
        return;
    }

    string prev_str = "    ";
    Trunk *trunk = new Trunk(prev, prev_str);

    printTree(root->right, trunk, true);

    if (!prev) {
        trunk->str = "———";
    }
    else if (isLeft)
    {
        trunk->str = ".———";
        prev_str = "   |";
    }
    else {
        trunk->str = "`———";
        prev->str = prev_str;
    }

    showTrunks(trunk);
    cout << root->key << endl;

    if (prev) {
        prev->str = prev_str;
    }
    trunk->str = "   |";

    printTree(root->left, trunk, false);
}
