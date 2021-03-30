#include<iostream>
#include<random>
#include<vector>
#include <bits/stdc++.h>
#include <fstream>
#include "random_va.h"
#define COUNT 10


using namespace std;

template<typename T>
vector<T> s(vector<T> const &v, int m, int n) {
   auto first = v.begin() + m;
   auto last = v.begin() + n + 1;
   vector<T> vector(first, last);
   return vector;
}

template<typename T>
void show(vector<T> const &v) {
   for (auto i: v) {
      cout << i << ' ';
   }
   cout << '\n';
}

template<typename T>
ostream & operator <<(ostream & os,const vector<T>& v){
  for (auto s: v){
    os << s << endl;
  }
  return os;
};


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

vector<double> time_scale(const path<double>& paths)
{
  int N = paths.size();
  vector<double> times(N);
  for (int i = 0; i< N; i++){
    times[i] = paths[i].time;
  };
  return times;
};

template <typename T>
std::ostream & operator<<(std::ostream & o, path<T> const & p)
{
    for (auto const & st : p)
        o << st << std::endl;
    return o << std::endl;
}

template<typename TDistrib1, typename TDistrib2, typename TGen >
path<double> path_sim(state<double> ini_position,int N,TDistrib1 X, TDistrib2 Life_time,TGen& gen){
  double life_time = Life_time(gen);
  return X(gen,ini_position,life_time,N);
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

template<typename TDistrib1,typename TDistrib2, typename TGen >
struct Node* create(state<double> init_position,double T,int N,TDistrib1 X,TDistrib2 Life_time,TGen& gen)
   {
    Node* p;
   	path<double> ini_recur = path_sim(init_position,N,X,Life_time,gen);
    if (ini_recur[0].time > T ){
      return NULL;
    }
    else
    {
    p= newNode(ini_recur);
   	p->left  = create(ini_recur[N],T,N,X,Life_time,gen);
   	p->right = create(ini_recur[N],T,N,X,Life_time,gen);
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

void LastLeaves(Node* root,ofstream& of){
  if(root == NULL)
  {
    return;
  }
  if (root-> left != NULL | root -> right != NULL )
  {
      LastLeaves(root->right,of);
      LastLeaves(root->left,of);
  }
  else
  {
    of << root -> key << std::endl;
    // std::cout << root-> key << std::endl;
    // std::cout << endl;
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

path<double> read_vect_particles_T (const char *Nomfich , int N,double T)
{
  ifstream fileIn(Nomfich);
  double data;
  vector<double> myVec ;
  while (fileIn >> data)
  {
    myVec.push_back(data);
  }
  int N_pat = myVec.size()/(2*N);
  path<double> pat(N_pat);
  for (int i = 0 ; i<N_pat; i++){
    int k = 0;
    double min_temps = 10.;
    for (int j = 0;j<N;j++){
      min_temps = min(min_temps,abs(myVec[2*j +i*2*N]-T));
    };
    while( min_temps != abs(myVec[2*k+i*2*N] - T) ){
      k+=1;
    };
    pat[i].value = 0.5*(myVec[2*k+i*2*N+1] + myVec[2*k+i*2*N+3]) ;
    pat[i].time  = 0.5*(myVec[2*k+i*2*N] + myVec[2*k+i*2*N+2] );
  };
  return pat;
};



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

void showTrunks(Trunk *p)
{
    if (p == nullptr) {
        return;
    }

    showTrunks(p->prev);
    cout << p->str;
}


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

struct Brownien_geo{

  Brownien_geo(double r, double sigma): r(r),sigma(sigma),G(0,1) {};

  template<typename TGen>
  path<double> operator()(TGen & gen, state<double> ini_position, double t, int N)
  {
    path<double> trajectories(N+1);
    trajectories[0] = ini_position;
    double mu = r - sigma * sigma * 0.5 ;
    double delta_t = t/N;
    for (int i = 1 ; i< N+1 ; i++){
      trajectories[i].time = trajectories[i-1].time + delta_t;
      trajectories[i].value= trajectories[i-1].value * exp( mu * delta_t + sigma * sqrt(delta_t) * G(gen) ) ;
    }
    return trajectories;
  };

private:
  double sigma;
  double r;
  std::normal_distribution<> G;
};



double prod(const path<double>& Times_traj_T,std::function<double(double const &)> payoff){
  int N = Times_traj_T.size();
  double product = 1.;
  for (int i = 0; i< N; i++){
    product *= payoff(Times_traj_T[i].value);
  };
  return product;
};


template<typename TDistrib1,typename TDistrib2>
struct Branch_diffusion_simple{ 
  Branch_diffusion_simple(TDistrib1 X, TDistrib2 Life_time, std::function<double(double const & )> payoff, double Maturity, double spot):
  X(X),Life_time(Life_time),payoff(payoff),Maturity(Maturity),spot(spot){};

  template<typename TGen>
  double operator()(int N_simulations,TGen& gen){
    double somme_esperance = 0;
    state<double> ini_position = {0.,spot };
    int N = 50;
    for (int i = 0; i< N_simulations;i++)
    {
      Node* root;
      root=create(ini_position,Maturity,N,X,Life_time,gen);
      ofstream of1("last_leaves.txt");
      LastLeaves(root,of1);
      of1.close();
      path<double> Times_traj_T = read_vect_particles_T("last_leaves.txt",N,Maturity);
      somme_esperance += prod(Times_traj_T,payoff);
    };
    return somme_esperance / N_simulations;
  };

private:
  TDistrib1 X;
  TDistrib2 Life_time;
  std::function<double(double const & )> payoff;
  double Maturity;
  double spot;

};
