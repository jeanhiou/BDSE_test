#include<iostream>
#include<random>
#include<vector>
#include <bits/stdc++.h>
#include <fstream>

using namespace std;

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

vector<double> cumul_proba(vector<double> proba){
  int N = proba.size();
  for (int i = 1 ; i< N ; i++){
    proba[i] += proba[i-1];
  };
  return proba;
};

template <typename T>
struct path : protected std::vector<state<T>> {
    using vec = std::vector<state<T>>;  // alias de nom
    using vec::vec;             // constructeur de la classe vector utilisable
    using vec::operator[];      // opérateur [] utilisable (public)
    using vec::begin;           // itérateurs utilisables (for-range loop)
    using vec::end;
    using vec::size;            // utile !
};

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
  double r;
  double sigma;
  std::normal_distribution<> G;

};

struct Brownien{

    Brownien(): G(0,1) {};

    template<typename TGen>
    path<double> operator()(TGen & gen, state<double> ini_position, double t, int N)
    {
      path<double> trajectories(N+1);
      trajectories[0] = ini_position;
      double delta_t = t/N;
      for (int i = 1 ; i< N+1 ; i++){
        trajectories[i].time = trajectories[i-1].time + delta_t ;
        trajectories[i].value= trajectories[i-1].value + sqrt(delta_t) * G(gen) ;
      }
      return trajectories;
    };

private:
  std::normal_distribution<> G;
};

struct Children{
  Children(vector<double> probs) : probs(probs), U(0,1) {};

  template<typename TGen>
  int operator()(TGen &gen){
    double u = U(gen);
    int k = 0;
    while (u>probs[k]){
      k+=1;
    }
    return k+1;
  };
  private:
    vector<double> probs;
    std::uniform_real_distribution<double> U;
};

template<typename TDistrib1, typename TDistrib2, typename TGen >
path<double> path_sim(state<double> ini_position,int N,TDistrib1 X, TDistrib2 Life_time,TGen& gen){
  double life_time = Life_time(gen);
  return X(gen,ini_position,life_time,N);
};

struct Node
  {
      path<double> key;
      int number_of_childs;
      vector<Node*> children;
  };


template<typename T>
vector<T> remplir(int N, T value)
  {
    vector<T> vect1(N);
    fill(vect1.begin(), vect1.end(), value);
    return vect1;
  };

struct Node* newNode(path<double> data,int child) {
    struct Node* tree = new (struct Node);
    tree->key = data;
    tree->number_of_childs = child;
    tree->children = remplir<Node*>(child,NULL);
    return tree;
  }

template<typename TDistrib1,typename TDistrib2, typename TGen >
struct Node* create(state<double> init_position,double T,int N,TDistrib1 X,TDistrib2 Life_time,Children& enfants,TGen& gen)
     {
      Node* p;
     	path<double> ini_recur = path_sim(init_position,N,X,Life_time,gen);
      if (ini_recur[0].time > T ){
        return NULL;
      }
      else
      {
      int N_enfants = enfants(gen);
      p= newNode(ini_recur,N_enfants);
      for (int i = 0 ; i<N_enfants ; i++){
        p -> children[i] = create(ini_recur[N],T,N,X,Life_time,enfants,gen);
      };
    };
    return p;
  };

template <typename T>
std::ostream & operator<<(std::ostream & o, path<T> const & p)
  {
      for (auto const & st : p)
          o << st << std::endl;
      return o << std::endl;
  }


void SaveNodes(Node* root,ofstream& of){
    if(root == NULL)
    {
      return;
    }
    else
    {   of << root-> key << std::endl;
        of << endl;
        int n = root-> number_of_childs;
        for (int i = 0;i<n;i++){
          SaveNodes((root->children)[i],of);
        };
    };
  };

vector<int> indices_null(Node* root){
  vector<int> indi_null;
  for (int i = 0 ; i< root -> number_of_childs;i++){
    if ((root -> children)[i] == NULL){
      indi_null.push_back(i);
    };
  };
  return indi_null;
};

vector<int> indices_no_null(Node* root){
  vector<int> indi_null;
  int N_childs = root->number_of_childs;
  for (int i = 0 ; i< N_childs ;i++){
    if ((root->children)[i] != NULL){
      indi_null.push_back(i);
    };
  };
  return indi_null;
};


void LastLeaves(Node* root, ofstream& of){
    if(root == NULL)
    {
      return;
    }
    int N_childs = root -> number_of_childs;
    if (root -> children == remplir<Node*>(N_childs,NULL)){
      of << root -> key << endl;
    }
    else
    {
      for (int i = 0; i< N_childs;i++)
      {
      if (root -> children[i] != NULL){
        LastLeaves(root -> children[i],of);
      };
    }
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

double prod(const path<double>& Times_traj_T,std::function<double(double const &)> payoff){
  int N = Times_traj_T.size();
  std::cout << "nombre de particules = " << N << std::endl;
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
  double operator()(TGen& gen)
  {
    state<double> ini_position = {0.,spot };
    int N = 20;
    Node* root;
    root=create(ini_position,Maturity,N,X,Life_time,gen);
    ofstream of("particles.txt");
    SaveNodes(root,of);
    of.close();
    ofstream of1("last_leaves.txt");
    LastLeaves(root,of1);
    of1.close();
    path<double> Times_traj_T = read_vect_particles_T("last_leaves.txt",N,Maturity);
    return prod(Times_traj_T,payoff);
  };

private:
  TDistrib1 X;
  TDistrib2 Life_time;
  std::function<double(double const & )> payoff;
  double Maturity;
  double spot;

};
