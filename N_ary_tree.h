#include<iostream>
#include<random>
#include<vector>
#include <bits/stdc++.h>
#include <fstream>

#define NEGATIVE_INFINITY - std::numeric_limits<double>::infinity();

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
    using vec::size;
    using vec::push_back;
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
      if (ini_recur[0].time > T){
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

void read_vect_particles_T (Node* root , int N,double T , path<double>& pat)
{
  if(root == NULL)
  {
    return;
  }
  int N_childs = root -> number_of_childs;
  if (root -> children == remplir<Node*>(N_childs,NULL)){
        path<double> data2 =  root -> key ;
        int k = 0;
        double min_temps = 100.;
          for (int j = 0;j<N;j++){
            min_temps = min(min_temps,abs(data2[j].time-T));
          };
          while( min_temps != abs(data2[k].time - T) ){
            k+=1;
          };
            pat.push_back({0.5*(data2[k].time + data2[k+1].time),0.5*(data2[k].value + data2[k+1].value)});
  }
  else
  {
    for (int i = 0; i< N_childs;i++)
    {
    if (root -> children[i] != NULL){
      read_vect_particles_T (root -> children[i],N,T,pat);
    };
  }
}
};

double prod(const path<double>& Times_traj_T,std::function<double(double const &)> payoff, const vector<int>& v, const vector<double>& ak, const vector<double>& pk){
  int N = Times_traj_T.size();
  int Ak = ak.size();
  double product_g = 1.;
  double product_ak = 1.;
  for (int i = 0; i< N; i++){
    product_g *= payoff(Times_traj_T[i].value);
  };
  for (int i = 0 ; i< Ak ; i++){
    if (ak[i] != 0){
    product_ak *= pow(ak[i]/pk[i],v[i]);
  }
  }
  return product_ak * product_g ;
};

void Number_of_k_descendants(Node* root,vector<int>& v){
  if (root == NULL){
    return;
  }
  else{
    v[root-> number_of_childs - 1 ]  = v[root -> number_of_childs - 1] + 1 ;
    for (int i = 0 ; i<root -> number_of_childs; i++){
      Number_of_k_descendants(( root -> children)[i],v);
    };
  }
}

double norme_infini(std::function<double(const double &)> f , double a, double b){
  double max_fx = NEGATIVE_INFINITY;
  double delta = 0.01;
  double x = a;
  for (double x = a ; x <= b; x += delta)
  {
    const float fx = abs(f(x)) ;
    if (fx > max_fx) // note any fx will be greater than NEGATIVE_INFINITY
      max_fx = fx;
  };
  return max_fx;
};

//cumulated_proba vector //
vector<double> proba_k(double norme_inf, vector<double> ak){
  double sum = 0;
  int N_k = ak.size();
  vector<double> p_k(N_k);
  for (int i = 0 ; i< N_k ; i++){
    sum +=   ak[i] * pow(norme_inf,i);}
  for (int i = 0 ; i< N_k ; i++){
    p_k[i] = ak[i] * pow(norme_inf,i)/sum;
  };
  return p_k ;
};

template<typename TGen>
vector<double> random_p_choice(TGen &gen,const vector<double>& ak){
  int N_taille = ak.size();
  vector<double> proba(N_taille);
  std::uniform_real_distribution<double> U;
  double somme = 0;
  for (int i = 0 ; i< N_taille-1 ; i++){
    if (ak [i] == 0.)
    {
      proba[i] = 0;
    }
    else
    {
      double alea = (1 - somme) * U(gen);
      proba[i] = alea;
      somme += alea;
    };
  };
  proba[N_taille-1] = 1 - somme;
  return proba;
};


template<typename TDistrib1,typename TDistrib2>
struct Branch_diffusion{
  Branch_diffusion(TDistrib1 X, TDistrib2 Life_time, std::function<double(double const & )> payoff, double Maturity, double spot,vector<double> ak,vector<double> proba, Children childs):
  X(X),Life_time(Life_time),payoff(payoff),Maturity(Maturity),spot(spot),ak(ak),proba(proba), childs(childs) {};

  template<typename TGen>
  double operator()(TGen& gen)
  {
    state<double> ini_position = {0.,spot };
    int N = 40;

    Node* root = create(ini_position,Maturity,N,X,Life_time,childs,gen);

    vector<int> num_desc(proba.size());
    Number_of_k_descendants(root,num_desc);

    path<double> Times_traj_T ;
    read_vect_particles_T (root ,N,Maturity,Times_traj_T);
    return prod(Times_traj_T,payoff,num_desc,ak,proba);
  };

private:
  TDistrib1 X;
  TDistrib2 Life_time;
  std::function<double(double const & )> payoff;
  double Maturity;
  double spot;
  vector<double> ak;
  Children childs;
  vector<double> proba;

};
