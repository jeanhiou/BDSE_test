#include<iostream>
#include<cmath>

using namespace std;


struct Bernouilli_distribution
{
  Bernouilli_distribution(double p): p(p),U(0,1){};

  template<typename TGen>
  double operator()(TGen &gen){
    if (U(gen) < p)
    {
      return 0;
    }
    else
    {
      return 1;
    }
  }

private:
  double p;
  uniform_real_distribution<> U;
};

struct Binomial_distribution
{
  Binomial_distribution(double p, int n):N(n),B(p){};

  template<typename TGen>
  double operator()(TGen &gen){
    double somme = 0;
    for (int i = 0 ; i <N ; i++){
      somme += B(gen);
    };
    return somme;
  }

private:
  int N ;
  Bernouilli_distribution B;
};

struct exponentiel_distribution
{
  exponentiel_distribution(double lambda): lambda(lambda), U(0,1) {};

  template<typename TGen>
  double operator()(TGen &gen){
    return -log (U(gen))/lambda;
  };

private:
  double lambda;
  uniform_real_distribution<> U;
};

struct pareto_distribution{
  pareto_distribution(double theta): theta(theta),U(0,1) {};

  template<typename TGen>
  double operator()(TGen &gen){
    return pow(U(gen),-1/theta);
  };
private:
  double theta;
  uniform_real_distribution<> U;
};


template <typename TDistrib1, typename TDistrib2>
struct Mixture_loi{
  Mixture_loi(TDistrib1 & distrib1, TDistrib2 & distrib2,double alpha1 = 0.5 ,double alpha2 = 0.5 ):alpha1(alpha1),alpha2(alpha2),distrib1(distrib1),distrib2(distrib2){};

  template<typename TGen>
  double operator()(TGen &gen){
    return alpha1 * distrib1(gen) + alpha2 * distrib2(gen);
  }

private:
  double alpha1;
  double alpha2;
  TDistrib1 distrib1;
  TDistrib2 distrib2;
};


struct Normalized_Poisson_process{
  Normalized_Poisson_process(double t):t(t), U(0,1){};

  template<typename TGen>
  int operator()(TGen & gen){
    double somme = U(gen);
    int k = 0;
    while (somme > exp(-t)){
      somme *= U(gen);
      k +=1;
    };
    return k;
  }

private:
  uniform_real_distribution<> U;
  double t;
};
