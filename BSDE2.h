#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath>
#include <iostream>
#include <random>
#include<functional>
#include <eigen3/Eigen/Dense>

using namespace std;

inline unsigned p2(unsigned n) { unsigned u = 1; return u <<= n; }

Eigen::MatrixXd base_phi( int m, Eigen::VectorXd X_t){
  int N = X_t.size();
  Eigen::MatrixXd base_phi_m(N,m);
  for (int i = 0;i<m;i++){
    for (int j = 0;j<N;j++){
      base_phi_m(j,i) = std::hermite(i,X_t(j));
    }
  }
  return base_phi_m;
};

Eigen::VectorXd conditional_expect(int m, Eigen::VectorXd X_t, Eigen::VectorXd X_t_plus){
  int N = X_t.size();
  Eigen::MatrixXd J = base_phi(m,X_t);
  Eigen::MatrixXd A = J.transpose() * J ;
  Eigen::VectorXd b = J.transpose() * X_t_plus;
  Eigen::VectorXd alpha =  A.ldlt().solve(b);
  Eigen::VectorXd Conditional = Eigen::VectorXd::Zero(N);
  for (int i = 0 ; i<N; i++){
      Conditional(i) = alpha.dot(J.row(i));
  };
  return Conditional;
};

template<int iter_max>
double point_fixe(const double & x, const double & z, const double & time, const double & delta_t, const double & conditional_expe , std::function<double(const double &,const double &, const double& , const double &)> generator){
  double y_0 = x  ;
  for (int i = 0; i<iter_max;i++){
    y_0 = conditional_expe + generator(time,x,y_0,z)*delta_t  ;
  };
  return y_0;
};

Eigen::VectorXd Vectorize_f( function<double(const double &,const double &,const double &,const double & )> f, Eigen::VectorXd Y, Eigen::VectorXd X, Eigen::VectorXd Z, double time, double delta_t){
    int N = X.size();
    Eigen::VectorXd Vector_f(N);
    for (int i = 0; i< N; i++){
      Vector_f(i) = f(time,X(i),Y(i),Z(i));
    };
    return Vector_f*delta_t;
  };


struct BS_traj{
  BS_traj(double s0, double sigma, double T,double mu): s0(s0),sigma(sigma),T(T),mu(mu),G(0,1),N_temps(50),N_traj(100000) {};

  template<typename TGen>
  Eigen::MatrixXd operator()(TGen & gen)
  {
    double delta_t=T/N_temps;
    Eigen::MatrixXd Result = Eigen::MatrixXd::Zero(N_traj,N_temps+1);
    for (int i = 0 ; i < N_traj ; i++){
      Result(i,0) = s0;
    }
    for (int i= 0 ; i < N_traj; i++){
        for (int j=1;j<N_temps+1;j++){
                Result(i,j)= Result(i,j-1) * exp( mu * delta_t + sigma * sqrt(delta_t) * G(gen)) ;}
    };
    return Result;
  };


private:
  double s0;
  double sigma;
  double T;
  int N_temps;
  int N_traj;
  double mu;
  std::normal_distribution<> G;
};

template<typename TDistrib>
struct BSDE{

  BSDE( TDistrib X_target, double T, std::function<double(const double &)> payoff,std::function<double(const double &,const double &, const double &, const double &)> generator):
  X_target(X_target),T(T),payoff(payoff),generator(generator) {};

  template<typename TGen>
  double resolution_implicite( TGen & gen){

    std::normal_distribution<double> G(0,1);

    Eigen::MatrixXd X_0_t = X_target(gen);

    int N_temps = X_0_t.row(0).size();
    int N_traj = X_0_t.col(0).size();
    double delta_t = T/N_temps;

    Eigen::VectorXd Y_T(N_traj);
    Eigen::VectorXd Y_new(N_traj);

    Eigen::VectorXd Z_T(N_traj);
    Eigen::VectorXd Z_new(N_traj);

    Eigen::VectorXd Y_delta_W(N_traj) ;

    Eigen::VectorXd Conditional(N_traj);

    for (int i = 0 ; i< N_traj; i++ ){
      Y_T(i) = payoff(X_0_t.col(N_temps-1)(i));
    };

    for (int i = 0;i<N_traj;i++){
      Y_delta_W(i) = sqrt(delta_t)*G(gen) * Y_T(i);
    };

    Z_T = conditional_expect(3,X_0_t.col(N_temps-2),Y_delta_W);
    Conditional = conditional_expect(3,X_0_t.col(N_temps-2),Y_T);

    for (int i = 0; i< N_traj;i++){
      Y_new(i) = point_fixe<30>(X_0_t.col(N_temps-2)(i),Z_T(i), (N_temps-1)*delta_t,delta_t,Conditional(i),generator);
    };

    /* demarrage de l'algo */

    for (int i = N_temps-2;i>0;i--){
      for (int k = 0;k < N_traj;k++){
        Y_delta_W(k) = sqrt(delta_t)*G(gen) * Y_new(k);
      };
      Z_new = conditional_expect(3,X_0_t.col(i),Y_delta_W);
      Conditional = conditional_expect(3,X_0_t.col(i),Y_new);

      for (int j = 0;j<N_traj;j++){
        Y_new(j) = point_fixe<30>(X_0_t.col(i)(j),Z_new(j), j *delta_t,delta_t,Conditional(j),generator);
      };
    }
    return Y_new.sum()/N_traj;
  };

  template<typename TGen>
  double resolution_explicite( TGen & gen){

    std::normal_distribution<double> G(0,1);

    Eigen::MatrixXd X_0_t = X_target(gen);

    int N_temps = X_0_t.row(0).size();
    int N_traj = X_0_t.col(0).size();

    double delta_t = T/N_temps;

    Eigen::VectorXd Y_new(N_traj);

    Eigen::VectorXd Z_new(N_traj);

    Eigen::VectorXd Y_delta_W(N_traj);

    for (int i = 0 ; i< N_traj; i++ ){
      Y_new(i) = payoff(X_0_t.col(N_temps-1)(i));
    };

    for (int i = 0;i<N_traj;i++){
      Y_delta_W(i) = sqrt(delta_t)*G(gen) * Y_new(i);
    };

    Z_new = conditional_expect(3,X_0_t.col(N_temps-2),Y_delta_W);

    Y_new = conditional_expect(3,X_0_t.col(N_temps-2), Y_new + Vectorize_f( generator, Y_new, X_0_t.col(N_temps-1), Z_new , N_temps * delta_t , delta_t));

    for (int i = N_temps-2;i>0;i--){
      for (int k = 0;k < N_traj;k++){
        Y_delta_W(k) = sqrt(delta_t)*G(gen) * Y_new(k);
      };
      Y_new = conditional_expect(3,X_0_t.col(i-1),Y_new + Vectorize_f(generator,Y_new,X_0_t.col(i),Z_new, (i+1) * delta_t,delta_t));
      Z_new = conditional_expect(3,X_0_t.col(i-1),Y_delta_W);
    };
    return Y_new.sum()/N_traj;
  };

  template<typename TGen>
  double resolution_explicite_2( TGen & gen){

    std::normal_distribution<double> G(0,1);

    Eigen::MatrixXd X_0_t = X_target(gen);

    int N_temps = X_0_t.row(0).size();
    int N_traj = X_0_t.col(0).size();

    double delta_t = T/N_temps;

    Eigen::VectorXd Y_new(N_traj);

    Eigen::VectorXd Z_new(N_traj);

    Eigen::VectorXd Y_delta_W(N_traj);

    Eigen::VectorXd Conditional(N_traj);

    for (int i = 0 ; i< N_traj; i++ ){
      Y_new(i) = payoff(X_0_t.col(N_temps-1)(i));
    };

    for (int i = 0;i<N_traj;i++){
      Y_delta_W(i) = sqrt(delta_t)*G(gen) * Y_new(i);
    };

    Z_new = conditional_expect(3,X_0_t.col(N_temps-2),Y_delta_W);
    Conditional = conditional_expect(3,X_0_t.col(N_temps-2),Y_new);
    Y_new = Conditional + conditional_expect(3,X_0_t.col(N_temps-2),Vectorize_f( generator,Conditional, X_0_t.col(N_temps-1), Z_new , (N_temps - 1) * delta_t, delta_t));

    for (int i = N_temps-2;i>0;i--){
        for (int k = 0;k < N_traj;k++){
        Y_delta_W(k) = sqrt(delta_t)*G(gen) * Y_new(k);
      };
      Z_new = conditional_expect(3,X_0_t.col(i-1),Y_delta_W);
      Conditional = conditional_expect(3,X_0_t.col(i-1),Y_new);
      Y_new = Conditional + conditional_expect(3,X_0_t.col(i-1),Vectorize_f(generator,Conditional,X_0_t.col(i),Z_new,i * delta_t,delta_t));
    };
    return Y_new.sum()/N_traj;
  };

private:
  TDistrib X_target;
  double T;
  std::function<double(const double &)> payoff;
  std::function<double(const double &,const double &, const double &, const double &)> generator;
};
