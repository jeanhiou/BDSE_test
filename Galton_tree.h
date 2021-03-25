#include<iostream>
#include<random>
#include<vector>
#include <bits/stdc++.h>
#include <fstream>

using namespace std;

struct lien{
  int parent;
  int fils;
};

std::ostream & operator<<(std::ostream & o, lien const & p) {
    o << p.parent << p.fils << std::endl;
    return o << std::endl;
};

std::ostream & operator<<(std::ostream & o, vector<lien> const & p){
  for (auto const & st: p){
        o << st << std::endl;
}
      return o << std::endl;
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

template<typename T>
struct State{
  double time;
  T value;
};

template<typename T>
std::ostream & operator<<(std::ostream & o, State<T> const & p) {
    o << p.time <<" " <<  p.value << std::endl;
    return o << std::endl;
};

template <typename T>
std::ostream & operator<<(std::ostream & o, vector<T> const & p) {
    for (auto const & st : p)
        o << st << std::endl;
    return o;
};

void save_fichier_path(string filename, State<double> & p)
{
    ofstream of(filename); // ouverture du fichier
    of << p.time << " " << p.value << endl;           // fermature du fichier
};



State<double> path_sim(State<double> ini_position,int N){
  random_device rd;
  mt19937_64 gen(rd());
  exponentiel_distribution E(1);
  std::normal_distribution<double> G(0,1);
  double left_time = E(gen);
  double delta_t = left_time/N ;
  for (int i = 0;i<N;i++){
    ini_position.value += sqrt(delta_t) *G(gen);
  };
  return {ini_position.time + left_time, ini_position.value };
};

struct Node
{
    State<double> key;
    Node* left;
    Node* right;
};


struct Node* newNode(State<double> value)
  {
  	Node* n = new Node;
  	n->key = value;
  	n->left = NULL;
  	n->right = NULL;
  	return n;
  }

struct Node* create(State<double> init_position,double T)
 {
 	Node* p;
  State<double> ini_recur = path_sim(init_position,10);

  if (ini_recur.time > T){
    return NULL;
  };
  p=(Node*)malloc(sizeof(Node));
  p->key=ini_recur;

 	p->left=create(ini_recur,T);

 	p->right=create(ini_recur,T);
  init_position = ini_recur;

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


unsigned int getLeafCount(struct Node* node)
{
    if(node == NULL)
        return 0;
    if(node->left == NULL && node->right == NULL)
        return 1;
    else
        return getLeafCount(node->left)+
            getLeafCount(node->right);
}


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

void printNodes(Node* root)
{
    // return if the tree is empty
    if (root == nullptr) {
        return;
    }

    // print the root node
    cout << root->key << " ";

    // create two empty queues and enqueue root's left and
    // right child, respectively
    queue<Node*> q1, q2;
    q1.push(root->left);
    q2.push(root->right);

    // loop till queue is empty
    while (!q1.empty())
    {
        // calculate the total number of nodes at the current level
        int n = q1.size();

        // process every node of the current level
        while (n--)
        {
            // dequeue front node from the first queue and print it
            Node* x = q1.front();
            q1.pop();

            cout << x->key << " ";

            // enqueue left and right child of `x` to the first queue
            if (x->left) {
                q1.push(x->left);
            }

            if (x->right) {
                q1.push(x->right);
            }

            // dequeue front node from the second queue and print it
            Node* y = q2.front();
            q2.pop();

            cout << y->key << " ";

            // enqueue right and left child of `y` to the second queue
            if (y->right) {
                q2.push(y->right);
            }

            if (y->left) {
                q2.push(y->left);
            }
        }
    }
}


void SaveNodes(Node* root,ofstream& of)
{
    // return if the tree is empty
    if (root == nullptr) {
        return;
    }

    // print the root node
    of << root->key << " ";

    // create two empty queues and enqueue root's left and
    // right child, respectively
    queue<Node*> q1, q2;
    q1.push(root->left);
    q2.push(root->right);

    // loop till queue is empty
    while (!q1.empty())
    {
        // calculate the total number of nodes at the current level
        int n = q1.size();

        // process every node of the current level
        while (n--)
        {
            // dequeue front node from the first queue and print it
            Node* x = q1.front();
            q1.pop();

            of << x->key << " ";

            // enqueue left and right child of `x` to the first queue
            if (x->left) {
                q1.push(x->left);
            }

            if (x->right) {
                q1.push(x->right);
            }

            // dequeue front node from the second queue and print it
            Node* y = q2.front();
            q2.pop();

            of << y->key << " ";

            // enqueue right and left child of `y` to the second queue
            if (y->right) {
                q2.push(y->right);
            }

            if (y->left) {
                q2.push(y->left);
            }
        }
    }
}



vector<State<double>> read_vect (const char *Nomfich)
{
  ifstream fileIn(Nomfich); // fileName is not a good name for a file!

  double data;
  vector<double> myVec ;
  while (fileIn >> data)
  {
    myVec.push_back(data);
  }
  int N = myVec.size()/2;
  vector<State<double>> mes_petites_particules(N);
  for (int i = 0;i<N;i++){
    mes_petites_particules[i] = {myVec[2*i],myVec[2*i+1]};
  };
  return mes_petites_particules;
};


std::ostream & operator<<(std::ostream &o,const Node& n){
  o << n.key << std::endl;
  return o << std::endl;
};
