//g++ main.cpp -o app

#include <math.h>
#include <vector>
#include <iostream>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
using namespace std;

void print(std::vector <double> const &a) {
   for(int i=0; i < a.size(); i++){
      std::cout << a.at(i) <<' ';
    }
    cout << endl;
}

double Average ( vector<double> v )
{
        double total;
        for ( int i=0; i < v.size(); i++)
        {
            total += v[i];
        }
        return ( total / v.size());
}

double Variance(std::vector<double> samples)
{
     double size = samples.size();

     double variance = 0;
     double t = samples[0];
     for ( int i = 1; i < size; i++)
     {
          t += samples[i];
          double diff = ((i + 1) * samples[i]) - t;
          variance += (diff * diff) / ((i + 1.0) *i);
     }

     return variance / (size - 1);
}

double StandardDeviation(std::vector<double> samples)
{
     return sqrt(Variance(samples));
}



void Standardization(std::vector<double> &v){
  double avg = Average(v);
  double dev = StandardDeviation(v);
  for ( int i=0; i < v.size(); i++)
  {
      v[i] = (v[i]-avg)/dev;
  }

}

double cov(std::vector<double> vA, std::vector<double> vB){
  double total;
  for ( int i=0; i < vA.size(); i++)
  {
      total += vA[i]*vB[i];
  }
  return ( total / vA.size());
}

std::vector<std::vector<double>> Covarience(std::vector<double> v1, std::vector<double> v2, std::vector<double> v3, std::vector<double> v4){
  std::vector<double> vA{cov(v1,v1),cov(v1,v2),cov(v1,v3),cov(v1,v4)};
  std::vector<double> vB{cov(v2,v1),cov(v2,v2),cov(v2,v3),cov(v2,v4)};
  std::vector<double> vC{cov(v3,v1),cov(v3,v2),cov(v3,v3),cov(v3,v4)};
  std::vector<double> vD{cov(v4,v1),cov(v4,v2),cov(v4,v3),cov(v4,v4)};
  std::vector<std::vector<double>> covMatrix{vA,vB,vC,vD};
  return covMatrix;
}

int main() {
  vector<double> f1{1,5,1,5,8};
  vector<double> f2{2,5,4,3,1};
  vector<double> f3{3,6,2,2,2};
  vector<double> f4{4,7,3,1,2};

  cout << Average(f1)<< " "<< Average(f2)<< " "<< Average(f3)<< " "<< Average(f4)<< " "<<endl;
  cout << StandardDeviation(f1)<< " "<<StandardDeviation(f2)<< " "<<StandardDeviation(f3)<< " "<<StandardDeviation(f4)<< " "<<endl;
  Standardization(f1);
  Standardization(f2);
  Standardization(f3);
  Standardization(f4);
  cout << "before"<<endl;
  print(f1);
  print(f2);
  print(f3);
  print(f4);
  cout << "before 2"<<endl;
  std::vector<std::vector<double>> covMat = Covarience(f1,f2,f3,f4);
  for (int i = 0; i<covMat.size(); i++){
    print(covMat[i]);
  }
  Eigen::Vector4d v1(covMat[0].data());
  cout << v1 << endl;
  Eigen::Vector4d v2(covMat[1].data());
  cout << v2 << endl;
  Eigen::Vector4d v3(covMat[2].data());
  cout << v3 << endl;
  Eigen::Vector4d v4(covMat[3].data());
  cout << v4 << endl;
  Eigen::Matrix4d m;

  m.col(0) << v1;
  m.col(1) << v2;
  m.col(2) << v3;
  m.col(3) << v4;
  cout << m << endl;
  Eigen::EigenSolver<Eigen::MatrixXd> eigensolver;
  eigensolver.compute(m);
  Eigen::VectorXd eigen_values = eigensolver.eigenvalues().real();
  Eigen::MatrixXd eigen_vectors = eigensolver.eigenvectors().real();

  std::cout<< "eigen_vectors" <<endl << eigen_vectors.real() << std::endl;
  std::cout<< "eigen_values" << endl << eigen_values.real() << std::endl;
  Eigen::MatrixXd toChange(5,4);
  Eigen::VectorXd o1(5);
  Eigen::VectorXd o2(5);
  Eigen::VectorXd o3(5);
  Eigen::VectorXd o4(5);
  for (int i = 0; i<f1.size(); i++){
    o1[i] = f1[i];
    o2[i] = f2[i];
    o3[i] = f3[i];
    o4[i] = f4[i];
  }
  toChange.col(0) << o1;
  toChange.col(1) << o2;
  toChange.col(2) << o3;
  toChange.col(3) << o4;
  cout<< toChange << endl;
  Eigen::MatrixXd trans(4,2);
  trans.col(0) <<eigen_vectors.col(0);
  trans.col(1)<< eigen_vectors.col(1);

  cout<< "start" <<endl<< trans <<endl << "end" << endl;

  cout<< "Final" <<endl<< toChange*trans <<endl << "Result" << endl;


  return 0;
}
