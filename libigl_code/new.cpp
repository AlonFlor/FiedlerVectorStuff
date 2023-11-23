#include <igl/heat_geodesics.h>
#include <igl/massmatrix.h>
#include <igl/writePLY.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Sparse>
#include <iostream>
#include <chrono>

using namespace std;

double get_geodesic_distances(string model_path){
  using namespace Eigen;
  using namespace igl;

  MatrixXd V;
  MatrixXi F;

  if(!read_triangle_mesh(model_path.data(),V,F))
  {
    cout<<"failed to load mesh"<<endl;
  }

  /*VectorXi gamma = {{1}};
  //VectorXi gamma(V.rows());//{{1, 2}}
  for(int i=0; i<V.rows(); ++i){
    gamma[i] = i;
  }*/
  VectorXd D;

  //time start
  const auto start_time {chrono::steady_clock::now()};
  
  HeatGeodesicsData<double> data;
  heat_geodesics_precompute(V,F,data);

  /*//for(int i=0; i<V.rows(); ++i){
    int i = 0;
    VectorXi gamma (1,1);
    gamma << i;
    D = VectorXd::Zero(data.Grad.cols(),1);
    heat_geodesics_solve(data,gamma,D);
  //  break;
  //}*/

  //time end
  const auto end_time {chrono::steady_clock::now()};
  const chrono::duration<double> elapsed_seconds{end_time - start_time};
  double time_elapsed = elapsed_seconds.count();

  //cout<<"results:"<<endl<<D<<endl;
  //cout<<endl<<endl<<time_elapsed<<endl;
  //cout<<"sizes:"<<endl<<data.Grad.cols()<<endl<<V.rows()<<endl;
  return time_elapsed;
}

int main(int argc, char * argv[])
{
  //cout << get_geodesic_distances("elephant.obj") << endl;
  //cout << get_geodesic_distances("reordered_elephant.obj") << endl;
  //cout << get_geodesic_distances("scrambled_elephant.obj") << endl;
  for(int i=0; i<5; ++i){
    string thing = "scrambled_xyzrgb_dragon_digital_connected.ply";//"scrambled_dragon_vrip_connected.ply";//"scrambled_dragon_vrip_res2_connected.ply";
    cout <<thing << endl << get_geodesic_distances(thing) << endl;
  }
  cout<<endl;
  for(int i=0; i<5; ++i){
    string thing = "xyzrgb_dragon_digital_connected.ply";//"dragon_vrip_connected.ply";//"dragon_vrip_res2_connected.ply";
    cout <<thing << endl << get_geodesic_distances(thing) << endl;
  }
  cout<<endl;
  for(int i=0; i<5; ++i){
    string thing = "reordered_xyzrgb_dragon_digital_connected.ply";//"reordered_dragon_vrip_connected.ply";//"reordered_dragon_vrip_res2_connected.ply";
    cout <<thing << endl << get_geodesic_distances(thing) << endl;
  }


  //twod = V.col(2).minCoeff()==V.col(2).maxCoeff();
  //bbd = (V.colwise().maxCoeff()-V.colwise().minCoeff()).norm();
  //SparseMatrix<double> L,M;
  
  //massmatrix(V,F,MASSMATRIX_TYPE_DEFAULT,M);
  
}
