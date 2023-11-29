#include <igl/heat_geodesics.h>
#include <igl/massmatrix.h>
#include <igl/writePLY.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Sparse>
#include <iostream>
#include <chrono>

using namespace std;

double get_geodesic_distances(string model_path, int attempt_number){
  using namespace Eigen;
  using namespace igl;

  //to be used for printing the output
  std::string output_file_prefix = model_path + "";
  output_file_prefix = output_file_prefix.substr(0,output_file_prefix.size()-4);
  output_file_prefix += "_" + std::to_string(attempt_number);

  MatrixXd V;
  MatrixXi F;

  if(!read_triangle_mesh(model_path.data(),V,F))
  {
    cout<<"failed to load mesh"<<endl;
  }

  VectorXd D;

  //time start
  const auto start_time {chrono::steady_clock::now()};
  
  HeatGeodesicsData<double> data;
  heat_geodesics_precompute(V,F,data);

  //not computing any particular distance since we want an equal comparison between different model orderings.
  /*for(int i=0; i<V.rows(); ++i){
    int i = 0;
    VectorXi gamma (1,1);
    gamma << i;
    D = VectorXd::Zero(data.Grad.cols(),1);
    heat_geodesics_solve(data,gamma,D);
    //break;
  }*/

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
  int num_times = 0;
  for(int i=1; i<argc; ++i){
    if(i%2==1){
        num_times=std::stoi(argv[i]);
    }else{
        printf("Simulating %s for %d times.\n", argv[i], num_times);

        //run the program and get the timings
        std::vector<double> timings;
        double average_time = 0.;
        for(int j=0; j<num_times; ++j){
          double time_elapsed_one_run = get_geodesic_distances(argv[i], j);
          timings.push_back(time_elapsed_one_run);
          average_time += time_elapsed_one_run;
        }
        average_time /= (double)num_times;

        //get standard deviation of timings
        double std_dev_time = 0.;
        for(int j=0; j<num_times; ++j){
          double diff = timings[j] - average_time;
          std_dev_time += diff*diff;
        }
        std_dev_time /= (double)num_times;
        std_dev_time = sqrt(std_dev_time);

        //print timings
        char file_path[70];
        std::string file_name = argv[i];
        file_name = file_name.substr(0,file_name.size()-4);
        sprintf(file_path, "heat_geodesics_output/%s_timings.txt", file_name.data());
        FILE* fp = fopen(file_path, "w");
        for(int j=0; j<num_times; ++j){
          fprintf(fp, "%lf\n", timings[j]);
        }
        fprintf(fp, "\n\nAverage: %lf\nStandard deviation: %lf\n", average_time, std_dev_time);
        fclose(fp);
    }
  }


  //twod = V.col(2).minCoeff()==V.col(2).maxCoeff();
  //bbd = (V.colwise().maxCoeff()-V.colwise().minCoeff()).norm();
  //SparseMatrix<double> L,M;
  
  //massmatrix(V,F,MASSMATRIX_TYPE_DEFAULT,M);
  
}
