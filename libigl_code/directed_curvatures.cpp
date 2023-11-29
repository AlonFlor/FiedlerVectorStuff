#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/parula.h>
#include <igl/per_corner_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Sparse>
#include <iostream>
#include <chrono>

using namespace std;
using namespace Eigen;
using namespace igl;

int write_color_PLY(string output_file_prefix, MatrixXd V, MatrixXi F, VectorXd K){
  int v_rows = V.rows();
  int f_rows = F.rows();

  int num_bins = 9;
  double K_sorted[K.size()];
  for(int i=0; i<v_rows; ++i){
    K_sorted[i] = K(i);
  }
  sort(K_sorted,K_sorted+K.size());
  //bin thresholds
  double max_bin_0 = K_sorted[K.size()/10];
  double max_bin_1 = K_sorted[2*K.size()/10];
  double max_bin_2 = K_sorted[3*K.size()/10];
  double max_bin_3 = K_sorted[4*K.size()/10];
  double max_bin_4 = K_sorted[5*K.size()/10];
  double max_bin_5 = K_sorted[6*K.size()/10];
  double max_bin_6 = K_sorted[7*K.size()/10];
  double max_bin_7 = K_sorted[8*K.size()/10];
  double max_bin_8 = K_sorted[9*K.size()/10];

  string header = "ply\nformat ascii 1.0\nelement vertex "+to_string(v_rows)+"\nproperty float x\nproperty float y\nproperty float z\n" +
             "property uchar red\nproperty uchar green\nproperty uchar blue\n" +
             "element face "+to_string(f_rows)+"\n" +
             "property list uchar int vertex_indices\nend_header\n";
  string vec_data_str = "";
  for(int i=0; i<v_rows; ++i){
      int color = 0.;//int(255 * K(i) * scale_factor);
      bool included = false;
      if(K(i) <= max_bin_0){
        color = 255;
        vec_data_str += to_string(V(i,0)) + " " + to_string(V(i,1)) + " " + to_string(V(i,2)) + " "+ to_string(color) +" 0 255\n";
        included=true;
      }
      if(K(i) > max_bin_0 && K(i) <= max_bin_1){
        color = 255*4/5;
        vec_data_str += to_string(V(i,0)) + " " + to_string(V(i,1)) + " " + to_string(V(i,2)) + " "+ to_string(color) +" 0 255\n";
        included=true;
      }
      if(K(i) > max_bin_1 && K(i) <= max_bin_2){
        color = 255*3/5;
        vec_data_str += to_string(V(i,0)) + " " + to_string(V(i,1)) + " " + to_string(V(i,2)) + " "+ to_string(color) +" 0 255\n";
        included=true;
      }
      if(K(i) > max_bin_2 && K(i) <= max_bin_3){
        color = 255*2/5;
        vec_data_str += to_string(V(i,0)) + " " + to_string(V(i,1)) + " " + to_string(V(i,2)) + " "+ to_string(color) +" 0 255\n";
        included=true;
      }
      if(K(i) > max_bin_3 && K(i) <= max_bin_4){
        color = 255*1/5;
        vec_data_str += to_string(V(i,0)) + " " + to_string(V(i,1)) + " " + to_string(V(i,2)) + " "+ to_string(color) +" 0 255\n";
        included=true;
      }

      if(K(i) > max_bin_4 && K(i) <= max_bin_5){
        color = 255*1/5;
        vec_data_str += to_string(V(i,0)) + " " + to_string(V(i,1)) + " " + to_string(V(i,2)) + " " +
                        to_string(color) + " " + to_string(color) + " " + to_string(255-color) + "\n";
        included=true;
      }
      if(K(i) > max_bin_5 && K(i) <= max_bin_6){
        color = 255*2/5;
        vec_data_str += to_string(V(i,0)) + " " + to_string(V(i,1)) + " " + to_string(V(i,2)) + " " +
                        to_string(color) + " " + to_string(color) + " " + to_string(255-color) + "\n";
        included=true;
      }
      if(K(i) > max_bin_6 && K(i) <= max_bin_7){
        color = 255*3/5;
        vec_data_str += to_string(V(i,0)) + " " + to_string(V(i,1)) + " " + to_string(V(i,2)) + " " +
                        to_string(color) + " " + to_string(color) + " " + to_string(255-color) + "\n";
        included=true;
      }
      if(K(i) > max_bin_7 && K(i) <= max_bin_8){
        color = 255*4/5;
        vec_data_str += to_string(V(i,0)) + " " + to_string(V(i,1)) + " " + to_string(V(i,2)) + " " +
                        to_string(color) + " " + to_string(color) + " " + to_string(255-color) + "\n";
        included=true;
      }
      if(K(i) > max_bin_8){
        color = 255;
        vec_data_str += to_string(V(i,0)) + " " + to_string(V(i,1)) + " " + to_string(V(i,2)) + " " +
                        to_string(color) + " " + to_string(color) + " " + to_string(255-color) + "\n";
        included=true;
      }

      if(!included){
        //value is nan
        vec_data_str += to_string(V(i,0)) + " " + to_string(V(i,1)) + " " + to_string(V(i,2)) + " 0 0 0\n";
      }
  }

  string face_data_str = "";
  for(int i=0; i<f_rows; ++i){
    face_data_str += "3 " + to_string(F(i,0)) + " " + to_string(F(i,1)) + " " + to_string(F(i,2));
    if(i<f_rows-1){
      face_data_str += "\n";
    }
  }

  string string_to_write = header + vec_data_str + face_data_str;

  char file_path[70];
  sprintf(file_path, "directed_curvatures_output/%s.ply",output_file_prefix.data());
  FILE* data_file = fopen(file_path, "w");
  fprintf(data_file, "%s", string_to_write.data());
  fclose(data_file);

  return 0;
}

double get_curvatures(string model_path, int attempt_number){

  //to be used for printing the output
  string output_file_prefix = model_path + "";
  output_file_prefix = output_file_prefix.substr(0,output_file_prefix.size()-4);
  output_file_prefix += "_" + to_string(attempt_number);

  MatrixXd V;
  MatrixXi F;

  if(!read_triangle_mesh(model_path.data(),V,F))
  {
    cout<<"failed to load mesh"<<endl;
  }

  //time start
  const auto start_time {chrono::steady_clock::now()};

  // Alternative discrete mean curvature
  MatrixXd HN;
  SparseMatrix<double> L,M,Minv;
  igl::cotmatrix(V,F,L);
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
  igl::invert_diag(M,Minv);
  // Laplace-Beltrami of position
  HN = -Minv*(L*V);
  // Extract magnitude as mean curvature
  VectorXd H = HN.rowwise().norm();

  // Compute curvature directions via quadric fitting
  MatrixXd PD1,PD2;
  VectorXd PV1,PV2;
  igl::principal_curvature(V,F,PD1,PD2,PV1,PV2);
  // mean curvature
  H = 0.5*(PV1+PV2);

  //time end
  const auto end_time {chrono::steady_clock::now()};
  const chrono::duration<double> elapsed_seconds{end_time - start_time};
  double time_elapsed = elapsed_seconds.count();

  //cout<<"Vector"<<endl<<K<<endl;
  //write_color_PLY(output_file_prefix, V, F, K);

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
        num_times=stoi(argv[i]);
    }else{
        printf("Simulating %s for %d times.\n", argv[i], num_times);

        //run the program and get the timings
        vector<double> timings;
        double average_time = 0.;
        for(int j=0; j<num_times; ++j){
          double time_elapsed_one_run = get_curvatures(argv[i], j);
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
        string file_name = argv[i];
        file_name = file_name.substr(0,file_name.size()-4);
        sprintf(file_path, "directed_curvatures_output/%s_timings.txt", file_name.data());
        FILE* fp = fopen(file_path, "w");
        for(int j=0; j<num_times; ++j){
          fprintf(fp, "%lf\n", timings[j]);
        }
        fprintf(fp, "\n\nAverage: %lf\nStandard deviation: %lf\n", average_time, std_dev_time);
        fclose(fp);
    }
  }
  
}
