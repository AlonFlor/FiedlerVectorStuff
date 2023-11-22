//#include <igl/cotmatrix.h>
//#include <igl/massmatrix.h>
//#include <igl/readPLY.h>
//#include <igl/eigs.h>
//#include <Eigen/Dense>
#include <cmath>
//#include <Eigen/Sparse>
#include <iostream>
#include <fstream>
#include <sstream>
#include "volumetricMeshLoader.h"
#include "corotationalLinearFEM.h"
#include "corotationalLinearFEMStencilForceModel.h"
#include "forceModelAssembler.h"
#include "generateMassMatrix.h"
#include "implicitBackwardEulerSparse.h"
#include <chrono>

int print_frame(std::string file_prefix, std::string header, std::string suffix, int i, double* u_original, double* u, int r, FILE* fp){
   char file_path[50];
   sprintf(file_path, "output/%s_%04d.ply",file_prefix.data(), i);
   printf("%s\n",file_path);
   fp = fopen(file_path, "w");
   fprintf(fp, "%s", header.data());
   for(int i=0; i<r; i++){
      fprintf(fp, "%f ",u_original[i] + u[i]);
      if(i>0 && (i+1)%3==0){
         fprintf(fp, "\n");
      }
   }
   fprintf(fp, "%s", suffix.data());
   fclose(fp);

   return 0;
}

std::vector<std::string> split(const std::string& s)
{
    std::stringstream str_stream(s);
    std::vector<std::string> words;
    for (std::string w; str_stream>>w; ) words.push_back(w);
    return words;
}

std::vector<std::string> get_header_and_suffix(int num_vertices, int num_faces, std::string header, std::string suffix, char*file_path){
   header += "ply\nformat ascii 1.0\nelement vertex "+std::to_string(num_vertices)+"\nproperty float x\nproperty float y\nproperty float z\n"+
             "element face "+std::to_string(num_faces)+"\nproperty list uchar int vertex_indices\nend_header\n";
   std::ifstream file(file_path);
   
   std::string line_std;
   while(std::getline(file, line_std)){
      if(line_std=="*ELEMENTS"){
         break;
      }
   }
   std::getline(file, line_std);
   std::getline(file, line_std);
   while(std::getline(file, line_std)){
      if(line_std.size()==0){
         break;
      }
      std::vector<std::string> face_elements = split(line_std);
      int size = face_elements.size();
      std::string to_add = std::to_string(size-1);
      for(int i=1; i<size; ++i){
         to_add += " " + std::to_string(std::stoi(face_elements[i])-1);
      }
      to_add += "\n";
      suffix += to_add;
   }

   std::vector<std::string> result;
   result.push_back(header);
   result.push_back(suffix);
   return result;
}

double run_one_sim(char* inputFilename, int attempt_number){
   std::string output_file_prefix = inputFilename;
   output_file_prefix = output_file_prefix.substr(0,output_file_prefix.size()-4);
   output_file_prefix += "_" + std::to_string(attempt_number);
   VolumetricMesh * volumetricMesh = VolumetricMeshLoader::load(inputFilename);
   if (volumetricMesh == NULL) printf("Error: failed to load mesh.\n");
   else printf("Success. Number of vertices: %d . Number of elements: %d .\n", volumetricMesh->getNumVertices(), volumetricMesh->getNumElements());

   CorotationalLinearFEM * deformableModel = new CorotationalLinearFEM(volumetricMesh);

   // Create the class to compute internal elastic forces and stiffness matrices
   // on each mesh element.
   StencilForceModel * stencilForceModel = new CorotationalLinearFEMStencilForceModel(deformableModel);
   // Create the class to connect the deformable model to the integrator.
   // ForceModelAssembler automatically computes in parallel if Intel TBB is enabled.
   ForceModel * forceModel = new ForceModelAssembler(stencilForceModel);
   int r = 3 * volumetricMesh->getNumVertices(); // total number of DOFs

   //get the original displacements
   std::string header = "";
   std::string suffix = "";
   std::vector<std::string> hs = get_header_and_suffix(volumetricMesh->getNumVertices(), volumetricMesh->getNumElements(), header, suffix, inputFilename);
   header = hs[0];
   suffix = hs[1];
   double * u_original = (double*) malloc (sizeof(double) * r);
   int count = 0;
   for(int i=0; i<r;i++){
      u_original[i]=volumetricMesh->getVertex(count)[i%3];
      if(i%3==2){
         count+=1;
      }
   }
   double timestep = 0.0333; // the timestep, in seconds
   SparseMatrix * massMatrix;
   // Create consistent (non-lumped) mass matrix.
   GenerateMassMatrix::computeMassMatrix(volumetricMesh, &massMatrix, true);
   // constraining vertices 4, 10, 14 (constrained DOFs are specified 0-indexed):
   int numConstrainedDOFs = 0;//9;
   int constrainedDOFs[0] = {};
   //int constrainedDOFs[9] = { 12, 13, 14, 30, 31, 32, 42, 43, 44 };
   // (tangential) Rayleigh damping
   double dampingMassCoef = 0.0; // "underwater"-like damping (here turned off)
   double dampingStiffnessCoef = 0.01; // (primarily) high-frequency damping
   // initialize the integrator
   ImplicitBackwardEulerSparse * implicitBackwardEulerSparse = new
   ImplicitBackwardEulerSparse(r, timestep, massMatrix, forceModel,
   numConstrainedDOFs, constrainedDOFs,
   dampingMassCoef, dampingStiffnessCoef);

   // allocate buffer for external forces
   double * f = (double*) malloc (sizeof(double) * r);

   // allocate buffer to read the resulting displacements
   double * u = (double*) malloc (sizeof(double) * r);
   for(int i=0; i<r; i++){
      u[i]=0.;
   }

   double time_elapsed = 0.;
   int numTimesteps = 10;
   FILE* fp;
   for(int i=0; i<numTimesteps; i++){
      //print the current time step
      print_frame(output_file_prefix, header, suffix, i, u_original, u, r, fp);

      //time start
      const auto start_time {std::chrono::steady_clock::now()};

      //printf("%d\n",i);
      // important: must always clear forces, as they remain in effect unless changed
      implicitBackwardEulerSparse->SetExternalForcesToZero();
      if (i == 0){ // set some force at the first timestep
         for(int j=0; j<r; j++)f[j] = -50.; // set forces for every node
         //f[37] = -500; // apply force of -500 N to vertex 12, in y-direction, 3*12+1 = 37
         implicitBackwardEulerSparse->SetExternalForces(f);
      }
      implicitBackwardEulerSparse->DoTimestep();

      //time end
      const auto end_time {std::chrono::steady_clock::now()};
      const std::chrono::duration<double> elapsed_seconds{end_time - start_time};
      time_elapsed += elapsed_seconds.count();

      //update the displacements
      implicitBackwardEulerSparse->GetqState(u);
   }
   //print the final time step
   print_frame(output_file_prefix, header, suffix, numTimesteps, u_original, u, r, fp);

   //print time
   std::cout << time_elapsed << std::endl;

   free(u_original);
   free(u);
   free(f);

   return time_elapsed;
}


int main(int argc, char *argv[])
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
            double time_elapsed_one_run = run_one_sim(argv[i], j);
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
         char file_path[50];
         std::string file_name = argv[i];
         file_name = file_name.substr(0,file_name.size()-4);
         sprintf(file_path, "output/%s_timings.txt", file_name.data());
         FILE* fp = fopen(file_path, "w");
         for(int j=0; j<num_times; ++j){
            fprintf(fp, "%lf\n", timings[j]);
         }
         fprintf(fp, "\n\nAverage: %lf\nStandard deviation: %lf\n", average_time, std_dev_time);
         fclose(fp);


      }
      
   }


   return 0;
}

