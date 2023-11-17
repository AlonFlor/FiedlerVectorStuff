//#include <igl/cotmatrix.h>
//#include <igl/massmatrix.h>
//#include <igl/readPLY.h>
//#include <igl/eigs.h>
//#include <Eigen/Dense>
#include <cmath>
//#include <Eigen/Sparse>
#include <iostream>
#include "volumetricMeshLoader.h"
#include "corotationalLinearFEM.h"
#include "corotationalLinearFEMStencilForceModel.h"
#include "forceModelAssembler.h"
#include "generateMassMatrix.h"
#include "implicitBackwardEulerSparse.h"
#include <chrono>


int main(int argc, char *argv[])
{
   char inputFilename[96] = "myInputMeshFile.veg";
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
   double timestep = 0.0333; // the timestep, in seconds
   SparseMatrix * massMatrix;
   // Create consistent (non-lumped) mass matrix.
   GenerateMassMatrix::computeMassMatrix(volumetricMesh, &massMatrix, true);
   // constraining vertices 4, 10, 14 (constrained DOFs are specified 0-indexed):
   int numConstrainedDOFs = 9;
   int constrainedDOFs[9] = { 12, 13, 14, 30, 31, 32, 42, 43, 44 };
   // (tangential) Rayleigh damping
   double dampingMassCoef = 0.0; // "underwater"-like damping (here turned off)
   double dampingStiffnessCoef = 0.01; // (primarily) high-frequency damping
   // initialize the integrator
   ImplicitBackwardEulerSparse * implicitBackwardEulerSparse = new
   ImplicitBackwardEulerSparse(r, timestep, massMatrix, forceModel,
   numConstrainedDOFs, constrainedDOFs,
   dampingMassCoef, dampingStiffnessCoef);

   // alocate buffer for external forces
   double * f = (double*) malloc (sizeof(double) * r);
   int numTimesteps = 10;
   const auto start_time {std::chrono::steady_clock::now()};
   for(int i=0; i<numTimesteps; i++){
      printf("%d\n",i);
      // important: must always clear forces, as they remain in effect unless changed
      implicitBackwardEulerSparse->SetExternalForcesToZero();
      if (i == 0){ // set some force at the first timestep
         for(int j=0; j<r; j++)f[j] = 0; // clear to 0
         f[37] = -500; // apply force of -500 N to vertex 12, in y-direction, 3*12+1 = 37
         implicitBackwardEulerSparse->SetExternalForces(f);
      }
      implicitBackwardEulerSparse->DoTimestep();
   }
   // alocate buffer to read the resulting displacements
   double * u = (double*) malloc (sizeof(double) * r);
   implicitBackwardEulerSparse->GetqState(u);

   const auto end_time {std::chrono::steady_clock::now()};
   const std::chrono::duration<double> elapsed_seconds{end_time - start_time};
   std::cout << elapsed_seconds.count() << std::endl;

   //TODO
   // initialize mesh
   // implement the force

   /*std::vector<std::vector<double>> vector_data = read_vectors("vectors.txt");
   int vectors_size = vector_data.size();
   Eigen::MatrixXd V = Eigen::MatrixXd::Zero(vectors_size,3);
   for(int i=0; i<vectors_size; ++i){
      for(int j=0; j<3; ++j)
      V(i,j) = vector_data[i][j];
   }
   //std::cout<<"V:"<<std::endl<<V<<std::endl;

   std::vector<std::vector<int>> faces_data = read_faces("faces.txt");
   int faces_size = faces_data.size();
   Eigen::MatrixXi F = Eigen::MatrixXi::Zero(faces_size,3);
   for(int i=0; i<faces_size; ++i){
      for(int j=0; j<3; ++j)
      F(i,j) = faces_data[i][j];
   }
   //std::cout<<"F:"<<std::endl<<F<<std::endl;

   Eigen::MatrixXd V;
   Eigen::MatrixXi F;
   igl::readPLY("../models/Armadillo_digital.ply", V, F);
   std::cout<<V<<std::endl;
   std::cout<<F<<std::endl;*/
   
   /*Eigen::SparseMatrix<double> L;
   igl::cotmatrix(V,F,L);
   L = (-L).eval();

   //print L
   FILE* cot_matrix_file = fopen("cotangent_Laplacian.csv", "w");
   fprintf(cot_matrix_file, "empty_header\n");
   //int Fiedler_vector_size = Fiedler_vector.size();
   for(int i=0; i<L.outerSize(); ++i){
      for(Eigen::SparseMatrix<double>::InnerIterator it(L,i); it; ++it){
         double val = it.value();
         if(std::isinf(val)){
            char sign[10];
            snprintf(sign, 10, "%lf", val);
            if(sign[0]=='-'){
               val = -1000000.;
            }else{
               val = 1000000.;
            }
         }else{
            if(std::isnan(val)){
               val=0.;
            }
         }
         fprintf(cot_matrix_file, "%ld,%ld,%lf\n",it.col(),it.row(),val);
      }
   }
   fclose(cot_matrix_file);*/


   return 0;
}

