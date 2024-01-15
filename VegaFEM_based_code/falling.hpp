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

#include "neoHookeanIsotropicMaterial.h"
//#include "MooneyRivlinIsotropicMaterial.h"
#include "isotropicHyperelasticFEM.h"
#include "isotropicHyperelasticFEMStencilForceModel.h"

#include "linearFEMStencilForceModel.h"

#include "massSpringSystem.h"
#include "massSpringSystemFromTetMesh.h"
#include "massSpringStencilForceModel.h"

#include "StVKElementABCDLoader.h"
#include "StVKFEM.h"
#include "StVKStencilForceModel.h"

#include "forceModelAssembler.h"
#include "generateMassMatrix.h"
#include "implicitBackwardEulerSparse.h"
#include <chrono>

#include "helper_methods.hpp"


double run_one_sim_falling(char* inputFilename, int attempt_number, std::string constitutive_model, int numTimesteps){
   std::string output_file_prefix = inputFilename;
   std::string input_file = "models/";
   input_file += inputFilename;
   output_file_prefix = output_file_prefix.substr(0,output_file_prefix.size()-4) + "_" + constitutive_model + "_falling";
   VolumetricMesh * volumetricMesh = VolumetricMeshLoader::load(input_file.data());
   if (volumetricMesh == NULL) printf("Error: failed to load mesh.\n");
   else printf("Success. Number of vertices: %d . Number of elements: %d .\n", volumetricMesh->getNumVertices(), volumetricMesh->getNumElements());

   // Create the class to compute internal elastic forces and stiffness matrices
   // on each mesh element.
   IsotropicHyperelasticFEM * IsotropicHyperelasticDeformableModel;
   IsotropicMaterial * material;
   MassSpringSystem * MassSpringDeformableModel;
   StVKElementABCD * precomputedIntegrals;
   StVKFEM * StVKDeformableModel;
   CorotationalLinearFEM * CorotationalLinearDeformableModel;
   TetMesh* tetMesh = (TetMesh*)volumetricMesh;

   StencilForceModel * stencilForceModel;
   StencilForceModel * stencilForceModel_1;
   std::cout<<constitutive_model.data()<<std::endl;
   if(constitutive_model.compare("neoHookean")==0){
      material = new NeoHookeanIsotropicMaterial(tetMesh);
      IsotropicHyperelasticDeformableModel = new IsotropicHyperelasticFEM(tetMesh, material);
      stencilForceModel = new IsotropicHyperelasticFEMStencilForceModel(IsotropicHyperelasticDeformableModel);
   }/*else if(constitutive_model.compare("MR")==0){
      material = new MooneyRivlinIsotropicMaterial(tetMesh);
      IsotropicHyperelasticDeformableModel = new IsotropicHyperelasticFEM(tetMesh, material);
      stencilForceModel = new IsotropicHyperelasticFEMStencilForceModel(IsotropicHyperelasticDeformableModel);
   }*/else if(constitutive_model.compare("massSpring")==0){
      MassSpringSystemFromTetMesh massSpringSystemFromTetMesh;
      massSpringSystemFromTetMesh.GenerateMassSpringSystem(tetMesh, &MassSpringDeformableModel, 1.0, 50000., 0.01);
      stencilForceModel = new MassSpringStencilForceModel(MassSpringDeformableModel);
   }else if(constitutive_model.compare("StVK")==0){
      unsigned int loadingFlag = 0; // 0 = use the low-memory version, 1 = use the high-memory version
      precomputedIntegrals = StVKElementABCDLoader::load(volumetricMesh, loadingFlag);
      StVKDeformableModel = new StVKFEM(volumetricMesh, precomputedIntegrals);
      stencilForceModel = new StVKStencilForceModel(StVKDeformableModel);
   }else if(constitutive_model.compare("corotationalLinear")==0){
      CorotationalLinearDeformableModel = new CorotationalLinearFEM(volumetricMesh);
      stencilForceModel = new CorotationalLinearFEMStencilForceModel(CorotationalLinearDeformableModel);
   }else if(constitutive_model.compare("linear")==0){
      CorotationalLinearDeformableModel = new CorotationalLinearFEM(volumetricMesh);
      stencilForceModel_1 = new CorotationalLinearFEMStencilForceModel(CorotationalLinearDeformableModel);
      stencilForceModel = new LinearFEMStencilForceModel(stencilForceModel_1);
   }


   // Create the class to connect the deformable model to the integrator.
   // ForceModelAssembler automatically computes in parallel if Intel TBB is enabled.
   ForceModel * forceModel = new ForceModelAssembler(stencilForceModel);
   int r = 3 * volumetricMesh->getNumVertices(); // total number of DOFs

   //get the original displacements
   std::string header = "";
   std::string suffix = "";
   std::string suffix_file_loc = "models/triangle_lists/";
   suffix_file_loc += inputFilename;
   suffix_file_loc += "_faces.txt";
   std::vector<std::string> hs = get_header_and_suffix(volumetricMesh->getNumVertices(), volumetricMesh->getNumElements(), header, suffix, suffix_file_loc.data());
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
   int numConstrainedDOFs = 0;//9;
   int constrainedDOFs[0] = {};
   /*// constraining vertices 4, 10, 14 (constrained DOFs are specified 0-indexed):
   int constrainedDOF_0 = get_start_index_of_vertex(1.22424, -2.25693, 1.43553, u_original, r);
   int constrainedDOF_1 = get_start_index_of_vertex(1.19529, -2.30061, 1.44579, u_original, r);
   int constrainedDOF_2 = get_start_index_of_vertex(1.15391, -2.35621, 1.40962, u_original, r);
   //int constrainedDOFs[9] = { 12, 13, 14, 30, 31, 32, 42, 43, 44 };
   int constrainedDOFs[9] = { constrainedDOF_0, constrainedDOF_0+1, constrainedDOF_0+2, constrainedDOF_1, constrainedDOF_1+1, constrainedDOF_1+2,
                              constrainedDOF_2, constrainedDOF_2+1, constrainedDOF_2+2 };
   //std::cout<<std::to_string(constrainedDOFs[0])+" "<<std::to_string(constrainedDOFs[1])+" "<<std::to_string(constrainedDOFs[2])+" ";
   //std::cout<<std::to_string(constrainedDOFs[3])+" "<<std::to_string(constrainedDOFs[4])+" "<<std::to_string(constrainedDOFs[5])+" ";
   //std::cout<<std::to_string(constrainedDOFs[6])+" "<<std::to_string(constrainedDOFs[7])+" "<<std::to_string(constrainedDOFs[8])<<std::endl;*/
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

   int force_index = get_start_index_of_vertex(1.3844, -2.18075, 1.35392, u_original, r) + 1;
   //std::cout<<"force_index "<<force_index<<std::endl;
   //std::cout<<"u[force_index] "<<u_original[force_index]<<std::endl;

   double time_elapsed = 0.;
   FILE* fp;
   for(int i=0; i<numTimesteps; i++){
      //print the current time step
      if(attempt_number==0)print_frame(output_file_prefix, header, suffix, i, u_original, u, r, fp);

      //time start
      const auto start_time {std::chrono::steady_clock::now()};

      //printf("%d\n",i);
      // important: must always clear forces, as they remain in effect unless changed
      implicitBackwardEulerSparse->SetExternalForcesToZero();
      if (i == 0){ // set some force at the first timestep
         for(int j=0; j<r; j++)f[j] = -50.; // set forces for every node
         //f[force_index] = -500; // apply force of -500 N to vertex 12, in y-direction, 3*12+1 = 37
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
   if(attempt_number==0)print_frame(output_file_prefix, header, suffix, numTimesteps, u_original, u, r, fp);

   //print time
   std::cout << time_elapsed << std::endl;

   //free memory
   free(u_original);
   free(u);
   free(f);
   if(constitutive_model.compare("neoHookean")==0){//} || constitutive_model.compare("MR")==0){
      delete material;
      delete IsotropicHyperelasticDeformableModel;
   }else if(constitutive_model.compare("massSpring")==0){
      delete MassSpringDeformableModel;
   }else if(constitutive_model.compare("StVK")==0){
      //delete precomputedIntegrals;      //Not sure why, but freeing the memory from precomputedIntegrals leads to a segfault.
      delete StVKDeformableModel;
   }else if(constitutive_model.compare("corotationalLinear")==0){
      delete CorotationalLinearDeformableModel;
   }else if(constitutive_model.compare("linear")==0){
      delete CorotationalLinearDeformableModel;
      delete stencilForceModel_1;
   }
   delete stencilForceModel;
   delete forceModel;
   delete implicitBackwardEulerSparse;
   delete volumetricMesh; //I think this should also delete tetMesh.

   return time_elapsed;
}




