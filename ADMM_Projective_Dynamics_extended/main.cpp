
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

/*#include "System.hpp"
#include "TriangleForce.hpp"
#include "BendForce.hpp"
#include "AnchorForce.hpp"
#include "TetForce.hpp"
#include "CollisionForce.hpp"
#include "CollisionFloor.hpp"*/
#include "Solver.hpp"
#include "AddMeshes.hpp"

#include <Eigen/Sparse>
#include <chrono>


int print_frame(std::string file_prefix, Eigen::MatrixXi render_faces, std::vector<int> mesh_lengths, int num, Eigen::VectorXd x, FILE* fp){
    char file_path[100];
    sprintf(file_path, "output/%s_%04d.obj",file_prefix.data(), num);
    std::cout<<file_path<<std::endl;
    fp = fopen(file_path, "w");

    int verts_mesh_lengths_index = 0;
    int next_object_index = 0;
    int faces_offset = 0;
    int num_faces = render_faces.rows(); //faces per mesh

    //print vertices and faces
    int x_len = x.rows();
    for(int i=0; i<x_len; ++i){
        //print object labels
        if(i == next_object_index){
            std::string object_label = "o object_";
            object_label += std::to_string(verts_mesh_lengths_index);
            fprintf(fp, "%s\n", object_label.data());

            //update number needed to reach before writing next object
            next_object_index += 3*mesh_lengths[verts_mesh_lengths_index];
            verts_mesh_lengths_index += 1;
        }
        //print vertices
        if(i%3==0)fprintf(fp, "v ");
        fprintf(fp, "%f",x[i]);
        if(i>0 && (i+1)%3==0){
            fprintf(fp, "\n");
        }else{
            fprintf(fp, " ");
        }
        //print faces
        if(i+1 == next_object_index){
            for(int j=0; j<num_faces; ++j){
                fprintf(fp, "f %d %d %d\n", render_faces(j,0)+faces_offset, render_faces(j,1)+faces_offset, render_faces(j,2)+faces_offset);
            }
            faces_offset += mesh_lengths[verts_mesh_lengths_index];
        }
    }

    fclose(fp);
    return 0;
}

Eigen::MatrixXi get_render_faces(std::string mesh_name){
    std::string mesh_path = "../models/triangle_lists/";
    mesh_path += mesh_name + "_faces.txt";
    std::fstream fp;
    fp.open(mesh_path.data(), std::ios::in);
    int num_faces;

    std::string line;
    std::getline(fp, line);
    std::stringstream line_stream(line);
    line_stream >> num_faces;
    Eigen::MatrixXi F = Eigen::MatrixXi::Zero(num_faces,3);
    for(int i=0; i<num_faces; ++i){
        int v1, v2, v3;
        std::getline(fp, line);
        std::stringstream line_stream(line);
        line_stream >> v1 >> v2 >> v3;
        F(i,0)=v1+1;
        F(i,1)=v2+1;
        F(i,2)=v3+1;
    }
    fp.close();

    return F;
}


mcl::TetMesh::Ptr create_mesh(int mesh_index, float base_height, std::string mesh_path, int number_of_meshes){
    mcl::TetMesh::Ptr mesh1 = mcl::TetMesh::create();
    mcl::meshio::load_elenode( mesh1.get(), mesh_path.data() );

    //declare mesh constitutive model. If only doing one constitutive model, then no need to read constitutive model off of args
    mesh1->flags |= binding::NOSELFCOLLISION | binding::LINEAR;//NEOHOOKEAN;//

    //rescale object to be 0.07 m height
    float min_y = 10000.;
    float max_y = -10000.;
    for(int i=0; i<mesh1->vertices.size(); ++i){
        float candidate_y = mesh1->vertices[i][1];
        if(candidate_y < min_y){
            min_y = candidate_y;
        }
        if(candidate_y > max_y){
            max_y = candidate_y;
        }
    }
    float object_height = max_y - min_y;
    float scale_factor = 0.07f/object_height;

    float shape_move_number = 2*3.14159265*mesh_index / number_of_meshes;
    mcl::XForm<float> scale = mcl::xform::make_scale<float>(scale_factor,scale_factor,scale_factor);
    mcl::XForm<float> rotate = mcl::xform::make_rot<float>(30.f*mesh_index,mcl::Vec3f(mesh_index%3==0,(mesh_index+1)%3==0,(mesh_index+2)%3==0));

    mesh1->apply_xform(rotate*scale);

    //set translations
    float x_trans = 0.07*sin(shape_move_number);
    float y_trans = base_height + 0.08*mesh_index;
    float z_trans = 0.07*cos(1.5*shape_move_number);

    //apply translations
    for(int i=0; i<mesh1->vertices.size(); ++i){
        mesh1->vertices[i][0] += x_trans;
        mesh1->vertices[i][1] += y_trans;
        mesh1->vertices[i][2] += z_trans;
    }

    return mesh1;
}



double run_one_sim(std::string mesh_name, std::string scenario_name, int number_of_time_steps, int attempt_number)
{
    std::string output_file_prefix = mesh_name;
    if(!mesh_name.compare(mesh_name.size()-4, mesh_name.size(), ".ply") || !mesh_name.compare(mesh_name.size()-4, mesh_name.size(), ".PLY")){
        output_file_prefix = output_file_prefix.substr(0,output_file_prefix.size()-4);
    }
    output_file_prefix += "_" + scenario_name;

    std::cout<<"mesh name: "<<mesh_name<<std::endl;
    //TODO: load ply objects, discriminate between ply and node/ele files based on file extension.
    std::string mesh_path = "../models/";
    mesh_path += mesh_name;

    std::shared_ptr<admm::Solver> solver = std::make_shared<admm::Solver>();
    admm::Solver::Settings settings;

    //process scenarios
    int number_of_meshes = 1;
    float base_height = 0.f; //note: floor is at -0.2 m to prevent some stupid inverted tetrahedra error involving small numbers in lines 41-44 of TetEnergyTerm.cpp
    admm::Lame stiffness_to_use;
    settings.admm_iters = 30;
    if(!scenario_name.compare(0, scenario_name.size(), "falling_in_bowl")){
        std::cout<<"scenario name: "<<scenario_name<<std::endl;
        //add bowl
        Eigen::Matrix<double,3,1> center;
        center << 0.,-0.1,0.;
        std::shared_ptr<admm::PassiveCollision> bowl = std::make_shared<admm::Bowl>( admm::Bowl(center, 0.089, 0.0935) );
        solver->add_obstacle(bowl);
        //set number of meshes
        number_of_meshes = 1;
        //set the base height
        base_height = 0.2f;
        //alter the settings
        settings.linsolver = 1; // NodalMultiColorGS //2;// UzawaCG
        settings.gravity = -2.;
        stiffness_to_use = admm::Lame::rubber();//admm::Lame(1000000,0.1);
    }else if(!scenario_name.compare(0, scenario_name.size(), "bounce")){
        std::cout<<"scenario name: "<<scenario_name.data()<<std::endl;
        //alter the settings
        settings.linsolver = 1; // NodalMultiColorGS
        settings.gravity = -3;
        stiffness_to_use = admm::Lame::soft_rubber();
        //set number of meshes
        number_of_meshes = 1;
        //set the base height
        base_height = 0.f;
    }else{
        std::cout<<"scenario not available"<<std::endl;
        exit(1);
    }

    //load the meshes
    std::vector<int> mesh_lengths;
    for(int mesh_index=0; mesh_index<number_of_meshes; ++mesh_index){
        mcl::TetMesh::Ptr new_mesh = create_mesh(mesh_index, base_height, mesh_path, number_of_meshes);
        mesh_lengths.push_back(new_mesh->vertices.size());
        binding::add_tetmesh( solver.get(), new_mesh, stiffness_to_use, settings.verbose );
    }

    //add floor
    std::shared_ptr<admm::PassiveCollision> floor_collider = std::make_shared<admm::Floor>( admm::Floor(-0.2f) );
    solver->add_obstacle(floor_collider);

    // Try to init the solver
    if( !solver->initialize(settings) ){ return EXIT_FAILURE; }
    
    //get render faces for the mesh
    Eigen::MatrixXi render_faces = get_render_faces(mesh_name);
    
    FILE* fp;
    double time_elapsed = 0.;
    
    for(int count = 0; count<number_of_time_steps; ++count){
        Eigen::VectorXd m_x = solver->m_x;
        if(attempt_number==0)print_frame(output_file_prefix, render_faces, mesh_lengths, count, m_x, fp);

        //time start
        const auto start_time {std::chrono::steady_clock::now()};

        //run time step
        solver->step();
        
        //time end
        const auto end_time {std::chrono::steady_clock::now()};
        const std::chrono::duration<double> elapsed_seconds{end_time - start_time};
        time_elapsed += elapsed_seconds.count();
        
        Eigen::VectorXd diff = m_x - solver->m_x;
        std::cout<<"step "<<count<<"\t\tstep diff: "<<diff.norm()<<std::endl;//<<solver->m_x<<std::endl;
    }
    Eigen::VectorXd m_x = solver->m_x;
    if(attempt_number==0)print_frame(output_file_prefix, render_faces, mesh_lengths, number_of_time_steps, m_x, fp);

    return time_elapsed;
}

int main(int argc, char *argv[])
{
    std::string scenario_name = argv[1];
    int number_of_time_steps = std::stoi(argv[2]);
    int attempt_number = 0;

    int num_times = 0;
    for(int i=3; i<argc; ++i){
        if(i%2!=0){
            num_times=std::stoi(argv[i]);
        }else{
            printf("Simulating %s for %d times.\n", argv[i], num_times);
            std::string mesh_name = argv[i];

            //run the program and get the timings
            std::vector<double> timings;
            double average_time = 0.;
            for(int j=0; j<num_times; ++j){
                double time_elapsed_one_run = run_one_sim(mesh_name, scenario_name, number_of_time_steps, j);
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
            if(!file_name.compare(mesh_name.size()-4, mesh_name.size(), ".ply") || !file_name.compare(mesh_name.size()-4, mesh_name.size(), ".PLY")){
                file_name = file_name.substr(0,file_name.size()-4);
            }
            file_name = file_name + "_" + scenario_name;
            sprintf(file_path, "timings/%s_timings.txt", file_name.data());
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
