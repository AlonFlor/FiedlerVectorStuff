
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


int print_frame(std::string file_prefix, Eigen::MatrixXi render_faces, int num, Eigen::VectorXd x, FILE* fp){
    char file_path[100];
    sprintf(file_path, "output/%s_%04d.obj",file_prefix.data(), num);
    std::cout<<file_path<<std::endl;
    fp = fopen(file_path, "w");

    //print vertices
    int x_len = x.rows();
    fprintf(fp, "v ");
    for(int i=0; i<x_len; ++i){
        fprintf(fp, "%f ",x[i]);
        if(i>0 && (i+1)%3==0){
            fprintf(fp, "\n");
            if(i<x_len-1)fprintf(fp, "v ");
        }
    }

    //print faces
    int num_faces = render_faces.rows();
    for(int i=0; i<num_faces; ++i){
        fprintf(fp, "f %d %d %d\n", render_faces(i,0), render_faces(i,1), render_faces(i,2));
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



double run_one_sim(std::string mesh_name, std::string scenario_name, int number_of_time_steps, int attempt_number)
{
    std::string output_file_prefix = mesh_name;
    if(!mesh_name.compare(mesh_name.size()-4, mesh_name.size(), ".ply") || !mesh_name.compare(mesh_name.size()-4, mesh_name.size(), ".PLY")){
        output_file_prefix = output_file_prefix.substr(0,output_file_prefix.size()-4);
    }
    output_file_prefix += "_" + scenario_name;

    std::cout<<"mesh name: "<<mesh_name<<std::endl;
    std::cout<<"scenario name: "<<scenario_name<<std::endl;
    //TODO: load ply objects, discriminate between ply and node/ele files based on file extension.
    std::string mesh_path = "../models/";
    mesh_path += mesh_name;

    std::shared_ptr<admm::Solver> solver = std::make_shared<admm::Solver>();
    admm::Solver::Settings settings;

    mcl::TetMesh::Ptr mesh1 = mcl::TetMesh::create();
    mcl::meshio::load_elenode( mesh1.get(), mesh_path.data() );

    //declare mesh constitutive model. TODO: if only doing one or a few scenarios with this, then no need to read constitutive model off of args
    mesh1->flags |= binding::NOSELFCOLLISION | binding::NEOHOOKEAN;

    //rescale object to be 0.07 m height
    double min_y = 1000000.;
    double max_y = -1000000.;
    for(int i=0; i<mesh1->vertices.size(); ++i){
        double candidate_y = mesh1->vertices[i][1];
        if(candidate_y < min_y){
            min_y = candidate_y;
        }
        if(candidate_y > max_y){
            max_y = candidate_y;
        }
    }
    double object_height = max_y - min_y;
    float scale_factor = 0.07f/object_height;

    mcl::XForm<float> scale = mcl::xform::make_scale<float>(scale_factor,scale_factor,scale_factor);
    mcl::XForm<float> rotate = mcl::xform::make_rot<float>(20.f,mcl::Vec3f(1,0,0));
    mesh1->apply_xform(rotate*scale);

    binding::add_tetmesh( solver.get(), mesh1, admm::Lame::soft_rubber(), settings.verbose );
    std::shared_ptr<admm::PassiveCollision> floor_collider = std::make_shared<admm::Floor>( admm::Floor(0.f) );
    solver->add_obstacle(floor_collider);

    //mcl::TetMesh::Ptr bowl_mesh = mcl::TetMesh::create();
    //mcl::meshio::load_elenode( bowl_mesh.get(), "glass_mixing_bowl" );
    //std::shared_ptr<admm::PassiveCollision> bowl = std::make_shared<admm::PassiveMesh>( admm::PassiveMesh(bowl_mesh) );
    Eigen::Matrix<double,3,1> center;
    center << 0.,0.1,0.;
    std::shared_ptr<admm::PassiveCollision> bowl = std::make_shared<admm::Bowl>( admm::Bowl(center, 0.089, 0.0935) );
    solver->add_obstacle(bowl);

    //TODO add other meshes and obstacles

    // Try to init the solver
    settings.linsolver = 1; // NodalMultiColorGS
    settings.gravity = -3.;//-9.8;//0;
    if( !solver->initialize(settings) ){ return EXIT_FAILURE; }
    
    //get faces
    Eigen::MatrixXi render_faces = get_render_faces(mesh_name);

    //set the object according to the given scenario. Only falling_in_bowl is available right now.
    if(!scenario_name.compare(0, scenario_name.size(), "falling_in_bowl")){
        std::cout<<scenario_name.data()<<std::endl;
    }else{
        std::cout<<"scenario not available"<<std::endl;
        exit(1);
    }
    //make the object be above the floor
    min_y = 10000.;
    for(int i=1; i<solver->m_x.size(); i+=3){
        double candidate_min_y = solver->m_x[i];
        if(candidate_min_y < min_y){
            min_y = candidate_min_y;
        }
    }
    double height = 0.3;
    if(min_y < height){
        double y_adjust = height - min_y;
        for(int i=1; i<solver->m_x.size(); i+=3){
            solver->m_x[i] += y_adjust;
        }
    }
    
    FILE* fp;
    double time_elapsed = 0.;
    
    for(int count = 0; count<number_of_time_steps; ++count){
        Eigen::VectorXd m_x = solver->m_x;
        if(attempt_number==0)print_frame(output_file_prefix, render_faces, count, m_x, fp);

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
    if(attempt_number==0)print_frame(output_file_prefix, render_faces, number_of_time_steps, m_x, fp);

    return 0.;//time_elapsed;
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
