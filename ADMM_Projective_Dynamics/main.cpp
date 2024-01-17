
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "System.hpp"
#include "TriangleForce.hpp"
#include "BendForce.hpp"
#include "AnchorForce.hpp"
#include "TetForce.hpp"
#include "CollisionForce.hpp"
#include "CollisionFloor.hpp"

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


struct mesh{
    Eigen::VectorXd V;
    Eigen::MatrixXi connections;
};


mesh read_mesh_file(std::string mesh_name){

    mesh results;

    std::string mesh_path = "../models/";
    mesh_path += mesh_name;
    if(!mesh_name.compare(mesh_name.size()-4, mesh_name.size(), ".ply") || !mesh_name.compare(mesh_name.size()-4, mesh_name.size(), ".PLY")){
        //Code for PLY file
        std::fstream fp;
        fp.open(mesh_path.data(), std::ios::in);
        int num_verts;
        int num_faces;
        std::string line;
        std::string verts_prefix("element vertex ");
        std::string faces_prefix("element face ");
        std::string end_header("end_header");
        while(true){
            std::getline(fp, line);
            if(!line.compare(0,verts_prefix.size(), verts_prefix)){ //line.starts_with(verts_prefix)
                num_verts = std::stoi(line.substr(verts_prefix.size()));
            }
            if(!line.compare(0,faces_prefix.size(), faces_prefix)){ //line.starts_with(faces_prefix)
                num_faces = std::stoi(line.substr(faces_prefix.size()));
            }
            if(!line.compare(0,faces_prefix.size(), end_header))break;
        }
        Eigen::VectorXd V = Eigen::VectorXd::Zero(num_verts*3);
        for(int i=0; i<num_verts; ++i){
            //get verts
            std::getline(fp, line);
            double x, y, z;
            std::stringstream line_stream(line);
            line_stream >> x >> y >> z;
            V(3*i+0)=x;
            V(3*i+1)=y;
            V(3*i+2)=z;
        }
        Eigen::MatrixXi F = Eigen::MatrixXi::Zero(num_faces,3);
        for(int i=0; i<num_faces; ++i){
            //get faces
            std::getline(fp, line);
            int filler, v1, v2, v3;
            std::stringstream line_stream(line);
            line_stream >> filler >> v1 >> v2 >> v3;
            F(i,0)=v1;
            F(i,1)=v2;
            F(i,2)=v3;
        }
        std::cout<<"num_verts "<<num_verts<<std::endl;
        std::cout<<"num_faces "<<num_faces<<std::endl;
        fp.close();
        results.V = V;
        results.connections = F;
    }else{
        //Code for .node and .ele files
        std::fstream fp;
        fp.open((mesh_path+".node").data(), std::ios::in);
        int num_verts;
        std::string header;
        getline(fp, header);
        std::stringstream header_stream(header);
        header_stream >> num_verts;
        std::cout<<"num_verts "<<num_verts<<std::endl;

        Eigen::VectorXd V = Eigen::VectorXd::Zero(num_verts*3);
        std::string line;
        for(int i=0; i<num_verts; ++i){
            //get verts
            std::getline(fp, line);
            double filler, x, y, z;
            std::stringstream line_stream(line);
            line_stream >> filler >> x >> y >> z;
            V(3*i+0)=x;
            V(3*i+1)=y;
            V(3*i+2)=z;
        }
        fp.close();
        results.V = V;
        
        fp.open((mesh_path+".ele").data(), std::ios::in);
        int num_tets;
        getline(fp, header);
        std::stringstream new_header_stream(header);
        new_header_stream >> num_tets;
        std::cout<<"num_tets "<<num_tets<<std::endl;

        Eigen::MatrixXi T = Eigen::MatrixXi::Zero(num_tets,4);
        for(int i=0; i<num_tets; ++i){
            //get tetrahedra
            std::getline(fp, line);
            int filler, v1, v2, v3, v4;
            std::stringstream line_stream(line);
            line_stream >> filler >> v1 >> v2 >> v3 >> v4;
            T(i,0)=v1;
            T(i,1)=v2;
            T(i,2)=v3;
            T(i,3)=v4;
        }
        fp.close();
        results.connections = T;
    }

    return results;
}


int add_forces(std::shared_ptr<admm::System> system, std::string constitutive_model, mesh input_mesh){
    int verts_per_connection = input_mesh.connections.cols(); //triangular or tet

    //force types:
    //  spring
    //  lineartrianglestrain
    //  bend
    //  lineartetstrain
    //  neohookeantet
    //  stvktet
    //  volumetet
    //loop along the faces/tets

    for(int i=0; i<input_mesh.connections.rows(); ++i){
        Eigen::VectorXi connection = input_mesh.connections.row(i);
        
        if(!constitutive_model.compare(0,constitutive_model.size(), "spring")){
            //spring
            double stiffness = 1000.;
            for(int j=0; j<verts_per_connection; ++j){
                for(int k=0; k<j; ++k){
                    std::shared_ptr<admm::Force> new_force(new admm::Spring(connection[j], connection[k], stiffness));
                    system->forces.push_back(new_force);
                }
            }
        }else if(!constitutive_model.compare(0,constitutive_model.size(), "lineartrianglestrain")){
            //linear triangle strain
            if(verts_per_connection==4){
                std::cout<<"Shape should be triangle mesh, not tet mesh"<<std::endl;
                return 1;
            }
            double stiffness = 100.;
            double limit_l = 0.95;
            double limit_u = 1.05;
            std::shared_ptr<admm::Force> new_force(new admm::LimitedTriangleStrain(connection[0], connection[1], connection[2], stiffness, limit_l, limit_u));
            system->forces.push_back(new_force);
        }/*else if(!constitutive_model.compare(0,constitutive_model.size(), "bend")){
            //bend
            if(verts_per_connection==4){
                std::cout<<"Shape should be triangle mesh, not tet mesh"<<std::endl;
                return 1;
            }
            double stiffness = 20.;
            std::shared_ptr<admm::Force> new_force(new admm::BendForce(connection[0], ??, connection[1], connection[2], stiffness));
            system->forces.push_back(new_force);
        }*/else if(!constitutive_model.compare(0,constitutive_model.size(), "lineartetstrain")){
            //lineartetstrain
            if(verts_per_connection==3){
                std::cout<<"Shape should be tet mesh, not triangle mesh"<<std::endl;
                return 1;
            }
            double stiffness = 100000.;
            double weight_scale = 1.;
            std::shared_ptr<admm::Force> new_force(new admm::LinearTetStrain(connection[0], connection[1], connection[2], connection[3], stiffness, weight_scale));
            system->forces.push_back(new_force);
        }else if(!constitutive_model.compare(0,constitutive_model.size(), "neohookeantet")){
            //neohookeantet
            if(verts_per_connection==3){
                std::cout<<"Shape should be tet mesh, not triangle mesh"<<std::endl;
                return 1;
            }
            double mu = 100000.;
            double lambda = 100000.;
            int max_iters = 50;
            std::shared_ptr<admm::Force> new_force(new admm::HyperElasticTet(connection[0], connection[1], connection[2], connection[3], mu, lambda, max_iters, "nh"));
            system->forces.push_back(new_force);
        }else if(!constitutive_model.compare(0,constitutive_model.size(), "stvktet")){
            //stvktet
            if(verts_per_connection==3){
                std::cout<<"Shape should be tet mesh, not triangle mesh"<<std::endl;
                return 1;
            }
            double mu = 100.;
            double lambda = 100.;
            int max_iters = 5;
            std::shared_ptr<admm::Force> new_force(new admm::HyperElasticTet(connection[0], connection[1], connection[2], connection[3], mu, lambda, max_iters, "stvk"));
            system->forces.push_back(new_force);
        }else if(!constitutive_model.compare(0,constitutive_model.size(), "volumetet")){
            //volumetet
            if(verts_per_connection==3){
                std::cout<<"Shape should be tet mesh, not triangle mesh"<<std::endl;
                return 1;
            }
            double stiffness = 100000.;
            double rangeMin = 0.95;
            double rangeMax = 1.05;
            std::shared_ptr<admm::Force> new_force(new admm::TetVolume(connection[0], connection[1], connection[2], connection[3], stiffness, rangeMin, rangeMax ));
            system->forces.push_back(new_force);
        }
            
    }
    
    //std::shared_ptr<admm::Force> af( new admm::StaticAnchor( idx0 ) );
    //system->forces.push_back( af );

    //floor collision
    std::vector<std::shared_ptr<admm::CollisionShape>> shapes;
    Eigen::Vector3d center(0,0,0), scale(1,1,1);
    std::shared_ptr<admm::CollisionShape> shape(new admm::CollisionFloor(center));
    shapes.push_back(shape);
    std::shared_ptr<admm::Force> collision_force(new admm::CollisionForce(shapes));
    system->forces.push_back(collision_force);

    return 0;
}

double run_one_sim(std::string mesh_name, std::string constitutive_model, std::string scenario_name, int number_of_time_steps, int attempt_number)
{
    std::string output_file_prefix = mesh_name;
    if(!mesh_name.compare(mesh_name.size()-4, mesh_name.size(), ".ply") || !mesh_name.compare(mesh_name.size()-4, mesh_name.size(), ".PLY")){
        output_file_prefix = output_file_prefix.substr(0,output_file_prefix.size()-4);
    }
    output_file_prefix += "_" + constitutive_model + "_" + scenario_name;

    std::cout<<"mesh name: "<<mesh_name<<std::endl;
    std::cout<<"constitutive model: "<<constitutive_model<<std::endl;
    std::cout<<"scenario name: "<<scenario_name<<std::endl;
    //Instead: determines the mesh type via whether or not it has a .ply extension. If not, it assumes that .node and .ele files are available.
    mesh mesh1 = read_mesh_file(mesh_name);
    
    std::shared_ptr<admm::System> system = std::shared_ptr<admm::System>(new admm::System());
    double mass = 14.;
    int V_length = mesh1.V.rows();
    int num_elements = V_length/3;
    Eigen::VectorXd m(V_length);
    m.fill(mass/num_elements); // set node masses
    system->add_nodes( mesh1.V, m );
    //std::cout<<system->m_x<<std::endl;

    //rescale object to be 1 m height
    double min_y = 1000000.;
    double max_y = -1000000.;
    for(int i=1; i<system->m_x.size(); i+=3){
        double candidate_y = system->m_x[i];
        if(candidate_y < min_y){
            min_y = candidate_y;
        }
        if(candidate_y > max_y){
            max_y = candidate_y;
        }
    }
    double object_height = max_y - min_y;
    for(int i=0; i<system->m_x.size(); ++i){
        system->m_x[i] /= object_height;
    }



    //add forces for the implicit solve
    if(add_forces(system, constitutive_model, mesh1))exit(1);

    //add gravity if applicable to the scenario
    if(!scenario_name.compare(0, scenario_name.size(), "bounce")){
        std::shared_ptr<admm::ExplicitForce> ef(new admm::ExplicitForce());
        ef->direction[1] = -9.8;
        system->explicit_forces.push_back(ef);
    }


    //initialize system
    double timestep_s = 0.04;
    system->settings.timestep_s = timestep_s;
    system->settings.admm_iters = 20;
    std::cout<<"system->settings.admm_iters: "<<system->settings.admm_iters<<std::endl;
    if( !system->initialize() ){
        std::cout<<"failed to initialize"<<std::endl;
    }

    Eigen::MatrixXi render_faces = get_render_faces(mesh_name);
    

    //set the objects according to the given scenario
    if(!scenario_name.compare(0, scenario_name.size(), "stretch")){
        //setup stretch
        double mid_x = 0.;
        for(int i=0; i<system->m_x.size(); i+=3){
            mid_x += system->m_x[i];
        }
        mid_x /= num_elements;
        for(int i=0; i<system->m_x.size(); i+=3){
            if(system->m_x[i] < mid_x){
                system->m_v[i] += -0.1;
            }else{
                system->m_v[i] = 0.1;
            }
        }
    }/*else if(!scenario_name.compare(0, scenario_name.size(), "bounce")){
        //no need to do anything, gravity was turned on before the system was initialized.
    }*/else if(!scenario_name.compare(0, scenario_name.size(), "unsquash")){
        //squash to a point
        for(int i=0; i<system->m_x.size(); ++i){
            system->m_x[i] = 0.;
        }
    }
    //make the object be above the floor
    min_y = 10000.;
    for(int i=1; i<system->m_x.size(); i+=3){
        double candidate_min_y = system->m_x[i];
        if(candidate_min_y < min_y){
            min_y = candidate_min_y;
        }
    }
    double height = 1.5;
    if(min_y < height){
        double y_adjust = height - min_y;
        for(int i=1; i<system->m_x.size(); i+=3){
            system->m_x[i] += y_adjust;
        }
    }
    
    FILE* fp;
    double time_elapsed = 0.;
    for(int count = 0; count<number_of_time_steps; ++count){
        Eigen::VectorXd m_x = system->m_x;
        if(attempt_number==0)print_frame(output_file_prefix, render_faces, count, m_x, fp);

        //time start
        const auto start_time {std::chrono::steady_clock::now()};

        //run time step
        if(!system->step())break;

        //time end
        const auto end_time {std::chrono::steady_clock::now()};
        const std::chrono::duration<double> elapsed_seconds{end_time - start_time};
        time_elapsed += elapsed_seconds.count();
        
        Eigen::VectorXd diff = m_x - system->m_x;
        std::cout<<"step "<<count<<"\t\tstep diff: "<<diff.norm()<<std::endl;//<<system->m_x<<std::endl;
    }
    Eigen::VectorXd m_x = system->m_x;
    if(attempt_number==0)print_frame(output_file_prefix, render_faces, number_of_time_steps, m_x, fp);

    return time_elapsed;
}

int main(int argc, char *argv[])
{
    std::string constitutive_model = argv[1];
    std::string scenario_name = argv[2];
    int number_of_time_steps = std::stoi(argv[3]);
    int attempt_number = 0;

    int num_times = 0;
    for(int i=4; i<argc; ++i){
        if(i%2==0){
            num_times=std::stoi(argv[i]);
        }else{
            printf("Simulating %s for %d times.\n", argv[i], num_times);
            std::string mesh_name = argv[i];

            //run the program and get the timings
            std::vector<double> timings;
            double average_time = 0.;
            for(int j=0; j<num_times; ++j){
                double time_elapsed_one_run = run_one_sim(mesh_name, constitutive_model, scenario_name, number_of_time_steps, j);
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
            file_name = file_name + "_" + scenario_name + "_" + constitutive_model;
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
