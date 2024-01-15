
#pragma once
#include <iostream>
#include <fstream>
#include <cmath>



int print_frame(std::string file_prefix, std::string header, std::string suffix, int i, double* u_original, double* u, int r, FILE* fp){
    char file_path[100];
    sprintf(file_path, "output/%s_%04d.ply",file_prefix.data(), i);
    printf("%s\n",file_path);
    fp = fopen(file_path, "w");
    fprintf(fp, "%s", header.data());
    for(int i=0; i<r; i++){
        fprintf(fp, "%f ",u_original[i] + u[i]);
        if(i%3==2){
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "%s", suffix.data());
    fclose(fp);

    return 0;
}

std::vector<std::string> split(const std::string& s){
    std::stringstream str_stream(s);
    std::vector<std::string> words;
    for (std::string w; str_stream>>w; ) words.push_back(w);
    return words;
}


std::vector<std::string> get_header_and_suffix(int num_vertices, int num_faces, std::string header, std::string suffix, const char*file_path){
    header += "ply\nformat ascii 1.0\nelement vertex "+std::to_string(num_vertices)+"\nproperty float x\nproperty float y\nproperty float z\n"+
              "element face "+std::to_string(num_faces)+"\nproperty list uchar int vertex_indices\nend_header\n";
    std::ifstream file(file_path);
    
    std::string line;
    std::getline(file, line);
    
    for(int i=0; i<num_faces; ++i){
        int v1, v2, v3;
        std::getline(file, line);
        std::stringstream line_stream(line);
        line_stream >> v1 >> v2 >> v3;
        std::string to_add = "3";
        to_add += " " + std::to_string(v1-1);
        to_add += " " + std::to_string(v2-1);
        to_add += " " + std::to_string(v3-1);
        to_add += "\n";
        suffix += to_add;
    }

    std::vector<std::string> result;
    result.push_back(header);
    result.push_back(suffix);
    return result;
}

int get_start_index_of_vertex(double x, double y, double z, double* u_original, int r){
    //get index of x-value of vertex from dragon model or its reorderings. Needed because the same vertex could be in different locations in different models.
    double dist_sq = 0.;
    double min_dist = 1000000.;
    int min_dist_start_index = -1;
    for(int i=0; i<r; ++i){
        if(i%3==0)dist_sq += (u_original[i] - x)*(u_original[i] - x);
        if(i%3==1)dist_sq += (u_original[i] - y)*(u_original[i] - y);
        if(i%3==2){
            dist_sq += (u_original[i] - z)*(u_original[i] - z);
            if(dist_sq < min_dist){
                min_dist = dist_sq;
                min_dist_start_index = i-2;
            }
            dist_sq = 0.;
        }
    }
    return min_dist_start_index;
}
