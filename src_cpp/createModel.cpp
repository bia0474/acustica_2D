#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>


//-------------------------------
// media speeds
//-------------------------------

float* velocity(int nx, int nz, float c1, float c2, int interface_Z){

    float *velocity = (float*) malloc(nx * nz * sizeof(float));

    for(int i = 0; i < nx; i++){

        for(int j = 0; j < nz; j++){

            if(j < interface_Z){
                velocity[i * nz + j] = c1;
            }
            else{
                velocity[i * nz + j] = c2;
            }
        }
    }

    return velocity;
}

//----------------------------------
// CFL condition
//----------------------------------


bool CFL(const float* c, float dt, float dx, float dz, int nx, int nz){ //function of the stability codition

    float cmax = 0.0f;

    for(int i = 0; i < nx; i++){

        for(int j = 0; j < nz; j++){

            cmax = std::max(cmax, c[i * nz + j]);
        }
    }

    float courant = cmax * dt / dx;

    if(courant > 0.7f){

        std::cout << "ERROR! NOT STABLE" << std::endl;
        return false;
    }

    return true;
}
int main(){

    float L = 5000.0; //model size

    float dx = 10.0; //space step in X
    float dz = 10.0; //space step in Z
    float dt = 0.0005; //time lapse

    int nx = int(L/dx) + 1; //number of spatial points in X
    int nz = int(L/dz) + 1; //number of spatial points in Z

    //speeds

    float c1 = 1500.0f;
    float c2 = 4000.0f;

    //interface in Z

    int interface_Z = nz/2;

    float *c = (float*) malloc(nx * nz * sizeof(float));

    c = velocity(nx, nz, c1, c2, interface_Z);
//----------------------------------
// CFL check
//----------------------------------

    if(CFL(c, dt, dx, dz, nx, nz)){

        std::cout << "Stable simulation" << std::endl;
    }
    else{

        return 1;
    }

//------------------------------------------
// Save the documento of the velocity model
//-----------------------------------------

    FILE *file_velocities = fopen("/home/processamento/acustica_2D/inputs/velocityModel.csv", "w"); //cria um ponteiro para um arquivo e abre um arquivo chamado "sources.csv" no modo write("w") (escrita)

    if(file_velocities == NULL){
        printf("Erro ao abrir sources.csv\n");
        return 1;
    }

    fprintf(file_velocities, "c\n"); //escrevendo o cabeçalho

    for(int i = 0; i < nx * nz; i++){

        fprintf(file_velocities, "%d,%d,%.1f,%.1f\n", c[i]);

    }

    fclose(file_velocities);

    std::cout << "Velocity csv file saved!" << std::endl;

    free(c);

    return 0;
}