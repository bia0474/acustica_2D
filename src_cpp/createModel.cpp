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

    //----------------------------------
    // open the document of PARAMETERS
    //----------------------------------

    FILE *file_parameters = fopen("/home/processamento/acustica_2D/inputs/parameters.txt", "r");

    if(file_parameters == NULL){
        printf("Erro ao abrir arquivo de parametros\n");
        return 1;
    }

    char linha[256];

    float dx = 0.0f;
    float dz = 0.0f;
    float dt = 0.0f;
    float c1 = 0.0f;
    float c2 = 0.0f;

    int nx = 0;
    int nz = 0;
    int interface_Z = 0;
    int Nboudary = 0;


    while(fgets(linha, sizeof(linha), file_parameters)){

        sscanf(linha, "nx = %d", &nx);
        sscanf(linha, "nz  = %d", &nz);
        sscanf(linha, "dx = %f", &dx);
        sscanf(linha, "dz = %f", &dz);
        sscanf(linha, "dt = %f", &dt);

        sscanf(linha, "c1 = %f", &c1);
        sscanf(linha, "c2 = %f", &c2);

        sscanf(linha, "interface_Z = %d", &interface_Z);
        sscanf(linha, "Nboudary = %d", &Nboudary);
    }

    fclose(file_parameters);

    float *c = velocity(nx, nz, c1, c2, interface_Z);

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
// Save the velocity model as binary (.bin)
//------------------------------------------

    FILE *file_velocities = fopen("/home/processamento/acustica_2D/inputs/velocityModel.bin", "wb"); // "wb" = write binary

    if(file_velocities == NULL){
        printf("Erro ao abrir velocityModel.bin\n");
        return 1;
    }

    // Grava o array inteiro de uma vez (mais rápido que escrever elemento a elemento)
    fwrite(c, sizeof(float), nx * nz, file_velocities);

    fclose(file_velocities);

    std::cout << "Velocity bin file saved!" << std::endl;

    free(c);

    return 0;
}