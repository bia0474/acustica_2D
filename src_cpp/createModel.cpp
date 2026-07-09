#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>


//-------------------------------
// media speeds
//-------------------------------

float* velocity(int nx_abc, int nz_abc, float c1, float c2, int interface_Z){

    float *velocity = (float*) malloc(nx_abc * nz_abc * sizeof(float));

    for(int i = 0; i < nx_abc; i++){

        for(int j = 0; j < nz_abc; j++){

            if(j < interface_Z){
                velocity[i * nz_abc + j] = c1;
            }
            else{
                velocity[i * nz_abc + j] = c2;
            }
        }
    }

    return velocity;
}

//----------------------------------
// CFL condition
//----------------------------------


bool CFL(const float* c, float dt, float dx, float dz, int nx_abc, int nz_abc){ //function of the stability codition

    float cmax = 0.0f;

    for(int i = 0; i < nx_abc; i++){

        for(int j = 0; j < nz_abc; j++){

            cmax = std::max(cmax, c[i * nz_abc + j]);
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

    int nx_abc = 0;
    int nz_abc = 0;
    int interface_Z = 0;


    while(fgets(linha, sizeof(linha), file_parameters)){

        sscanf(linha, "nx_abc = %d", &nx_abc);
        sscanf(linha, "nz_abc = %d", &nz_abc);
        sscanf(linha, "dx = %f", &dx);
        sscanf(linha, "dz = %f", &dz);
        sscanf(linha, "dt = %f", &dt);

        sscanf(linha, "c1 = %f", &c1);
        sscanf(linha, "c2 = %f", &c2);

        sscanf(linha, "interface_Z = %d", &interface_Z);
    }

    fclose(file_parameters);

    float *c = velocity(nx_abc, nz_abc, c1, c2, interface_Z);

//----------------------------------
// CFL check
//----------------------------------

    if(CFL(c, dt, dx, dz, nx_abc, nz_abc)){

        std::cout << "Stable simulation" << std::endl;
    }
    else{

        return 1;
    }

//------------------------------------------
// Save the documento of the velocity model
//-----------------------------------------

    FILE *file_velocities = fopen("/home/processamento/acustica_2D/inputs/velocityModel.csv", "w"); //Creates a pointer to a file and opens a file named "sources.csv" in write mode (w)

    if(file_velocities == NULL){
        printf("Erro ao abrir sources.csv\n");
        return 1;
    }

    fprintf(file_velocities, "c\n"); //writing the header

    for(int i = 0; i < nx_abc; i++){
        for(int j = 0; j < nz_abc; j++){
            
            fprintf(file_velocities, "%.1f", c[i * nz_abc + j]);

            if(j < nz_abc - 1){
                fprintf(file_velocities, ",");
            }
        }

        fprintf(file_velocities, "\n");
    }

    fclose(file_velocities);

    std::cout << "Velocity csv file saved!" << std::endl;

    free(c);

    return 0;
}