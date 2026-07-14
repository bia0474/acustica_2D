#include <stdlib.h>
#include <stdio.h>
#include <cstdio> 

//----------------------------------
// linscpace function
//----------------------------------

float* linspace(float start, int end, int quantity){//this function calculates the step between the points and fills the vector accordingly

    float *number = (float*) malloc(quantity * sizeof(float));

    float dx = (end - start) / (quantity - 1);

    for(int i = 0; i < quantity; i++){

        number[i] = start + i * dx;
    }

    return number;
}

int main(){

    //----------------------------------
    // model parameters
    //----------------------------------

    float L = 5000.0; //model size
    int T = 2.0; //total simulation time

    float dx = 10.0; //space step in X
    float dz = 10.0; //space step in Z
    float dt = 0.0005; //time lapse 

    float c1 = 1500.0f;
    float c2 = 4000.0f;

    float f0 = 15.0; //dominant frequency

    //Cerjan boudary

    int Nboudary = 100; //number of the edge points
        
    int nx = int(L/dx) + 1; //number of spatial points in X
    int nz = int(L/dz) + 1; //number of spatial points in Z
    int nt = int(T/dt) + 1; //number of temporal steps

    int nx_abc = nx + 2 * Nboudary;
    int nz_abc = nz + 2 * Nboudary; //expeding the model

    int interface_Z = nz_abc/2;

    float *x = linspace(0.0, nx, nx);
    float *z = linspace(0.0, nz, nz);
    float *t = linspace(0.0, (nt - 1) * dt, nt);

    int nrec = 381;
    int Nsource = 1;

    //----------------------------------
    // PARAMETERS DATA in txt
    //----------------------------------

    FILE *file_parameters = fopen("/home/processamento/acustica_2D/inputs/parameters.txt", "w");

    if(file_parameters == NULL){
        printf("Erro ao criar arquivo de parametros\n");
        return 1;
    }

    fprintf(file_parameters, "T = %d\n", T);
    
    fprintf(file_parameters, "nx = %d\n", nx);
    fprintf(file_parameters, "nz = %d\n", nz);
    fprintf(file_parameters, "nx_abc = %d\n", nx_abc);
    fprintf(file_parameters, "nz_abc = %d\n", nz_abc);
    fprintf(file_parameters, "nt = %d\n", nt);

    fprintf(file_parameters, "dx = %.1f\n", dx);
    fprintf(file_parameters, "dz = %.1f\n", dz);
    fprintf(file_parameters, "dt = %.6f\n", dt);

    fprintf(file_parameters, "c1 = %.1f\n", c1);
    fprintf(file_parameters, "c2 = %.1f\n", c2);

    fprintf(file_parameters, "interface_Z = %d\n", interface_Z);

    fprintf(file_parameters, "f0 = %.1f\n", f0);

    fprintf(file_parameters, "Nboudary = %d\n", Nboudary);
    fprintf(file_parameters, "nrec = %d\n", nrec);
    fprintf(file_parameters, "Nsource = %d\n", Nsource);

    fprintf(file_parameters, "sources_file = /home/processamento/acustica_2D/inputs/sources.csv\n");
    fprintf(file_parameters, "receivers_file = /home/processamento/acustica_2D/inputs/receivers.csv\n");
    fprintf(file_parameters, "velocity_file = /home/processamento/acustica_2D/inputs/velocityModel.csv\n");

    fclose(file_parameters);

    return 0;
}
