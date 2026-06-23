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
    float f0 = 15.0; //dominant frequency

    //Cerjan boudary

    int Nboudary = 60; //number of the edge points
        
    int nx = int(L/dx) + 1; //number of spatial points in X
    int nz = int(L/dz) + 1; //number of spatial points in Z
    int nt = int(T/dt) + 1; //number of temporal steps

    //fountain position

    int sx = nx/2;
    int sz = 100;

    //----------------------------------
    // PARAMETERS DATA in json
    //----------------------------------

    FILE *file = fopen("/home/processamento/acustica_2D/inputs/parameters.json", "w");

    fprintf(file, "{\n");
    fprintf(file, "    \"L\": %.1f,\n", L);
    fprintf(file, "    \"T\": %d,\n", T);
    fprintf(file, "    \"nx\": %d,\n", nx);
    fprintf(file, "    \"nz\": %d,\n", nz);
    fprintf(file, "    \"nt\": %d,\n", nt);
    fprintf(file, "    \"dx\": %.1f,\n", dx);
    fprintf(file, "    \"dz\": %.1f,\n", dz);
    fprintf(file, "    \"dt\": %.1f\n", dt);
    fprintf(file, "    \"f0\": %.1f\n", f0);
    fprintf(file, "    \"Nboudary\": %d,\n", Nboudary);
    fprintf(file, "    \"sx\": %d,\n", sx);
    fprintf(file, "    \"sz\": %d,\n", sz);

    fprintf(file, "    \"sources_file\": \"sources.csv\",\n");
    fprintf(file, "    \"receivers_file\": \"receivers.csv\"\n");

    fprintf(file, "}\n");

    fclose(file);

    return 0;
}
