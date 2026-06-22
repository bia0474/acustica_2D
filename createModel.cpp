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
// Save the documento of the velocity model
//-----------------------------------------

    std::ofstream file_vel("velocity.bin", std::ios::binary);

    file_vel.write(reinterpret_cast<char*>(c), nx * nz * sizeof(float)); //"Pegue esse endereço e trate-o como um ponteiro para bytes."

    file_vel.close();

    std::cout << "Velocity Binary file saved!" << std::endl;

    free(c);

    return 0;
}