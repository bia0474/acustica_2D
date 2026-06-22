#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>

/* 

bool checkSources(float *sx, float *sz, int Nsource, int nx, int nz, int Nboundary){

    for(int i = 0; i < Nsource; i++){

        if(sx[i] < Nboundary || sx[i] >= nx - Nboundary || sz[i] < Nboundary || sz[i] >= nz - Nboundary){

            std::cout << "Erro na posição da fonte " << i << std::endl;

            return false;
        }
    }

    return true;
}

bool checkReceivers(float *rx, float *rz, int nx, int nz, int Nboundary){

    for(int i = 0; i < Nsource; i++){

        if(rx[i] < Nboundary || rx[i] >= nx - Nboundary || rz[i] < Nboundary || rz[i] >= nz - Nboundary){

            std::cout << "Erro na posição do receptor " << i << std::endl;

            return false;
        }
    }

    return true;
}

*/

//-------------------------------
// Struct of the receivers
//-------------------------------

    typedef struct{
        int x;
        int z;
    } Receiver;

//-------------------------------
// Cerjan - absorving boudanry
//-------------------------------

float* AbsorbingBoudanry(int Nboudary, int nx, int nz, const float* A){ //it returns a vector f that will be injected into future and present fields to decrease the energy and consequently the amplitude of the wave
    
    float *f = (float*) malloc(nx * nz * sizeof(float)); //matriz of size nx with all values equal to 1 (f(x) = 1 -> withuot cushioning)

    if(f == NULL){ //checks if memory has been allocated
        return NULL; //null means it's not pointing anywhere
    }

    for(int i = 0; i < nx * nz; i++){

        f[i] = 1.0f;
    }

    for(int x = 0; x < Nboudary; x++){//Left

        for(int z = 0; z < nz; z++){

            f[x * nz + z] *= A[x];
        }
    }

    for(int x = nx - Nboudary; x < nx; x++){//right

        int k = nx - 1 - x;

        for(int z = 0; z < nz; z++){

            f[x * nz + z] *= A[k];
        }
    }   

    for(int z = 0; z < Nboudary; z++){ //Top

        for(int x = 0; x < nx; x++){

            f[x * nz + z] *= A[z];
        }
    }

    for(int z = nz - Nboudary; z < nz; z++){

        int k = nz - 1 - z;

        for(int x = 0; x < nx; x++){ //Base

            f[x * nz + z] *= A[k];
        }
    }

    return f;
}

//-------------------------------
// Cerjan Vector
//-------------------------------

float* createCerjanVector(int Nboudary){ //generates the damping coefficients

    float Sb = 6.0f * Nboudary; //parameter that controls the width of the damping

    float *A = (float*) malloc(Nboudary * sizeof(float)); //stores the coefficients

    if(A == NULL){ //checks if memory has been allocated
        return NULL; //null means it's not pointing anywhere
    }

    for(int i = 0; i < Nboudary; i++){

        float fb = (float) (Nboudary - i) / (1.4142f * Sb); //for each position of the absorbent layer, a normalized distance is calculated

        A[i] = std::exp(-fb * fb); //the coefficients follow a Gaussian curve, where a smooth transition occurs
    }
    
    return A;   
}

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

//----------------------------------
// Ricker source
//----------------------------------

float* source(float f0, const float* t, int nt){

    float *s = (float*) malloc(nt * sizeof(float)); 

    float t0 = 1.0 / f0; //wavelet time delay

    for(int n = 0; n < nt; n++){

        float a = M_PI * M_PI * f0 * f0 * pow(t[n] - t0, 2);

        s[n] = (1.0 - 2.0 * a) * std::exp(-a); //wavelet equation 
    }

    return s;
}

//----------------------------------
// Wave equation
//----------------------------------


float* derivates(const float* c, float dt, float dx, float dz, const float* fonte, int nx, int nz, int nt, const float* f, int Nboudary, int sx, int sz, Receiver *receivers, int nrec){

    float *u_old = (float*) malloc(nx * nz * sizeof(float)); //passed field
    float *u_curr = (float*) malloc(nx * nz * sizeof(float)); //present field
    float *u_next = (float*) malloc(nx * nz * sizeof(float)); //future field

    for(int i = 0; i < nx * nz; i++){

        u_old[i] = 0.0f;
        u_curr[i] = 0.0f;
        u_next[i] = 0.0f;
    }

//----------------------------------
// Courant number for speeds
//----------------------------------

   float *e = (float*) malloc(nx * nz * sizeof(float));

    for(int i = 0; i < nx; i++){ 

        for(int j = 0; j < nz; j++){

            e[i * nz + j] = c[i * nz + j] * dt / dx; //because dx == dz !!
        }
    }


//----------------------------------
// SEISMOGRAM
//----------------------------------
    //stores seismic traces (nrec x nt)

    float *seismogram = (float*) malloc(nrec * nt * sizeof(float));

//----------------------------------
// time loop - 2nd order
//----------------------------------

    for(int n = 1; n < nt; n++){ //each iteration calculates the wave at the next instant

        std::fill(u_next, u_next + nx * nz, 0.0f); //resets the futures field before the next calculation

        //----------------------------------
        // space loop - 4nd order
        //----------------------------------

        for(int j = 2; j < nx - 2; j++){ //traverses all points of the grid in X

            for(int i = 2; i < nz - 2; i++){ //traverses all points of the grid in Z

                float d2x = (-u_curr[(j + 2) * nz + i] + 16 * u_curr[(j + 1) * nz + i] - 30 * u_curr[j * nz + i] + 16 * u_curr[(j - 1) * nz + i] -u_curr[(j - 2) * nz + i])/(12 * dx * dx);

                float d2z = (-u_curr[j * nz + (i + 2)] + 16 * u_curr[j * nz + (i + 1)] - 30 * u_curr[j * nz + i] + 16 * u_curr[j * nz + (i-1)] - u_curr[j * nz + (i - 2)])/(12 * dz * dz);

                u_next[j * nz + i] = 2 * u_curr[j * nz + i] - u_old[j * nz + i] + c[j * nz + i] * c[j * nz + i] * dt * dt * (d2x + d2z);
            }
        }

        //----------------------------------
        // source injection
        //----------------------------------

        //adds energy to the grid
        u_next[sx * nz + sz] += fonte[n]; 

        //----------------------------------
        // CERJAN
        //----------------------------------
        
        for(int j = 0; j < nx; j++){ //The amplitude of each field at each point gradually decreases

            for(int i = 0; i < nz; i++){

                u_next[j * nz + i] *= f[j * nz + i];
                u_curr[j * nz + i] *= f[j * nz + i];
            }
        }
        
        //----------------------------------
        // save the receiver
        //----------------------------------

        for(int i = 0; i < nrec; i++){

            int xr = receivers[i].x;
            int zr = receivers[i].z;

            seismogram[i * nt + n] = u_next[xr * nz + zr];

        }

        //-------------------------------------
        // SAVE SNAPSHOT HERE (binary document)
        //------------------------------------
         
         if(n % 100 == 0){

            std::ofstream file("snapshot_" + std::to_string(n) + ".bin", std::ios::binary);

            file.write(reinterpret_cast<char*>(u_next), nx * nz * sizeof(float));

            file.close();
        }
        //----------------------------------
        // advance in time
        //----------------------------------

        float *tmp = u_old; //temporary pointer to save the memory adress of u_old

        u_old = u_curr;
        u_curr = u_next;
        u_next = tmp;
    }

    //-----------------------------------
    // SAVE THE DOCUMENT OF THE SISMOGRAM
    //-----------------------------------
            
    std::ofstream file("seismogram.bin", std::ios::binary);

    file.write(reinterpret_cast<char*>(seismogram), nrec * nt * sizeof(float));

    file.close();

    std::cout << "Seismogram binary file saved!" << std::endl;

    //-----------------------------------
    // SAVE the copy of the final field
    //-----------------------------------
                
    float *result = (float*) malloc(nx * nz *sizeof(float));

    for(int i = 0; i < nx * nz; i++){

        result[i] = u_curr[i];
    }

    free(u_old);
    free(u_curr);
    free(u_next);
    free(seismogram);

    return result;
}

//----------------------------------
// MAIN
//----------------------------------

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

    //speeds

    float c1 = 1500.0f;
    float c2 = 4000.0f;

    //Cerjan boudary

    int Nboudary = 60; //number of the edge points
    
    int nx = int(L/dx) + 1; //number of spatial points in X
    int nz = int(L/dz) + 1; //number of spatial points in Z
    int nt = int(T/dt) + 1; //number of temporal steps

    //interface em Z

    int interface_Z = nz/2;
    
    float *x = linspace(0.0, nx, nx);
    float *z = linspace(0.0, nz, nz);
    float *t = linspace(0.0, (nt - 1) * dt, nt);

    //fountain position

    int sx = nx/2;
    int sz = 100;

//----------------------------------
// RECEIVERS
//----------------------------------
    //the quantity of the receivers
    int nrec = nx - 2 * Nboudary;

    Receiver *receivers = (Receiver*) malloc(nrec * sizeof(Receiver));

    for(int i = 0; i < nrec; i++){
        receivers[i].x = Nboudary + i;
        receivers[i].z = 60;
    }

//user can't change

//velocity

float *c = velocity(nx, nz, c1, c2, interface_Z);

//source

float *fonte = source(f0, t, nt);  
//------------------------------------
//  CERJAN
//------------------------------------

float *A = createCerjanVector(Nboudary);

float *f = AbsorbingBoudanry(Nboudary, nx, nz, A);

//----------------------------------
// CFL check
//----------------------------------

    if(CFL(c, dt, dx, dz, nx, nz)){

        std::cout << "Stable simulation" << std::endl;
    }
    else{

        return 1;
    }

//----------------------------------
// simulation 
//----------------------------------

    float *wavefield = derivates(c, dt, dx, dz, fonte, nx, nz, nt, f, Nboudary, sx, sz, receivers, nrec); 

    
//------------------------------------------
// Save the documento of the velocity model
//-----------------------------------------

    std::ofstream file_vel("velocity.bin", std::ios::binary);

    file_vel.write(reinterpret_cast<char*>(c), nx * nz * sizeof(float)); //"Pegue esse endereço e trate-o como um ponteiro para bytes."

    file_vel.close();

    std::cout << "Velocity Binary file saved!" << std::endl;

//---------------------------------------
// Save binary document of the simulation
//---------------------------------------

    std::ofstream file("wave.bin", std::ios::binary);

    file.write(reinterpret_cast<char*>(wavefield), nx * nz * sizeof(float)); //"Pegue os bytes que formam a matriz wavefield e grave-os no arquivo exatamente como estão na memória."

    file.close();

    std::cout << "Wavefield binary file saved!" << std::endl;

    free(wavefield);
    free(f);
    free(A);
    free(fonte);
    free(c);
    free(x);
    free(z);
    free(t);
    free(receivers);

    return 0;

}