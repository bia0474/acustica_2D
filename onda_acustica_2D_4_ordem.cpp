#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

//-------------------------------
// Cerjan - absorving boudanry
//-------------------------------

std::vector<std::vector<float>> AbsorbingBoudanry(int Nboudary, int nx, int nz, const std::vector<float>& A){ //it returns a vector f that will be injected into future and present fields to decrease the energy and consequently the amplitude of the wave
    
    std::vector<std::vector<float>> f(nx, std::vector<float>(nz, 1.0f)); //matriz of size nx with all values equal to 1 (f(x) = 1 -> withuot cushioning)

    for(int x = 0; x < Nboudary; x++){//Esquerda
        for(int z = 0; z < nz; z++){
            f[x][z] *= A[x];
        }
    }
    for(int x = nx - Nboudary; x < nx; x++){//Direita

        int k = nx - 1 - x;

        for(int z = 0; z < nz; z++){
            f[x][z] *= A[k];
        }
    }   
    for(int z = 0; z < Nboudary; z++){ //Topo
        for(int x = 0; x < nx; x++){
            f[x][z] *= A[z];
        }
    }
    for(int z = nz - Nboudary; z < nz; z++){

        int k = nz - 1 - z;

        for(int x = 0; x < nx; x++){ //Base
            f[x][z] *= A[k];
        }
    }

    return f;
}

//-------------------------------
// Cerjan Vector
//-------------------------------

std::vector<float> createCerjanVector(int Nboudary){ //generates the damping coefficients

    float Sb = 6.0f * Nboudary; //parameter that controls the width of the damping

    std::vector<float> A(Nboudary); //stores the coefficients

    float fb;

    for(int i = 0; i < Nboudary; i++){
        fb = (Nboudary - i) / (1.4142f * Sb); //for each position of the absorbent layer, a normalized distance is calculated
        A[i] = std::exp(-fb * fb); //the coefficients follow a Gaussian curve, where a smooth transition occurs
    }
    
    return A;   
}

//-------------------------------
// media speeds
//-------------------------------

std::vector<std::vector<float>> velocity(int nx, int nz){

    float c1 = 1500.0f;
    float c2 = 2000.0f;

    int interface_Z = nz/2;

    std::vector<std::vector<float>> velocity(nx, std::vector<float>(nz));

    for(int i = 0; i < nx; i++){
        for(int j = 0; j < nz; j++){

            if(j < interface_Z){
                velocity[i][j] = c1;
            }
            else{
                velocity[i][j] = c2;
            }
        }
    }

    return velocity;
}

//----------------------------------
// linscpace function
//----------------------------------

std::vector<float> linspace(float start, int end, int quantity){//this function calculates the step between the points and fills the vector accordingly

    std::vector<float> number(quantity);

    float dx = (end - start) / (quantity - 1);

    for(int i = 0; i < quantity; i++){
        number[i] = start + i * dx;
    }

    return number;
}

//----------------------------------
// CFL condition
//----------------------------------


bool CFL(const std::vector<std::vector<float>>& c, float dt, float dx, float dz, int nx, int nz){ //function of the stability codition

    float cmax = 0.0f;

    for(int i = 0; i < nx; i++){
        for(int j = 0; j < nz; j++){
            cmax = std::max(cmax, c[i][j]);
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

std::vector<float> source(float f0, const std::vector<float>& t){

    int nt = t.size(); //total number of samples over time

    std::vector<float> s(nt);

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


std::vector<std::vector<float>> derivates(std::vector<std::vector<float>>& c, float dt, float dx, float dz, const std::vector<float>& fonte, int nx, int nz, int nt, const std::vector<std::vector<float>>& f, int Nboudary){

    
     std::vector<std::vector<float>> u_old(nx, std::vector<float>(nz, 0.0)); //passed field
     std::vector<std::vector<float>> u_curr(nx, std::vector<float>(nz, 0.0)); //present field
     std::vector<std::vector<float>> u_next(nx, std::vector<float>(nz, 0.0)); //future field

//----------------------------------
// Courant number for speeds
//----------------------------------

    std::vector<std::vector<float>> e(nx, std::vector<float>(nz));

    for(int i = 0; i < nx; i++){ 
        for(int j = 0; j < nz; j++){
            e[i][j] = c[i][j] * dt / dx; //because dx == dz !!
        }
    }

//----------------------------------
// fountain position
//----------------------------------

    int sx = nx/2;
    int sz = 100;

//----------------------------------
// RECEIVERS
//----------------------------------
    struct Receiver{
    int x;
    int z;
    };

    std::vector<Receiver> receivers;

    for(int x = Nboudary; x < nx - Nboudary; x++){
        receivers.push_back({x, 60});
    }
//----------------------------------
// SEISMOGRAM
//----------------------------------
    //the quantity of the receivers
    int nrec = receivers.size();    

    //stores seismic traces (nrec x nt)
    std::vector<std::vector<float>> seismogram(nrec, std::vector<float>(nt, 0.0));

//----------------------------------
// time loop - 2nd order
//----------------------------------

    for(int n = 1; n < nt; n++){ //each iteration calculates the wave at the next instant

        for(int k = 0; k < nx; k++){
            std::fill(u_next[k].begin(), u_next[k].end(), 0.0f); //resets the futures field before the next calculation
        } 

        //----------------------------------
        // space loop - 4nd order
        //----------------------------------

        for(int j = 2; j < nx - 2; j++){ //traverses all points of the grid in X
            for(int i = 2; i < nz - 2; i++){ //traverses all points of the grid in Z
                float d2x = (-u_curr[j+2][i] + 16 * u_curr[j+1][i] - 30 * u_curr[j][i] + 16 * u_curr[j-1][i] -u_curr[j-2][i])/(12 * dx * dx);
                float d2z = (-u_curr[j][i+2] + 16 * u_curr[j][i+1] - 30 * u_curr[j][i] + 16 * u_curr[j][i-1] - u_curr[j][i-2])/(12 * dz * dz);

                u_next[j][i] = 2 * u_curr[j][i] - u_old[j][i] + c[j][i] * c[j][i] * dt * dt * (d2x + d2z);
            }
        }

        //----------------------------------
        // source injection
        //----------------------------------

        //adds energy to the grid
        u_next[sx][sz] += fonte[n]; 

        //----------------------------------
        // CERJAN
        //----------------------------------
        
        for(int j = 0; j < nx; j++){ //The amplitude of each field at each point gradually decreases
            for(int i = 0; i < nz; i++){
                u_next[j][i] *= f[j][i];
                u_curr[j][i] *= f[j][i];
            }
        }
        
        //----------------------------------
        // save the receiver
        //----------------------------------

        for(int ir = 0; ir < nrec; ir++){

            int xr = receivers[ir].x;
            int zr = receivers[ir].z;

            seismogram[ir][n] = u_next[xr][zr];
        }

        //----------------------------------
        // SAVE SNAPSHOT HERE (binary document)
        //----------------------------------

        if(n % 100 == 0){
            std::ofstream file("snapshot_" + std::to_string(n) + ".bin", std::ios::binary);

            for(int j = 0; j < nx; j++){
                file.write(
                 reinterpret_cast<char*>(u_next[j].data()), nz * sizeof(float));
            }

            file.close();
        }

        //----------------------------------
        // advance in time
        //----------------------------------

        u_old = u_curr;
        u_curr = u_next;
    }

    //-----------------------------------
    // SAVE THE DOCUMENT OF THE SISMOGRAM
    //-----------------------------------

    std::ofstream file("seismogram.bin", std::ios::binary);

    for(int ir = 0; ir < nrec; ir++){

        file.write(reinterpret_cast<char*>(seismogram[ir].data()), nt * sizeof(float));
    }

    file.close();

    return u_curr;
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
    float f0 = 30.0; //dominant frequency

    int nx = int(L/dx) + 1; //number of spatial points in X
    int nz = int(L/dz) + 1; //number of spatial points in Z
    int nt = int(T/dt) + 1; //number of temporal steps
    
    std::vector<float> x = linspace(0.0, nx, nx);
    std::vector<float> z = linspace(0.0, nz, nz);
    std::vector<float> t = linspace(0.0, (nt - 1) * dt, nt);

    //velocity

    std::vector<std::vector<float>> c = velocity(nx, nz);

    //Save the documento of the velocity model

    std::ofstream file_vel("velocity.bin", std::ios::binary);

    for(int i = 0; i < nx; i++){
        file_vel.write(reinterpret_cast<char*>(c[i].data()), nz * sizeof(float));
    }

    file_vel.close();

    //source

    std::vector<float> fonte;

    fonte = source(f0, t);

//------------------------------------
//  CERJAN
//------------------------------------

int Nboudary = 60; //number of the edge points

std::vector<float> A = createCerjanVector(Nboudary);

std::vector<std::vector<float>> f = AbsorbingBoudanry(Nboudary, nx, nz, A);

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

    std::vector<std::vector<float>> wavefield;

    wavefield = derivates(c, dt, dx, dz, fonte, nx, nz, nt, f, Nboudary); 

//---------------------------------------
// save binary document of the simulation
//---------------------------------------

    std::ofstream file("wave.bin", std::ios::binary);

    for(int j = 0; j < nx; j++){
        file.write(reinterpret_cast<char*>(wavefield[j].data()), nz * sizeof(float));
    }

    file.close();

    std::cout << "Binary file saved!" << std::endl;

    return 0;

}