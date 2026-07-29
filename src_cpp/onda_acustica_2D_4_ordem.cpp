#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

//-------------------------------
// Struct of the receivers
//-------------------------------

typedef struct
{
    int x;
    int z;
} Receiver;

//----------------------------------
// linscpace function
//----------------------------------

float *linspace(float start, int end, int quantity)
{ // this function calculates the step between the points and fills the vector accordingly

    float *number = (float *)malloc(quantity * sizeof(float));

    float dx = (end - start) / (quantity - 1);

    for (int i = 0; i < quantity; i++)
    {

        number[i] = start + i * dx;
    }

    return number;
}

//----------------------------------
// read parameters function
//----------------------------------

void readParameters(const char *filename, int *T, int *nx, int *nz, int *nx_abc, int *nz_abc, int *nt, float *dx, float *dz, float *dt, float *f0, int *Nboudary, int *Nsource, int *nrec, char receivers_file[], char sources_file[], char velocity_file[], float **x, float **z, float **t)
{

    FILE *file_parameters = fopen(filename, "r");

    if (file_parameters == NULL)
    {
        printf("Erro ao abrir arquivo de parametros\n");
        exit(1);
    }

    char linha[256];

    while (fgets(linha, sizeof(linha), file_parameters))
    {

        if (sscanf(linha, "T = %d", T) == 1)
            continue;

        if (sscanf(linha, "nx = %d", nx) == 1)
            continue;

        if (sscanf(linha, "nz = %d", nz) == 1)
            continue;

        if (sscanf(linha, "nx_abc = %d", nx_abc) == 1)
            continue;

        if (sscanf(linha, "nz_abc = %d", nz_abc) == 1)
            continue;

        if (sscanf(linha, "nt = %d", nt) == 1)
            continue;

        if (sscanf(linha, "dx = %f", dx) == 1)
            continue;

        if (sscanf(linha, "dz = %f", dz) == 1)
            continue;

        if (sscanf(linha, "dt = %f", dt) == 1)
            continue;

        if (sscanf(linha, "f0 = %f", f0) == 1)
            continue;

        if (sscanf(linha, "Nboudary = %d", Nboudary) == 1)
            continue;

        if (sscanf(linha, "nrec = %d", nrec) == 1)
            continue;

        if (sscanf(linha, "Nsource = %d", Nsource) == 1)
            continue;

        if (sscanf(linha, "receivers_file = %255s", receivers_file) == 1)
            continue;

        if (sscanf(linha, "sources_file = %255s", sources_file) == 1)
            continue;

        if (sscanf(linha, "velocity_file = %255s", velocity_file) == 1)
            continue;
    }

    fclose(file_parameters);

    *x = linspace(0.0f, *nx_abc, *nx_abc);
    *z = linspace(0.0f, *nz_abc, *nz_abc);
    *t = linspace(0.0f, (*nt - 1) * (*dt), *nt);
}

//----------------------------------
// read receivers function
//----------------------------------

Receiver *readReceivers(const char *receivers_file, int nrec, int Nboudary)
{

    Receiver *receivers = (Receiver *)malloc(nrec * sizeof(Receiver));

    if (receivers == NULL)
    {
        std::cout << "Erro ao alocar memoria para os receptores.\n";
        exit(1);
    }

    std::ifstream file(receivers_file);

    if (!file.is_open())
    {
        std::cout << "Erro ao abrir receivers.csv\n";
        free(receivers);
        exit(1);
    }

    std::string linha;

    std::getline(file, linha);

    int i = 0;

    while (std::getline(file, linha) && i < nrec)
    {

        std::stringstream ss(linha);

        std::string index, rx, rz;

        std::getline(ss, index, ',');
        std::getline(ss, rx, ',');
        std::getline(ss, rz, ',');

        receivers[i].x = std::stoi(rx) + Nboudary;
        receivers[i].z = std::stoi(rz) + Nboudary;

        i++;
    }

    file.close();

    return receivers;
}

//----------------------------------
// read sources function
//----------------------------------

void readSources(const char *sources_file, int Nsource, int **sx, int **sz, int Nboudary)
{

    *sx = (int *)malloc(Nsource * sizeof(int));
    *sz = (int *)malloc(Nsource * sizeof(int));

    if (*sx == NULL || *sz == NULL)
    {
        std::cout << "Erro ao alocar memoria para as fontes.\n";
        exit(1);
    }

    std::ifstream file(sources_file);

    if (!file.is_open())
    {
        std::cout << "Erro ao abrir sources.csv\n";
        free(*sx);
        free(*sz);
        exit(1);
    }

    std::string linha;

    std::getline(file, linha);

    int i = 0;

    while (std::getline(file, linha) && i < Nsource)
    {

        std::stringstream ss(linha);

        std::string index_str;
        std::string sx_str;
        std::string sz_str;

        std::getline(ss, index_str, ',');
        std::getline(ss, sx_str, ',');
        std::getline(ss, sz_str, ',');

        (*sx)[i] = std::stoi(sx_str) + Nboudary;
        (*sz)[i] = std::stoi(sz_str) + Nboudary;

        i++;
    }

    file.close();
}

//---------------------------------------
// read velocity model (function)
//---------------------------------------

float *readVelocity(const char *velocity_file, int nx, int nz, int nx_abc, int nz_abc, int Nboudary)
{

    FILE *file = fopen(velocity_file, "rb");

    if (file == NULL)
    {
        printf("Erro ao abrir o arquivo do modelo de velocidade.\n");
        exit(1);
    }

    // Lê no formato como está salvo: nz x nx (linha = z fixo, todas as posições x)
    float *c_raw = (float *)malloc(nx * nz * sizeof(float));

    fread(c_raw, sizeof(float), nx * nz, file);

    fclose(file);

    // Transpõe para o formato que o resto do código espera: nx x nz
    float *c = (float *)malloc(nx * nz * sizeof(float));

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            // c_raw está em ordem (nz, nx): índice = j * nx + i
            // c deve ficar em ordem (nx, nz): índice = i * nz + j
            c[i * nz + j] = c_raw[j * nx + i];
        }
    }

    free(c_raw);

    float *c_exp = (float *)calloc(nx_abc * nz_abc, sizeof(float));

    //----------------------------------
    // Centro
    //----------------------------------

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < nz; j++)
        {

            c_exp[(i + Nboudary) * nz_abc + (j + Nboudary)] = c[i * nz + j];
        }
    }

    //----------------------------------
    // Borda superior
    //----------------------------------

    for (int i = 0; i < Nboudary; i++)
    {
        for (int j = Nboudary; j < nz_abc - Nboudary; j++)
        {

            c_exp[i * nz_abc + j] = c_exp[Nboudary * nz_abc + j];
        }
    }

    //----------------------------------
    // Borda inferior
    //----------------------------------

    for (int i = nx_abc - Nboudary; i < nx_abc; i++)
    {
        for (int j = Nboudary; j < nz_abc - Nboudary; j++)
        {

            c_exp[i * nz_abc + j] = c_exp[(nx_abc - Nboudary - 1) * nz_abc + j];
        }
    }

    //----------------------------------
    // Borda esquerda
    //----------------------------------

    for (int i = Nboudary; i < nx_abc - Nboudary; i++)
    {
        for (int j = 0; j < Nboudary; j++)
        {

            c_exp[i * nz_abc + j] = c_exp[i * nz_abc + Nboudary];
        }
    }

    //----------------------------------
    // Borda direita
    //----------------------------------

    for (int i = Nboudary; i < nx_abc - Nboudary; i++)
    {
        for (int j = nz_abc - Nboudary; j < nz_abc; j++)
        {

            c_exp[i * nz_abc + j] = c_exp[i * nz_abc + (nz_abc - Nboudary - 1)];
        }
    }

    //----------------------------------
    // Canto superior esquerdo
    //----------------------------------

    for (int i = 0; i < Nboudary; i++)
    {
        for (int j = 0; j < Nboudary; j++)
        {

            c_exp[i * nz_abc + j] = c[0];
        }
    }

    //----------------------------------
    // Canto superior direito
    //----------------------------------

    for (int i = 0; i < Nboudary; i++)
    {
        for (int j = nz_abc - Nboudary; j < nz_abc; j++)
        {

            c_exp[i * nz_abc + j] = c[nz - 1];
        }
    }

    //----------------------------------
    // Canto inferior esquerdo
    //----------------------------------

    for (int i = nx_abc - Nboudary; i < nx_abc; i++)
    {
        for (int j = 0; j < Nboudary; j++)
        {

            c_exp[i * nz_abc + j] = c[(nx - 1) * nz];
        }
    }

    //----------------------------------
    // Canto inferior direito
    //----------------------------------

    for (int i = nx_abc - Nboudary; i < nx_abc; i++)
    {
        for (int j = nz_abc - Nboudary; j < nz_abc; j++)
        {

            c_exp[i * nz_abc + j] = c[(nx - 1) * nz + (nz - 1)];
        }
    }

    free(c);

    return c_exp;
}

//-------------------------------------
// read velocity model by me function
//-------------------------------------
/*
float* readVelocity(const char *velocity_file, int nx_abc, int nz_abc){

    float *c = (float*) malloc(nx_abc * nz_abc * sizeof(float));

    if(c == NULL){
        std::cout << "Erro ao alocar memoria para o modelo de velocidade.\n";
        exit(1);
    }

    std::ifstream file(velocity_file);

    if(!file.is_open()){
        std::cout << "Erro ao abrir o arquivo do modelo\n";
        free(c);
        exit(1);
    }

    std::string linha;

    std::getline(file, linha);

    int i = 0;

    while(std::getline(file, linha) && i < nx_abc){

        std::stringstream ss(linha);

        std::string velocity;

        int j = 0;

        while(std::getline(ss, velocity, ',') && j < nz_abc){

            c[i * nz_abc + j] = std::stof(velocity);

            j++;
        }

        i++;
    }

    file.close();

    return c;
}
*/

//----------------------------------
// check geometry function
//----------------------------------

bool checkGeometry(const int *sx, const int *sz, int Nsource, Receiver *receivers, int nrec, int nx, int nz, int Nboudary)
{

    for (int i = 0; i < Nsource; i++)
    {

        if (sx[i] < Nboudary || sx[i] >= nx + Nboudary || sz[i] < Nboudary || sz[i] >= nz + Nboudary)
        {

            std::cout << "Erro: Fonte " << i << " esta dentro da borda de absorcao.\n";
            return false;
        }
    }

    for (int j = 0; j < nrec; j++)
    {

        if (receivers[j].x < Nboudary || receivers[j].x >= nx + Nboudary || receivers[j].z < Nboudary || receivers[j].z >= nz + Nboudary)
        {

            std::cout << "Erro: Receptor " << j << " esta dentro da borda de absorcao.\n";
            return false;
        }
    }

    return true;
}

//-------------------------------
// Cerjan Vector
//-------------------------------

float *createCerjanVector(int Nboudary)
{ // generates the damping coefficients

    float Sb = 6.0f * Nboudary; // parameter that controls the width of the damping

    float *A = (float *)malloc(Nboudary * sizeof(float)); // stores the coefficients

    if (A == NULL)
    {                // checks if memory has been allocated
        return NULL; // null means it's not pointing anywhere
    }

#pragma omp parallel for
    for (int i = 0; i < Nboudary; i++)
    {

        float fb = (float)(Nboudary - i) / (1.4142f * Sb); // for each position of the absorbent layer, a normalized distance is calculated

        A[i] = std::exp(-fb * fb); // the coefficients follow a Gaussian curve, where a smooth transition occurs
    }

    return A;
}

//----------------------------------
// Ricker source
//----------------------------------

float *source(float f0, const float *t, int nt)
{

    float *s = (float *)malloc(nt * sizeof(float));

    float t0 = 1.0 / f0; // wavelet time delay

#pragma omp parallel for
    for (int n = 0; n < nt; n++)
    {

        float a = M_PI * M_PI * f0 * f0 * pow(t[n] - t0, 2);

        s[n] = (1.0 - 2.0 * a) * std::exp(-a); // wavelet equation
    }

    return s;
}

//----------------------------------
// Wave equation
//----------------------------------

float *derivates(float *c, float dt, float dx, float dz, const float *fonte, int nx, int nz, int nx_abc, int nz_abc, int nt, const float *A, int Nboudary, int *sx, int *sz, int Nsource, Receiver *receivers, int nrec)
{

    float *u_curr = (float *)calloc(nx_abc * nz_abc, sizeof(float)); // present field
    float *u_next = (float *)calloc(nx_abc * nz_abc, sizeof(float)); // future field

    //----------------------------------
    // vector of the PV
    //----------------------------------

    float *PVx = (float *)calloc(nx_abc * nz_abc, sizeof(float)); // vector d2P/dx2
    float *PVz = (float *)calloc(nx_abc * nz_abc, sizeof(float)); // vector d2P/dz2

    //----------------------------------
    // Optical Flow
    //----------------------------------

    float *px = (float *)calloc(nx_abc * nz_abc, sizeof(float));
    float *pz = (float *)calloc(nx_abc * nz_abc, sizeof(float));
    float *pt = (float *)calloc(nx_abc * nz_abc, sizeof(float));

    float *ux = (float *)calloc(nx_abc * nz_abc, sizeof(float));
    float *uz = (float *)calloc(nx_abc * nz_abc, sizeof(float));

    //----------------------------------
    // Courant number for speeds
    //----------------------------------

    std::cout << "Making the calculos of Courant number!" << std::endl;

    float *e = (float *)malloc(nx_abc * nz_abc * sizeof(float));

    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < nx_abc; i++)
    {

        for (int j = 0; j < nz_abc; j++)
        {

            e[i * nz_abc + j] = c[i * nz_abc + j] * dt / dx; // because dx == dz !!
        }
    }

    //----------------------------------
    // SEISMOGRAM
    //----------------------------------
    // stores seismic traces (nrec x nt)

    float *seismogram = (float *)calloc(nrec * nt, sizeof(float));

    //----------------------------------
    // time loop - 2nd order
    //----------------------------------

    std::cout << "Starting the temporal and spacial loops!" << std::endl;

    for (int n = 1; n < nt; n++)
    { // each iteration calculates the wave at the next instant

        //----------------------------------
        // space loop - 4nd order
        //----------------------------------

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 2; j < nx_abc - 2; j++)
        { // traverses all points of the grid in X

            for (int i = 2; i < nz_abc - 2; i++)
            { // traverses all points of the grid in Z

                //----------------------------------
                // finite dofferences
                //----------------------------------

                float d2x = (-u_curr[(j + 2) * nz_abc + i] + 16 * u_curr[(j + 1) * nz_abc + i] - 30 * u_curr[j * nz_abc + i] + 16 * u_curr[(j - 1) * nz_abc + i] - u_curr[(j - 2) * nz_abc + i]) / (12 * dx * dx);

                float d2z = (-u_curr[j * nz_abc + (i + 2)] + 16 * u_curr[j * nz_abc + (i + 1)] - 30 * u_curr[j * nz_abc + i] + 16 * u_curr[j * nz_abc + (i - 1)] - u_curr[j * nz_abc + (i - 2)]) / (12 * dz * dz);

                u_next[j * nz_abc + i] = 2 * u_curr[j * nz_abc + i] - u_next[j * nz_abc + i] + c[j * nz_abc + i] * c[j * nz_abc + i] * dt * dt * (d2x + d2z);

                if (n % 500 == 0)
                {
                    //----------------------------------
                    // Poynting vectors
                    //----------------------------------

                    float dUdt = (u_next[j * nz_abc + i] - u_curr[j * nz_abc + i]) / dt;

                    float Ux = (u_curr[(j - 2) * nz_abc + i] - 8 * u_curr[(j - 1) * nz_abc + i] + 8 * u_curr[(j + 1) * nz_abc + i] - u_curr[(j + 2) * nz_abc + i]) / (12 * dx);

                    float Uz = (u_curr[j * nz_abc + (i - 2)] - 8 * u_curr[j * nz_abc + (i - 1)] + 8 * u_curr[j * nz_abc + (i + 1)] - u_curr[j * nz_abc + (i + 2)]) / (12 * dz);

                    PVx[j * nz_abc + i] = -dUdt * Ux;
                    PVz[j * nz_abc + i] = -dUdt * Uz;

                    //--------------------------------------
                    // saves the derivatives to Optical Flow
                    //--------------------------------------

                    px[j * nz_abc + i] = Ux;
                    pz[j * nz_abc + i] = Uz;
                    pt[j * nz_abc + i] = dUdt;
                }

            }
        }

        //---------------------------------------------------------------
        // Using Optical Flow method (20 iterations over the entire mesh)
        //---------------------------------------------------------------

        if (n % 500 == 0)
        {

            std::fill(ux, ux + nx_abc * nz_abc, 0.0f);
            std::fill(uz, uz + nx_abc * nz_abc, 0.0f);

            int n_iter = 20;
            float alpha = 1.0f;

            for (int iter = 0; iter < n_iter; iter++)
            {
                for (int j = 2; j < nx_abc - 2; j++)
                {
                    for (int i = 2; i < nz_abc - 2; i++)
                    {
                        float somaux = 0.0f;
                        float somauz = 0.0f;

                        for (int a = -1; a <= 1; a++)
                        {
                            for (int b = -1; b <= 1; b++)
                            {
                                somaux += ux[(j + a) * nz_abc + (i + b)];
                                somauz += uz[(j + a) * nz_abc + (i + b)];
                            }
                        }

                        float ux_average = (1.0f / 12.0f) * (ux[j * nz_abc + i - nz_abc] + ux[j * nz_abc + i + nz_abc] + ux[j * nz_abc + i - 1] + ux[j * nz_abc + i + 1] - ux[j * nz_abc + i] + somaux);
                        float uz_average = (1.0f / 12.0f) * (uz[j * nz_abc + i - nz_abc] + uz[j * nz_abc + i + nz_abc] + uz[j * nz_abc + i - 1] + uz[j * nz_abc + i + 1] - uz[j * nz_abc + i] + somauz);

                        float denominator = alpha * alpha + px[j * nz_abc + i] * px[j * nz_abc + i] + pz[j * nz_abc + i] * pz[j * nz_abc + i];

                        ux[j * nz_abc + i] = ux_average - (px[j * nz_abc + i] * (px[j * nz_abc + i] * ux_average + pz[j * nz_abc + i] * uz_average + pt[j * nz_abc + i]) / denominator);
                        uz[j * nz_abc + i] = uz_average - (pz[j * nz_abc + i] * (px[j * nz_abc + i] * ux_average + pz[j * nz_abc + i] * uz_average + pt[j * nz_abc + i]) / denominator);
                    }
                }
            }
        }

        //----------------------------------
        // source injection
        //----------------------------------

        for (int k = 0; k < Nsource; k++)
        {
            // adds energy to the grid
            u_next[sx[k] * nz_abc + sz[k]] += fonte[n];
        }

        //-----------------------------------
        // CERJAN 
        //-----------------------------------

        if (n == 1)
        {
            std::cout << "Making the CERJAN boudary" << std::endl;
        }

        #pragma omp parallel for collapse(2)
        for (int x = 0; x < Nboudary; x++)
        { // Left

            for (int z = 0; z < nz_abc; z++)
            {

                u_next[x * nz_abc + z] *= A[x];
                u_curr[x * nz_abc + z] *= A[x];
            }
        }

        for (int x = nx_abc - Nboudary; x < nx_abc; x++)
        { // right

            int k = nx_abc - 1 - x;

            for (int z = 0; z < nz_abc; z++)
            {

                u_next[x * nz_abc + z] *= A[k];
                u_curr[x * nz_abc + z] *= A[k];
            }
        }

        #pragma omp parallel for collapse(2)
        for (int z = 0; z < Nboudary; z++)
        { // Top

            for (int x = 0; x < nx_abc; x++)
            {

                u_next[x * nz_abc + z] *= A[z];
                u_curr[x * nz_abc + z] *= A[z];
            }
        }

        for (int z = nz_abc - Nboudary; z < nz_abc; z++)
        {

            int k = nz_abc - 1 - z;

            for (int x = 0; x < nx_abc; x++)
            { // Base

                u_next[x * nz_abc + z] *= A[k];
                u_curr[x * nz_abc + z] *= A[k];
            }
        }

        //----------------------------------
        // save the receiver
        //----------------------------------

        #pragma omp parallel for
        for (int i = 0; i < nrec; i++)
        {

            int xr = receivers[i].x;
            int zr = receivers[i].z;

            seismogram[i * nt + n] = u_next[xr * nz_abc + zr];
        }

        //-------------------------------------------------------
        // SAVE ALL THE SNAPSHOT HERE WITH PVxz (binary document)
        //-------------------------------------------------------

        if (n == 1)
        {
            std::cout << "Saving the file of the snapshots, PVs and PVOFs!" << std::endl;
        }

        if (n % 500 == 0)
        {

            std::ofstream file("/home/processamento/acustica_2D/outputs/snapshot_" + std::to_string(n) + ".bin", std::ios::binary);

            std::ofstream file_PVx("/home/processamento/acustica_2D/outputs/PoyntingVectorDirectionX" + std::to_string(n) + ".bin", std::ios::binary);
            std::ofstream file_PVz("/home/processamento/acustica_2D/outputs/PoyntingVectorDirectionZ" + std::to_string(n) + ".bin", std::ios::binary);

            std::ofstream file_PVxOF("/home/processamento/acustica_2D/outputs/PoyntingVectorOFx" + std::to_string(n) + ".bin", std::ios::binary);
            std::ofstream file_PVzOF("/home/processamento/acustica_2D/outputs/PoyntingVectorOFz" + std::to_string(n) + ".bin", std::ios::binary);

            for (int x = Nboudary; x < nx_abc - Nboudary; x++)
            {

                file.write(reinterpret_cast<char *>(&u_next[x * nz_abc + Nboudary]), (nz_abc - 2 * Nboudary) * sizeof(float)); // saves snaps without the absorbent border

                file_PVx.write(reinterpret_cast<char *>(&PVx[x * nz_abc + Nboudary]), (nz_abc - 2 * Nboudary) * sizeof(float)); // saves PV values

                file_PVz.write(reinterpret_cast<char *>(&PVz[x * nz_abc + Nboudary]), (nz_abc - 2 * Nboudary) * sizeof(float)); // saves PV values

                file_PVxOF.write(reinterpret_cast<char *>(&ux[x * nz_abc + Nboudary]), (nz_abc - 2 * Nboudary) * sizeof(float)); // saves PV values

                file_PVzOF.write(reinterpret_cast<char *>(&uz[x * nz_abc + Nboudary]), (nz_abc - 2 * Nboudary) * sizeof(float)); // saves PV values
            }

            file.close();
            file_PVx.close();
            file_PVz.close();
            file_PVxOF.close();
            file_PVzOF.close();
        }

        //----------------------------------
        // advance in time
        //----------------------------------

        std::swap(u_curr, u_next);
    }

    //-----------------------------------
    // SAVE THE DOCUMENT OF THE SISMOGRAM
    //-----------------------------------

    std::cout << "Saving the seismogram binary file!" << std::endl;

    std::ofstream file("/home/processamento/acustica_2D/outputs/seismogram.bin", std::ios::binary);

    file.write(reinterpret_cast<char *>(seismogram), nrec * nt * sizeof(float));

    file.close();

    std::cout << "Seismogram binary file saved!" << std::endl;

    //-----------------------------------
    // SAVE the copy of the final field
    //-----------------------------------

    float *result = (float *) malloc(nx_abc * nz_abc * sizeof(float));

    for (int i = 0; i < nx_abc * nz_abc; i++)
    {

        result[i] = u_curr[i];
    }

    free(e);
    free(u_curr);
    free(u_next);
    free(seismogram);
    free(PVx);
    free(PVz);
    free(px);
    free(pz);
    free(pt);
    free(ux);
    free(uz);

    return result;
}

//----------------------------------
// MAIN
//----------------------------------

int main()
{

    //----------------------------------
    // open the document of PARAMETERS
    //----------------------------------

    int T;
    int nx;
    int nz;
    int nx_abc;
    int nz_abc;
    int nt;
    int Nboudary;
    int Nsource;
    int nrec;

    float dx;
    float dz;
    float dt;
    float f0;

    char receivers_file[256];
    char sources_file[256];
    char velocity_file[256];

    float *x = NULL;
    float *z = NULL;
    float *t = NULL;

    std::cout << "Reading the document of the parameters!" << std::endl;

    readParameters("/home/processamento/acustica_2D/inputs/parameters.txt", &T, &nx, &nz, &nx_abc, &nz_abc, &nt, &dx, &dz, &dt, &f0, &Nboudary, &Nsource, &nrec, receivers_file, sources_file, velocity_file, &x, &z, &t);

    //----------------------------------
    // open the document of RECEIVERS
    //----------------------------------

    std::cout << "Reading the document of the receivers!" << std::endl;

    Receiver *receivers = readReceivers(receivers_file, nrec, Nboudary);

    //-----------------------------------------
    // open the document of the VELOCITY MODEL
    //-----------------------------------------

    std::cout << "Reading the document of the velocity model!" << std::endl;

    float *c = readVelocity("/home/processamento/acustica_2D/src_cpp/vp_marmousi-ii_shape_(2801, 13601)_dh10m_Nz351_Nx851.bin", nx, nz, nx_abc, nz_abc, Nboudary);

    //------------------------------------------
    // open the document of the SOURCE
    //-----------------------------------------

    int *sx = NULL;
    int *sz = NULL;

    std::cout << "Reading the document of the sources!" << std::endl;

    readSources(sources_file, Nsource, &sx, &sz, Nboudary);

    //----------------------------------
    // check geometry
    //----------------------------------

    std::cout << "Cheking the geometry!" << std::endl;

    if (!checkGeometry(sx, sz, Nsource, receivers, nrec, nx, nz, Nboudary))
    {
        return 1;
    }

    //-----------------------------------------
    // call the source fuction
    //-----------------------------------------

    std::cout << "Starting to source function!" << std::endl;

    float *fonte = source(f0, t, nt);

    //------------------------------------
    //  CERJAN
    //------------------------------------

    std::cout << "Starting to CERJAN function!" << std::endl;

    float *A = createCerjanVector(Nboudary);

    //----------------------------------
    // wavefield
    //----------------------------------

    std::cout << "Starting to wavefield function!" << std::endl;

    float *wavefield = derivates(c, dt, dx, dz, fonte, nx, nz, nx_abc, nz_abc, nt, A, Nboudary, sx, sz, Nsource, receivers, nrec);

    //---------------------------------------
    // Save binary document of the simulation
    //---------------------------------------

    std::ofstream file("/home/processamento/acustica_2D/outputs/wave.bin", std::ios::binary);

    file.write(reinterpret_cast<char *>(wavefield), nx_abc * nz_abc * sizeof(float));

    file.close();

    std::cout << "Wavefield binary file saved!" << std::endl;

    free(wavefield);
    free(A);
    free(fonte);
    free(c);
    free(sx);
    free(sz);
    free(x);
    free(z);
    free(t);
    free(receivers);

    return 0;
}
