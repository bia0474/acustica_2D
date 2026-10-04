#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <stdlib.h>

//-------------------------------
// Struct of the receivers
//-------------------------------

typedef struct
{
    int x;
    int z;
} Receiver;


//-------------------------------------------------------
// Horn-Schunck optical flow for any wavefield
//-------------------------------------------------------

void opticalFlow(float *ux, float *uz, const float *px, const float *pz, const float *pt, int nx_abc, int nz_abc, float alpha, int n_iter)
{
    float *ux_temp = (float *)malloc(nx_abc * nz_abc * sizeof(float));
    float *uz_temp = (float *)malloc(nx_abc * nz_abc * sizeof(float));

    if (ux_temp == NULL || uz_temp == NULL)
    {
        free(ux_temp);
        free(uz_temp);
        return;
    }

    std::copy(ux, ux + nx_abc * nz_abc, ux_temp);
    std::copy(uz, uz + nx_abc * nz_abc, uz_temp);

    float *ux_current = ux;
    float *uz_current = uz;
    float *ux_next = ux_temp;
    float *uz_next = uz_temp;

    const int N = nx_abc * nz_abc;
#pragma acc data copyin(px[0:N], pz[0:N], pt[0:N]) copy(ux[0:N], uz[0:N], ux_temp[0:N], uz_temp[0:N])
    for (int iter = 0; iter < n_iter; iter++)
    {
#pragma acc parallel loop collapse(2)
        for (int j = 2; j < nx_abc - 2; j++)
        {
            for (int i = 2; i < nz_abc - 2; i++)
            {
                float sum_ux = 0.0f;
                float sum_uz = 0.0f;

                for (int a = -1; a <= 1; a++)
                {
                    for (int b = -1; b <= 1; b++)
                    {
                        sum_ux += ux_current[(j + a) * nz_abc + (i + b)];
                        sum_uz += uz_current[(j + a) * nz_abc + (i + b)];
                    }
                }

                const float ux_average = (1.0f / 12.0f) * (ux_current[j * nz_abc + i - nz_abc] + ux_current[j * nz_abc + i + nz_abc] + ux_current[j * nz_abc + i - 1] + ux_current[j * nz_abc + i + 1] - ux_current[j * nz_abc + i] + sum_ux);
                const float uz_average = (1.0f / 12.0f) * (uz_current[j * nz_abc + i - nz_abc] + uz_current[j * nz_abc + i + nz_abc] + uz_current[j * nz_abc + i - 1] + uz_current[j * nz_abc + i + 1] - uz_current[j * nz_abc + i] + sum_uz);
                const float denominator = alpha * alpha + px[j * nz_abc + i] * px[j * nz_abc + i] + pz[j * nz_abc + i] * pz[j * nz_abc + i];

                if (denominator > 0.0f)
                {
                    const float correction = (px[j * nz_abc + i] * ux_average + pz[j * nz_abc + i] * uz_average + pt[j * nz_abc + i]) / denominator;
                    ux_next[j * nz_abc + i] = ux_average - px[j * nz_abc + i] * correction;
                    uz_next[j * nz_abc + i] = uz_average - pz[j * nz_abc + i] * correction;
                }
                else
                {
                    ux_next[j * nz_abc + i] = ux_average;
                    uz_next[j * nz_abc + i] = uz_average;
                }
            }
        }

        std::swap(ux_current, ux_next);
        std::swap(uz_current, uz_next);
    }

    if (ux_current != ux)
    {
        std::copy(ux_current, ux_current + nx_abc * nz_abc, ux);
        std::copy(uz_current, uz_current + nx_abc * nz_abc, uz);
    }

    free(ux_temp);
    free(uz_temp);
}

//--------------------------------------------------
// Calculate wavefield derivatives for Optical Flow
//--------------------------------------------------

void calculateWavefieldDerivatives(const float *u_curr, const float *u_next, float *px, float *pz, float *pt, float dt, float dx, float dz, int nx_abc, int nz_abc)
{
    const int N = nx_abc * nz_abc;
#pragma acc data copyin(u_curr[0:N], u_next[0:N]) copy(px[0:N], pz[0:N], pt[0:N])
{
#pragma acc parallel loop collapse(2)
    for (int j = 2; j < nx_abc - 2; j++)
    { // traverses all points of the grid in X

        for (int i = 2; i < nz_abc - 2; i++)
        {
            pt[j * nz_abc + i] = (u_next[j * nz_abc + i] - u_curr[j * nz_abc + i]) / dt;
            px[j * nz_abc + i] = (u_curr[(j - 2) * nz_abc + i] - 8 * u_curr[(j - 1) * nz_abc + i] + 8 * u_curr[(j + 1) * nz_abc + i] - u_curr[(j + 2) * nz_abc + i]) / (12 * dx);
            pz[j * nz_abc + i] = (u_curr[j * nz_abc + (i - 2)] - 8 * u_curr[j * nz_abc + (i - 1)] + 8 * u_curr[j * nz_abc + (i + 1)] - u_curr[j * nz_abc + (i + 2)]) / (12 * dz);
        }
    }
}
}

//-------------------------------------------------------
// function to save all the PVs+OF here (binary document)
//-------------------------------------------------------

void saveSnapshots(std::string filenamex, std::string filenamez, int nx_abc, int nz_abc, int Nboudary, int n, int shot, const float *ux_fwd, const float *uz_fwd)
{
    std::ofstream file_PVxOF_fwd(filenamex, std::ios::binary);
    std::ofstream file_PVzOF_fwd(filenamez, std::ios::binary);

    for (int x = Nboudary; x < nx_abc - Nboudary; x++)
    {

        file_PVxOF_fwd.write(reinterpret_cast<const char *>(&ux_fwd[x * nz_abc + Nboudary]), (nz_abc - 2 * Nboudary) * sizeof(float)); // saves PV values

        file_PVzOF_fwd.write(reinterpret_cast<const char *>(&uz_fwd[x * nz_abc + Nboudary]), (nz_abc - 2 * Nboudary) * sizeof(float)); // saves PV values
    }

    file_PVxOF_fwd.close();
    file_PVzOF_fwd.close();

}

void saveWavefieldAndOpticalFlowSnapshots(const std::string &filename_wavefield, const std::string &filename_ux, const std::string &filename_uz, int nx_abc, int nz_abc, int Nboudary, const float *wavefield, const float *ux, const float *uz)
{
    std::ofstream file_wavefield(filename_wavefield, std::ios::binary);
    std::ofstream file_ux(filename_ux, std::ios::binary);
    std::ofstream file_uz(filename_uz, std::ios::binary);

    for (int x = Nboudary; x < nx_abc - Nboudary; x++)
    {
        const int row_size = (nz_abc - 2 * Nboudary) * sizeof(float);

        file_wavefield.write(reinterpret_cast<const char *>(&wavefield[x * nz_abc + Nboudary]), row_size);
        file_ux.write(reinterpret_cast<const char *>(&ux[x * nz_abc + Nboudary]), row_size);
        file_uz.write(reinterpret_cast<const char *>(&uz[x * nz_abc + Nboudary]), row_size);
    }
}

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
// read velocity model
//---------------------------------------

float *readVelocity(const char *velocity_file, int nx, int nz, int nx_abc, int nz_abc, int Nboudary)
{

    FILE *file = fopen(velocity_file, "rb");

    if (file == NULL)
    {
        printf("Erro ao abrir o arquivo do modelo de velocidade: %s\n", velocity_file);
        exit(1);
    }

    float *c = (float *)malloc(nx * nz * sizeof(float));

    fread(c, sizeof(float), nx * nz, file);

    fclose(file);

    // Reads in the format as saved: nz x nx (row = fixed z, all x positions)
    // float *c_raw = (float *)malloc(nx * nz * sizeof(float));

    // fread(c_raw, sizeof(float), nx * nz, file);

    // fclose(file);

    // Converts it to the format expected by the rest of the code: nx x nz.
    // float *c = (float *)malloc(nx * nz * sizeof(float));

    // for (int i = 0; i < nx; i++)
    //{
    // for (int j = 0; j < nz; j++)
    //{
    // c_raw is in (nz, nx) order: index = j * nx + i
    // c must be in (nx, nz) order: index = i * nz + j
    // c[i * nz + j] = c_raw[j * nx + i];
    //}
    //}

    // free(c_raw);

    float *c_exp = (float *)calloc(nx_abc * nz_abc, sizeof(float));

    //----------------------------------
    // Center
    //----------------------------------

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < nz; j++)
        {

            c_exp[(i + Nboudary) * nz_abc + (j + Nboudary)] = c[i * nz + j];
        }
    }

    //----------------------------------
    // Top edge
    //----------------------------------

    for (int i = 0; i < Nboudary; i++)
    {
        for (int j = Nboudary; j < nz_abc - Nboudary; j++)
        {

            c_exp[i * nz_abc + j] = c_exp[Nboudary * nz_abc + j];
        }
    }

    //----------------------------------
    // Bottom edge
    //----------------------------------

    for (int i = nx_abc - Nboudary; i < nx_abc; i++)
    {
        for (int j = Nboudary; j < nz_abc - Nboudary; j++)
        {

            c_exp[i * nz_abc + j] = c_exp[(nx_abc - Nboudary - 1) * nz_abc + j];
        }
    }

    //----------------------------------
    // Left edge
    //----------------------------------

    for (int i = Nboudary; i < nx_abc - Nboudary; i++)
    {
        for (int j = 0; j < Nboudary; j++)
        {

            c_exp[i * nz_abc + j] = c_exp[i * nz_abc + Nboudary];
        }
    }

    //----------------------------------
    // Right edge
    //----------------------------------

    for (int i = Nboudary; i < nx_abc - Nboudary; i++)
    {
        for (int j = nz_abc - Nboudary; j < nz_abc; j++)
        {

            c_exp[i * nz_abc + j] = c_exp[i * nz_abc + (nz_abc - Nboudary - 1)];
        }
    }

    //----------------------------------
    // Top-left corner
    //----------------------------------

    for (int i = 0; i < Nboudary; i++)
    {
        for (int j = 0; j < Nboudary; j++)
        {

            c_exp[i * nz_abc + j] = c[0];
        }
    }

    //----------------------------------
    // Top right corner
    //----------------------------------

    for (int i = 0; i < Nboudary; i++)
    {
        for (int j = nz_abc - Nboudary; j < nz_abc; j++)
        {

            c_exp[i * nz_abc + j] = c[nz - 1];
        }
    }

    //----------------------------------
    // Bottom-left corner
    //----------------------------------

    for (int i = nx_abc - Nboudary; i < nx_abc; i++)
    {
        for (int j = 0; j < Nboudary; j++)
        {

            c_exp[i * nz_abc + j] = c[(nx - 1) * nz];
        }
    }

    //----------------------------------
    // Bottom right corner
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

//----------------------------------
// CFL condition
//----------------------------------

bool CFL(const float* c, float dt, float dx, int nx, int nz){ //function of the stability codition

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

#pragma acc parallel loop copyout(A[0:Nboudary])
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

#pragma acc parallel loop copyin(t[0:nt]) copyout(s[0:nt])
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

    float *u_curr = (float *)calloc(nx_abc * nz_abc, sizeof(float)); // present field of the fwd
    float *u_next = (float *)calloc(nx_abc * nz_abc, sizeof(float)); // future field of the fwd

    float *u_back_curr = (float *)calloc(nx_abc * nz_abc, sizeof(float)); // present field of the back
    float *u_back_next = (float *)calloc(nx_abc * nz_abc, sizeof(float)); // future field of the back

    //----------------------------------
    // Poynting vector + Optical Flow
    //----------------------------------

    float *px_fwd = (float *)calloc(nx_abc * nz_abc, sizeof(float));
    float *pz_fwd = (float *)calloc(nx_abc * nz_abc, sizeof(float));
    float *pt_fwd = (float *)calloc(nx_abc * nz_abc, sizeof(float));

    float *ux_fwd = (float *)calloc(nx_abc * nz_abc, sizeof(float));
    float *uz_fwd = (float *)calloc(nx_abc * nz_abc, sizeof(float));

    float *px_back = (float *)calloc(nx_abc * nz_abc, sizeof(float));
    float *pz_back = (float *)calloc(nx_abc * nz_abc, sizeof(float));
    float *pt_back = (float *)calloc(nx_abc * nz_abc, sizeof(float));

    float *ux_back = (float *)calloc(nx_abc * nz_abc, sizeof(float));
    float *uz_back = (float *)calloc(nx_abc * nz_abc, sizeof(float));

    //----------------------------------
    // SEISMOGRAM
    //----------------------------------
    // Nshots pares, cada um com nt amostras

    // Nsource == Nrec (mesmo número de pares fonte-receptor)
    if (Nsource != nrec)
    {
        std::cerr << "ERRO: a geometria CMP exige Nsource == nrec (" << Nsource << " != " << nrec << ")\n";
        return nullptr;
    }

    int Nshots = Nsource;
    const int N = nx_abc * nz_abc;
    const int image_size = nx * nz;

    float *seismogram_shot = (float*) malloc(nt * sizeof(float));

    //----------------------------------
    // IMAGE CONDITION and ADCIGs
    //----------------------------------

    const float angle_step = 5.0f;
    const int max_angle = 90;
    const int n_angles = max_angle / angle_step;
    const int n_gathers = 1;

    float *image = (float *)calloc(nx * nz, sizeof(float));
    float *image_ADCIGs = (float *)calloc(n_angles * n_gathers * nz, sizeof(float));
    int *adcig_fold = (int *)calloc(n_angles * n_gathers * nz, sizeof(int));
    float *u_fwd_n = (float *)calloc(nx * nz, sizeof(float));
    float *ux_fwd_n = (float *)calloc(nx * nz, sizeof(float));
    float *uz_fwd_n = (float *)calloc(nx * nz, sizeof(float));

    float *theta = (float *)calloc(nx * nz, sizeof(float));

    //-----------------------------------
    // FORWARD FIELD
    //-----------------------------------

    std::cout << "Starting the temporal and spacial loops of the forward!" << std::endl;

#pragma acc data copyin(c[0:N], A[0:Nboudary], fonte[0:nt]) copy(u_curr[0:N], u_next[0:N], u_back_curr[0:N], u_back_next[0:N], px_fwd[0:N], pz_fwd[0:N], pt_fwd[0:N], ux_fwd[0:N], uz_fwd[0:N], px_back[0:N], pz_back[0:N], pt_back[0:N], ux_back[0:N], uz_back[0:N], image[0:image_size], u_fwd_n[0:image_size])

    for (int shot = 0; shot < Nshots; shot++)
    {
        std::cout << "Tiro " << shot << "\n";

        std::fill(u_curr, u_curr + nx_abc * nz_abc, 0.0f);
        std::fill(u_next, u_next + nx_abc * nz_abc, 0.0f);
        std::fill(seismogram_shot, seismogram_shot + nt, 0.0f);
        std::fill(px_fwd, px_fwd + nx_abc * nz_abc, 0.0f);
        std::fill(pz_fwd, pz_fwd + nx_abc * nz_abc, 0.0f);
        std::fill(pt_fwd, pt_fwd + nx_abc * nz_abc, 0.0f);
        std::fill(ux_fwd, ux_fwd + nx_abc * nz_abc, 0.0f);
        std::fill(uz_fwd, uz_fwd + nx_abc * nz_abc, 0.0f);  
        std::fill(u_back_curr, u_back_curr + nx_abc * nz_abc, 0.0f);
        std::fill(u_back_next, u_back_next + nx_abc * nz_abc, 0.0f);
        std::fill(px_back, px_back + nx_abc * nz_abc, 0.0f);
        std::fill(pz_back, pz_back + nx_abc * nz_abc, 0.0f);
        std::fill(pt_back, pt_back + nx_abc * nz_abc, 0.0f);
        std::fill(ux_back, ux_back + nx_abc * nz_abc, 0.0f);
        std::fill(uz_back, uz_back + nx_abc * nz_abc, 0.0f);

        //-----------------------------------
        // trace with mute 
        //-----------------------------------

        std::string filename = "../outputs/seismogram_shot" + std::to_string(shot) + "_mute.bin";

        std::ifstream file_shot(filename, std::ios::binary);

        if(!file_shot.is_open()){
            std::cout << "Erro ao abrir " << filename << "\n";
            exit(1);
        }

        file_shot.read(reinterpret_cast<char*>(seismogram_shot), nt * sizeof(float));

        if(!file_shot){
            std::cout << "ERRO: leitura incompleta em " << filename << ", leu " << file_shot.gcount() << " bytes (esperado " << nt * sizeof(float) << ")\n";
            exit(1);
        }

        file_shot.close();


        for (int n = 1; n < nt; n++)
        { // each iteration calculates the wave at the next instant

            //----------------------------------
            // space loop - 4nd order
            //----------------------------------

#pragma acc parallel loop collapse(2)
            for (int j = 2; j < nx_abc - 2; j++)
            { // traverses all points of the grid in X

                for (int i = 2; i < nz_abc - 2; i++)
                { // traverses all points of the grid in Z

                    //----------------------------------
                    // finite differences
                    //----------------------------------

                    float d2x = (-u_curr[(j + 2) * nz_abc + i] + 16 * u_curr[(j + 1) * nz_abc + i] - 30 * u_curr[j * nz_abc + i] + 16 * u_curr[(j - 1) * nz_abc + i] - u_curr[(j - 2) * nz_abc + i]) / (12 * dx * dx);

                    float d2z = (-u_curr[j * nz_abc + (i + 2)] + 16 * u_curr[j * nz_abc + (i + 1)] - 30 * u_curr[j * nz_abc + i] + 16 * u_curr[j * nz_abc + (i - 1)] - u_curr[j * nz_abc + (i - 2)]) / (12 * dz * dz);

                    u_next[j * nz_abc + i] = 2 * u_curr[j * nz_abc + i] - u_next[j * nz_abc + i] + c[j * nz_abc + i] * c[j * nz_abc + i] * dt * dt * (d2x + d2z);
                }
            }

            //----------------------------------
            // source injection (just one source)
            //----------------------------------

            u_next[sx[shot] * nz_abc + sz[shot]] += (fonte[n]) / (dx * dz);

            //-----------------------------------
            // CERJAN
            //-----------------------------------

            if (n == 1)
            {
                std::cout << "Making the CERJAN boudary of the forward" << std::endl;
            }

#pragma acc parallel loop collapse(2)
            for (int x = 0; x < Nboudary; x++)
            { // Left
                for (int z = 0; z < nz_abc; z++)
                {
                    u_next[x * nz_abc + z] *= A[x];
                    u_curr[x * nz_abc + z] *= A[x];
                }
            }

#pragma acc parallel loop collapse(2)
            for (int x = nx_abc - Nboudary; x < nx_abc; x++)
            { // right
                for (int z = 0; z < nz_abc; z++)
                {
                    const int k = nx_abc - 1 - x;

                    u_next[x * nz_abc + z] *= A[k];
                    u_curr[x * nz_abc + z] *= A[k];
                }
            }

#pragma acc parallel loop collapse(2)
            for (int z = 0; z < Nboudary; z++)
            { // Top
                for (int x = 0; x < nx_abc; x++)
                {
                    u_next[x * nz_abc + z] *= A[z];
                    u_curr[x * nz_abc + z] *= A[z];
                }
            }

#pragma acc parallel loop collapse(2)
            for (int z = nz_abc - Nboudary; z < nz_abc; z++)
            {
                for (int x = 0; x < nx_abc; x++)
                { // Base
                    const int k = nz_abc - 1 - z;

                    u_next[x * nz_abc + z] *= A[k];
                    u_curr[x * nz_abc + z] *= A[k];
                }
            }

            //--------------------------------------------------
            // Saving the derivatives to Optical Flow - FWD
            //--------------------------------------------------

            if(n == 1){
                std::cout << "Saving the derivatives to Optical Flow of the forward!" << std::endl;
            }

            if (n % 10 == 0)
            {
                calculateWavefieldDerivatives(u_curr, u_next, px_fwd, pz_fwd, pt_fwd, dt, dx, dz, nx_abc, nz_abc);
            }

            //---------------------------------------------------------------
            // Using Optical Flow method (20 iterations over the entire mesh)
            //---------------------------------------------------------------

            int n_iter = 20;
            float alpha = 1.0f;

            if(n == 1)
            {
                std::cout << "Calculating the Optical Flow of the forward!" << std::endl;
            }

            if (n % 10 == 0)
            {
                std::fill(ux_fwd, ux_fwd + nx_abc * nz_abc, 0.0f);
                std::fill(uz_fwd, uz_fwd + nx_abc * nz_abc, 0.0f);

                opticalFlow(ux_fwd, uz_fwd, px_fwd, pz_fwd, pt_fwd, nx_abc, nz_abc, alpha, n_iter);
            }

            //-------------------------------------------------------
            // SAVE ALL THE PVxz + OF (binary document)
            //-------------------------------------------------------

            if (n == 1)
            {
                std::cout << "Saving the file of the PVs + OF of the forward!" << std::endl;
            }

            if (n % 10 == 0)
            {
                std::string filenamex = "../outputs/PV+OF_fwd_x" + std::to_string(n) + "_shot " + std::to_string(shot) + ".bin";
                std::string filenamez = "../outputs/PV+OF_fwd_z" + std::to_string(n) + "_shot " + std::to_string(shot) + ".bin";

                saveSnapshots(filenamex, filenamez, nx_abc, nz_abc, Nboudary, n, shot, ux_fwd, uz_fwd);

            }

            //----------------------------------
            // advance in time
            //----------------------------------

            std::swap(u_curr, u_next);
        }

        //-----------------------------------
        // BACKWARD FIELD
        //-----------------------------------

        for (int n = 1; n < nt; n++)
        { // Each iteration calculates the backward field at the next instant (which, physically, is an earlier time).

            //----------------------------------
            // space loop - 4th order
            //----------------------------------

#pragma acc parallel loop collapse(2)
            for (int j = 2; j < nx_abc - 2; j++)
            {
                for (int i = 2; i < nz_abc - 2; i++)
                {
                    float d2x = (-u_back_curr[(j + 2) * nz_abc + i] + 16 * u_back_curr[(j + 1) * nz_abc + i] - 30 * u_back_curr[j * nz_abc + i] + 16 * u_back_curr[(j - 1) * nz_abc + i] - u_back_curr[(j - 2) * nz_abc + i]) / (12 * dx * dx);

                    float d2z = (-u_back_curr[j * nz_abc + (i + 2)] + 16 * u_back_curr[j * nz_abc + (i + 1)] - 30 * u_back_curr[j * nz_abc + i] + 16 * u_back_curr[j * nz_abc + (i - 1)] - u_back_curr[j * nz_abc + (i - 2)]) / (12 * dz * dz);

                    u_back_next[j * nz_abc + i] = 2 * u_back_curr[j * nz_abc + i] - u_back_next[j * nz_abc + i] + c[j * nz_abc + i] * c[j * nz_abc + i] * dt * dt * (d2x + d2z);

                }
            }

            //----------------------------------
            // injection energy to the grid
            //----------------------------------
            // The values ​​recorded on the seismogram are read at the receiver positions in reverse order and injected.

            int xr = receivers[shot].x;
            int zr = receivers[shot].z;

            int time_index = nt - 1 - n;

            u_back_next[xr * nz_abc + zr] += seismogram_shot[time_index] / (dx * dz);

            //-----------------------------------
            // CERJAN
            //-----------------------------------

            if (n == 1)
            {
                std::cout << "Making the CERJAN boudary of the backward" << std::endl;
            }

#pragma acc parallel loop collapse(2)
            for (int x = 0; x < Nboudary; x++)
            { // Left
                for (int z = 0; z < nz_abc; z++)
                {
                    u_back_next[x * nz_abc + z] *= A[x];
                    u_back_curr[x * nz_abc + z] *= A[x];
                }
            }

#pragma acc parallel loop collapse(2)
            for (int x = nx_abc - Nboudary; x < nx_abc; x++)
            { // right
                for (int z = 0; z < nz_abc; z++)
                {
                    const int k = nx_abc - 1 - x;

                    u_back_next[x * nz_abc + z] *= A[k];
                    u_back_curr[x * nz_abc + z] *= A[k];
                }
            }

#pragma acc parallel loop collapse(2)
            for (int z = 0; z < Nboudary; z++)
            { // Top
                for (int x = 0; x < nx_abc; x++)
                {
                    u_back_next[x * nz_abc + z] *= A[z];
                    u_back_curr[x * nz_abc + z] *= A[z];
                }
            }

#pragma acc parallel loop collapse(2)
            for (int z = nz_abc - Nboudary; z < nz_abc; z++)
            {
                for (int x = 0; x < nx_abc; x++)
                { // Base

                    const int k = nz_abc - 1 - z;

                    u_back_next[x * nz_abc + z] *= A[k];
                    u_back_curr[x * nz_abc + z] *= A[k];
                }
            }

            //--------------------------------------------------
            // Saving the derivatives to Optical Flow - BACK
            //--------------------------------------------------

            if(n == 1)
            {
                std::cout << "Saving the derivatives to Optical Flow of the backward!" << std::endl;
            }
            
            if (n % 10 == 0)
            {
                calculateWavefieldDerivatives(u_back_curr, u_back_next, px_back, pz_back, pt_back, dt, dx, dz, nx_abc, nz_abc);
            }

            //---------------------------------------------------------------
            // Using Optical Flow method (20 iterations over the entire mesh)
            //---------------------------------------------------------------

            if (n % 10 == 0)
            {

                std::fill(ux_back, ux_back + nx_abc * nz_abc, 0.0f);
                std::fill(uz_back, uz_back + nx_abc * nz_abc, 0.0f);

                int n_iter = 20;
                float alpha = 1.0f;

                opticalFlow(ux_back, uz_back, px_back, pz_back, pt_back, nx_abc, nz_abc, alpha, n_iter);
            }

            //-------------------------------------------------------
            // SAVE ALL THE SNAPSHOT HERE WITH PVs+OF (binary document)
            //-------------------------------------------------------

            if (n == 1)
            {
                std::cout << "Saving the file of the snapshots and PVs + OF of the backward!" << std::endl;
            }

            if (n % 10 == 0)
            {
                std::string filename_snapshot = "../outputs/snapshot_back_" + std::to_string(n) + "_shot " + std::to_string(shot) + ".bin";
                std::string filename_pv_of_x = "../outputs/PV+OF_back_x" + std::to_string(n) + "_shot " + std::to_string(shot) + ".bin";
                std::string filename_pv_of_z = "../outputs/PV+OF_back_z" + std::to_string(n) + "_shot " + std::to_string(shot) + ".bin";

                saveWavefieldAndOpticalFlowSnapshots(filename_snapshot, filename_pv_of_x, filename_pv_of_z, nx_abc, nz_abc, Nboudary, u_back_next, ux_back, uz_back);
            }

            //----------------------------------
            // imaging condition
            //----------------------------------

            int fwd_index = nt - 1 - n;

            if (fwd_index != 0 && fwd_index % 10 == 0) // At t=0, the forward field is zero by definition.
            {
                std::ifstream fwd_file("../outputs/snapshot_fwd_" + std::to_string(fwd_index) + "_shot_" + std::to_string(shot) + ".bin", std::ios::binary);

                if (!fwd_file.is_open()) // check when opening file
                {
                    std::cerr << "ERRO: nao abriu snapshot_fwd_" << fwd_index << ".bin" << std::endl;
                }

                fwd_file.read(reinterpret_cast<char *>(u_fwd_n), nx * nz * sizeof(float));

                if (!fwd_file) // verification during file reading
                {
                    std::cerr << "ERRO: leitura incompleta em snapshot_fwd_" << fwd_index << ".bin, leu " << fwd_file.gcount() << " bytes" << std::endl;
                }

                fwd_file.close();

#pragma acc parallel loop collapse(2)
                for (int x = 0; x < nx; x++)
                {
                    for (int z = 0; z < nz; z++)
                    {
                        image[x * nz + z] += u_fwd_n[x * nz + z] * u_back_next[(x + Nboudary) * nz_abc + (z + Nboudary)];
                    }
                }
            }

            //----------------------------
            // imaging condition + ADCIGs
            //----------------------------

            if (fwd_index != 0 && fwd_index % 10 == 0) // At t=0, the forward field is zero by definition.
            {

                // --------------------------------------------------------
                // PV+OF_fwd_x and PV+OF_fwd_z
                // --------------------------------------------------------

                std::ifstream fwd_ux_file("../outputs/PV+OF_fwd_x" + std::to_string(fwd_index) + "_shot " + std::to_string(shot) + ".bin", std::ios::binary);

                if (!fwd_ux_file.is_open())
                {
                    std::cerr << "ERRO: nao abriu PV+OF_fwd_x" << fwd_index << ".bin" << std::endl;
                }

                fwd_ux_file.read(reinterpret_cast<char *>(ux_fwd_n), nx * nz * sizeof(float));

                if (!fwd_ux_file)
                {
                    std::cerr << "ERRO: leitura incompleta em PV+OF_fwd_x" << fwd_index << ".bin, leu " << fwd_ux_file.gcount() << " bytes" << std::endl;
                }

                fwd_ux_file.close();

                std::ifstream fwd_uz_file("../outputs/PV+OF_fwd_z" + std::to_string(fwd_index) + "_shot " + std::to_string(shot) + ".bin", std::ios::binary);

                if (!fwd_uz_file.is_open())
                {
                    std::cerr << "ERRO: nao abriu PV+OF_fwd_z" << fwd_index << ".bin" << std::endl;
                }

                fwd_uz_file.read(reinterpret_cast<char *>(uz_fwd_n), nx * nz * sizeof(float));

                if (!fwd_uz_file)
                {
                    std::cerr << "ERRO: leitura incompleta em PV+OF_fwd_z" << fwd_index << ".bin, leu " << fwd_uz_file.gcount() << " bytes" << std::endl;
                }

                fwd_uz_file.close();

                //calcula as amplitudes usadas para decidir se um ponto da malha entra ou não no cálculo do gather
                float max_amp_fwd = 0.0f, max_amp_back = 0.0f;

                //encontrar a maior amplitude do campo forward e do campo backward no instante em questão
                for (int j = 0; j < nx; j++)
                {
                    for (int i = 0; i < nz; i++)
                    {
                        float m_f = sqrt(ux_fwd_n[j * nz + i] * ux_fwd_n[j * nz + i] + uz_fwd_n[j * nz + i] * uz_fwd_n[j * nz + i]);
                        if (m_f > max_amp_fwd) max_amp_fwd = m_f;

                        float m_b = sqrt(ux_back[(j + Nboudary) * nz_abc + (i + Nboudary)] * ux_back[(j + Nboudary) * nz_abc + (i + Nboudary)] + uz_back[(j + Nboudary) * nz_abc + (i + Nboudary)] * uz_back[(j + Nboudary) * nz_abc + (i + Nboudary)]);
                        if (m_b > max_amp_back) max_amp_back = m_b;

                    }
                }
                
                //define que qualquer ponto com amplitude menor que 0,001% é considerado "não confiável" e é descartado no cálculo do gather
                float thresh_fwd  = 1e-5f * max_amp_fwd;
                float thresh_back = 1e-5f * max_amp_back;

                for (int j = 0; j < nx; j++)
                {
                    for (int i = 0; i < nz; i++)
                    {
                        // -------------------
                        // opening angle
                        // -------------------

                        float modulo_fwd  = sqrt(ux_fwd_n[j * nz + i] * ux_fwd_n[j * nz + i] + uz_fwd_n[j * nz + i] * uz_fwd_n[j * nz + i]);

                        float modulo_back = sqrt(ux_back[(j + Nboudary) * nz_abc + (i + Nboudary)] * ux_back[(j + Nboudary) * nz_abc + (i + Nboudary)] + uz_back[(j + Nboudary) * nz_abc + (i + Nboudary)] * uz_back[(j + Nboudary) * nz_abc + (i + Nboudary)]);

                        int bin_angle = -1; //guarda o índice do bin de ângulo onde -1 signiifca "não classificado"
                        
                        if (modulo_fwd > thresh_fwd && modulo_back > thresh_back) //só deixa passar pontos com amplitude genuína o suficiente para confiar na direção estimada
                        {
                            //definição de cosseno do ângulo entre dois vetores:
                            float cos_2theta = (ux_fwd_n[j * nz + i] * ux_back[(j + Nboudary) * nz_abc + (i + Nboudary)] + uz_fwd_n[j * nz + i] * uz_back[(j + Nboudary) * nz_abc + (i + Nboudary)]) / (modulo_fwd * modulo_back);

                            //garante que o valor fique dentro do domínio válido da função aeco-cosseno [-1,1]
                            if (cos_2theta >  1.0f) cos_2theta =  1.0f;
                            if (cos_2theta < -1.0f) cos_2theta = -1.0f;

                            theta[j * nz + i] = 0.5f * acos(cos_2theta) * 180.0f / M_PI; //angulo para graus

                            //Se o ângulo estiver no intervalo físico esperado [0°, 90°), calcula em qual bin discreto ele cai, dividindo pelo tamanho do passo de ângulo (angle_step) 
                            // e truncando para inteiro. Por exemplo, se angle_step = 5° e theta = 23°, então gather = 4 (bin de 20°–25°)
                            
                            if (theta[j * nz + i] >= 0.0f && theta[j * nz + i] < max_angle)
                            {
                                bin_angle = (int)(theta[j * nz + i] / angle_step);
                            }
                        }

                        if (bin_angle >= 0)
                        {
                            image_ADCIGs[bin_angle * n_gathers * nz + i] += u_fwd_n[j * nz + i] * u_back_next[(j + Nboudary) * nz_abc + (i + Nboudary)];

                            adcig_fold[bin_angle * n_gathers * nz + i]++; //incrementa o contador de contribuições para cada bin de ângulo e profundidade
                        }
                    }
                }
            }

            //----------------------------------
            // advance in time
            //----------------------------------

            std::swap(u_back_curr, u_back_next);
        }
    }

    //--------------------------------------
    // SAVE the copy of the final field fwd
    //--------------------------------------

    float *result = (float *)malloc(nx_abc * nz_abc * sizeof(float));

    for (int i = 0; i < nx_abc * nz_abc; i++)
    {

        result[i] = u_curr[i];
    }

    //-----------------------------------------
    // SAVE THE DOCUMENT OF THE MIGRATED IMAGE
    //-----------------------------------------

    std::cout << "Saving the migrated image!" << std::endl;

    for (int i = 0; i < nx * nz; i++)
    {
        image[i] /= static_cast<float>(Nshots);
    }

    std::ofstream img_file("../outputs/image.bin", std::ios::binary);

    img_file.write(reinterpret_cast<char *>(image), nx * nz * sizeof(float));

    img_file.close();

    std::cout << "Migrated image binary file saved!" << std::endl;

    //------------------------------------------
    // SAVE THE DOCUMENT OF THE IMAGE_ADICIG
    //------------------------------------------

    //one gather of the CMP
    std::string filename = "../outputs/image_migrated_ADCIG.bin";

    std::ofstream out_file(filename, std::ios::binary);

    if (!out_file.is_open())
    {
        std::cerr << "ERRO: nao conseguiu salvar " << filename << std::endl;
        return nullptr;
    }

    for (int i = 0; i < n_angles * n_gathers * nz; i++)
    {
        if (adcig_fold[i] > 0)
        {
            image_ADCIGs[i] /= static_cast<float>(adcig_fold[i]); //normaliza o gather pelo número de contribuições
        }
    }

    out_file.write(reinterpret_cast<char *>(image_ADCIGs), n_angles * n_gathers * nz * sizeof(float));

    out_file.close();

    std::cout << "Salvo: " << filename << std::endl;

    free(u_curr);
    free(u_next);
    free(seismogram_shot);
    free(px_back);
    free(pz_back);
    free(pt_back);
    free(ux_back);
    free(uz_back);
    free(u_back_curr);
    free(u_back_next);
    free(image);
    free(image_ADCIGs);
    free(u_fwd_n);
    free(ux_fwd_n);
    free(uz_fwd_n);
    free(theta);
    free(adcig_fold);

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

    readParameters("../inputs/parameters.txt", &T, &nx, &nz, &nx_abc, &nz_abc, &nt, &dx, &dz, &dt, &f0, &Nboudary, &Nsource, &nrec, receivers_file, sources_file, velocity_file, &x, &z, &t);

    //----------------------------------
    // open the document of RECEIVERS
    //----------------------------------

    std::cout << "Reading the document of the receivers!" << std::endl;

    Receiver *receivers = readReceivers(receivers_file, nrec, Nboudary);

    //-----------------------------------------
    // open the document of the VELOCITY MODEL
    //-----------------------------------------

    std::cout << "Reading the document of the velocity model!" << std::endl;

    float *c = readVelocity("../inputs/velocityModel2.bin", nx, nz, nx_abc, nz_abc, Nboudary);

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

    //----------------------------------
    // CFL check
    //----------------------------------

    if(CFL(c, dt, dx, nx_abc, nz_abc)){

        std::cout << "Stable simulation" << std::endl;
    }
    else{

        std::cout << "unstable simulation" << std::endl;
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

    std::ofstream file("../outputs/wave.bin", std::ios::binary);

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