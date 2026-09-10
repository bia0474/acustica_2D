#include <stdlib.h>
#include <stdio.h>

//----------------------------------
// linscpace function
//----------------------------------

int *linspace(int start, int end, int quantity, int endpoint)
{

    int *number = (int *)malloc(quantity * sizeof(int));

    int dx;

    if (quantity == 1)
    {

        number[0] = end - start;

        return number;
    }

    if (endpoint)
        dx = (end - start) / (quantity - 1);
    else
        dx = (end - start) / quantity;

    for (int i = 0; i < quantity; i++)
    {

        number[i] = start + i * dx;
    }

    return number;
}

//----------------------------------
// Gera array de pontos de GRID, a partir de posicoes em metros
// start_m, end_m, step_m: em metros
// dx: espacamento do grid em metros/ponto
// count: recebe a quantidade de pontos gerados
//----------------------------------

int *arange_gridpoints(int start_m, int end_m, int step_m, int dx, int *count)
{

    int start_grid = start_m / dx;
    int end_grid = end_m / dx;
    int step_grid = step_m / dx;

    *count = (end_grid - start_grid) / step_grid + 1;

    int *arr = (int *)malloc((*count) * sizeof(int));
    if (arr == NULL)
        return NULL;

    for (int i = 0; i < *count; i++)
    {
        arr[i] = start_grid + i * step_grid;
    }

    return arr;
}

int main()
{

    //-------------------------------
    // Parametros do modelo (mesmos do solver)
    //-------------------------------

    int dx = 10; // m/ponto de grid em x
    int dz = 10; // m/ponto de grid em z

    //-------------------------------
    // Sources - de 5000m a 10000m, passo 50m
    //-------------------------------

    int sx_init_m = 5000;
    int sx_end_m = 10000;
    int step_m = 50;
    int depth_source_m = 100; // profundidade da fonte em metros

    int Nsource;
    int *sx = arange_gridpoints(sx_init_m, sx_end_m, step_m, dx, &Nsource);

    if (sx == NULL)
    {
        printf("Erro ao alocar memoria (sx)\n");
        return 1;
    }

    int *sz = (int *)malloc(Nsource * sizeof(int));
    if (sz == NULL)
    {
        printf("Erro ao alocar memoria (sz)\n");
        free(sx);
        return 1;
    }

    int sz_grid = depth_source_m / dz;
    for (int i = 0; i < Nsource; i++)
    {
        sz[i] = sz_grid;
    }

    //----------------------------------
    // RECEIVERS - de 0m a 5000m, passo 50m
    //----------------------------------

    int rx_init_m = 0;
    int rx_end_m = 5000;
    int depth_rec_m = 100; // profundidade do receptor em metros

    int Nrec;
    int *rx = arange_gridpoints(rx_init_m, rx_end_m, step_m, dx, &Nrec);

    if (rx == NULL)
    {
        printf("Erro ao alocar memoria (rx)\n");
        free(sx);
        free(sz);
        return 1;
    }

    int *rz = (int *)malloc(Nrec * sizeof(int));
    if (rz == NULL)
    {
        printf("Erro ao alocar memoria (rz)\n");
        free(sx);
        free(sz);
        free(rx);
        return 1;
    }

    int rz_grid = depth_rec_m / dz;
    for (int i = 0; i < Nrec; i++)
    {
        rz[i] = rz_grid;
    }

    //----------------------------------
    // Checagem de consistencia (pareamento 1 a 1)
    //----------------------------------

    if (Nsource != Nrec)
    {
        printf("Erro: Nsource (%d) e Nrec (%d) precisam ser iguais para o pareamento 1 a 1.\n", Nsource, Nrec);
        free(sx);
        free(sz);
        free(rx);
        free(rz);
        return 1;
    }

    printf("Nsource = %d, Nrec = %d (pontos de grid, SEM borda de absorcao)\n", Nsource, Nrec);

    //----------------------------------
    // SOURCES DATA in csv (indices de GRID, sem Nboudary)
    //----------------------------------

    FILE *file_sources = fopen("/home/processamento/acustica_2D/inputs/sources.csv", "w");
    if (file_sources == NULL)
    {
        printf("Erro ao abrir sources.csv\n");
        free(sx);
        free(sz);
        free(rx);
        free(rz);
        return 1;
    }

    fprintf(file_sources, "index,coordx,coordz\n");
    for (int i = 0; i < Nsource; i++)
    {
        fprintf(file_sources, "%d,%d,%d\n", i, sx[i], sz[i]);
    }
    fclose(file_sources);

    //----------------------------------
    // RECEIVERS DATA in csv (indices de GRID, sem Nboudary)
    //----------------------------------

    FILE *file_receivers = fopen("/home/processamento/acustica_2D/inputs/receivers.csv", "w");
    if (file_receivers == NULL)
    {
        printf("Erro ao abrir receivers.csv\n");
        free(sx);
        free(sz);
        free(rx);
        free(rz);
        return 1;
    }

    fprintf(file_receivers, "index,coordx,coordz\n");
    for (int i = 0; i < Nrec; i++)
    {
        fprintf(file_receivers, "%d,%d,%d\n", i, rx[i], rz[i]);
    }
    fclose(file_receivers);

    printf("Arquivos sources.csv e receivers.csv gerados com sucesso.\n");

    free(sx);
    free(sz);
    free(rx);
    free(rz);

    return 0;
}