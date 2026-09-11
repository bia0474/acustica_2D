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

int main()
{

    //-------------------------------
    // Parametros do modelo (mesmos do solver)
    //-------------------------------

    int dx = 10; // m/ponto de grid em x
    int dz = 10; // m/ponto de grid em z

    //-------------------------------
    // Parametros do CMP: midpoint fixo, offset variando
    //-------------------------------

    int M_m = 5000;          // midpoint em metros (centro do modelo)
    int offset_min_m = 200;  // menor afastamento fonte-receptor, em metros
    int offset_max_m = 7600; // maior afastamento fonte-receptor, em metros
    int Ncmp = 20;           // quantidade de pares (fold do CMP)

    int depth_source_m = 100; // profundidade da fonte em metros
    int depth_rec_m = 100;    // profundidade do receptor em metros

    // passo de offset, forcado a multiplo de 2*dx para sx/rx carem em ponto de grid exato
    int step_m = (offset_max_m - offset_min_m) / (Ncmp - 1);
    step_m = (step_m / (2 * dx)) * (2 * dx); // arredonda para baixo, multiplo de 2*dx

    if (step_m <= 0)
    {
        printf("Erro: intervalo de offset muito pequeno para Ncmp=%d pontos com dx=%d\n", Ncmp, dx);
        return 1;
    }

    printf("Midpoint fixo em %d m | offset de %d a %d m | passo de offset = %d m\n", M_m, offset_min_m, offset_min_m + (Ncmp - 1) * step_m, step_m);

    //-------------------------------
    // Aloca arrays de fontes e receptores (indices de GRID, sem Nboudary)
    //-------------------------------

    int *sx = (int *)malloc(Ncmp * sizeof(int));
    int *sz = (int *)malloc(Ncmp * sizeof(int));
    int *rx = (int *)malloc(Ncmp * sizeof(int));
    int *rz = (int *)malloc(Ncmp * sizeof(int));

    if (sx == NULL || sz == NULL || rx == NULL || rz == NULL)
    {
        printf("Erro ao alocar memoria\n");
        free(sx);
        free(sz);
        free(rx);
        free(rz);
        return 1;
    }

    int sz_grid = depth_source_m / dz;
    int rz_grid = depth_rec_m / dz;

    for (int i = 0; i < Ncmp; i++)
    {

        int offset_m = offset_min_m + i * step_m;

        int sx_m = M_m + offset_m / 2; // fonte fica a direita do midpoint
        int rx_m = M_m - offset_m / 2; // receptor fica a esquerda do midpoint

        sx[i] = sx_m / dx; // metros -> ponto de grid
        rx[i] = rx_m / dx;

        sz[i] = sz_grid;
        rz[i] = rz_grid;
    }

    int Nsource = Ncmp;
    int Nrec = Ncmp;

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

    //----------------------------------
    // OFFSET/MIDPOINT DATA in csv (para conferencia e para usar no mute depois)
    //----------------------------------

    FILE *file_cmp = fopen("/home/processamento/acustica_2D/inputs/cmp_geometry.csv", "w");
    if (file_cmp == NULL)
    {
        printf("Erro ao abrir cmp_geometry.csv\n");
        free(sx);
        free(sz);
        free(rx);
        free(rz);
        return 1;
    }

    fprintf(file_cmp, "shot,sx_m,rx_m,offset_m,midpoint_m\n");
    for (int i = 0; i < Ncmp; i++)
    {
        int sx_m = sx[i] * dx;
        int rx_m = rx[i] * dx;
        int offset_m = sx_m - rx_m;
        int midpoint_m = (sx_m + rx_m) / 2;
        fprintf(file_cmp, "%d,%d,%d,%d,%d\n", i, sx_m, rx_m, offset_m, midpoint_m);
    }
    fclose(file_cmp);

    printf("Arquivos sources.csv, receivers.csv e cmp_geometry.csv gerados com sucesso.\n");

    free(sx);
    free(sz);
    free(rx);
    free(rz);

    return 0;
}