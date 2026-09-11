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
    // Parametros do modelo 
    //-------------------------------
    int nx = 401;
    int nz = 161;

    double dx = 12.5;
    double dz = 12.5;

    //-------------------------------
    // Parametros da geometria fonte/receptor
    //-------------------------------

    int Nsource = 25;
    int Nrec = 25;

    double spacing_m = 50.0; // distancia entre fontes consecutivas e entre receptores consecutivos
    int spacing_grid = (int)(spacing_m / dx + 0.5);

    // Profundidade de fonte e receptor (em metros) - AJUSTE CONFORME NECESSARIO
    double depth_source_m = 25.0;
    double depth_rec_m = 25.0;

    int sz_grid = (int)(depth_source_m / dz + 0.5);
    int rz_grid = (int)(depth_rec_m / dz + 0.5);

    int mid_x = nx / 2; // indice de grid do meio do modelo em x

    if (spacing_grid <= 0)
    {
        printf("Erro: espacamento de %.1f m e menor que dx=%.2f m\n", spacing_m, dx);
        return 1;
    }

    // Meio espacamento: usado para afastar o par 0 do centro, para que o
    // primeiro offset ja nasca em spacing_m (ex: 50 m) em vez de zero.
    int half_spacing_grid = spacing_grid / 2;

    if (spacing_grid % 2 != 0)
    {
        printf("Aviso: spacing_grid=%d e impar, half_spacing_grid=%d nao reproduz exatamente %.1f m de offset inicial\n",
               spacing_grid, half_spacing_grid, spacing_m / 2.0);
    }

    //-------------------------------
    // Aloca arrays de fontes e receptores (indices de GRID, sem Nboudary)
    //-------------------------------
    
    int *sx = (int *)malloc(Nsource * sizeof(int));
    int *sz = (int *)malloc(Nsource * sizeof(int));
    int *rx = (int *)malloc(Nrec * sizeof(int));
    int *rz = (int *)malloc(Nrec * sizeof(int));

    if (sx == NULL || sz == NULL || rx == NULL || rz == NULL)
    {
        printf("Erro ao alocar memoria\n");
        free(sx);
        free(sz);
        free(rx);
        free(rz);
        return 1;
    }

    // Fontes: partindo de meio-espacamento a esquerda do centro, espalhando para a ESQUERDA
    for (int i = 0; i < Nsource; i++)
    {
        sx[i] = mid_x - half_spacing_grid - i * spacing_grid;
        sz[i] = sz_grid;
    }

    // Receptores: partindo de meio-espacamento a direita do centro, espalhando para a DIREITA
    for (int i = 0; i < Nrec; i++)
    {
        rx[i] = mid_x + half_spacing_grid + i * spacing_grid;
        rz[i] = rz_grid;
    }

    // Checagem de limites do modelo
    for (int i = 0; i < Nsource; i++)
    {
        if (sx[i] < 0 || sx[i] >= nx)
        {
            printf("Erro: fonte %d fora do modelo (sx=%d, nx=%d)\n", i, sx[i], nx);
            free(sx);
            free(sz);
            free(rx);
            free(rz);
            return 1;
        }
    }
    for (int i = 0; i < Nrec; i++)
    {
        if (rx[i] < 0 || rx[i] >= nx)
        {
            printf("Erro: receptor %d fora do modelo (rx=%d, nx=%d)\n", i, rx[i], nx);
            free(sx);
            free(sz);
            free(rx);
            free(rz);
            return 1;
        }
    }

    printf("Meio do modelo em x = %d (grid) | dx = %.2f m\n", mid_x, dx);
    printf("Espacamento fonte/receptor = %.1f m (%d pontos de grid)\n", spacing_m, spacing_grid);
    printf("Fontes: x=%d ate x=%d (indo para a esquerda), profundidade z=%d\n", sx[0], sx[Nsource - 1], sz_grid);
    printf("Receptores: x=%d ate x=%d (indo para a direita), profundidade z=%d\n", rx[0], rx[Nrec - 1], rz_grid);

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
        fprintf(file_sources, "%d,%d,%d\n", i, sx[i], sz[i]);
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
        fprintf(file_receivers, "%d,%d,%d\n", i, rx[i], rz[i]);
    fclose(file_receivers);

    //----------------------------------
    // PARES FONTE-RECEPTOR (offset/midpoint) DATA in csv
    //----------------------------------
    FILE *file_pairs = fopen("/home/processamento/acustica_2D/inputs/pares_geometry.csv", "w");
    if (file_pairs == NULL)
    {
        printf("Erro ao abrir pares_geometry.csv\n");
        free(sx);
        free(sz);
        free(rx);
        free(rz);
        return 1;
    }
    fprintf(file_pairs, "par,sx_m,rx_m,offset_m,midpoint_m\n");
    int Npairs = (Nsource < Nrec) ? Nsource : Nrec;
    for (int i = 0; i < Npairs; i++)
    {
        double sx_m = sx[i] * dx;
        double rx_m = rx[i] * dx;
        double offset_m = rx_m - sx_m;
        double midpoint_m = (sx_m + rx_m) / 2.0;
        fprintf(file_pairs, "%d,%.1f,%.1f,%.1f,%.1f\n", i, sx_m, rx_m, offset_m, midpoint_m);
    }
    fclose(file_pairs);

    printf("Arquivos sources.csv, receivers.csv e pares_geometry.csv gerados com sucesso.\n");

    free(sx);
    free(sz);
    free(rx);
    free(rz);

    return 0;
}