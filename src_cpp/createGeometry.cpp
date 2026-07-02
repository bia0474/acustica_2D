#include <stdlib.h>
#include <stdio.h>

//----------------------------------
// linscpace function
//----------------------------------

int* linspace(int start, int end, int quantity, int endpoint){

    int *number = (int*) malloc(quantity * sizeof(int));

    int dx;

    if(quantity == 1){
        
        number[0] = end - start;

        return number;
    }

    if(endpoint)
        dx = (end - start)/(quantity - 1);
    else
        dx = (end - start) / quantity;

    for(int i = 0; i < quantity; i++){

        number[i] = start + i * dx;
    }

    return number;
}

int main(){

    //----------------------------------
    // open the document of PARAMETERS
    //----------------------------------

    FILE *file = fopen("/home/processamento/acustica_2D/inputs/parameters.txt", "r");

    if(file == NULL){
        printf("Erro ao abrir arquivo de parametros\n");
        return 1;
    }

    char linha[256];

    int Nrec = 0;
    int Nsource = 0;

    while(fgets(linha, sizeof(linha), file)){

        if(sscanf(linha, "nrec = %d", &Nrec) == 1){
            continue;
        }

        if(sscanf(linha, "Nsource = %d", &Nsource) == 1){
            continue;
        }
    }

    fclose(file);

    //-------------------------------
    // Sources
    //-------------------------------

    int sx_init = 0;
    int sx_end = 250;

    int *sx = linspace(sx_init, sx_end, Nsource, 0);

    int *sz = (int*) malloc(Nsource * sizeof(int));

    if(sz == NULL){
        printf("Erro ao alocar memoria\n");
        return 1;
    }

    for(int i = 0; i < Nsource; i++){

        sz[i] = 100.0f;
    }

    int *sIdx = (int*) malloc(Nsource * sizeof(int));

    if(sIdx == NULL){
        printf("Erro ao alocar memoria\n");
        return 1;
    }

    for(int i = 0; i < Nsource; i++){

        sIdx[i] = i;
    }

    //----------------------------------
    // RECEIVERS
    //----------------------------------

    int rx_init = 60;
    int rx_end = 440;

    int *rx = linspace(rx_init, rx_end, Nrec, 1);

    int *rz = (int*) malloc(Nrec * sizeof(int));

    if(rz == NULL){
        printf("Erro ao alocar memoria\n");
        return 1;
    }

    for(int i = 0; i < Nrec; i++){

        rz[i] = 60.0f;
    }

    int *rIdx = (int*) malloc(Nrec * sizeof(int));

    if(rIdx == NULL){
        printf("Erro ao alocar memoria\n");
        return 1;
    }

    for(int i = 0; i < Nrec; i++){

        rIdx[i] = i;
    }

    //----------------------------------
    // SOURCES DATA in csv
    //----------------------------------

    FILE *file_sources = fopen("/home/processamento/acustica_2D/inputs/sources.csv", "w"); //Creates a pointer to a file and opens a file named "sources.csv" in write mode (w)

    if(file_sources == NULL){
        printf("Erro ao abrir sources.csv\n");
        return 1;
    }

    fprintf(file_sources, "index,coordx,coordz\n"); //writing the header

    for(int i = 0; i < Nsource; i++){

        fprintf(file_sources, "%d,%d,%d\n", i, sx[i], sz[i]);
    }

    fclose(file_sources);

    //----------------------------------
    // RECEIVERS DATA in csv
    //----------------------------------

    FILE *file_receivers = fopen("/home/processamento/acustica_2D/inputs/receivers.csv", "w"); //Creates a pointer to a file and opens a file named "receivers.csv" in write mode (w)


    fprintf(file_receivers, "index,coordx,coordz\n"); //writing the header

    for(int i = 0; i < Nrec; i++){

        fprintf(file_receivers, "%d,%d,%d\n", i, rx[i], rz[i]);
    }

    fclose(file_receivers);

    return 0;

}



