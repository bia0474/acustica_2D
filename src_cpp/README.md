## Programa de modelagem sísmica 2D

### Executar `teste_PVs_OF.cpp` com OpenACC

Compile com o compilador NVIDIA HPC SDK e execute a partir deste diretório (`src_cpp`), para que os caminhos relativos de `inputs/` e `outputs/` sejam resolvidos corretamente:

```bash
nvc++ -O3 -acc -Minfo=accel teste_PVs_OF.cpp -o teste_PVs_OF
./teste_PVs_OF
```

O cálculo da propagação, do amortecimento, dos vetores de Poynting e do fluxo óptico é executado com OpenACC na GPU NVIDIA.