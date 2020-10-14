/*==============================================
LNCC Laboratorio Nacional de Computacao Cientifica 
Relatorio de Interacoes Intermoleculares 
Raquel lopes Costa
data: 04/11/10
==============================================*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h> 
#include <time.h>

#include "interacoes.h"

int main(int argc, char **argv) {
  
  int qtd = 0;
  int atom = 0;
  int hetatom = 0;
  char line[1024];
  char arqPDB[20];
  char str[100];
  int vetQtd[3];
  
  /* Um vetor estruturado para cada tipo */
  PDB* vetProt;
  PDB* vetLig;
  PDB* vetSol;
  PDB* vetProt2;
   
  IMOL* parametros;
  
  vetQtd[0] = 0; /* Informacoes da quantidade de proteinas */
  vetQtd[1] = 0; /* Informacoes da quantidade de ligantes */
  vetQtd[2] = 0; /* Informacoes da quantidade de solvente */
  
  printf("Entre com o nome do arquivo PDB\n");
  scanf("%s",arqPDB);
  /* Arbetura de arquivos e tratamento dos erros */
  FILE* pdb = fopen(arqPDB, "r");
    if (!pdb) {
     printf("O Arquivo nao pode ser aberto.\n");
     exit(1);
  }
  
  FILE* pdb_out = fopen("saida.txt", "wb");
    if(!pdb_out){
       printf("Erro ao abrir arquivo de saida\n");
       exit(1);
    }    
  FILE* parFile = fopen("parametros.txt", "r");
    if(!parFile){
       printf("O Arquivo nao pode ser aberto.\n");
       exit(1);
    } 
  FILE* aaCarga = fopen("aaCarga.txt", "r");
   if(!aaCarga){
       printf("O Arquivo nao pode ser aberto.\n");
       exit(1);
   } 
  
  /* Le arquivo de parametros e armazena no struct */
  int t = 1;
  parametros = (IMOL*)malloc(t*sizeof(IMOL));
  preencheParametros(parametros, parFile);
  
  
  preencheCarga(aaCarga,parametros);
//   print_parametros(parametros);  
  
   /* Vê a quantidade de dados que serao usados */
   filtraArqPDB(pdb, vetQtd);
   /* Declaracoes das matrizes de distancia */
   float** matProtLig;
   float** matProtSol;
   float** matLigSol;
   float** matProt;
   
   /* Alocacao de espaco para matrizes*/
   matProtLig = mat_alocMatDist(vetQtd[0],vetQtd[1]);
   matProtSol = mat_alocMatDist(vetQtd[0], vetQtd[2]);
   matLigSol =  mat_alocMatDist(vetQtd[1], vetQtd[2]);
   matProt = mat_alocMatDist(vetQtd[0],vetQtd[0]);
printf("Quantidade de aminacidos:%d\nQuatidade de ligantes:%d\nQuantidade de solventes:%d\n\n", vetQtd[0],vetQtd[1],vetQtd[2]);

    /* Retorna ao inicio do arquivo */    
    fseek(pdb, 0, SEEK_SET); 
    
    vetProt = (PDB *) malloc(vetQtd[0] * sizeof(PDB));
    vetLig  = (PDB *) malloc(vetQtd[1] * sizeof(PDB));
    vetSol  = (PDB *) malloc(vetQtd[2] * sizeof(PDB));
    vetProt2 = (PDB *) malloc(vetQtd[0] * sizeof(PDB));
    
    /* Preenche os dados */
    fill_vetProt(pdb, vetProt);
     fseek(pdb, 0, SEEK_SET); 
    fill_vetSol(pdb, vetSol);
     fseek(pdb, 0, SEEK_SET); 
    fill_vetLig(pdb, vetLig);
     fseek(pdb, 0, SEEK_SET); 
    fill_vetProt(pdb, vetProt2);
    
//     print_vetStruct(vetProt, vetQtd[0]);
  //   print_vetStruct(vetLig, vetQtd[1]);
    // print_vetStruct(vetSol, vetQtd[2]);
     
     /* Preenche matriz de distancias */
     preencheMat(matProtLig, vetProt, vetLig, vetQtd[0], vetQtd[1]);
     preencheMat(matProtSol, vetProt, vetSol, vetQtd[0], vetQtd[2]);
     preencheMat(matLigSol, vetLig, vetSol, vetQtd[1], vetQtd[2]);
     preencheMat(matProt, vetProt, vetProt2, vetQtd[0], vetQtd[0]);

//     print_matDist(matProtLig, vetQtd[0], vetQtd[1]);

     /* tipos  1 ponte de hidrogenio; 2 coulombiana; 3 covalente 4 hidrofobico 5 van der waals */ 

   /* Procura dados e gera Relatorio */
  fprintf(pdb_out,"%s\n", "Laboratorio Nacional de Computacao Cientifica\nAluna: Raquel Lopes Costa\n\n*****************************************************************************\n************ Relatório das Iteracoes Inter Moleculares **********************\n*****************************************************************************\n");                
  sprintf(str, "Arquivo do PDB: ==================> %s <===============================\n\n", arqPDB);
  fprintf(pdb_out,"%s", str);
  
  fprintf(pdb_out,"%s", "\n=============================================================================\n");
  fprintf(pdb_out,"%s", "================ Entre a Proteína [1] e o Ligante [2] =======================\n");
  fprintf(pdb_out,"%s", "=============================================================================\n\n");
     procura(pdb_out, parametros, matProtLig, vetQtd[0],vetQtd[1], vetProt, vetLig);   

  fprintf(pdb_out,"%s", "\n=============================================================================\n");
  fprintf(pdb_out,"%s", "================ Entre o Ligante [1] e o Solvente [2] =======================\n");
  fprintf(pdb_out,"%s", "=============================================================================\n\n");
    procura(pdb_out, parametros, matLigSol, vetQtd[1],vetQtd[2], vetLig, vetSol);     

  fprintf(pdb_out,"%s", "\n=============================================================================\n");
  fprintf(pdb_out,"%s", "============== Entre Ligante [1] Solvente [2] e Proteina [3] ================\n");
  fprintf(pdb_out,"%s", "=============================================================================\n\n");
    procuraLigSolProt(pdb_out, parametros, matLigSol, matProtSol, vetQtd[1],vetQtd[2], vetQtd[0], vetLig, vetSol, vetProt);     
  
  fprintf(pdb_out,"%s", "\n=============================================================================\n");
  fprintf(pdb_out,"%s", "=========== Pontes Dissulfeto Proteína [1] Proteína [1] =====================\n");
  fprintf(pdb_out,"%s", "=============================================================================\n\n");
   verificaPSulfeto(pdb_out, matProt, vetProt, vetQtd[0]);

fclose(pdb);
fclose(pdb_out);
return (EXIT_SUCCESS);
}