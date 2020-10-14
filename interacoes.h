#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include <limits.h> 
#include <time.h>

#ifndef _PDB_H
#define _PDB_H

typedef struct pdb_dados{
   char tipo[8];
   float x;
   float y;
   float z;
   int num;
   char aa[4];
   char atm;   
   int id;
   char atm_tipo[4];
   int mol;
   char cadeia;
} PDB;

typedef struct imol{
   double ptH;
   double ptH_var;
   double coval;
   double coval_var;
   double hidrofo;
   double hidrofo_var;
   double eletro;
   double eletro_var;
   double vdw;
   double vdw_var;
   int His; /* 0 nao carregado, 1 carregado */
   int Arg;
   int Lys;
   int Glu;
   int Asp;

} IMOL;

/* Funcoes Print */
void print_vetStruct(PDB* vet, int qtd);
void print_matDist(float** matrix, int qt1, int qt2);
void print_parametros(IMOL* parametros);

/* Funcoes de alocacao */
float** mat_alocMatDist(int lines, int cols);

/* Funcoes de preenchimento*/
void fill_vetProt(FILE* Arqpdb, PDB* vet);
void fill_vetSol(FILE* Arqpdb, PDB* vet);
void fill_vetLig(FILE* Arqpdb, PDB* vet);
void preencheMat(float** matDist, PDB* vet1, PDB* vet2, int qt1, int qt2);
void preencheParametros(IMOL* parametros, FILE* parFile);

/* Funcoes busca */
int match(char* line);
void filtraArqPDB(FILE* pdb, int* vet);
void procuraLigSolProt(FILE* saida, IMOL* parametros, float** matDist, float** matDist2, int lines, int cols, int linesProt, PDB* vet1, PDB* vet2, PDB* vet3);   
void procura(FILE* saida, IMOL* parametros, float** matDist, int lines, int cols, PDB* vet1, PDB* vet2);
void verificaPSulfeto(FILE* saida, float** matProt, PDB* vetProt, int qtd);

/* Demais funcoes */
PDB preencheAux(char* line, PDB* vet, int v);


#endif