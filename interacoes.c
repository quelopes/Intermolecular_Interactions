#include "interacoes.h"

// FUNCOES PRINT
/* ========================================================== 
   ================= PRINT VETOR ESTRUTURADO================= 
   ========================================================== */
void print_vetStruct(PDB* vet, int qtd) {
  int i;
  for(i=0; i<qtd;i++){
     printf("Aminoacido: %c%c%c\nAtomo: %c\nAtomo Tipo: %s\nCoordenadas\n X: %f\n Y: %f\n Z: %f\nIdentificador: %d\nTipo %s\n\n", vet[i].aa[0], vet[i].aa[1], vet[i].aa[2], vet[i].atm, vet[i].atm_tipo, vet[i].x, vet[i].y, vet[i].z, vet[i].id, vet[i].tipo);
  }
}
/* ========================================================== 
   ============ PRINT MATRIZ DE DISTANCIAS ================== 
   ========================================================== */   
void print_matDist(float** matrix, int qt1, int qt2){
    int m,n;
    for(m = 0;m < qt1;m++){
     printf("[%d]\n", m);
     for(n = 0;n < qt2;n++){
       printf("%f ",matrix[m][n]);
     }
     printf("\n");
    }
}
/* ========================================================== 
   =============== PRINT PARAMETROS ========================= 
   ========================================================== */
void print_parametros(IMOL* parametros){
   printf("Parametros\nPt Hidrogenio: %f Var: %f\nhidrofobicos: %f var: %f\nCovalente: %f var: %f\nEletrostatica: %f var: %f\nHIS: %d, LYS: %d, GLU: %d, ASP: %d, ARG: %d\n\n", parametros->ptH, parametros->ptH_var, parametros->hidrofo, parametros->hidrofo_var, parametros->coval, parametros->coval_var, parametros->eletro, parametros->eletro_var, parametros->His, parametros->Lys, parametros->Glu, parametros->Asp, parametros->Arg); 
} 

// FUNCOES DE ALOCACAO
/* ========================================================== 
   ============ ALOCA MATRIZ DE DISTANCIAS ==================
   ========================================================== */   
float** mat_alocMatDist(int lines, int cols){
  int i;
  float **matrix;
  matrix = (float**)malloc(lines*sizeof(float*));
  if (!matrix) {
    printf("Falta memória para alocar a matriz de ponteiros\n");
    exit(1);
  }
  for (i = 0; i < lines; i++) {
    matrix[i] = (float*)malloc(cols*sizeof(float));
    if (!matrix[i]){
      printf("Falta memória para alocar a matriz de ponteiros.\n");
      exit(1);
    }
  }
  return matrix;
}

//FUNCOES BUSCA
/* ========================================================== 
   ===================== MATCH ============================== 
   ========================================================== */
int match(char* line) {
    if( (strncmp(line, "ATOM", 4) != 0) && (strncmp(line, "HETATM", 6) != 0))  
       return 0;
    int j = strlen(line) - 1;
    while(line[j] == ' ' || line[j] == '\n')
      j--;
    switch(line[j]) {
        case 'C':
        case 'O':
        case 'N':        
        case 'S':
// 	case 'H':
            if( line[j-1] == ' ' )
                return 1;
            break;
    }  
    return 0;    
}

/* ========================================================== 
   ================= FUNCAO FILTRA ========================== 
   ========================================================== */
void filtraArqPDB(FILE* pdb, int* vet){
  char line[1024];
  int a = 0; int b = 0; int c = 0;
  int qtd = 0;
  int j;
  char aa[4];
   while(fgets(line, 1024, pdb) !=  0) {
     if(match(line) == 1){
       if(strncmp(line, "A", 1) == 0){
         vet[0]++; 
       }
       else{
         strncpy(aa, line + 17, 3);
	 aa[3] = '\0';
         if((strncmp(aa, "HOH", 3) != 0 || (strncmp(aa, "HOH", 3) != 0 ))){
	    vet[1]++;
	 }else
	   vet[2]++;
       }
    }
   }
}

/* ========================================================== 
   ============= FUNCAO PREENCHE PROTEINA =================== 
   ========================================================== */
void fill_vetProt(FILE* Arqpdb, PDB* vet){
char line[1024];
 int i = 0;
 int n = 0;
 while(fgets(line, 1024, Arqpdb) !=  0) {
      if(match(line) == 1){
         if((strncmp(line, "A", 1) != 0)) 
           n++;
         else{	   
           sscanf(line, "%s %d %s %s %c %d %f %f %f %*f %*f %c", vet[i].tipo, &vet[i].mol, vet[i].atm_tipo, vet[i].aa, &vet[i].cadeia, &vet[i].mol, &vet[i].x, &vet[i].y, &vet[i].z, &vet[i].atm);
	   	   i++;
	   vet[i].id = i;	   
         } 
      }     
  }	  
}

/* ========================================================== 
   ============ FUNCAO PREENCHE SOLVENTE ==================== 
   ========================================================== */
void fill_vetSol(FILE* Arqpdb, PDB* vet){
char line[1024];
 int i = 0;
 int n = 0;
 char aa[4];
 while(fgets(line, 1024, Arqpdb) !=  0) {
      if(match(line) == 1){
         if((strncmp(line, "A", 1) == 0)) 
           n++;
	else{	
	  strncpy(aa, line + 17, 3);
	  aa[3] = '\0';
          if((strncmp(aa, "HOH", 3) != 0 || strncmp(aa, "HOH", 3) != 0 ))
	    n++;
	  else{
	     sscanf(line, "%s %d %s %s %s %d %f %f %f %*f %*f %c", vet[i].tipo, &vet[i].mol, vet[i].atm_tipo, vet[i].aa, &vet[i].cadeia, &vet[i].mol, &vet[i].x, &vet[i].y, &vet[i].z, &vet[i].atm);
	     i++;
	     vet[i].id = i;	   
           }
	} 
      }     
  }	  
}
  
/* ========================================================== 
   ============= FUNCAO PREENCHE LIGANTE ==================== 
   ========================================================== */   
void fill_vetLig(FILE* Arqpdb, PDB* vet){
char line[1024];
 int i = 0;
 int n = 0;
 char aa[4];
 while(fgets(line, 1024, Arqpdb) !=  0) {
      if(match(line) == 1){
         if((strncmp(line, "A", 1) == 0)) 
           n++;
	else{	
	  strncpy(aa, line + 17, 3);
	  aa[3] = '\0';
          if((strncmp(aa, "HOH", 3) == 0 || strncmp(aa, "HOH", 3) == 0 ))
	    n++;
	  else{
	     sscanf(line, "%s %d %s %s %s %d %f %f %f %*f %*f %c", vet[i].tipo, &vet[i].mol, vet[i].atm_tipo, vet[i].aa, &vet[i].cadeia, &vet[i].mol, &vet[i].x, &vet[i].y, &vet[i].z, &vet[i].atm);
	     i++;
	     vet[i].id = i;	   
           }
	} 
      }     
  }	  
}
 
/* ========================================================== 
   ============ FUNCAO CALCULA DISTANCIAS =================== 
   ========================================================== */   
void preencheMat(float** matDist, PDB* vet1, PDB* vet2, int qt1, int qt2){
  int i,j;
  float valor = 0.0; float v = 0.0;
  for(i = 0; i < qt1; i++){
     for(j = 0; j < qt2; j++){
         //printf("Cordenada x em prot: %f\nCoordenada x em lig: %f\n\n", vet1[i].x,vet2[j].x);
        valor = ((float)(pow(vet1[i].x - vet2[j].x,2)) + (float)(pow(vet1[i].y - vet2[j].y,2)) + (float)(pow(vet1[i].z - vet2[j].z,2))); 
        matDist[i][j] = sqrt(valor);
	valor = 0.0;
     }
  }
}

/* ========================================================== 
   ==================FUNCAO PREENCHE========================= 
   ========================================================== */
void preencheParametros(IMOL* parametros, FILE* parFile){
float valor, var;
int i = 1;
int n = 0;
  while (!feof(parFile)){
     fscanf(parFile, "%f %f", &valor, &var);
     switch(i){
       case 1:
         parametros->ptH = valor;
	 parametros->ptH_var = var;
	 break;
       case 2:
         parametros->eletro = valor;
	 parametros->eletro_var = var;
	 break;
       case 3:
         parametros->coval = valor;
	 parametros->coval_var = var;
	 break;
       case 4:
         parametros->hidrofo = valor;
	 parametros->hidrofo_var = var;     
	 break;
       case 5:
         parametros->vdw = valor;
	 parametros->vdw_var = var;     
	 break;
     }
     i++;
  }
}   

/* ========================================================== 
   ===========FUNCAO PREENCHE CARGA ========================= 
   ========================================================== */
void preencheCarga(FILE* arqCarga, IMOL* parametros){
 char aa[4];
 int n = 0;
 while (!feof(arqCarga)){
  fscanf(arqCarga, "%s", aa);
    printf("%s\n", aa);
    if((strncmp(aa, "HIS", 3) == 0))
      parametros->His = 1;
    else if((strncmp(aa, "ARG", 3) != 0))
      parametros->Arg = 1;  
    else if((strncmp(aa, "ASP", 3) != 0))
      parametros->Asp = 1;
    if((strncmp(aa, "GLU", 3) != 0))
      parametros->Glu = 1;  
    else if((strncmp(aa, "LYS", 3) != 0))
      parametros->Lys = 1;  
    else
      n++;
 }
} 

/* ========================================================== 
   ==================FUNCAO PREENCHE========================= 
   ========================================================== */
void procuraLigSolProt(FILE* saida, IMOL* parametros, float** matDist, float** matDist2, int lines, int cols, int linesProt, PDB* vet1, PDB* vet2, PDB* vet3){
 int qtd_pontes = 0;
 float valor = 0.0;
 float valor2 = 0.0;
 int i, j, k;
 char atm1, atm2;
 char str [300];
 int n; 
                           
 if((strncmp(vet1[1].tipo, "HETATM", 6) == 0) && ((strncmp(vet2[1].aa, "HOH", 4) == 0)|| strncmp(vet2[1].aa, "HOH", 3) == 0 )){   
   fprintf(saida,"%s\n", "                    ***  LIGACOES DE HIDROGENIO *** \n");
   sprintf(str,"-----------------------------------------------------------------------------");
    fprintf(saida, "%s\n", str);
   sprintf(str, "[1--2--3]  d[1-2] d[2-3] r[1] r[2]  r[3]  a[1]  a[2] a[3]  n[1]  n[2]   n[3]");
    fprintf(saida, "%s\n", str);
   sprintf(str,"-----------------------------------------------------------------------------");
    fprintf(saida, "%s\n", str);
  for(i = 0; i < lines; i++){
     for(j = 0; j < cols; j++){
         valor = matDist[i][j];
         if(valor < (parametros->ptH - parametros->ptH_var)){ //&& valor < (parametros->ptH + parametros->ptH_var)){
	   atm1 = vet1[i].atm;
	   atm2 = vet2[j].atm;
	   if(atm1 == 'C' || atm2 == 'C')
             n++;
	   else{ 
	     for(k = 0; k < linesProt; k++){
	       valor2 = matDist2[k][j];
	       if(valor2 < (parametros->ptH - parametros->ptH_var)){ //&& valor2 < (parametros->ptH + parametros->ptH_var)){
	         atm1 = vet2[j].atm;
	         atm2 = vet3[k].atm;
	         if(atm1 == 'C' || atm2 == 'C')
                   n++; 
	         else{
	           sprintf(str, "%c---%c---%c   %0.2f   %0.2f  %s   %s  %s_%c %4s  %2s %4s   %d   %d     %d", vet1[i].atm, vet2[j].atm, vet3[k].atm, valor, valor2, vet1[i].aa, vet2[j].aa, vet3[k].aa, vet3[i].cadeia, vet1[i].atm_tipo, vet2[j].atm_tipo, vet3[k].atm_tipo, vet1[i].mol, vet2[j].mol, vet3[k].mol);
	           fprintf(saida, "%s\n", str);
		  qtd_pontes = qtd_pontes + 1; 
	         }  
	      }
	     }
           }
        }
    }
  }  
   fprintf(saida,"Quantidade de Pontes de Hidrogenio: %d\n\n", qtd_pontes); 	 
 } 
}
/* ========================================================== 
   ==================FUNCAO PREENCHE========================= 
   ========================================================== */
void procura(FILE* saida, IMOL* parametros, float** matDist, int lines, int cols, PDB* vet1, PDB* vet2){
/* Relembrando: vetqtd[0]: qtd de proteinas lines
  vetqtd[1]: qtd ligante
  vetqtd[2]: qtd solvente
*/
  int qtd_ptH = 0;
  int qtd_eletro = 0;
  int qtd_cova = 0;
  int qtd_hidrofo = 0;
  int qtd_vdw = 0;
  int i,j;
  float valor = 0.0;
  int n = 0;
  char str [150];
  char atm1, atm2;
 
if((strncmp(vet1[1].tipo, "ATOM", 4) == 0) && ((strncmp(vet2[1].aa, "HOH", 4) != 0) || strncmp(vet2[1].aa, "HOH", 3) != 0)){  
fprintf(saida,"%s\n", "                    ***  LIGAÇOES COVALENTES *** ");  
   sprintf(str,"-----------------------------------------------------------------------------");
    fprintf(saida, "%s\n", str);
   sprintf(str, "[1]-[2]  Distância   Res[1]  Res[2]     Atm[1]   Atm[2]     nPDB[1]   nPDB[2]");
    fprintf(saida, "%s\n", str);
   sprintf(str,"-----------------------------------------------------------------------------");
    fprintf(saida, "%s\n", str);
  for(i = 0; i < lines; i++){
     for(j = 0; j < cols; j++){
         valor = matDist[i][j];
         if(valor < (parametros->coval)){
	   sprintf(str, " %c---%c   %f    %s_%c     %s       %4s      %3s       %5d       %d", vet1[i].atm, vet2[j].atm, valor, vet1[i].aa, vet1[i].cadeia, vet2[j].aa, vet1[i].atm_tipo, vet2[j].atm_tipo, vet1[i].mol, vet2[j].mol);
	   fprintf(saida, "%s\n", str);
           qtd_cova = qtd_cova + 1;
	  }
      }
  } 
  fprintf(saida,"Quantidade de Ligacoes Covalentes: %d\n\n", qtd_cova); 	 
}

if((strncmp(vet1[1].tipo, "ATOM", 4) == 0) && ((strncmp(vet2[1].aa, "HOH", 4) != 0) || strncmp(vet2[1].aa, "HOH", 3) != 0)){  
fprintf(saida,"%s\n", "                    ***  LIGACOES HIDROFOBICAS *** ");  
   sprintf(str,"-----------------------------------------------------------------------------");
    fprintf(saida, "%s\n", str);
   sprintf(str, "[1]-[2]  Distância   Res[1]  Res[2]     Atm[1]   Atm[2]     nPDB[1]   nPDB[2]");
    fprintf(saida, "%s\n", str);
   sprintf(str,"-----------------------------------------------------------------------------");
    fprintf(saida, "%s\n", str);
  for(i = 0; i < lines; i++){
     for(j = 0; j < cols; j++){
         valor = matDist[i][j];
         if(valor < (parametros->hidrofo - parametros->hidrofo_var)){ //&& valor > (parametros->hidrofo + parametros->hidrofo_var)){
	  atm1 = vet1[i].atm;
	  atm2 = vet2[j].atm;
	  if((atm1 == 'C' || atm1 == 'S') && (atm2 == 'C' || atm2 == 'S')){
             sprintf(str, " %c---%c   %f    %s_%c     %s       %4s      %3s       %5d       %d", vet1[i].atm, vet2[j].atm, valor, vet1[i].aa, vet1[i].cadeia, vet2[j].aa, vet1[i].atm_tipo, vet2[j].atm_tipo, vet1[i].mol, vet2[j].mol);
	     fprintf(saida, "%s\n", str);
	     qtd_hidrofo = qtd_hidrofo + 1;
	  } 
	}
      }
  } 
  fprintf(saida,"Quantidade de Ligacoes Hidrofobicas: %d\n\n", qtd_hidrofo); 	 
}

if((strncmp(vet1[1].tipo, "ATOM", 4) == 0) && ((strncmp(vet2[1].aa, "HOH", 4) != 0) || strncmp(vet2[1].aa, "HOH", 3) != 0 )){  
  fprintf(saida,"%s\n", "                    ***  LIGACOES ELETROSTATICAS *** ");  
   sprintf(str,"-----------------------------------------------------------------------------");
    fprintf(saida, "%s\n", str);
   sprintf(str, "[1]-[2]  Distância   Res[1]  Res[2]     Atm[1]   Atm[2]     nPDB[1]   nPDB[2]");
    fprintf(saida, "%s\n", str);
   sprintf(str,"-----------------------------------------------------------------------------");
    fprintf(saida, "%s\n", str);
   for(i = 0; i < lines; i++){
     for(j = 0; j < cols; j++){
         valor = matDist[i][j];
         atm1 = vet1[i].atm;
	 atm2 = vet2[j].atm;
	 if((atm1 == 'C' || atm2 == 'C') || (atm1 == atm2))
	   n++;
	 else{  
          if(valor < (parametros->eletro - parametros->eletro_var)){// && valor > (parametros->eletro + parametros->eletro_var)){
            if((strncmp(vet1[i].aa, "HIS",3) == 0) && (parametros->His == 1) && ((strncmp(vet1[i].atm_tipo, "ND1", 3) == 0) || (strncmp(vet1[i].atm_tipo, "NE2", 3) == 0))){
	       sprintf(str, " %c---%c   %f    %s_%c     %s       %4s      %3s       %5d       %d", vet1[i].atm, vet2[j].atm, valor, vet1[i].aa, vet1[i].cadeia, vet2[j].aa, vet1[i].atm_tipo, vet2[j].atm_tipo, vet1[i].mol, vet2[j].mol);
	       fprintf(saida, "%s\n", str);
	       qtd_eletro = qtd_eletro + 1;
	    }
	    else if((strncmp(vet1[i].aa, "ARG",3) == 0) && (parametros->Arg == 1) && ((strncmp(vet1[i].atm_tipo, "NH1", 3) == 0) || (strncmp(vet1[i].atm_tipo, "NH2", 3) == 0))){
	       sprintf(str, " %c---%c   %f    %s_%c     %s       %4s      %3s       %5d       %d", vet1[i].atm, vet2[j].atm, valor, vet1[i].aa, vet1[i].cadeia, vet2[j].aa, vet1[i].atm_tipo, vet2[j].atm_tipo, vet1[i].mol, vet2[j].mol);
	       fprintf(saida, "%s\n", str);
	       qtd_eletro = qtd_eletro + 1;
	    }
	    else if((strncmp(vet1[i].aa, "LYS",3) == 0) && (parametros->Lys == 1) && ((strncmp(vet1[i].atm_tipo, "NZ", 2) == 0))){
	       sprintf(str, " %c---%c   %f    %s_%c     %s       %4s      %3s       %5d       %d", vet1[i].atm, vet2[j].atm, valor, vet1[i].aa, vet1[i].cadeia, vet2[j].aa, vet1[i].atm_tipo, vet2[j].atm_tipo, vet1[i].mol, vet2[j].mol);
	       fprintf(saida, "%s\n", str);
	       qtd_eletro = qtd_eletro + 1;
	    }
            else if((strncmp(vet1[i].aa, "GLU",3) == 0) && (parametros->Glu == 1) && ((strncmp(vet1[i].atm_tipo, "OE1", 3) == 0) || (strncmp(vet1[i].atm_tipo, "OE2", 3) == 0))){
	       sprintf(str, " %c---%c   %f    %s_%c     %s       %4s      %3s       %5d       %d", vet1[i].atm, vet2[j].atm, valor, vet1[i].aa, vet1[i].cadeia, vet2[j].aa, vet1[i].atm_tipo, vet2[j].atm_tipo, vet1[i].mol, vet2[j].mol);
	       fprintf(saida, "%s\n", str);
	       qtd_eletro = qtd_eletro + 1;
	    }	    
	    else if((strncmp(vet1[i].aa, "ASP",3) == 0) && (parametros->Asp == 1) && ((strncmp(vet1[i].atm_tipo, "OD1", 3) == 0) || (strncmp(vet1[i].atm_tipo, "OD2", 3) == 0))){
	       sprintf(str, " %c---%c   %f    %s_%c     %s       %4s      %3s       %5d       %d", vet1[i].atm, vet2[j].atm, valor, vet1[i].aa, vet1[i].cadeia, vet2[j].aa, vet1[i].atm_tipo, vet2[j].atm_tipo, vet1[i].mol, vet2[j].mol);
	       fprintf(saida, "%s\n", str);
	       qtd_eletro = qtd_eletro + 1;
	    }	
	  }  
	}
     }
   }
   fprintf(saida,"Quantidade de Ligacoes eletrostaticas: %d\n\n", qtd_eletro); 	 
}   

if((strncmp(vet1[1].tipo, "ATOM", 4) == 0) && ((strncmp(vet2[1].aa, "HOH", 4) != 0) || strncmp(vet2[1].aa, "HOH", 3) != 0)){  
 fprintf(saida,"%s\n", "                    ***  LIGACOES DE VAN DER WAALS *** ");
   sprintf(str,"-----------------------------------------------------------------------------");
    fprintf(saida, "%s\n", str);
   sprintf(str, "[1]-[2]  Distância   Res[1]  Res[2]     Atm[1]   Atm[2]     nPDB[1]   nPDB[2]");
    fprintf(saida, "%s\n", str);
   sprintf(str,"-----------------------------------------------------------------------------");
    fprintf(saida, "%s\n", str);
  for(i = 0; i < lines; i++){
     for(j = 0; j < cols; j++){
         valor = matDist[i][j];
         if(valor < (parametros->vdw - parametros->vdw_var)){// && valor > (parametros->vdw + parametros->vdw_var)){
	  atm1 = vet1[i].atm;
	  atm2 = vet2[j].atm;
	  if((atm1 == 'N' || atm1 == 'C' || atm1 == 'S' || atm1 == 'O') && (atm2 == 'C' || atm2 == 'S' || atm2 == 'N' || atm2 == 'O')){
            sprintf(str, " %c---%c   %f    %s_%c     %s       %4s      %3s       %5d       %d", vet1[i].atm, vet2[j].atm, valor, vet1[i].aa, vet1[i].cadeia, vet2[j].aa, vet1[i].atm_tipo, vet2[j].atm_tipo, vet1[i].mol, vet2[j].mol);
	    fprintf(saida, "%s\n", str);
	    qtd_vdw = qtd_vdw + 1;
	  } 
	}
      }
  } 
   fprintf(saida,"Quantidade de Van der Waals: %d\n\n", qtd_vdw); 	 
}

if((strncmp(vet1[1].tipo, "ATOM", 4) == 0) && ((strncmp(vet2[1].aa, "HOH", 4) == 0) || strncmp(vet2[1].aa, "HOH", 3) == 0))
 n++;
else{ 
fprintf(saida,"%s\n", "                    ***  LIGACOES DE HIDROGENIO *** ");
   sprintf(str,"-----------------------------------------------------------------------------");
    fprintf(saida, "%s\n", str);
   sprintf(str, "[1]-[2]  Distância   Res[1]  Res[2]     Atm[1]   Atm[2]     nPDB[1]   nPDB[2]");
    fprintf(saida, "%s\n", str);
   sprintf(str,"-----------------------------------------------------------------------------");
    fprintf(saida, "%s\n", str);
  for(i = 0; i < lines; i++){
     for(j = 0; j < cols; j++){
         valor = matDist[i][j];
         if(valor < (parametros->ptH - parametros->ptH_var)){// && valor > (parametros->ptH + parametros->ptH_var)){
	  atm1 = vet1[i].atm;
	  atm2 = vet2[j].atm;
	  if(atm1 == 'C' || atm2 == 'C')
            n++;
	  else{ 
	    sprintf(str, " %c---%c   %f    %s_%c     %s       %4s      %3s       %5d       %d", vet1[i].atm, vet2[j].atm, valor, vet1[i].aa, vet1[i].cadeia, vet2[j].aa, vet1[i].atm_tipo, vet2[j].atm_tipo, vet1[i].mol, vet2[j].mol);
	    fprintf(saida, "%s\n", str);
	    qtd_ptH = qtd_ptH + 1;
	  }
	}
    }
  } 
   fprintf(saida,"Quantidade de Pontes de Hidrogenio: %d\n\n", qtd_ptH); 	 
} 

  fprintf(saida,"--- RESUMO ---\nLigacoes Eletrostaticas: %d\nLigacoes Hidrofobicas: %d\nLigacoes Covalentes: %d\nLigacoes de Hidrogenio: %d\nLigacoes de Van der Waals: %d\n\n", qtd_eletro, qtd_hidrofo, qtd_cova, qtd_ptH, qtd_vdw); 	 

} 

/* ========================================================== 
   ============ FUNCAO PT DISSULFETO ======================== 
   ========================================================== */
void verificaPSulfeto(FILE* saida, float** matProt, PDB* vetProt, int qtd){
 int i,j, n;
 int qtd_ptS = 0;
 float valor = 0.0;
 char atm1, atm2;
 char str[200];
 fprintf(saida,"%s\n", "                    ***  PONTES DISSULFETO *** ");
         sprintf(str,"-----------------------------------------------------------------------------");
        fprintf(saida, "%s\n", str);
         sprintf(str, "[1]-[1]  Distância   Res[1]  Res[1]     Atm[1]   Atm[1]     nPDB[1]   nPDB[1]");
        fprintf(saida, "%s\n", str);
         sprintf(str,"-----------------------------------------------------------------------------");
	fprintf(saida, "%s\n", str); 
 for(i = 0; i < qtd; i++){
 for(j = 0; j < qtd; j++){
   if(i == j || i > j)
    n++;
   else{
    valor = matProt[i][j];
    if(valor < 5){
      atm1 = vetProt[i].atm;
      atm2 = vetProt[j].atm;
      if((atm1 == 'S') && (atm2 == 'S')){
        sprintf(str, " %c---%c   %f     %s     %s       %4s      %3s       %5d       %d", vetProt[i].atm, vetProt[j].atm, valor, vetProt[i].aa, vetProt[j].aa, vetProt[i].atm_tipo, vetProt[j].atm_tipo, vetProt[i].mol, vetProt[j].mol);
	fprintf(saida, "%s\n", str);
        qtd_ptS = qtd_ptS + 1;
      }
    }
  }
 }
}
fprintf(saida,"Quantidade de Pontes Dissulfeto: %d\n\n", qtd_ptS);


}