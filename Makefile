PDB: interacoes.o interacoes_main.o 
	gcc -lm -o Interacoes interacoes.o interacoes_main.o 
interacoes.o: interacoes.c interacoes.h
	gcc -lm -c interacoes.c
interacoes_main.o: interacoes_main.c interacoes.h
	gcc -lm -c interacoes_main.c
