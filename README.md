---
output:
  pdf_document: default
  html_document: default
---

# Intermolecular Interactions in C

Trabalho final do curso de Introdução à Estrutura de Proteínas e Modelagem Molecular

Laboratório Nacional de Computação Científica.

Professor: Laurent Dardene

Data: Dezembro de 2010


## Description

Intermolecular interactions program. From the spatial coordinates of the atoms, it is possible to know how many and which molecular interactions are taking place between proteins, ligands and solvents. The coordinate data were obtained from molecules determined experimentally and deposited in the PDB.


## Program structure

O programa foi construído utilizando a linguagem de programação C. Contém os seguintes arquivos:

* interacoes.c
* interacoes.h
* interacoes main.c
* parametros.txt
* aaCarga.txt

Os arquivos (parametros.txt e aaCarga.txt) são arquivos de configuração dos parâmetros de distâncias intermoleculares e dos aminoácidos carregados respectivamente.
Para compliar o programa digite make no terminal. Para rodar digite ./Interacoes, em seguinda entrar com o nome do arquivo PDB. Por exemplo, 1POP.pdb

## How to execute -- Executar o código
    make
    ./Interacoes


