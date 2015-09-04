%{
#include <stdio.h>
#include <stdlib.h>
#include <igraph/igraph.h>
#include "util.h"

int yylex(void);
void yyerror(igraph_vector_t *tree, int *size, double *branch_length, const char *str);
int yywrap(void);

int i, error, node = 0, cur = 0;
double tmp;

%}

%union
{
    double number;
    char *string;
}
%token COLON SEMICOLON LPAREN RPAREN COMMA
%token <string> STRING
%token <number> NUMBER

%parse-param {igraph_vector_t *edge} {int *size} {double *branch_length}

%start tree

%%

tree:
    subtree SEMICOLON;

subtree:
    node
    {
        size[node++] = 1;
    }
    |
    LPAREN subtree COMMA subtree RPAREN node
    {
        VECTOR(*edge)[cur++] = node;
        VECTOR(*edge)[cur++] = node - 1;
        VECTOR(*edge)[cur++] = node;
        VECTOR(*edge)[cur++] = node - 1 - size[node-1];
        size[node] = size[node-1] + size[node-1-size[node-1]] + 1;
        ++node;
    }
    ;

node:
    label length;

label:
    |
    NUMBER
    |
    STRING;

length:
    {
        branch_length[node] = 0.0;
    }
    |
    COLON NUMBER
    {
        branch_length[node] = yylval.number;
    }
    ;

%%

void yyerror(igraph_vector_t *edge, int *size, double *branch_length, const char *str)
{
    fprintf(stderr, "invalid Newick format or non-binary tree: %s\n", str);
    exit(EXIT_FAILURE);
}
 
int yywrap(void)
{
    return 1;
} 
