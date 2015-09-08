%{
#include <stdio.h>
#include <stdlib.h>
#include <igraph/igraph.h>
#include "util.h"

int yylex(void);
void yyerror(igraph_vector_t *tree, igraph_vector_t *size, igraph_vector_t *branch_length, const char *str);
int yywrap(void);

int yynode = 0;

%}

%union
{
    double number;
    char *string;
}
%token COLON SEMICOLON LPAREN RPAREN COMMA
%token <string> STRING
%token <number> NUMBER

%parse-param {igraph_vector_t *edge} {igraph_vector_t *size} {igraph_vector_t *branch_length}

%start tree

%%

tree:
    subtree SEMICOLON
    {
        YYACCEPT;
    }
    ;

subtree:
    node
    {
        igraph_vector_push_back(size, 1.0);
        ++yynode;
    }
    |
    LPAREN subtree COMMA subtree RPAREN node
    {
        igraph_vector_push_back(edge, yynode);
        igraph_vector_push_back(edge, yynode-1);
        igraph_vector_push_back(edge, yynode);
        igraph_vector_push_back(edge, yynode-1-VECTOR(*size)[yynode-1]);
        igraph_vector_push_back(size, VECTOR(*size)[yynode-1] + VECTOR(*size)[yynode-1-(int) VECTOR(*size)[yynode-1]] + 1);
        ++yynode;
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
        igraph_vector_push_back(branch_length, 0.0);
    }
    |
    COLON NUMBER
    {
        igraph_vector_push_back(branch_length, yylval.number);
    }
    ;

%%

void yyerror(igraph_vector_t *tree, igraph_vector_t *size, igraph_vector_t *branch_length, const char *str)
{
    fprintf(stderr, "invalid Newick format or non-binary tree: %s\n", str);
    exit(EXIT_FAILURE);
}
 
int yywrap(void)
{
    return 1;
} 
