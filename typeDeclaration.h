#ifndef TYPEDECLARATION_H_INCLUDED
#define TYPEDECLARATION_H_INCLUDED
#include "gmp.h"
typedef unsigned long ul;
///=============================================================================================
typedef struct
{
    int id;
    char *genotype;
    int phenotype[10000];
    ///************************************
    int contrainteViole;
    int nbrCluster;
    int clustersSize[1000];
    int coutCoupeNormalise;
    int constraintVector[10];
    ///************************************
    int coutCoupe;
    float fitness;
    /// 26/11/2015 : ajouter pour le calcul de fitness scaling
    float expectedValue;

} partitionEA;
///=============================================================================================
typedef struct
{
    int id;
    ul genotypeEntier[10000]; ///je vais utiliser ce vecteur pour générer des solution sans redandance.
    int *cocycles; /// cette matrice contient la solution générer aléatoirement et sur sa base on va créer les partition
    ///int cocycles[100000];
    ///int genotype[100000];
    int *genotype;
    int phenotype[10000];
    ///************************************
    int contrainteViole;
    int nbrCluster;
    int clustersSize[1000];
    int coutCoupeNormalise;
    int constraintVector[10];
    ///************************************
    int coutCoupe;
    float fitness;
    ///************************************
    float expectedValue;
} partitionDC;
///=============================================================================================

typedef struct
{
    int id;
    ///ul genotypeEntier[1000]; ///je vais utiliser ce vecteur pour générer des solution sans redandance.
    mpz_t genotypeEntier[10000]; ///je vais utiliser ce vecteur pour générer des solution sans redandance.
    ///ul genotypeEcart[1000];
    mpz_t genotypeEcart[10000];
    ///int *cocycles; /// cette matrice contient la solution générer aléatoirement et sur sa base on va créer les partition
    int cocycles[100000];
    ///int *genotype;
    int genotype[100000];
    int phenotype[10000];
    ///************************************
    int contrainteViole;
    int nbrCluster;
    int clustersSize[1000];
    int coutCoupeNormalise;
    int constraintVector[10];
    ///************************************
    int coutCoupe;
    float fitness;
    ///************************************
    float expectedValue;
} partitionCGE;
///=============================================================================================

typedef struct
{
    int id;
    float genotype[10000];
    int phenotype[10000];
///************************************
    int contrainteViole;
    int nbrCluster;
    int clustersSize[10000];
    int coutCoupeNormalise;
    int constraintVector[10];
///************************************
    int coutCoupe;
    float fitness;
///************************************
    float expectedValue;
} partitionFVTC;
///=============================================================================================

typedef struct
{
    int id;
    int genotype[10000];
    int phenotype[10000];
    int mediansVertices[10000];
    int medians;
    int nonMediansVertices[10000];
    int nonMedians;
///************************************
    int contrainteViole;
    int clustersSize[1000];
    int coutCoupeNormalise;
    int constraintVector[10];
///************************************
    int coutCoupe;
    float fitness;
///************************************
    float expectedValue;
} partitionPMP;
///=============================================================================================
typedef struct
{
    int id;
    int genotype[10000];
///************************************
    int contrainteViole;
    int nbrCluster;
    int clustersSize[1000];
    int coutCoupeNormalise;
    int constraintVector[10];
///************************************
    int coutCoupe;
    float fitness;
///************************************
    float expectedValue;
} partitionVTC;
///=============================================================================================
typedef struct
{
    int id;
    mpz_t genotypeEntier[10000]; ///je vais utiliser ce vecteur pour générer des solution sans redandance.
    ///int cocycles[1000]; /// cette matrice contient la solution générer aléatoirement et sur sa base on va créer les partition
    int *cocycles; /// cette matrice contient la solution générer aléatoirement et sur sa base on va créer les partition
    ///int genotype[1000];
    int *genotype;
    int phenotype[10000];
    ///************************************
    int contrainteViole;
    int nbrCluster;
    int clustersSize[1000];
    int coutCoupeNormalise;
    int constraintVector[10];
    ///************************************
    int coutCoupe;
    float fitness;
    ///************************************
    float expectedValue;
} partitionIC;
///=============================================================================================

typedef struct
{
    int id;
    int genotypeEntier[10000]; ///je vais utiliser ce vecteur pour générer des solution sans redandance.
    char *genotype;     /// c'est le vecteur qui contient les solutions
    int phenotype[10000];
    ///************************************
    int contrainteViole;
    int nbrCluster;
    int clustersSize[1000];
    int coutCoupeNormalise;
    int constraintVector[10];
    ///************************************
    int coutCoupe;
    float fitness;
    ///************************************
    float expectedValue;
} partitionVAE;
///==================================================================================================
typedef struct
{
    int id;
    int genotype[10000];
///************************************
    int contrainteViole;
    int nbrCluster;
    int clustersSize[1000];
    int coutCoupeNormalise;
    int constraintVector[10];
///************************************
    int coutCoupe;
    float fitness;
///************************************
    float expectedValue;
} partitionSVTC;
///=============================================================================================
typedef struct
{
    int id;
    int genotype[10000];
///************************************
    int contrainteViole;
    int nbrCluster;
    int clustersSize[1000];
    int coutCoupeNormalise;
    int constraintVector[10];
///************************************
    int coutCoupe;
    float fitness;
///************************************
    float expectedValue;
///************************************
    int delimiterVector[100];
    int phenotype[1000];

} partitionGCE;
///=============================================================================================


#endif // TYPEDECLARATION_H_INCLUDED
