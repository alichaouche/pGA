#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "typeDeclaration.h"
#include "./communesFunctions.h"
#include "./vertex_to_cluster.h"
#include "./fractionalEncoding.h"
#include "./sortedCutBasedEncoding.h"
#include "./cutBasedEncoding.h"
#include "./binaryEncoding.h"
#include "gmp.h"

#define  RAND_INV_RANGE(r) ((int) ((RAND_MAX + 1 )/ (r)))


extern char *cocyclesDeBase;
extern int taillePopulation,tauxReproduction,sommeTotalFlux,nbrNoeuds,nbr_constraint,nbrArretes,nbrParties,nbrNoeudsInCluster, sommeTotalFluxReal;
extern double tauxMutation ,averageTraffic;
///extern ul cocyclesDeBaseEntier[1000];
extern int fluxMatrix[1000][1000],*fluxVector, maximumDistanceMatrix[1000][1000],shortestPath[1000][1000];
extern edge *edgeVector,cohabitationConstraintes[100],nonCohabitationConstraintes[100];
extern int max_clusters , min_clusters,min_sizeCluster ,max_sizeCluster,nbrCohabitationConstraintes,nbrNonCohabitationConstraintes;
extern int bestSolution, nbrApparition,nbrGeneration, iteration,nbrCluster;
extern mpz_t gmpCocyclesDeBaseEntier[1000];
extern edge *edgeVectorIntra;


///*******************************************************************
///*************** Les Corps des fonctions ************************///
///*******************************************************************
///**************************************************************************************
int rnd(int v1, int v2)
{
    ///srand(time(NULL));
    ///int returnedValue = (rand()%(v2 - v1) + v1);
    if(v1 == v2){ printf("la fonction rnd : v1 == v2 impossible de générer un nombre aléatoire");}
    int returnedValue = (rand()%(v2 - v1)+v1);
    return returnedValue;
}

int randBis(int RANGE){

    int returnedValue;
    do{
        returnedValue = rand();
    }while(returnedValue >= RANGE * RAND_INV_RANGE(RANGE) );
    returnedValue /= RAND_INV_RANGE(RANGE);
}

///**************************************************************************************
int rndUl(ul v1, ul v2)
{
    return (rand()%(v2 - v1)+v1);
}
///**************************************************************************************
double drand48ForWindows(int v1, int v2)
{
    ///return (double)(rand()) / (double)(RAND_MAX);
    return (double)(rnd(v1,v2)) /1000.0;
}

///**************************************************************************************
void trieGenotypeEntier(mpz_t A[])
{
    int i,j;
    mpz_t tmp; mpz_init(tmp);
    for(i=0; i<nbrParties-2; i++)
    {
        for(j=i+1; j<nbrParties-1; j++)
        {
            if(mpz_cmp(A[i],A[j])<0)
            {
                mpz_set(tmp , A[i]);
                mpz_set(A[i] ,A[j]);
                mpz_set(A[j] ,tmp);
            }
        }
    }
}

///************************************************************************************************************
void afficherMatriceFlux()
{

    printf("afficherMatriceFlux ...\n");
    int i,j;
    for (i = 0; i < nbrNoeuds; i++)
    {
        ///for (j = i+1; j < n; j++) {
        for (j = 0; j < nbrNoeuds; j++)
        {
            if(i==j)
            {
                printf("%d|\t",fluxMatrix[i][j]);
            }
            else
            {
                printf("%d \t",fluxMatrix[i][j]);
            }
        }
        printf("\n");
    }
}
///************************************************************************************************************
void remplissageMatriceFlux()
{

    int i,j,vi;
    for (i = 0; i < nbrNoeuds; i++)
    {
        ///for (j = i+1; j < n; j++) {
        for (j = 0; j < nbrNoeuds; j++)
        {
            if(i==j)
            {
                fluxMatrix[i][j] = 0;
            }
            else
            {
                vi = rnd(0,6);
                fluxMatrix[i][j] = vi;
                fluxMatrix[j][i] = vi; /// pour ne pas calcule deux fois le même poids
            }
        }
    }
}

///************************************************************************************************************
void connecterGraphe()
{

    int i;
    for (i = 0; i < nbrNoeuds-1; i++)
    {
        if(fluxMatrix[i][i+1] == -1)
        {
            fluxMatrix[i][i+1] = 0;
            fluxMatrix[i+1][i] = 0;
        }
    }
}
///**************************************************************************************
int calculerNombreArreteEtFluxVectorEtEdgeVector()
{
    int i,j;
    int nbrArretes = 0,nbrArretesTmp=0; /// puique la boucle commence à partir de 0 donc je devrais rajouter un 1 au nbrArretes

    for(i=0; i<nbrNoeuds; i++)
    {
        for(j=i+1; j<nbrNoeuds; j++)
        {
            if(fluxMatrix[i][j]!= -1 && i!=j)
            {
                nbrArretes++;
            }
        }
    }
    fluxVector = (int*)malloc(nbrArretes*sizeof(int));
    if(fluxVector == NULL) {fprintf(stderr,"Memory allocation for fluxVector failled\n"); (EXIT_FAILURE);    }

    edgeVector = (edge*)malloc(nbrArretes*sizeof(edge));
    if(edgeVector == NULL) {
            fprintf(stderr,"Memory allocation for edgeVector failled\n");  exit(EXIT_FAILURE);
    }
    nbrArretesTmp = 0;
    ///nbrArretesTmp = nbrArretes-1;
    for(i=0; i<nbrNoeuds; i++)
    {
        for(j=i+1; j<nbrNoeuds; j++)
        {
            if(fluxMatrix[i][j]!= -1 && i!=j)
            {
                edgeVector[nbrArretesTmp].numeroEdge = nbrArretesTmp;
                edgeVector[nbrArretesTmp].nouedDepart = i;
                edgeVector[nbrArretesTmp].nouedArrive = j;
                fluxVector[nbrArretesTmp] = fluxMatrix[i][j];
                nbrArretesTmp++;
                ///nbrArretesTmp--;
            }
        }
    }
    /// allocation de la mémoire pour edgeVectorIntra

    edgeVectorIntra = (edge*)malloc(nbrArretes*sizeof(edge));
    if(edgeVectorIntra == NULL) {printf("erreur d'allocation de la memoir pour edgeVectorIntra \n"); exit(-1) ;}

    return nbrArretes;
}

void DisplayFlowVector()
{

    int i;
    ///printf("DisplayFlowVector...\n");
    for(i=0; i<nbrArretes; i++)
    {
        printf("%d ",fluxVector[i]);
    }
    printf("\n");

}
///************************************************************************************************************
int calculSommeTotalFlux()
{


    int i,j,sommeFlux=0;
    for(i=0; i<nbrNoeuds; i++)
    {
        for(j=i+1; j<nbrNoeuds; j++)
        {
            if(fluxMatrix[i][j] > 0 )
            {
                sommeFlux+=fluxMatrix[i][j];
            }
        }
    }
    return sommeFlux;
}

void mon_sleep( int nbr_seconds )
{
    clock_t goal;
    goal = ( nbr_seconds * CLOCKS_PER_SEC ) +  clock();
    while( goal > clock() )
    {
       ; /* loop */
    }

}

///**************************************************************************************
void creeCocycleDeBase()
{
    int i,j;
    ///int base = 2, exp;
    mpz_t puis2; mpz_init(puis2);
    cocyclesDeBase = (char*)malloc(((nbrNoeuds-1)*nbrArretes)*sizeof(char));
    if(cocyclesDeBase == NULL){
            fprintf(stderr,"Memory allocation for cocyclesDeBase failled\n");
            exit(EXIT_FAILURE);
    }
    ///for(i=nbrNoeuds-2; i>=0; i--)
    for(i=0; i<=nbrNoeuds-2; i++)
    {
        ///cocyclesDeBaseEntier[i]=0;
        mpz_init(gmpCocyclesDeBaseEntier[i]);
        ///for(j=nbrArretes-1; j>=0; j--)
        for(j=0; j<=nbrArretes-1; j++)
        {
            ///if (edgeVector[j].nouedDepart == nbrNoeuds-2-i || edgeVector[j].nouedArrive == nbrNoeuds-2-i)
            if (edgeVector[j].nouedDepart == i || edgeVector[j].nouedArrive == i)
            {
                cocyclesDeBase[(i*nbrArretes)+j] = 1; /// y a pas d'inversion l'ordre des cocycles de base
                ///cocyclesDeBaseEntier[i] = cocyclesDeBaseEntier[i] + pow(2,(nbrArretes-1-j));
                ///exp = nbrArretes-j-1;
                mpz_ui_pow_ui(puis2,2,(nbrArretes-j-1));
                ///gmp_printf("puis2 = %Zd  |",puis2);
                mpz_add(gmpCocyclesDeBaseEntier[i],puis2,gmpCocyclesDeBaseEntier[i]);
            }
            else
            {
                cocyclesDeBase[(i)*nbrArretes+j] = 0;
            }
        }
    }

}
///**************************************************************************************
void afficherCocyclesDeBase()
{
    int i,j;

    for(i=0; i<nbrNoeuds-1; i++)
    {
        for(j=(i*nbrArretes); j<(i+1)*nbrArretes; j++)
        {
            printf("%d\t",cocyclesDeBase[j]);
        }
        printf("\n");
    }

}

///**************************************************************************************
/// agencement des gènes
void agencementDeGene(int t[1000])
{

    int i,j,m=nbrNoeuds, t2[m],tmp;

    for(i=0; i<nbrNoeuds; i++)
    {
        t2[i]=t[i];
    }

    for(i=0; i<m-1; i++)
    {
        for(j=i+1; j<m; j++)
        {
            if (t2[j] < t2[i])
            {
                tmp = t2[i];
                t2[i] = t2[j];
                t2[j] = tmp;
            }
        }
    }
    for(i=0; i<m; i++)
    {
        while(t2[i] == t2[i+1])
        {
            for(j=i; j<nbrNoeuds; j++)
            {
                t2[j] = t2[j+1];
            }
            m--;
        }

    }
///*****************************************************

    for(i=0; i<m; i++)
    {
        for(j=0; j<nbrNoeuds; j++)
        {
            if (t[j] == t2[i])
            {
                t[j] = i;
            }
        }
    }
}
///**************************************************************************************
void writeCocyclesDeBaseInFile(FILE *outputFileCocyclesDeBase)
{
    int i,j;
    fprintf(outputFileCocyclesDeBase,"les cocycle de base \n");

    for(i=0; i<nbrNoeuds-1; i++)
    {
        for(j=(i*nbrArretes); j<(i+1)*nbrArretes; j++)
        {
            fprintf(outputFileCocyclesDeBase,"%d\t",cocyclesDeBase[j]);
        }
        fprintf(outputFileCocyclesDeBase,"\n\n");
    }
    fprintf(outputFileCocyclesDeBase,"\n\n");
    for(j=0; j<nbrNoeuds-1; j++)
    {
        ///fprintf(outputFileCocyclesDeBase,"%d\t",cocyclesDeBaseEntier[j]);
        gmp_fprintf(outputFileCocyclesDeBase,"%Zd\t",gmpCocyclesDeBaseEntier[j]);

    }

}

///**************************************************************************************************
/// fonction déterminant la matrice des distances minimales entre tous les paires de noeuds du graphe
/// All Pairs Shortest Path (APSP) Problem
void themaximumDistanceMatrix()
{
    int i,j,k;

    for (i = 0; i <nbrNoeuds; i++)
    {
        for (j = 0; j < nbrNoeuds; j++)
        {
            maximumDistanceMatrix[i][j] = fluxMatrix[i][j];
            shortestPath[i][j] = -1;
        }
        maximumDistanceMatrix[i][i] = 0;
    }

    for (k = 0; k < nbrNoeuds; k++)
    {
        for (i = 0; i < nbrNoeuds; i++)
        {
            for (j = 0; j < nbrNoeuds; j++)
            {

                if ((maximumDistanceMatrix[i][k] + maximumDistanceMatrix[k][j]) > maximumDistanceMatrix[i][j])
                {
                    maximumDistanceMatrix[i][j] = maximumDistanceMatrix[i][k] + maximumDistanceMatrix[k][j];
                    shortestPath[i][j] = k;
                }
            }
        }
    }
} /* maximumDistanceMatrix */


///******************************************************************************************
int compareCroissant (void const *a, void const *b)
{
    /* definir des pointeurs type's et initialise's
       avec les parametres */
    int const *pa = a;
    int const *pb = b;

    /* evaluer et retourner l'etat de l'evaluation (tri déscroissant) */
    return  *pa-*pb;
}

int compareTrieDecroissant (void const *a, void const *b)
{
    /* definir des pointeurs type's et initialise's
       avec les parametres */
    int const *pa = a;
    int const *pb = b;

    /* evaluer et retourner l'etat de l'evaluation (tri déscroissant) */
    return  *pb-*pa;
}

///********************************************************************************************
void writeDetailsProblemInSolutionsFile(FILE* solutionsFile)
{
    int i,j;
    fprintf(solutionsFile,"We are partitioning the following graph : \n");
    for(i=0; i<nbrNoeuds; i++)
    {
        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(solutionsFile,"%d\t",fluxMatrix[i][j]);
        }
        fprintf(solutionsFile,"\n");
    }
    fprintf(solutionsFile,"nbrNoeuds = %d , nbrArretes = %d , somme total des flux = %d \n", nbrNoeuds, nbrArretes, sommeTotalFlux);
    fprintf(solutionsFile,"Constraints: the max cluster size = %d, the max number of clusters = %d \n", max_sizeCluster, max_clusters);
    if(nbrCohabitationConstraintes != 0)
    {
        fprintf(solutionsFile,"Cohabitation constraints : \n");
        for(i=0; i<nbrCohabitationConstraintes; i++) fprintf(solutionsFile,"%d\t%d\n",cohabitationConstraintes[i].nouedDepart,cohabitationConstraintes[i].nouedArrive);
    }
    if(nbrNonCohabitationConstraintes != 0)
    {
        fprintf(solutionsFile,"Non cohabitation constraints : \n");
        for(i=0; i<nbrNonCohabitationConstraintes; i++) fprintf(solutionsFile,"%d\t%d\n",nonCohabitationConstraintes[i].nouedDepart,nonCohabitationConstraintes[i].nouedArrive);
    }

}


///***********************************************************************************************
void openingFile(int numeroGraphe){

///================================================================================

    char cheminBestSolutionSBGNE[100],cheminAllPopulationSBGNE[100],cheminOptimalSolutionSBGNE[150];

    sprintf(cheminBestSolutionSBGNE,"../benchmark/sortedBinaryGroupNumberEncoding/bestSolutionSBGNE%d.txt",numeroGraphe);
    sprintf(cheminAllPopulationSBGNE,"../benchmark/sortedBinaryGroupNumberEncoding/allPopulationSBGNE%d.txt",numeroGraphe);
    sprintf(cheminOptimalSolutionSBGNE,"../benchmark/sortedBinaryGroupNumberEncoding/OptimalSolutionSBGNE%d.txt",numeroGraphe);

    FILE *outputFileSBGNE,*outputFilePopSBGNE,*outputOptimalSolutionFileSBGNE;

    outputFileSBGNE=fopen(cheminBestSolutionSBGNE,"w");
    outputFilePopSBGNE=fopen(cheminAllPopulationSBGNE,"w");
    outputOptimalSolutionFileSBGNE=fopen(cheminOptimalSolutionSBGNE,"w");
    if(outputFilePopSBGNE== NULL)
    {
        printf("Impossible d'ouvrir le fichier allPopulationSBGNE.txt en ecriture\n");
        ///exit(EXIT_FAILURE);
    }
    if(outputFileSBGNE == NULL )
    {
        printf("Impossible d'ouvrir le fichier bestSolutionSBGNE.txt en ecriture\n");
        ///exit(EXIT_FAILURE);
    }

    if(outputOptimalSolutionFileSBGNE== NULL)
    {
        printf("Impossible d'ouvrir le fichier OptimalSolutionSBGNE.txt en ecriture\n");
        ///exit(EXIT_FAILURE);
    }
///================================================================================

    char cheminBestSolutionFC[100],cheminAllPopulationFC[100],cheminOptimalSolutionFC[150];

    sprintf(cheminBestSolutionFC,"../benchmark/fractionnalEncoding/bestSolutionFE%d.txt",numeroGraphe);
    sprintf(cheminAllPopulationFC,"../benchmark/fractionnalEncoding/allPopulationFE%d.txt",numeroGraphe);
    sprintf(cheminOptimalSolutionFC,"../benchmark/fractionnalEncoding/OptimalSolutionFC%d.txt",numeroGraphe);

    FILE *outputFileFC,*outputFilePopFC,*outputOptimalSolutionFileFC;
    outputFileFC=fopen(cheminBestSolutionFC,"w");
    outputFilePopFC=fopen(cheminAllPopulationFC,"w");
    outputOptimalSolutionFileFC=fopen(cheminOptimalSolutionFC,"w");
    if(outputFileFC == NULL)
    {
        printf("Impossible d'ouvrir le fichier bestSolution.txt en ecriture\n");
        ///exit(EXIT_FAILURE);
    }
    if(outputFilePopFC== NULL)
    {
        printf("Impossible d'ouvrir le fichier allPopulation.txt en ecriture\n");
        ///exit(EXIT_FAILURE);
    }
    if(outputOptimalSolutionFileFC== NULL)
    {
        printf("Impossible d'ouvrir le fichier OptimalSolutionFC.txt en ecriture\n");
        ///exit(EXIT_FAILURE);
    }
///================================================================================

    char cheminBestSolutionBE[100],cheminAllPopulationBE[100],cheminOptimalSolutionBE[150];

    sprintf(cheminBestSolutionBE,"../benchmark/binaryEncoding/bestSolutionBE%d.txt",numeroGraphe);
    sprintf(cheminAllPopulationBE,"../benchmark/binaryEncoding/allPopulationBE%d.txt",numeroGraphe);
    sprintf(cheminOptimalSolutionBE,"../benchmark/binaryEncoding/OptimalSolutionBE%d.txt",numeroGraphe);

    FILE *outputFileBE,*outputFilePopBE,*outputOptimalSolutionFileBE;
    outputFileBE=fopen(cheminBestSolutionBE,"w");
    outputFilePopBE=fopen(cheminAllPopulationBE,"w");
    outputOptimalSolutionFileBE=fopen(cheminOptimalSolutionBE,"w");

    if(outputFileBE == NULL)
    {
        printf("Impossible d'ouvrir le fichier outputFile en ecriture\n");
        ///exit(EXIT_FAILURE);
    }
    if(outputFilePopBE== NULL)
    {
        printf("Impossible d'ouvrir le fichier outputFilePop en ecriture\n");
        ///exit(EXIT_FAILURE);
    }
    if(outputOptimalSolutionFileBE== NULL)
    {
        printf("Impossible d'ouvrir le fichier OptimalSolutionBE en ecriture\n");
        ///exit(EXIT_FAILURE);
    }
///================================================================================
    char cheminBestSolutionCBE[100],cheminAllPopulationCBE[100],cheminOptimalSolutionCBE[150];

    sprintf(cheminBestSolutionCBE,"../benchmark/cutBasedEncoding/bestSolutionCBE%d.txt",numeroGraphe);
    sprintf(cheminAllPopulationCBE,"../benchmark/cutBasedEncoding/allPopulationCBE%d.txt",numeroGraphe);
    sprintf(cheminOptimalSolutionCBE,"../benchmark/cutBasedEncoding/OptimalSolutionCBE%d.txt",numeroGraphe);

    FILE *outputFileCBE,*outputFilePopCBE,*outputOptimalSolutionFileCBE;
    outputFileCBE=fopen(cheminBestSolutionCBE,"w");
    outputFilePopCBE=fopen(cheminAllPopulationCBE,"w");
    outputOptimalSolutionFileCBE=fopen(cheminOptimalSolutionCBE,"w");
    if(outputFileCBE == NULL)
    {
        printf("Impossible d'ouvrir le fichier outputFile en ecriture\n");
        ///exit(EXIT_FAILURE);
    }
    if(outputFilePopCBE== NULL)
    {
        printf("Impossible d'ouvrir le fichier outputFilePop en ecriture\n");
        ///exit(EXIT_FAILURE);
    }
    if(outputOptimalSolutionFileCBE== NULL)
    {
        printf("Impossible d'ouvrir le fichier OptimalSolutionCBE en ecriture\n");
        ///exit(EXIT_FAILURE);
    }
///================================================================================
    char cheminBestSolutionSCBE[100],cheminAllPopulationSCBE[100],cheminOptimalSolutionSCBE[150];

    sprintf(cheminBestSolutionSCBE,"../benchmark/sortedCutBasedEncoding/bestSolutionSCBE%d.txt",numeroGraphe);
    sprintf(cheminAllPopulationSCBE,"../benchmark/sortedCutBasedEncoding/allPopulationSCBE%d.txt",numeroGraphe);
    sprintf(cheminOptimalSolutionSCBE,"../benchmark/sortedCutBasedEncoding/OptimalSolutionSCBE%d.txt",numeroGraphe);

    FILE *outputFileSCBE,*outputFilePopSCBE,*outputOptimalSolutionFileSCBE;
    outputFileSCBE=fopen(cheminBestSolutionSCBE,"w");
    outputFilePopSCBE=fopen(cheminAllPopulationSCBE,"w");
    outputOptimalSolutionFileSCBE=fopen(cheminOptimalSolutionSCBE,"w");
    if(outputFileSCBE == NULL )
    {
        printf("Impossible d'ouvrir le fichier bestSolutionSCBE en ecriture\n");
        ///exit(EXIT_FAILURE);
    }
    if(outputFilePopSCBE== NULL)
    {
        printf("Impossible d'ouvrir le fichier allPopulationSCBE en ecriture\n");
        exit(EXIT_FAILURE);
    }
    if(outputOptimalSolutionFileSCBE== NULL)
    {
        printf("Impossible d'ouvrir le fichier OptimalSolutionSCBE en ecriture\n");
        exit(EXIT_FAILURE);
    }
///================================================================================
    char cheminBestSolutionPMPEdgeBased[150],cheminAllPopulationPMPEdgeBased[150],cheminOptimalSolutionPMPEdgeBased[150];

    sprintf(cheminBestSolutionPMPEdgeBased,"../benchmark/pMedianEncoding/affectationParBaseSurArretes/bestSolutionPMP%d.txt",numeroGraphe);
    sprintf(cheminAllPopulationPMPEdgeBased,"../benchmark/pMedianEncoding/affectationParBaseSurArretes/allPopulationPMP%d.txt",numeroGraphe);
    sprintf(cheminOptimalSolutionPMPEdgeBased,"../benchmark/pMedianEncoding/affectationParBaseSurArretes/OptimalSolutionPMP%d.txt",numeroGraphe);

    FILE *outputFilePMPEdgeBased,*outputFilePopPMPEdgeBased,*outputOptimalSolutionFilePMPEdgeBased;
    outputFilePMPEdgeBased=fopen(cheminBestSolutionPMPEdgeBased,"w");
    outputFilePopPMPEdgeBased=fopen(cheminAllPopulationPMPEdgeBased,"w");
    outputOptimalSolutionFilePMPEdgeBased=fopen(cheminOptimalSolutionPMPEdgeBased,"w");
    if(outputFilePMPEdgeBased == NULL)
    {
        printf("Impossible d'ouvrir le fichier données bestSolutionPMP affectation base sur arretes en ecriture\n");
        exit(EXIT_FAILURE);
    }
    if(outputFilePopPMPEdgeBased== NULL)
    {
        printf("Impossible d'ouvrir le fichier données en ecriture\n");
        exit(EXIT_FAILURE);
    }
    if(outputOptimalSolutionFilePMPEdgeBased== NULL)
    {
        printf("Impossible d'ouvrir le fichier OptimalSolutionPMP.txt en ecriture\n");
        exit(EXIT_FAILURE);
    }

///================================================================================
    char cheminBestSolutionPMPClusterBased[150],cheminAllPopulationPMPClusterBased[150],cheminOptimalSolutionPMPClusterBased[150];


    sprintf(cheminBestSolutionPMPClusterBased,"../benchmark/pMedianEncoding/affectationBaseSurCluster/bestSolutionPMP%d.txt",numeroGraphe);
    sprintf(cheminAllPopulationPMPClusterBased,"../benchmark/pMedianEncoding/affectationBaseSurCluster/allPopulationPMP%d.txt",numeroGraphe);
    sprintf(cheminOptimalSolutionPMPClusterBased,"../benchmark/pMedianEncoding/affectationBaseSurCluster/OptimalSolutionPMP%d.txt",numeroGraphe);

    FILE *outputFilePMPClusterBased,*outputFilePopPMPClusterBased,*outputOptimalSolutionFilePMPClusterBased;
    outputFilePMPClusterBased=fopen(cheminBestSolutionPMPClusterBased,"w");
    outputFilePopPMPClusterBased=fopen(cheminAllPopulationPMPClusterBased,"w");
    outputOptimalSolutionFilePMPClusterBased=fopen(cheminOptimalSolutionPMPClusterBased,"w");
    if(outputFilePMPClusterBased == NULL)
    {
        printf("Impossible d'ouvrir le fichier données en ecriture\n");
        exit(EXIT_FAILURE);
    }
    if(outputFilePopPMPClusterBased == NULL)
    {
        printf("Impossible d'ouvrir le fichier données en ecriture\n");
        exit(EXIT_FAILURE);
    }
    if(outputOptimalSolutionFilePMPClusterBased == NULL)
    {
        printf("Impossible d'ouvrir le fichier OptimalSolutionPMP.txt en ecriture\n");
        exit(EXIT_FAILURE);
    }
///================================================================================
    char cheminBestSolutionVTC[100],cheminAllPopulationVTC[100],cheminOptimalSolutionVTC[150];

    sprintf(cheminBestSolutionVTC,"../benchmark/vertexToClusterEncoding/bestSolutionVTC%d.txt",numeroGraphe);
    sprintf(cheminAllPopulationVTC,"../benchmark/vertexToClusterEncoding/allPopulationVTC%d.txt",numeroGraphe);
    sprintf(cheminOptimalSolutionVTC,"../benchmark/vertexToClusterEncoding/OptimalSolutionVTC%d.txt",numeroGraphe);



    FILE *outputFileVTC,*outputFilePopVTC,*outputOptimalSolutionFileVTC;
    outputFileVTC=fopen(cheminBestSolutionVTC,"w");
    outputFilePopVTC=fopen(cheminAllPopulationVTC,"w");
    outputOptimalSolutionFileVTC=fopen(cheminOptimalSolutionVTC,"w");

    if(outputFileVTC == NULL)
    {
        printf("Impossible d'ouvrir le fichier bestSolutionVTC en ecriture\n");
        exit(EXIT_FAILURE);
    }
    if(outputFilePopVTC== NULL)
    {
        printf("Impossible d'ouvrir le fichier allPopulationVTC en ecriture\n");
        exit(EXIT_FAILURE);
    }
    if(outputOptimalSolutionFileVTC== NULL)
    {
        printf("Impossible d'ouvrir le fichier OptimalSolutionVTC en ecriture\n");
        exit(EXIT_FAILURE);
    }
///================================================================================
    char cheminBestSolutionEE[100],cheminAllPopulationEE[100],cheminOptimalSolutionEE[150];

    sprintf(cheminBestSolutionEE,"../benchmark/ecartEncoding/bestSolutionEE%d.txt",numeroGraphe);
    sprintf(cheminAllPopulationEE,"../benchmark/ecartEncoding/allPopulationEE%d.txt",numeroGraphe);
    sprintf(cheminOptimalSolutionEE,"../benchmark/ecartEncoding/OptimalSolutionEE%d.txt",numeroGraphe);


    FILE *outputFileEE,*outputFilePopEE, *outputOptimalSolutionFileEE;;
    outputFileEE=fopen(cheminBestSolutionEE,"w");
    outputFilePopEE=fopen(cheminAllPopulationEE,"w");
    outputOptimalSolutionFileEE=fopen(cheminOptimalSolutionEE,"w");

    if(outputFileEE == NULL)
    {
        printf("Impossible d'ouvrir le fichier bestSolutionEE en ecriture\n");
        exit(EXIT_FAILURE);
    }
    if(outputFilePopEE== NULL)
    {
        printf("Impossible d'ouvrir le fichier allPopulationEE en ecriture\n");
        exit(EXIT_FAILURE);
    }
    if(outputOptimalSolutionFileEE== NULL)
    {
        printf("Impossible d'ouvrir le fichier OptimalSolutionEE en ecriture\n");
        exit(EXIT_FAILURE);
    }
///================================================================================
    char cheminBestSolutionSVTC[100]="../benchmark/SortedVTC/bestSolutionSVTC.txt";
    char cheminAllPopulationSVTC[100]="../benchmark/SortedVTC/allPopulationSVTC.txt";
    char cheminOptimalSolutionSVTC[150]="../benchmark/SortedVTC/OptimalSolutionSVTC.txt";

    sprintf(cheminBestSolutionSVTC,"../benchmark/SortedVTC/bestSolutionSVTC%d.txt",numeroGraphe);
    sprintf(cheminAllPopulationSVTC,"../benchmark/SortedVTC/allPopulationSVTC%d.txt",numeroGraphe);
    sprintf(cheminOptimalSolutionSVTC,"../benchmark/SortedVTC/OptimalSolutionSVTC%d.txt",numeroGraphe);

    FILE *outputFileSVTC,*outputFilePopSVTC,*outputOptimalSolutionFileSVTC;
    outputFileSVTC=fopen(cheminBestSolutionSVTC,"w");
    outputFilePopSVTC=fopen(cheminAllPopulationSVTC,"w");
    outputOptimalSolutionFileSVTC=fopen(cheminOptimalSolutionSVTC,"w");

    if(outputFileSVTC == NULL)
    {
        printf("Impossible d'ouvrir le fichier bestSolutionSVTC en ecriture\n");
        exit(EXIT_FAILURE);
    }
    if(outputFilePopSVTC== NULL)
    {
        printf("Impossible d'ouvrir le fichier allPopulationSVTC en ecriture\n");
        exit(EXIT_FAILURE);
    }
    if(outputOptimalSolutionFileSVTC== NULL)
    {
        printf("Impossible d'ouvrir le fichier OptimalSolutionSVTC en ecriture\n");
        exit(EXIT_FAILURE);
    }
///================================================================================
}

///================================================================================

int isConnex(){

    int i,j,k;
    int voisin[2000],nonVoisin[2000],v=0,nv=0;
    voisinEtNonVoisinArray(voisin,nonVoisin,&v,&nv);

    if(nv!=0)
    {
        for(i=0; i<v; i++)
        {
            for(j=0; j<nv; j++)
            {
                if(fluxMatrix[voisin[i]][nonVoisin[j]]>=0)
                {
                    voisin[v]=nonVoisin[j];
                    v++;
                    /// rajouter le non voisin à la liste des voisins puis le supprimer de la liste des voisins
                    for(k=j; k<nv; k++)
                    {
                        nonVoisin[k] = nonVoisin[k+1];
                    }
                    nv--;
                }

            }
        }
    }
    if(nv == 0) return 1;
    else return 0;
}

///================================================================================
void voisinEtNonVoisinArray(int voisin[],int nonVoisin[],int *v,int *nv)
{

    int i,j,k;
    int vLocal=0,nvLocal=0;
    for(i=1; i<nbrNoeuds; i++)
    {
        if(fluxMatrix[0][i]>=0)
        {
            voisin[vLocal] = i;
            vLocal++;
        }
        else
        {
            nonVoisin[nvLocal] = i;
            nvLocal++;
        }
    }
    *v = vLocal;
    *nv = nvLocal;

}
///================================================================================
/// on fait appel à cette fonction uniquement dans le cas ou le graphe n'est pas connexe
void graphConnexite(int voisin[],int nonVoisin[],int v,int nv)
{
    int i,j,k;
    for(i=0; i<v; i++)
    {
        for(j=0; j<nv; j++)
        {
            if(fluxMatrix[voisin[i]][nonVoisin[j]]>=0)
            {
                voisin[v]=nonVoisin[j];
                v++;
                for(k=j; k<nv; k++)
                {
                    nonVoisin[k] = nonVoisin[k+1];
                }
                nv--;
            }
        }
    }
    if(nv!=0)
    {
        if(v==0){
            fluxMatrix[0][nonVoisin[0]]=0; /// une arrêtes fictive
            fluxMatrix[nonVoisin[0]][0]=0;
            voisin[0]=nonVoisin[0];
            v++;
            for(k=0; k<nv; k++)
            {
                nonVoisin[k] = nonVoisin[k+1];
            }
            nv--;
            if(!isConnex())
            {
                graphConnexite(voisin,nonVoisin,v,nv);
            }
        }
        else if(v > 0){
            fluxMatrix[voisin[v-1]][nonVoisin[0]]=0; /// une arrêtes fictive
            fluxMatrix[nonVoisin[0]][voisin[v-1]]=0;
            voisin[v]=nonVoisin[0];
            v++;
            for(k=0; k<nv; k++)
            {
                nonVoisin[k] = nonVoisin[k+1];
            }
            nv--;
            if(!isConnex())
            {
                graphConnexite(voisin,nonVoisin,v,nv);
            }
        }
    }
}


void cohabitationFictivesEdges(){

    int i,na,nd;
    for(i=0; i<nbrCohabitationConstraintes; i++)
    {
        nd = cohabitationConstraintes[i].nouedDepart;
        na = cohabitationConstraintes[i].nouedArrive;
        if( fluxMatrix[nd][na] < 0  ||  fluxMatrix[na][nd] <0){
            fluxMatrix[nd][na] = 0; fluxMatrix[na][nd] = 0;
        }
    }
}


///************************************************************************************
int maxFlux(){
    int i,j,maxFlow;
    maxFlow = fluxMatrix[0][0];
    for(i=0;i<nbrNoeuds;i++){
        for(j=0;j<nbrNoeuds;j++){
            if(i!=j && maxFlow < fluxMatrix[i][j]){
                maxFlow = fluxMatrix[i][j];
            }
        }
    }
    return maxFlow;
}

/**
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 cette fonction devrait avoir en entrée le vecteur de la meilleur solution puis dans
 un premier temps la convertir en SVTC puis caluler le nombre de clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

///**********************************************************************************
int SVTC_Representation(int solution[])
{
    int i,j,k,m=nbrNoeuds, t1[m],t[m],sortIndex,tmp,exist,min;

    for(j=0; j<nbrNoeuds; j++)
    {
        t1[j] = -1;
    }

    t[0]=solution[0]; sortIndex=1;
    for(j=1;j<nbrNoeuds;j++){
        exist = 0;
        for(k=0;k<sortIndex;k++){
            if(solution[j]==t[k]){
                exist = 1;
                break;
            }
        }
        if(!exist){
            t[sortIndex]=solution[j];
            sortIndex++;
        }
    }

    for(j=0; j<sortIndex; j++)
    {
        for(k=0; k<nbrNoeuds; k++)
        {
            if (solution[k] == t[j] && t1[k] <0)
            {
                solution[k] = j;
                t1[k] = 1;
            }
        }
    }
    return sortIndex; /// ça represente le nombre de cluster de la partition
}

/// conversion phenotype => genotype des codage choisis
/// l'élement à inseré va prendre la place d'un autre individu au sein de la population
/// comment choisir cet element ???

void conversionVersFVTC(int solution[nbrNoeuds], partitionFVTC* populationFVTC,int indiceFVTC){
    int i;
    (populationFVTC+indiceFVTC)->genotype[nbrNoeuds] = nbrCluster/nbrNoeuds;

    for(i=0;i<nbrNoeuds;i++){
        (populationFVTC+indiceFVTC)->genotype[i] = solution[i]/nbrCluster;
    }
}

void conversionVersEA(int solution[nbrNoeuds], partitionEA* populationEA, int indiceEA){
    int i,j,nd,na;

    for(i=0;i<nbrArretes;i++){
        nd = edgeVector[i].nouedDepart;
        na = edgeVector[i].nouedArrive;
        if(solution[nd] == solution[na]){
            (populationEA+indiceEA)->genotype[i] = 0;
        }
        else {
            (populationEA+indiceEA)->genotype[i] = 1;
        }
    }

  }
///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void conversionVersPMEB(int solution[nbrNoeuds], partitionPMP* populationPMP, int indicePMEB){
    int i,j,nd,na,firstIndex,tailleCluster,indiceMedian;
    int cluster[nbrNoeuds];
    for(i=0;i<nbrNoeuds;i++)
        (populationPMP+indicePMEB)->genotype[i] = 0;
    ///Déterminer les liens entre les sommets pour choisir le médians
    for(i=0;i<nbrCluster;i++){
            tailleCluster = 0;
        for(j=0;j<nbrNoeuds;j++){
            if(solution[j]==i){
                cluster[tailleCluster]= solution[j];
                tailleCluster++;
            }
        }
        indiceMedian = choixDuSommetMedian(cluster, tailleCluster);

        (populationPMP+indicePMEB)->genotype[indiceMedian] = 1;
    }

}
///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/// compter le nombre de liens avec l'ensemble de sommets pour choisi le median
int choixDuSommetMedian(int cluster[], int tailleCluster){
    int i,j,maxLink=0, nbrLink = 0,indiceM = cluster[0];
    for(i=0;i<tailleCluster;i++){
        nbrLink = 0;
        for(j=0;j<tailleCluster;j++){
            if(fluxMatrix[cluster[i]][cluster[j]]>=0 && i!=j){
                nbrLink++;
            }
        }
        if(maxLink < nbrLink){
            maxLink = nbrLink;
            indiceM = i;
        }
    }
    return indiceM;
}
///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void conversionVersPMCA(int solution[nbrNoeuds], partitionPMP* populationPMP, int indicePMCA){
    int i,j,nd,na,firstIndex,tailleCluster,indiceMedian;
    int cluster[nbrNoeuds]; ///nbrNoeuds c'est le nombre max des clusters qu'on peut avoir.
    for(i=0;i<nbrNoeuds;i++)
        (populationPMP+indicePMCA)->genotype[i] = 0;
    ///Déterminer les liens entre les sommets pour choisir le médians
    for(i=0;i<nbrCluster;i++){
            tailleCluster = 0;
        for(j=0;j<nbrNoeuds;j++){
            if(solution[j]==i){
                cluster[tailleCluster]= solution[j]; /// c'est pour choisir le median au sien d'un cluster.
                tailleCluster++;
            }
        }
        indiceMedian = choixDuSommetMedian(cluster, tailleCluster);
        (populationPMP+indicePMCA)->genotype[indiceMedian] = 1;
    }

}
///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void conversionVersDC(int solution[nbrNoeuds], partitionDC* populationDC, int indiceDC){

    int i,j,nd,na,firstIndex,selectedCut;
    for(i=0;i<nbrParties*(nbrNoeuds-1);i++)
        (populationDC+indiceDC)->genotype[i] = 0;
    if(nbrCluster > (nbrParties+1)){
        /// choisir nbrParties parties clusters.
        for(i=0;i<nbrParties;i++){
            selectedCut=rnd(0,nbrCluster);
            for(j=0;j<nbrNoeuds;j++){
                if(solution[j]==selectedCut){
                    (populationDC+indiceDC)->genotype[i*nbrParties+j] = 1;
                }
            }

        }

    }
    else {
        for(i=0;i<nbrCluster;i++){
            for(j=0;j<nbrNoeuds;j++){
                if(solution[j]==i){
                    (populationDC+indiceDC)->genotype[i*nbrParties+j] = 1;
                }
            }

        }
    }
}
///*********************************************************************************************
void writeBestSolutionInFilePGA(partitionFVTC *solutionDominanteFVTC, partitionEA *solutionDominanteEA,
                                partitionDC *solutionDominanteDC, partitionPMP *solutionDominantePMEB , partitionPMP *solutionDominantePMCA , FILE *outputFilePGA ){
    int i,PMEB,PMCA, EA, FVTC, DC;
    PMEB = solutionDominantePMEB->coutCoupeNormalise;
    PMCA = solutionDominantePMCA->coutCoupeNormalise;
    FVTC  = solutionDominanteFVTC->coutCoupeNormalise;
    EA = solutionDominanteEA->coutCoupeNormalise;
    DC = solutionDominanteDC->coutCoupeNormalise;

    if(FVTC > PMEB && FVTC > PMCA && FVTC > EA && FVTC > DC){

        fprintf(outputFilePGA,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,solutionDominanteFVTC->id,
                solutionDominanteFVTC->coutCoupe,solutionDominanteFVTC->fitness,solutionDominanteFVTC->coutCoupeNormalise,
                solutionDominanteFVTC->contrainteViole,solutionDominanteFVTC->nbrCluster);
        for(i=0; i<=nbrNoeuds; i++)
        {
            fprintf(outputFilePGA,"%0.2f ",solutionDominanteFVTC->genotype[i]);
        }
        fprintf(outputFilePGA,"\t\t\t\t");
        for(i=0; i<=nbrNoeuds; i++)
        {
            fprintf(outputFilePGA,"%d ",solutionDominanteFVTC->phenotype[i]);
        }
        fprintf(outputFilePGA,"\n");

   }
   else if(EA > PMEB && EA > PMCA && EA > FVTC && EA > DC ){

        fprintf(outputFilePGA,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,solutionDominanteEA->id,
                solutionDominanteEA->coutCoupe,solutionDominanteEA->fitness,solutionDominanteEA->coutCoupeNormalise,
                solutionDominanteEA->contrainteViole,solutionDominanteEA->nbrCluster);
        for(i=0; i<=nbrNoeuds; i++)
        {
            fprintf(outputFilePGA,"%d ",solutionDominanteEA->genotype[i]);
        }
        fprintf(outputFilePGA,"\t\t\t\t");
        for(i=0; i<=nbrNoeuds; i++)
        {
            fprintf(outputFilePGA,"%d ",solutionDominanteEA->phenotype[i]);
        }
        fprintf(outputFilePGA,"\n");

   }
   else if(DC > PMEB && DC > PMCA && DC > FVTC && DC > EA){

        fprintf(outputFilePGA,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,solutionDominanteDC->id,
                solutionDominanteDC->coutCoupe,solutionDominanteDC->fitness,solutionDominanteDC->coutCoupeNormalise,
                solutionDominanteDC->contrainteViole,solutionDominanteDC->nbrCluster);
        for(i=0; i<=nbrNoeuds; i++)
        {
            fprintf(outputFilePGA,"%d ",solutionDominanteDC->genotype[i]);
        }
        fprintf(outputFilePGA,"\t\t\t\t");
        for(i=0; i<=nbrNoeuds; i++)
        {
            fprintf(outputFilePGA,"%d ",solutionDominanteDC->phenotype[i]);
        }
        fprintf(outputFilePGA,"\n");

   }
   else if(PMCA > PMEB && PMCA > FVTC && PMCA > EA && PMCA > DC){

        fprintf(outputFilePGA,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,solutionDominantePMCA->id,
                solutionDominantePMCA->coutCoupe,solutionDominantePMCA->fitness,solutionDominantePMCA->coutCoupeNormalise,
                solutionDominantePMCA->contrainteViole,solutionDominantePMCA->medians);
        for(i=0; i<=nbrNoeuds; i++)
        {
            fprintf(outputFilePGA,"%d ",solutionDominantePMCA->genotype[i]);
        }
        fprintf(outputFilePGA,"\t\t\t\t");
        for(i=0; i<=nbrNoeuds; i++)
        {
            fprintf(outputFilePGA,"%d ",solutionDominantePMCA->phenotype[i]);
        }
        fprintf(outputFilePGA,"\n");

   }
   else {

        fprintf(outputFilePGA,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,solutionDominantePMEB->id,
                solutionDominantePMEB->coutCoupe,solutionDominantePMEB->fitness,solutionDominantePMEB->coutCoupeNormalise,
                solutionDominantePMEB->contrainteViole,solutionDominantePMEB->medians);
        for(i=0; i<=nbrNoeuds; i++)
        {
            fprintf(outputFilePGA,"%d ",solutionDominantePMEB->genotype[i]);
        }
        fprintf(outputFilePGA,"\t\t\t\t");
        for(i=0; i<=nbrNoeuds; i++)
        {
            fprintf(outputFilePGA,"%d ",solutionDominantePMEB->phenotype[i]);
        }
        fprintf(outputFilePGA,"\n");

   }

}
