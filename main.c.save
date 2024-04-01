#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "typeDeclaration.h"
#include "./communesFunctions.h"
#include "./compilationConditionnelle.h"
#include "./ecartEncoding.h"
#include "./sortedBinaryGroupNumberEncoding.h"
#include "./sortedCutBasedEncoding.h"
#include "./cutBasedEncoding.h"
#include "./pMediansEncoding.h"
#include "./fractionalEncoding.h"
#include "./vertex_to_cluster.h"
#include "./binaryEncoding.h"

#define max_runs 30
#define NBR_INSTANCES 27
#define SUB_POP_SIZE 100


char *cocyclesDeBase;
int taillePopulation=100,tauxReproduction=10,sommeTotalFlux,nbrNoeuds,nbr_constraint=2,nbrArretes,
    nbrParties,nbrNoeudsInCluster, regeneration = 10;
double tauxMutation = 0.002,averageTraffic;
int fluxMatrix[1000][1000],fluxVector[1000], maximumDistanceMatrix[1000][1000],shortestPath[1000][1000];
edge *edgeVector,cohabitationConstraintes[100],nonCohabitationConstraintes[100], *edgeVectorIntra;
int max_clusters, min_clusters= 2,min_sizeCluster = 1,max_sizeCluster, nbrCluster = 0;
int bestSolution, nbrApparition = 0, iteration = 1,nbrCohabitationConstraintes,nbrNonCohabitationConstraintes;
int nbrRun = 1;
mpz_t gmpCocyclesDeBaseEntier[1000];


int main()
{

    printf("===================================================================================\n");
    printf("          application de l'algorithme genetique pour le probleme GPP               \n");
    printf("    utilisation de différents types de codage pour la representation des solution  \n");
    printf("===================================================================================\n");

    char choix = 'Y';
    int i,j,nbrGeneration=100;
    int numeroGraphe;
    int nbrSolution, indicePMEB, indiceFVTC, indiceEA, indiceDC, indicePMCA;

    printf("Veuillez introduire le numero de graphe \n");
    scanf("%d",&numeroGraphe);
    char cheminGraph[100],cheminCocyclesDeBase[100],cheminInputFileConstrainte[100];
    FILE *inputFile,*outputFileCocyclesDeBase,*inputFileConstrainte;
    srand(time(NULL));

///=============================================================================================================
#if EA
    partitionEA *populationEA1,*populationEA2, *solutionDominanteEA;
    populationEA1 = (partitionEA*)malloc(taillePopulation*sizeof (partitionEA));
    if(populationEA1==NULL) printf("memory allocation failed for the populationEA1 \n");
    populationEA2 = (partitionEA*)malloc(taillePopulation*sizeof(partitionEA));
    if(populationEA2==NULL) printf("memory allocation failed for the populationEA2 \n");
    solutionDominanteEA = (partitionEA*)malloc(sizeof(partitionEA));
    if(solutionDominanteEA==NULL) printf("memory allocation failed for the solutionDominanteEA \n");

    char cheminBestSolutionEA[200],cheminAllPopulationEA[200],cheminOptimalSolutionEA[200];
    FILE *outputFileEA,*outputFilePopEA,*outputOptimalSolutionFileEA;

#endif // EA
///=============================================================================================================
#if DC
    partitionDC *populationDC1,*populationDC2, *solutionDominanteDC;
    populationDC1 = (partitionDC*)malloc(taillePopulation*sizeof (partitionDC));
    if(populationDC1==NULL) printf("memory allocation failed for the populationDC1 \n");
    populationDC2 = (partitionDC*)malloc(taillePopulation*sizeof(partitionDC));
    if(populationDC2==NULL) printf("memory allocation failed for the populationDC2 \n");
    solutionDominanteDC = (partitionDC*)malloc(sizeof(partitionDC));
    if(solutionDominanteDC==NULL) printf("memory allocation failed for the solutionDominanteDC \n");

    char cheminBestSolutionDC[200],cheminAllPopulationDC[200],cheminOptimalSolutionDC[200];
    FILE *outputFileDC,*outputFilePopDC,*outputOptimalSolutionFileDC;

#endif // DC
///=============================================================================================================
#if CGE
    partitionCGE *populationCGE1,*populationCGE2, *solutionDominanteCGE;
    populationCGE1 = (partitionCGE*)malloc(taillePopulation*sizeof (partitionCGE));
    if(populationCGE1==NULL) printf("memory allocation failed for the populationCGE1 \n");
    populationCGE2 = (partitionCGE*)malloc(taillePopulation*sizeof(partitionCGE));
    if(populationCGE2==NULL) printf("memory allocation failed for the populationCGE2 \n");
    solutionDominanteCGE = (partitionCGE*)malloc(sizeof(partitionCGE));
    if(solutionDominanteCGE==NULL) printf("memory allocation failed for the solutionDominanteCGE \n");

    char cheminBestSolutionCGE[200],cheminAllpopulationCGE[200],cheminOptimalSolutionCGE[200];
    FILE *outputFileCGE,*outputFilePopCGE, *outputOptimalSolutionFileCGE;

#endif // CGE
///=============================================================================================================
#if FVTC
    partitionFVTC *populationFVTC1,*populationFVTC2, *solutionDominanteFVTC;
    populationFVTC1 = (partitionFVTC*)malloc(taillePopulation*sizeof (partitionFVTC));
    if(populationFVTC1==NULL) printf("memory allocation failed for the populationFVTC1 \n");
    populationFVTC2 = (partitionFVTC*)malloc(taillePopulation*sizeof(partitionFVTC));
    if(populationFVTC2==NULL) printf("memory allocation failed for the populationFVTC2 \n");
    solutionDominanteFVTC = (partitionFVTC*)malloc(sizeof(partitionFVTC));
    if(solutionDominanteFVTC==NULL) printf("memory allocation failed for the solutionDominanteFVTC \n");

    char cheminBestSolutionFVTC[200],cheminAllPopulationFVTC[200],cheminOptimalSolutionFVTC[200];
    FILE *outputFileFVTC,*outputFilePopFVTC,*outputOptimalSolutionFileFVTC;

#endif // FVTC
///=============================================================================================================

#if PMEB
    partitionPMP *populationPMEB1,*populationPMEB2, *solutionDominantePMEB;
    populationPMEB1 = (partitionPMP*)malloc(taillePopulation*sizeof (partitionPMP));
    if(populationPMEB1==NULL) printf("memory allocation failed for the populationPMP1 \n");
    populationPMEB2 = (partitionPMP*)malloc(taillePopulation*sizeof(partitionPMP));
    if(populationPMEB2==NULL) printf("memory allocation failed for the populationPMP2 \n");
    solutionDominantePMEB = (partitionPMP*)malloc(sizeof(partitionPMP));
    if(solutionDominantePMEB==NULL) printf("memory allocation failed for the solutionDominantePMP \n");

    char cheminBestSolutionPMPEdgeBased[200],cheminAllPopulationPMPEdgeBased[200],cheminOptimalSolutionPMPEdgeBased[200];
    FILE *outputFilePMPEdgeBased,*outputFilePopPMPEdgeBased,*outputOptimalSolutionFilePMPEdgeBased;
#endif // PMEB
#if PMCA
    partitionPMP *populationPMCA1,*populationPMCA2, *solutionDominantePMCA;
    populationPMCA1 = (partitionPMP*)malloc(taillePopulation*sizeof (partitionPMP));
    if(populationPMCA1==NULL) printf("memory allocation failed for the populationPMP1 \n");
    populationPMCA2 = (partitionPMP*)malloc(taillePopulation*sizeof(partitionPMP));
    if(populationPMCA2==NULL) printf("memory allocation failed for the populationPMP2 \n");
    solutionDominantePMCA = (partitionPMP*)malloc(sizeof(partitionPMP));
    if(solutionDominantePMCA==NULL) printf("memory allocation failed for the solutionDominantePMP \n");
    char cheminBestSolutionPMPClusterBased[200],cheminAllPopulationPMPClusterBased[200],cheminOptimalSolutionPMPClusterBased[200];
    FILE *outputFilePMPClusterBased,*outputFilePopPMPClusterBased,*outputOptimalSolutionFilePMPClusterBased;
#endif // PMCA

///=============================================================================================================
#if VAE
    partitionVAE *populationVAE1,*populationVAE2, *solutionDominanteVAE;
    populationVAE1 = (partitionVAE*)malloc(taillePopulation*sizeof (partitionVAE));
    if(populationVAE1==NULL) printf("memory allocation failed for the populationVAE1 \n");
    populationVAE2 = (partitionVAE*)malloc(taillePopulation*sizeof(partitionVAE));
    if(populationVAE2==NULL) printf("memory allocation failed for the populationVAE2 \n");
    solutionDominanteVAE = (partitionVAE*)malloc(sizeof(partitionVAE));
    if(solutionDominanteVAE==NULL) printf("memory allocation failed for the solutionDominanteVAE \n");

    char cheminBestSolutionVAE[200],cheminAllPopulationVAE[200],cheminOptimalSolutionVAE[200];
    FILE *outputFileVAE,*outputFilePopVAE,*outputOptimalSolutionFileVAE;
#endif // VAE
///=============================================================================================================
#if IC
    partitionIC *populationIC1,*populationIC2, *solutionDominanteIC;
    populationIC1 = (partitionIC*)malloc(taillePopulation*sizeof (partitionIC));
    if(populationIC1==NULL) printf("memory allocation failed for the populationIC1 \n");
    populationIC2 = (partitionIC*)malloc(taillePopulation*sizeof(partitionIC));
    if(populationIC2==NULL) printf("memory allocation failed for the populationIC2 \n");
    solutionDominanteIC = (partitionIC*)malloc(sizeof(partitionIC));
    if(solutionDominanteIC==NULL) printf("memory allocation failed for the solutionDominanteIC \n");

    char cheminBestSolutionIC[200],cheminAllPopulationIC[200],cheminOptimalSolutionIC[200];
    FILE *outputFileIC,*outputFilePopIC,*outputOptimalSolutionFileIC;

#endif // IC
///=============================================================================================================
#if DVTC
    partitionVTC *populationVTC1,*populationVTC2, *solutionDominanteVTC;
    populationVTC1 = (partitionVTC*)malloc(taillePopulation*sizeof (partitionVTC));
    if(populationVTC1==NULL) printf("memory allocation failed for the populationVTC1 \n");
    populationVTC2 = (partitionVTC*)malloc(taillePopulation*sizeof(partitionVTC));
    if(populationVTC2==NULL) printf("memory allocation failed for the populationVTC2 \n");
    solutionDominanteVTC = (partitionVTC*)malloc(sizeof(partitionVTC));
    if(solutionDominanteVTC==NULL) printf("memory allocation failed for the solutionDominanteVTC \n");

    char cheminBestSolutionVTC[200],cheminAllPopulationVTC[200],cheminOptimalSolutionVTC[200];
    FILE *outputFileVTC,*outputFilePopVTC,*outputOptimalSolutionFileVTC;

#endif // DVTC
///=============================================================================================================
#if SVTC
    partitionSVTC *populationSVTC1,*populationSVTC2, *solutionDominanteSVTC;
    populationSVTC1 = (partitionSVTC*)malloc(taillePopulation*sizeof (partitionSVTC));
    if(populationSVTC1==NULL) printf("memory allocation failed for the populationSVTC1 \n");
    populationSVTC2 = (partitionSVTC*)malloc(taillePopulation*sizeof(partitionSVTC));
    if(populationSVTC2==NULL) printf("memory allocation failed for the populationSVTC2 \n");
    solutionDominanteSVTC = (partitionSVTC*)malloc(sizeof(partitionSVTC));
    if(solutionDominanteSVTC==NULL) printf("memory allocation failed for the solutionDominanteSVTC \n");

    char cheminBestSolutionSVTC[200],cheminAllPopulationSVTC[200],cheminOptimalSolutionSVTC[200];
    FILE *outputFileSVTC,*outputFilePopSVTC,*outputOptimalSolutionFileSVTC;

#endif // SVTC
///=============================================================================================================
#if PGA
    char cheminBestSolutionPGA[200],cheminAllPopulationPGA[200],cheminOptimalSolutionPGA[200];
    FILE *outputFilePGA,*outputFilePopPGA,*outputOptimalSolutionFilePGA;
#endif // PGA



    ///int iteration=1 ;
    int bestSolutionVTC,bestSolutionIterationVTC,ES_VTC = 0,nbrApparitionVTC;
    int bestSolutionSVTC,bestSolutionIterationSVTC,ES_SVTC = 0,nbrApparitionSVTC;
    int bestSolutionFVTC,bestSolutionIterationFVTC,ES_FVTC = 0,nbrApparitionFVTC;
    int bestSolutionVAE,bestSolutionIterationVAE,ES_VAE = 0,nbrApparitionVAE;
    int bestSolutionEA,bestSolutionIterationEA,ES_EA = 0,nbrApparitionEA;
    int bestSolutionDC,bestSolutionIterationDC,ES_DC = 0,nbrApparitionDC;
    int bestSolutionIC,bestSolutionIterationIC,ES_IC = 0,nbrApparitionIC;
    int bestSolutionCGE,bestSolutionIterationCGE,ES_CGE = 0,nbrApparitionCGE;
    int bestSolutionPMEB,bestSolutionIterationPMEB,ES_PMEB = 0,nbrApparitionPMEB;
    int bestSolutionPMCA,bestSolutionIterationPMCA,ES_PMCA = 0,nbrApparitionPMCA;

    clock_t t1,t2;
    double temps;

    for(; numeroGraphe <=NBR_INSTANCES; numeroGraphe++)
    {
#include "opningFiles.txt"
        ///lecture de la matrice des flux
        fscanf(inputFile,"%d",&nbrNoeuds);
        for( i=0; i<nbrNoeuds; i++)
        {
            for( j=0; j<nbrNoeuds; j++)
            {
                fscanf(inputFile,"%d",&fluxMatrix[i][j]);
                /// cette partie est utile pour connecter le graphe
                if(i!=j && fluxMatrix[i][j]==0 )fluxMatrix[i][j] = -1;
            }
        }

        if(!isConnex())
        {
            int voisin[2000],nonVoisin[2000],v,nv;
            voisinEtNonVoisinArray(voisin,nonVoisin,&v,&nv);
            /**
            printf("le graphe N %d avec %d sommets est non connexe v = %d \n",numeroGraphe,nbrNoeuds,v);
            numeroGraphe++;
            mon_sleep(pause);
            goto x;
            */
            graphConnexite(voisin,nonVoisin,v,nv);
        }
        /**
        else{
            printf("le graphe N %d avec %d sommets est connexe\n",numeroGraphe,nbrNoeuds);
            numeroGraphe++;
            mon_sleep(pause);
            goto x;
        }*/
        ///lecture de fichier des contraintes
        nbr_constraint=2;
        fscanf(inputFileConstrainte,"%d",&max_sizeCluster);
        fscanf(inputFileConstrainte,"%d",&nbrCohabitationConstraintes);
        if(nbrCohabitationConstraintes!=0)
        {
            nbr_constraint++;
            for( i=0; i<nbrCohabitationConstraintes; i++)
            {
                fscanf(inputFileConstrainte,"%d %d",&cohabitationConstraintes[i].nouedDepart,&cohabitationConstraintes[i].nouedArrive);
                /// Ajouter le 12-10-2016 pour connecter les noeuds à cohabiter et ceux qui ne doivent pas cohabiter
                if (fluxMatrix[cohabitationConstraintes[i].nouedDepart][cohabitationConstraintes[i].nouedArrive] < 0)
                {
                    fluxMatrix[cohabitationConstraintes[i].nouedDepart][cohabitationConstraintes[i].nouedArrive] = 0;
                    fluxMatrix[cohabitationConstraintes[i].nouedArrive][cohabitationConstraintes[i].nouedDepart] = 0;
                }
            }
        }
        fscanf(inputFileConstrainte,"%d",&nbrNonCohabitationConstraintes);
        if(nbrNonCohabitationConstraintes!=0)
        {
            nbr_constraint++;
            for( i=0; i<nbrNonCohabitationConstraintes; i++)
            {
                fscanf(inputFileConstrainte,"%d %d",&nonCohabitationConstraintes[i].nouedDepart,&nonCohabitationConstraintes[i].nouedArrive);
                /// Ajouter le 12-10-2016 pour connecter les noeuds à cohabiter et ceux qui ne doivent pas cohabiter
                /**
                if (fluxMatrix[nonCohabitationConstraintes[i].nouedDepart][cohabitationConstraintes[i].nouedArrive] < 0) {
                    fluxMatrix[nonCohabitationConstraintes[i].nouedDepart][cohabitationConstraintes[i].nouedArrive] = 0;
                    fluxMatrix[nonCohabitationConstraintes[i].nouedArrive][cohabitationConstraintes[i].nouedDepart] = 0;
                }
                */
            }
        }
        /// 02/03/2017 Cette appel me permet de relier par une arête fictive les sommet àcohabiter (si l'arête n'existe pas)
        //cohabitationFictivesEdges(nbrCohabitationConstraintes,cohabitationConstraintes);
        /// calculer le vecteur des flux tout en retournant le nombre d'arrêtes dans le graphe
        nbrArretes = calculerNombreArreteEtFluxVectorEtEdgeVector();

        creeCocycleDeBase();
        writeCocyclesDeBaseInFile(outputFileCocyclesDeBase);

        nbrParties = ceil((float)nbrNoeuds / max_sizeCluster);
        ///max_clusters = nbrParties;
        /// Rajouter par CHAOUCHE ALI le : 25.01.2017
        max_clusters = nbrNoeuds;
        sommeTotalFlux = calculSommeTotalFlux();
///^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        printf("Parallel genetic algorithm : FVTC, EA, DC ,PMEB,PMCA \n");
        printf("\nle graphe utilise est %s  \n",cheminGraph);
///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        writeDetailsProblemInSolutionsFile(outputFilePGA);
        fprintf(outputFilePGA,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
        fprintf(outputFilePGA,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");
///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        while(nbrRun<=max_runs)
        {
            printf("\n===================================================\n");
            printf("\nThe number of the RUN est %d \n",nbrRun);
            printf("\n===================================================\n");
///************************ Vertex based catégorie *********************************************
            t1=clock();
#if DVTC
            writeDetailsProblemInSolutionsFile(outputFileVTC);
            fprintf(outputFileVTC,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFileVTC,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");

            writeDetailsProblemInSolutionsFile(outputFilePopVTC);
            fprintf(outputFilePopVTC,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFilePopVTC,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");

            generatePopulationVTC(populationVTC1,0);
            checkContrainstAndFitnessPenalizationVTC(populationVTC1);
            calculCoutCoupeEtFitnessVTC(populationVTC1);
#if writingPopulationInFile
            writeSolutionInFileVTC(populationVTC1,outputFilePopVTC,iteration);
#endif // writingPopulationInFile

            bestSolutionVTC=findTheBestSolutionVTC(populationVTC1);
            *solutionDominanteVTC=*(populationVTC1+bestSolutionVTC);
            bestSolutionIterationVTC= 1;
#if writingPopulationInFile
            writeBestSolutionInFileVTC(solutionDominanteVTC,outputFileVTC,iteration);
#endif // writingPopulationInFile
            nbrApparitionVTC=1;
#endif // DVTC
///==============================================================================================

#if SVTC
            writeDetailsProblemInSolutionsFile(outputFileSVTC);
            fprintf(outputFileSVTC,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFileSVTC,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");
            writeDetailsProblemInSolutionsFile(outputFilePopSVTC);
            fprintf(outputFilePopSVTC,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFilePopSVTC,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");

            generatePopulationSVTC(populationSVTC1,0);
            agencementDesGenotype(populationSVTC1);
            checkContrainstAndFitnessPenalizationSVTC(populationSVTC1);
            calculCoutCoupeEtFitnessSVTC(populationSVTC1);
#if writingPopulationInFile
            writeSolutionInFileSVTC(populationSVTC1,outputFilePopSVTC,iteration);
#endif // writingPopulationInFile

            bestSolutionSVTC=findTheBestSolutionSVTC(populationSVTC1);
            *solutionDominanteSVTC=*(populationSVTC1+bestSolutionSVTC);
            bestSolutionIterationSVTC= 1;
#if writingPopulationInFile
            writeBestSolutionInFileSVTC(solutionDominanteSVTC,outputFileSVTC,iteration);
#endif // writingPopulationInFile
            nbrApparitionSVTC=1;
#endif // SVTC
///=====================================================================================================
#if FVTC
            writeDetailsProblemInSolutionsFile(outputFileFVTC);
            fprintf(outputFileFVTC,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFileFVTC,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");
            writeDetailsProblemInSolutionsFile(outputFilePopFVTC);
            fprintf(outputFilePopFVTC,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFilePopFVTC,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");

            /// génération de la population initiale
            generatePopulationFVTC(populationFVTC1,0);
            calculPhenotypeFVTC(populationFVTC1);
            checkContrainstAndFitnessPenalizationFVTC(populationFVTC1);
            calculCoutCoupeEtFitnessFVTC(populationFVTC1,0,taillePopulation);
#if writingPopulationInFile
            writeSolutionInFileFVTC(populationFVTC1,outputFilePopFVTC,iteration);
#endif // writingPopulationInFile
            bestSolutionFVTC=findTheBestSolutionFVTC(populationFVTC1);
            *solutionDominanteFVTC=*(populationFVTC1+bestSolutionFVTC);
            bestSolutionIterationFVTC=1;
            nbrApparitionFVTC=1;
#if writingPopulationInFile
            writeBestSolutionInFileFVTC(solutionDominanteFVTC,outputFileFVTC,iteration);
#endif // writingPopulationInFile

#endif // FVTC

#if VAE
            writeDetailsProblemInSolutionsFile(outputFileVAE);
            fprintf(outputFileVAE,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFileVAE,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");

            writeDetailsProblemInSolutionsFile(outputFilePopVAE);
            fprintf(outputFilePopVAE,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFilePopVAE,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");

            ///genererPopulationInitialeVAE(populationVAE1);
            genererPopulationInitialeRandomlyVAE(populationVAE1,0);
            getPartitionFromSolutionVAE(populationVAE1);
            checkContrainstAndFitnessPenalizationVAE(populationVAE1);
            calculCoutCoupeEtFitnessVAE(populationVAE1);
#if writingPopulationInFile
            writeSolutionInFileVAE(populationVAE1,outputFilePopVAE,iteration); /// le fichier contenant toute la population
#endif // writingPopulationInFile

            ///initialisation de la solution dominante
            bestSolutionVAE=findTheBestSolutionVAE(populationVAE1);
            *solutionDominanteVAE=*(populationVAE1+bestSolutionVAE);
            bestSolutionIterationVAE = 1;

#if writingPopulationInFile
            writeBestSolutionInFileVAE(solutionDominanteVAE,outputFileVAE,iteration);/// fichier contenat que les meilleurs solutions de chaque itération
#endif // writingPopulationInFile
            nbrApparitionVAE=1;
#endif // VAE

///************************ Edge based catégorie *********************************************
#if EA
            writeDetailsProblemInSolutionsFile(outputFileEA);
            fprintf(outputFileEA,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFileEA,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");
            writeDetailsProblemInSolutionsFile(outputFilePopEA);
            fprintf(outputFilePopEA,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFilePopEA,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");

            generatePopulationEA(populationEA1,0);
            ///generatePopulationRandomlyEA(populationEA1);
            getPartitionFromSolutionEA(populationEA1);
            ///getPartitionFromSolutionWithoutRepetitionEA(populationEA1);
            checkContrainstAndFitnessPenalizationEA(populationEA1);
            calculCoutCoupeEtFitnessEA(populationEA1);
            ///calculCoutCoupeEtFitnessWithFlowVectorEA(populationEA1);
#if writingPopulationInFile
            writeSolutionInFileEA(populationEA1,outputFilePopEA,iteration);
#endif // writingPopulationInFile
            ///initialisation de la solution dominante
            bestSolutionEA=findTheBestSolutionEA(populationEA1);
            *solutionDominanteEA=*(populationEA1+bestSolutionEA);
            bestSolutionIterationEA = 1;
            nbrApparitionEA=1;
#if writingPopulationInFile
            writeBestSolutionInFileEA(solutionDominanteEA,outputFileEA,iteration);
#endif // writingPopulationInFile

#endif // EA

#if DC

            writeDetailsProblemInSolutionsFile(outputFileDC);
            fprintf(outputFileDC,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFileDC,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");
            writeDetailsProblemInSolutionsFile(outputFilePopDC);
            fprintf(outputFilePopDC,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFilePopDC,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");

            generatePopulationDC(populationDC1,0);
            ///calculerGenotypeDC(populationDC1);  /// l'appel à cette fonction est dans generatePopulationDC
            getPartitionFromSolutionDC(populationDC1);
            checkContrainstAndFitnessPenalizationDC(populationDC1);
            calculCoutCoupeEtFitnessDC(populationDC1);
#if writingPopulationInFile
            writeSolutionInFileDC(populationDC1,outputFilePopDC,iteration);
#endif // writingPopulationInFile
            bestSolutionDC=findTheBestSolutionDC(populationDC1);
            *solutionDominanteDC=*(populationDC1+bestSolution);
            bestSolutionIterationDC=1;
#if writingPopulationInFile
            writeBestSolutionInFileDC(solutionDominanteDC,outputFileDC,iteration);
#endif // writingPopulationInFile

            nbrApparitionDC=1;
#endif // DC

#if IC

            writeDetailsProblemInSolutionsFile(outputFileIC);
            fprintf(outputFileIC,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFileIC,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");
            writeDetailsProblemInSolutionsFile(outputFilePopIC);
            fprintf(outputFilePopIC,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFilePopIC,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");

            generatePopulationWithoutRedandancyIC(populationIC1,0);
            ///calculerGenotypeIC(populationIC1); l'appel à cette fonction est dans la fonction generatePopulationWithoutRedundancy
            getPartitionFromSolutionIC(populationIC1);
            checkContrainstAndFitnessPenalizationIC(populationIC1);
            calculCoutCoupeEtFitnessIC(populationIC1);
#if writingPopulationInFile
            writeSolutionInFileIC(populationIC1,outputFilePopIC,iteration);
#endif // writingPopulationInFile

            ///initialisation de la solution dominante
            bestSolutionIC=findTheBestSolutionIC(populationIC1);
            *solutionDominanteIC=*(populationIC1+bestSolutionIC);
            bestSolutionIterationIC =1;
#if writingPopulationInFile
            writeBestSolutionInFileIC(solutionDominanteIC,outputFileIC,iteration);
#endif // writingPopulationInFile
            nbrApparitionIC=1;

#endif // IC


#if CGE
            writeDetailsProblemInSolutionsFile(outputFileCGE);
            fprintf(outputFileCGE,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFileCGE,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");
            writeDetailsProblemInSolutionsFile(outputFilePopCGE);
            fprintf(outputFilePopCGE,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFilePopCGE,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");

            generatePopulationWithoutRedandancyCGE(populationCGE1,0);

            calculerGenotypeCGE(populationCGE1);
            getPartitionFromSolutionCGE(populationCGE1);
            checkContrainstAndFitnessPenalizationCGE(populationCGE1);
            calculCoutCoupeEtFitnessCGE(populationCGE1);

#if writingPopulationInFile
            writeSolutionInFileCGE(populationCGE1,outputFilePopCGE,iteration);
#endif // writingPopulationInFile

            ///initialisation de la solution dominante
            bestSolutionCGE=findTheBestSolutionCGE(populationCGE1);

            *solutionDominanteCGE=*(populationCGE1+bestSolutionCGE);
            bestSolutionIterationCGE=1;
#if writingPopulationInFile
            ///writeBestSolutionInFileCGE(solutionDominanteCGE,outputFileCGE,iteration);
            writeBestSolutionInFileCGE((populationCGE1+bestSolutionCGE),outputFileCGE,iteration);
#endif // writingPopulationInFile
            nbrApparitionCGE=1;
#endif // CGE
///********************************P-Median based encoding ******************************************
#if PMEB
            writeDetailsProblemInSolutionsFile(outputFilePMPEdgeBased);
            fprintf(outputFilePMPEdgeBased,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFilePMPEdgeBased,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");
            writeDetailsProblemInSolutionsFile(outputFilePopPMPEdgeBased);
            fprintf(outputFilePopPMPEdgeBased,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFilePopPMPEdgeBased,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");

            generatePopulationPMP(populationPMEB1,0);
            calculerMediansNonMediansVertices(populationPMEB1);
            affectationDesNonMediansParArretes(populationPMEB1);
            checkContrainstAndFitnessPenalizationPMP(populationPMEB1);
            calculcoutCoupeeEtFitnessPMP(populationPMEB1);
#if writingPopulationInFile
            writeSolutionInFilePMP(populationPMEB1,outputFilePopPMPEdgeBased,iteration);
#endif // writingPopulationInFile

            ///initialisation de la solution dominante
            bestSolutionPMEB=findTheBestSolutionPMP(populationPMEB1);
            *solutionDominantePMEB=*(populationPMEB1+bestSolutionPMEB);
            bestSolutionIterationPMEB = 1;
#if writingPopulationInFile
            writeBestSolutionInFilePMP(solutionDominantePMEB,outputFilePMPEdgeBased,iteration);
#endif // writingPopulationInFile
            nbrApparitionPMEB=1;


#endif // PMEB

#if PMCA
            writeDetailsProblemInSolutionsFile(outputFilePMPClusterBased);
            fprintf(outputFilePMPClusterBased,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFilePMPClusterBased,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");
            writeDetailsProblemInSolutionsFile(outputFilePopPMPClusterBased);
            fprintf(outputFilePopPMPClusterBased,"\nLe numero de graphe = %d ,Le nombre de generation = %d , la taille de la population = %d \n\n",numeroGraphe,nbrGeneration,taillePopulation);
            fprintf(outputFilePopPMPClusterBased,"itr\tid\tct\tft\t\tctN\tcv\tnbC\tsol\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpar\n");

            generatePopulationPMP(populationPMCA1,0);
            calculerMediansNonMediansVertices(populationPMCA1);
            affectationDesNonMediansBaseSurClusterSize(populationPMCA1);
            checkContrainstAndFitnessPenalizationPMP(populationPMCA1);
            calculcoutCoupeeEtFitnessPMP(populationPMCA1);
#if writingPopulationInFile
            writeSolutionInFilePMP(populationPMCA1,outputFilePopPMPClusterBased,iteration);
#endif // writingPopulationInFile


            ///initialisation de la solution dominante
            bestSolutionPMCA=findTheBestSolutionPMP(populationPMCA1);
            *solutionDominantePMCA=*(populationPMCA1+bestSolutionPMCA);
            bestSolutionIterationPMCA = 1;

#if writingPopulationInFile
            writeBestSolutionInFilePMP(solutionDominantePMCA,outputFilePMPClusterBased,iteration);
#endif // writingPopulationInFile
            nbrApparitionPMCA=1;

#endif // PMCA
///*************************************************************************************************************
            writeBestSolutionInFilePGA(solutionDominanteFVTC, solutionDominanteEA,solutionDominanteDC, solutionDominantePMEB, solutionDominantePMCA,outputFilePGA);
///*************************************************************************************************************

            for(iteration =2 ; iteration <= nbrGeneration; iteration++)
            {
                /// immigration des solutions entre les population
                if(iteration % 10 == 0)
                {
                    printf("le nombre d iteration = %d \n",iteration);
                    nbrSolution =0;

                    /// difusion des solutions PMEB
                    for(indicePMEB = 0; indicePMEB < 10; indicePMEB++)
                    {

                        nbrCluster = SVTC_Representation((populationPMEB1+indicePMEB)->phenotype);
                        indiceFVTC = rnd(10,SUB_POP_SIZE);
                        if((populationPMEB1+indicePMEB)->coutCoupeNormalise < (populationFVTC1+indiceFVTC)->coutCoupeNormalise)
                            conversionVersFVTC((populationPMEB1+indicePMEB)->phenotype,populationFVTC1, indiceFVTC);

                        indiceEA = rnd(10,SUB_POP_SIZE);
                        if((populationPMEB1+indicePMEB)->coutCoupeNormalise < (populationEA1+indiceEA)->coutCoupeNormalise)
                            conversionVersEA((populationPMEB1+indicePMEB)->phenotype,populationEA1,indiceEA);

                        indiceDC = rnd(10,SUB_POP_SIZE);
                        if((populationPMEB1+indicePMEB)->coutCoupeNormalise < (populationDC1+indiceDC)->coutCoupeNormalise)
                            conversionVersDC((populationPMEB1+indicePMEB)->phenotype,populationDC1,indiceDC);

                        indicePMCA = rnd(10,SUB_POP_SIZE);
                        if((populationPMEB1+indicePMEB)->coutCoupeNormalise < (populationPMCA1+indicePMCA)->coutCoupeNormalise)
                            conversionVersPMCA((populationPMEB1+indicePMEB)->phenotype,populationPMCA1,indicePMCA);
                    }

                    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    for(indiceFVTC = 0; indiceFVTC < 10; indiceFVTC++)
                    {

                        nbrCluster = SVTC_Representation((populationFVTC1+indiceFVTC)->phenotype);
                        indicePMEB = rnd(10,SUB_POP_SIZE);
                        if((populationFVTC1+indiceFVTC)->coutCoupeNormalise < (populationPMEB1+indicePMEB)->coutCoupeNormalise)
                            conversionVersPMEB((populationFVTC1+indiceFVTC)->phenotype,populationPMEB1, indicePMEB);

                        indiceEA = rnd(10,SUB_POP_SIZE);
                        if((populationFVTC1+indiceFVTC)->coutCoupeNormalise < (populationEA1+indiceEA)->coutCoupeNormalise)
                            conversionVersEA((populationFVTC1+indiceFVTC)->phenotype,populationEA1,indiceEA);

                        indiceDC = rnd(10,SUB_POP_SIZE);
                        if((populationFVTC1+indiceFVTC)->coutCoupeNormalise < (populationDC1+indiceDC)->coutCoupeNormalise)
                            conversionVersDC((populationFVTC1+indiceFVTC)->phenotype,populationDC1,indiceDC);

                        indicePMCA = rnd(10,SUB_POP_SIZE);
                        if((populationFVTC1+indicePMEB)->coutCoupeNormalise < (populationPMCA1+indicePMCA)->coutCoupeNormalise)
                            conversionVersPMCA((populationFVTC1+indiceFVTC)->phenotype,populationPMCA1,indicePMCA);
                    }

                    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    for(indiceEA = 0; indiceEA < 10; indiceEA++)
                    {

                        nbrCluster = SVTC_Representation((populationEA1+indiceEA)->phenotype);
                        indicePMEB = rnd(10,SUB_POP_SIZE);
                        if((populationEA1+indiceEA)->coutCoupeNormalise < (populationPMEB1+indicePMEB)->coutCoupeNormalise)
                            conversionVersPMEB((populationEA1+indiceEA)->phenotype,populationPMEB1, indicePMEB);

                        indiceFVTC = rnd(10,SUB_POP_SIZE);
                        if((populationEA1+indiceEA)->coutCoupeNormalise < (populationEA1+indiceEA)->coutCoupeNormalise)
                            conversionVersFVTC((populationEA1+indiceEA)->phenotype,populationFVTC1,indiceFVTC);

                        indiceDC = rnd(10,SUB_POP_SIZE);
                        if((populationEA1+indiceEA)->coutCoupeNormalise < (populationDC1+indiceDC)->coutCoupeNormalise)
                            conversionVersDC((populationEA1+indiceEA)->phenotype,populationDC1,indiceDC);

                        indicePMCA = rnd(10,SUB_POP_SIZE);
                        if((populationEA1+indiceEA)->coutCoupeNormalise < (populationPMCA1+indicePMCA)->coutCoupeNormalise)
                            conversionVersPMCA((populationEA1+indiceEA)->phenotype,populationPMCA1,indicePMCA);
                    }

                    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    for(indiceDC = 0; indiceDC < 10; indiceDC++)
                    {

                        nbrCluster = SVTC_Representation((populationDC1+indiceDC)->phenotype);
                        indicePMEB = rnd(10,SUB_POP_SIZE);
                        if((populationDC1+indiceDC)->coutCoupeNormalise < (populationPMEB1+indicePMEB)->coutCoupeNormalise)
                            conversionVersPMEB((populationDC1+indiceDC)->phenotype,populationPMEB1, indicePMEB);

                        indiceFVTC = rnd(10,SUB_POP_SIZE);
                        if((populationDC1+indiceDC)->coutCoupeNormalise < (populationDC1+indiceDC)->coutCoupeNormalise)
                            conversionVersFVTC((populationDC1+indiceDC)->phenotype,populationFVTC1,indiceFVTC);

                        indiceEA = rnd(10,SUB_POP_SIZE);
                        if((populationDC1+indiceDC)->coutCoupeNormalise < (populationEA1+indiceEA)->coutCoupeNormalise)
                            conversionVersEA((populationDC1+indiceDC)->phenotype,populationEA1,indiceEA);

                        indicePMCA = rnd(10,SUB_POP_SIZE);
                        if((populationDC1+indiceDC)->coutCoupeNormalise < (populationPMCA1+indicePMCA)->coutCoupeNormalise)
                            conversionVersPMCA((populationDC1+indiceDC)->phenotype,populationPMCA1,indicePMCA);
                    }

                    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    for(indicePMCA = 0; indicePMCA < 10; indicePMCA++)
                    {

                        nbrCluster = SVTC_Representation((populationPMCA1+indicePMCA)->phenotype);

                        indicePMEB = rnd(10,SUB_POP_SIZE);
                        if((populationPMCA1+indicePMCA)->coutCoupeNormalise < (populationPMEB1+indicePMEB)->coutCoupeNormalise)
                            conversionVersPMEB((populationPMCA1+indicePMCA)->phenotype,populationPMEB1, indicePMEB);

                        indiceFVTC = rnd(10,SUB_POP_SIZE);
                        if((populationPMCA1+indicePMCA)->coutCoupeNormalise < (populationDC1+indiceDC)->coutCoupeNormalise)
                            conversionVersFVTC((populationPMCA1+indicePMCA)->phenotype,populationFVTC1,indiceFVTC);

                        indiceEA = rnd(10,SUB_POP_SIZE);
                        if((populationPMCA1+indicePMCA)->coutCoupeNormalise < (populationEA1+indiceEA)->coutCoupeNormalise)
                            conversionVersEA((populationPMCA1+indicePMCA)->phenotype,populationEA1,indiceEA);

                        indiceDC = rnd(10,SUB_POP_SIZE);
                        if((populationPMCA1+indicePMCA)->coutCoupeNormalise < (populationDC1+indiceDC)->coutCoupeNormalise)
                            conversionVersPMCA((populationPMCA1+indicePMCA)->phenotype,populationDC1,indiceDC);
                    }

                }
#if DVTC
                vertexToClusterEncoding(nbrGeneration,outputFileVTC,outputFilePopVTC,outputOptimalSolutionFileVTC,
                                        populationVTC1,populationVTC2, solutionDominanteVTC, iteration,&bestSolutionIterationVTC,
                                        &nbrApparitionVTC);
#endif // DVTC

#if SVTC
                sortedVertexToClusterEncoding(nbrGeneration,outputFileSVTC,outputFilePopSVTC,outputOptimalSolutionFileSVTC,
                                              populationSVTC1,populationSVTC2, solutionDominanteSVTC, iteration,&bestSolutionIterationSVTC, &nbrApparitionSVTC);
#endif // SVTC

#if FVTC
                fractionalEncoding(nbrGeneration,outputFileFVTC,outputFilePopFVTC,outputOptimalSolutionFileFVTC,populationFVTC1,populationFVTC2, solutionDominanteFVTC,
                                   iteration,&bestSolutionIterationFVTC, &nbrApparitionFVTC);
#endif // FVTC

#if VAE
                sortedBinaryGroupNumberEncoding(nbrGeneration,outputFileVAE,outputFilePopVAE,outputOptimalSolutionFileVAE,populationVAE1,populationVAE2,
                                                solutionDominanteVAE,iteration,&bestSolutionIterationVAE, &nbrApparitionVAE);
#endif // VAE
#if EA
                binaryEncoding(nbrGeneration,outputFileEA,outputFilePopEA,outputOptimalSolutionFileEA,populationEA1,
                               populationEA2, solutionDominanteEA,iteration,&bestSolutionIterationEA, &nbrApparitionEA);
#endif // EA
#if DC
                cutBasedEncoding(nbrGeneration,outputFileDC,outputFilePopDC,outputOptimalSolutionFileDC,populationDC1,
                                 populationDC2, solutionDominanteDC,iteration,&bestSolutionIterationDC, &nbrApparitionDC);
#endif // DC
#if IC
                cutBasedEncodingWithoutRedandancyIC(nbrGeneration,outputFileIC,outputFilePopIC,outputOptimalSolutionFileIC,
                                                    populationIC1,populationIC2, solutionDominanteIC,iteration,&bestSolutionIterationDC, &nbrApparitionDC);
#endif // IC
#if CGE
                ecartEncoding(nbrGeneration,outputFileCGE,outputFilePopCGE,outputOptimalSolutionFileCGE,populationCGE1,populationCGE2,
                              solutionDominanteCGE,iteration,&bestSolutionIterationCGE, &nbrApparitionCGE);
#endif // CGE
#if PMEB
                pMedianEncodingAffectationParArretes(nbrGeneration,outputFilePMPEdgeBased,outputFilePopPMPEdgeBased,outputOptimalSolutionFilePMPEdgeBased,
                                                     populationPMEB1,populationPMEB2, solutionDominantePMEB,iteration,&bestSolutionIterationPMEB, &nbrApparitionPMEB);
#endif // PMEB
#if PMCA
                pMedianEncodingBaseSurClusterSize(nbrGeneration,outputFilePMPClusterBased,outputFilePopPMPClusterBased,outputOptimalSolutionFilePMPClusterBased,
                                                  populationPMCA1,populationPMCA2, solutionDominantePMCA,iteration,&bestSolutionIterationPMCA, &nbrApparitionPMCA);
#endif // PMCA
                /**
                après chaque itération je compare les résultats des codages (solutionDominante) pour inscrire la meilleur
                dans un fichier que j'appelle bestSolutionParallelGA / optimalSolutionParallelGA en précisant le codage qui
                a abouti à cette meilleur solution
                */

                writeBestSolutionInFilePGA(solutionDominanteFVTC, solutionDominanteEA,solutionDominanteDC, solutionDominantePMEB, solutionDominantePMCA,outputFilePGA);
            }
            fprintf(outputFilePGA,"\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");

            t2=clock();
            temps = (double)(t2-t1)/CLOCKS_PER_SEC;

#if DVTC
            ///displayTheBestSolutionVTC(solutionDominanteVTC);
            fprintf(outputFileVTC,"\n\n");
            writeBestSolutionInFileVTC(solutionDominanteVTC,outputFileVTC,bestSolutionIterationVTC);
            fprintf(outputFileVTC,"\n\nLe nombre d'apparition de la meilleur solution est  = %d\n",nbrApparitionVTC);
            fprintf(outputFileVTC,"\n\n==============================================================================\n");
            fprintf(outputFilePopVTC,"\n\n==============================================================================\n");

            if(bestSolutionIterationVTC>=2)
            {
                ES_VTC = ((taillePopulation-tauxReproduction)*(bestSolutionIterationVTC-2))+(taillePopulation+ solutionDominanteVTC->id -tauxReproduction+1);
            }
            else
            {
                ES_VTC = solutionDominanteVTC->id +1;
            }
            writeOptimalSolutionInFileVTC(solutionDominanteVTC,outputOptimalSolutionFileVTC,nbrRun,bestSolutionIterationVTC,temps, ES_VTC);

#endif // DVTC

#if activatePause
            mon_sleep(pause);
#endif // activatePause

#if SVTC
            ///displayTheBestSolutionSVTC(solutionDominante);
            fprintf(outputFileSVTC,"\n\n");
            writeBestSolutionInFileSVTC(solutionDominanteSVTC,outputFileSVTC,bestSolutionIterationSVTC);
            fprintf(outputFileSVTC,"\n\nLe nombre d'apparition de la meilleur solution est  = %d\n",nbrApparitionSVTC);
            fprintf(outputFileSVTC,"\n\n==============================================================================\n");
            fprintf(outputFilePopSVTC,"\n\n==============================================================================\n");
            if(bestSolutionIterationSVTC>=2)
            {
                ES_SVTC = ((taillePopulation-tauxReproduction)*(bestSolutionIterationSVTC-2))+(taillePopulation+ solutionDominanteSVTC->id -tauxReproduction+1);
            }
            else
            {
                ES_SVTC = solutionDominanteSVTC->id +1;
            }

            writeOptimalSolutionInFileSVTC(solutionDominanteSVTC,outputOptimalSolutionFileSVTC,nbrRun,bestSolutionIterationSVTC,temps, ES_SVTC);

#endif // SVTC

#if activatePause
            mon_sleep(pause);
#endif // activatePause

#if FVTC
            ///displayTheBestSolutionFVTC(solutionDominante);
            fprintf(outputFileFVTC,"\n\n");
            writeBestSolutionInFileFVTC(solutionDominanteFVTC,outputFileFVTC,bestSolutionIterationFVTC);
            fprintf(outputFileFVTC,"\n\nLe nombre d'apparition de la meilleur solution est  = %d\n",nbrApparitionFVTC);
            fprintf(outputFileFVTC,"\n\n==============================================================================\n");
            fprintf(outputFilePopFVTC,"\n\n==============================================================================\n");

            if(bestSolutionIterationFVTC>=2)
            {
                ES_FVTC = ((taillePopulation-tauxReproduction)*(bestSolutionIterationFVTC-2))+(taillePopulation+ solutionDominanteFVTC->id -tauxReproduction+1);
            }
            else
            {
                ES_FVTC = solutionDominanteFVTC->id +1;
            }
            writeOptimalSolutionInFileFVTC(solutionDominanteFVTC,outputOptimalSolutionFileFVTC,nbrRun,bestSolutionIterationFVTC,temps,ES_FVTC);

#endif // FVTC

#if activatePause
            mon_sleep(pause);
#endif // activatePause

#if VAE

            ///displayTheBestSolutionVAE(solutionDominante);

            if(bestSolutionIterationVAE>=2)
            {
                ES_VAE = ((taillePopulation-tauxReproduction)*(bestSolutionIterationVAE-2))+(taillePopulation+ solutionDominanteVAE->id -tauxReproduction+1);
            }
            else
            {
                ES_VAE = solutionDominanteVAE->id +1;
            }
            writeOptimalSolutionInFileVAE(solutionDominanteVAE,outputOptimalSolutionFileVAE,nbrRun,bestSolutionIterationVAE,temps,ES_VAE);
            fprintf(outputFileVAE,"\n\n");
            writeBestSolutionInFileVAE(solutionDominanteVAE,outputFileVAE,bestSolutionIterationVAE);
            fprintf(outputFileVAE,"\n\nLe nombre d'apparition de la meilleur solution est  = %d\n",nbrApparitionVAE);
            fprintf(outputFileVAE,"\n\n==============================================================================\n");
            fprintf(outputFilePopVAE,"\n\n==============================================================================\n");

#endif // VAE

#if activatePause
            mon_sleep(pause);
#endif // activatePause

#if PMEB
            ///displayTheBestSolutionPMP(solutionDominante);

            if(bestSolutionIterationPMEB>=2)
            {
                ES_PMEB = ((taillePopulation-tauxReproduction)*(bestSolutionIterationPMEB-2))+(taillePopulation+ solutionDominantePMEB->id -tauxReproduction+1);
            }
            else
            {
                ES_PMEB = solutionDominantePMEB->id +1;
            }

            writeOptimalSolutionInFilePMP(solutionDominantePMEB,outputOptimalSolutionFilePMPEdgeBased,nbrRun, bestSolutionIterationPMEB, temps,ES_PMEB);
            fprintf(outputFilePMPEdgeBased,"\n\n");
            writeBestSolutionInFilePMP(solutionDominantePMEB,outputFilePMPEdgeBased,bestSolutionIterationPMEB);
            fprintf(outputFilePMPEdgeBased,"\n\nLe nombre d'apparition de la meilleur solution est  = %d\n",nbrApparitionPMEB);
            fprintf(outputFilePMPEdgeBased,"\n\n==============================================================================\n");
            fprintf(outputFilePopPMPEdgeBased,"\n\n==============================================================================\n");


#endif // PMEB

#if activatePause
            mon_sleep(pause);
#endif // activatePause

#if PMCA

            if(bestSolutionIterationPMCA>=2)
            {
                ES_PMCA = ((taillePopulation-tauxReproduction)*(bestSolutionIterationPMCA-2))+(taillePopulation+ solutionDominantePMCA->id -tauxReproduction+1);
            }
            else
            {
                ES_PMCA = solutionDominantePMCA->id +1;
            }
            writeOptimalSolutionInFilePMP(solutionDominantePMCA,outputOptimalSolutionFilePMPClusterBased,nbrRun, bestSolutionIterationPMCA, temps,ES_PMCA);

            fprintf(outputFilePMPClusterBased,"\n\n");

            writeBestSolutionInFilePMP(solutionDominantePMCA,outputFilePMPClusterBased,bestSolutionIterationPMCA);
            fprintf(outputFilePMPClusterBased,"\n\nLe nombre d'apparition de la meilleur solution est  = %d\n",nbrApparitionPMCA);
            fprintf(outputFilePMPClusterBased,"\n\n==============================================================================\n");
            fprintf(outputFilePopPMPClusterBased,"\n\n==============================================================================\n");

#endif // PMCA

#if activatePause
            mon_sleep(pause);
#endif // activatePause

#if EA
            if(bestSolutionIterationEA>=2)
            {
                ES_EA = ((taillePopulation-tauxReproduction)*(bestSolutionIterationEA-2))+(taillePopulation+ solutionDominanteEA->id -tauxReproduction+1);
            }
            else
            {
                ES_EA = solutionDominanteEA->id +1;
            }

            writeOptimalSolutionInFileEA(solutionDominanteEA,outputOptimalSolutionFileEA,nbrRun, bestSolutionIterationEA,temps,ES_EA);
            fprintf(outputFileEA,"\n\n");
            writeBestSolutionInFileEA(solutionDominanteEA,outputFileEA,bestSolutionIterationEA);
            fprintf(outputFileEA,"\n\nLe nombre d'apparition de la meilleur solution est  = %d\n",nbrApparitionEA);

            fprintf(outputFileEA,"\n\n==============================================================================\n");
            fprintf(outputFilePopEA,"\n\n==============================================================================\n");

#endif // EA
#if activatePause
            mon_sleep(pause);
#endif // activatePause

#if DC

            ///displayTheBestSolutionDC(solutionDominante);

            if(bestSolutionIterationDC>=2)
            {
                ES_DC = ((taillePopulation-tauxReproduction)*(bestSolutionIterationDC-2))+(taillePopulation+ solutionDominanteDC->id -tauxReproduction+1);
            }
            else
            {
                ES_DC = solutionDominanteDC->id +1;
            }

            writeOptimalSolutionInFileDC(solutionDominanteDC,outputOptimalSolutionFileDC,nbrRun,bestSolutionIterationDC,temps,ES_DC);

            fprintf(outputFileDC,"\n\n");
            writeBestSolutionInFileDC(solutionDominanteDC,outputFileDC,bestSolutionIterationDC);
            fprintf(outputFileDC,"\n\nLe nombre d'apparition de la meilleur solution est  = %d\n",nbrApparitionDC);
            fprintf(outputFileDC,"\n\n==============================================================================\n");
            fprintf(outputFilePopDC,"\n\n==============================================================================\n");

#endif // DC

#if activatePause
            mon_sleep(pause);
#endif // activatePause

#if IC

            if(bestSolutionIterationIC>=2)
            {
                ES_IC = ((taillePopulation-tauxReproduction)*(bestSolutionIterationIC-2))+(taillePopulation+ solutionDominanteIC->id -tauxReproduction+1);
            }
            else
            {
                ES_IC = solutionDominanteIC->id +1;
            }

            writeOptimalSolutionInFileIC(solutionDominanteIC,outputOptimalSolutionFileIC,nbrRun,bestSolutionIterationIC,temps,ES_IC);

            fprintf(outputFileIC,"\n\n");
            writeBestSolutionInFileIC(solutionDominanteIC,outputFileIC,bestSolutionIterationIC);
            fprintf(outputFileIC,"\n\nLe nombre d'apparition de la meilleur solution est  = %d\n",nbrApparitionIC);
            fprintf(outputFileIC,"\n\n==============================================================================\n");
            fprintf(outputFilePopIC,"\n\n==============================================================================\n");

#endif // IC

#if activatePause
            mon_sleep(pause);
#endif // activatePause
#if CGE
            ///displayTheBestSolutionCGE(solutionDominante);

            if(bestSolutionIterationCGE>=2)
            {
                ES_CGE = ((taillePopulation-tauxReproduction)*(bestSolutionIterationCGE-2))+(taillePopulation+ solutionDominanteCGE->id -tauxReproduction+1);
            }
            else
            {
                ES_CGE = solutionDominanteCGE->id +1;
            }

            ///writeOptimalSolutionInFileCGE(solutionDominante,outputOptimalSolutionFileCGE,nbrRun,bestSolutionIteration,temps,ES);
            writeOptimalSolutionInFileCGE((populationCGE1+bestSolutionCGE),outputOptimalSolutionFileCGE,nbrRun,bestSolutionIterationCGE,temps,ES_CGE);
            fprintf(outputFileCGE,"\n\n");
            ///writeBestSolutionInFileCGE(solutionDominante,outputFileCGE,bestSolutionIteration);
            writeBestSolutionInFileCGE((populationCGE1+bestSolutionCGE),outputFileCGE,bestSolutionIterationCGE);
            fprintf(outputFileCGE,"\n\nLe nombre d'apparition de la meilleur solution est  = %d\n",nbrApparitionCGE);
            fprintf(outputFileCGE,"\n\n==============================================================================\n");
            fprintf(outputFilePopCGE,"\n\n==============================================================================\n");

#endif // CGE

#if activatePause
            mon_sleep(pause);
#endif // activatePause


            ///optimalSolutionPGA
            if(solutionDominanteFVTC->coutCoupeNormalise >= solutionDominanteDC->coutCoupeNormalise &&
                    solutionDominanteFVTC->coutCoupeNormalise >= solutionDominanteEA->coutCoupeNormalise &&
                    solutionDominanteFVTC->coutCoupeNormalise >= solutionDominantePMEB->coutCoupeNormalise &&
                    solutionDominanteFVTC->coutCoupeNormalise >= solutionDominantePMCA->coutCoupeNormalise)
            {

                writeOptimalSolutionInFileFVTC(solutionDominanteFVTC,outputOptimalSolutionFilePGA,nbrRun, bestSolutionIterationFVTC, temps,ES_FVTC);
            }
            else if(solutionDominanteEA->coutCoupeNormalise >= solutionDominanteFVTC->coutCoupeNormalise &&
                    solutionDominanteEA->coutCoupeNormalise >= solutionDominanteDC->coutCoupeNormalise &&
                    solutionDominanteEA->coutCoupeNormalise >= solutionDominantePMEB->coutCoupeNormalise &&
                    solutionDominantePMCA->coutCoupeNormalise >= solutionDominantePMCA->coutCoupeNormalise)
            {

                writeOptimalSolutionInFileEA(solutionDominanteEA,outputOptimalSolutionFilePGA,nbrRun, bestSolutionIterationEA, temps,ES_EA);
            }
            else if(solutionDominanteDC->coutCoupeNormalise >= solutionDominanteFVTC->coutCoupeNormalise &&
                    solutionDominanteDC->coutCoupeNormalise >= solutionDominanteEA->coutCoupeNormalise &&
                    solutionDominanteDC->coutCoupeNormalise >= solutionDominantePMEB->coutCoupeNormalise &&
                    solutionDominantePMCA->coutCoupeNormalise >= solutionDominantePMCA->coutCoupeNormalise)
            {

                writeOptimalSolutionInFileDC(solutionDominanteDC,outputOptimalSolutionFilePGA,nbrRun, bestSolutionIterationDC, temps,ES_DC);

            }
            else if(solutionDominantePMEB->coutCoupeNormalise >= solutionDominanteFVTC->coutCoupeNormalise &&
                    solutionDominantePMEB->coutCoupeNormalise >= solutionDominanteEA->coutCoupeNormalise &&
                    solutionDominantePMEB->coutCoupeNormalise >= solutionDominanteDC->coutCoupeNormalise &&
                    solutionDominantePMCA->coutCoupeNormalise >= solutionDominantePMCA->coutCoupeNormalise)
            {

                writeOptimalSolutionInFilePMP(solutionDominantePMEB,outputOptimalSolutionFilePGA,nbrRun, bestSolutionIterationPMEB, temps,ES_PMEB);
            }
            else if(
                solutionDominantePMCA->coutCoupeNormalise >= solutionDominanteFVTC->coutCoupeNormalise &&
                solutionDominantePMCA->coutCoupeNormalise >= solutionDominanteEA->coutCoupeNormalise &&
                solutionDominantePMCA->coutCoupeNormalise >= solutionDominanteDC->coutCoupeNormalise &&
                solutionDominantePMCA->coutCoupeNormalise >= solutionDominantePMCA->coutCoupeNormalise)
            {

                writeOptimalSolutionInFilePMP(solutionDominantePMCA,outputOptimalSolutionFilePGA,nbrRun, bestSolutionIterationPMCA, temps,ES_PMCA);
            }

            nbrRun++;
        }

        fclose(inputFile);
        fclose(outputFileCocyclesDeBase);
        fclose(inputFileConstrainte);
#if DVTC
        fclose(outputFileVTC);
        fclose(outputFilePopVTC);
        fclose(outputOptimalSolutionFileVTC);
#endif // DVTC
#if SVTC
        fclose(outputFileSVTC);
        fclose(outputFilePopSVTC);
        fclose(outputOptimalSolutionFileSVTC);
#endif // SVTC
#if FVTC
        fclose(outputFileFVTC);
        fclose(outputFilePopFVTC);
        fclose(outputOptimalSolutionFileFVTC);
#endif // FVTC
#if VAE
        fclose(outputFileVAE);
        fclose(outputFilePopVAE);
        fclose(outputOptimalSolutionFileVAE);
#endif // VAE
#if PMEB
        fclose(outputFilePMPEdgeBased);
        fclose(outputFilePopPMPEdgeBased);
        fclose(outputOptimalSolutionFilePMPEdgeBased);
#endif // PMEB
#if PMCA
        fclose(outputFilePMPClusterBased);
        fclose(outputFilePopPMPClusterBased);
        fclose(outputOptimalSolutionFilePMPClusterBased);
#endif // PMCA
#if EA
        fclose(outputFileEA);
        fclose(outputFilePopEA);
        fclose(outputOptimalSolutionFileEA);
#endif // EA
#if DC
        fclose(outputFileDC);
        fclose(outputFilePopDC);
        fclose(outputOptimalSolutionFileDC);
#endif // DC
#if IC
        fclose(outputFileIC);
        fclose(outputFilePopIC);
        fclose(outputOptimalSolutionFileIC);
#endif // IC
#if CGE
        fclose(outputFileCGE);
        fclose(outputFilePopCGE);
        fclose(outputOptimalSolutionFileCGE);
#endif // CGE
        nbrRun = 1;
    }

#define libre 0
#if libre
///=======================================================================================================
    if(inputFile!=NULL)free(inputFile);

    printf("inputFile = %P\n",inputFile);

    if(outputFileCocyclesDeBase!=NULL)free(outputFileCocyclesDeBase);
    if(inputFileConstrainte!=NULL)free(inputFileConstrainte);
///=======================================================================================================
    fclose(outputFileVAE);
    if(outputFileVAE!=NULL)free(outputFileVAE);
    fclose(outputFilePopVAE);
    if(outputFilePopVAE!=NULL)free(outputFilePopVAE);
    fclose(outputOptimalSolutionFileVAE);
    if(outputOptimalSolutionFileVAE!=NULL)free(outputOptimalSolutionFileVAE);

    if(populationVAE1!=NULL)free(populationVAE1);
    printf("populationVAE1 = %p\n",populationVAE1);
    if(populationVAE2!=NULL)free(populationVAE2);
    if(solutionDominanteVAE!=NULL)free(solutionDominanteVAE);
///=======================================================================================================

    fclose(outputFileFC);
    if(outputFileFC!=NULL)free(outputFileFC);
    fclose(outputFilePopFC);
    if(outputFilePopFC!=NULL)free(outputFilePopFC);
    fclose(outputOptimalSolutionFileFC);
    if(outputOptimalSolutionFileFC!=NULL)free(outputOptimalSolutionFileFC);

    if(populationFC1!=NULL)free(populationFC1);
    if(populationFC2!=NULL)free(populationFC2);
    if(solutionDominanteFC!=NULL)free(solutionDominanteFC);

///=======================================================================================================
    fclose(outputFileEA);
    if(outputFileEA!=NULL)free(outputFileEA);
    fclose(outputFilePopEA);
    if(outputFilePopEA!=NULL)free(outputFilePopEA);
    fclose(outputOptimalSolutionFileEA);
    if(outputOptimalSolutionFileEA!=NULL)free(outputOptimalSolutionFileEA);

    if(populationEA1!=NULL)free(populationEA1);
    if(populationEA2!=NULL)free(populationEA2);
    if(solutionDominanteEA!=NULL)free(solutionDominanteEA);

///=======================================================================================================

    fclose(outputFileDC);
    if(outputFileDC!=NULL)free(outputFileDC);
    fclose(outputFilePopDC);
    if(outputFilePopDC!=NULL)free(outputFilePopDC);
    fclose(outputOptimalSolutionFileDC);
    if(outputOptimalSolutionFileDC!=NULL)free(outputOptimalSolutionFileDC);

    if(populationDC1!=NULL)free(populationDC1);
    if(populationDC2!=NULL)free(populationDC2);
    if(solutionDominanteDC!=NULL)free(solutionDominanteDC);
///=======================================================================================================

    fclose(outputFileIC);
    if(outputFileIC!=NULL)free(outputFileIC);
    fclose(outputFilePopIC);
    if(outputFilePopIC!=NULL)free(outputFilePopIC);
    fclose(outputOptimalSolutionFileIC);
    if(outputOptimalSolutionFileIC!=NULL)free(outputOptimalSolutionFileIC);

    if(populationIC1!=NULL)free(populationIC1);
    if(populationIC2!=NULL)free(populationIC2);
    if(solutionDominanteIC!=NULL)free(solutionDominanteIC);
///=======================================================================================================

    fclose(outputFilePMPEdgeBased);
    if(outputFilePMPEdgeBased!=NULL)free(outputFilePMPEdgeBased);
    fclose(outputFilePopPMPEdgeBased);
    if(outputFilePopPMPEdgeBased!=NULL)free(outputFilePopPMPEdgeBased);
    fclose(outputOptimalSolutionFilePMPEdgeBased);
    if(outputOptimalSolutionFilePMPEdgeBased!=NULL)free(outputOptimalSolutionFilePMPEdgeBased);

    fclose(outputFilePMPClusterBased);
    if(outputFilePMPClusterBased!=NULL)free(outputFilePMPClusterBased);
    fclose(outputFilePopPMPClusterBased);
    if(outputFilePopPMPClusterBased!=NULL)free(outputFilePopPMPClusterBased);
    fclose(outputOptimalSolutionFilePMPClusterBased);
    if(outputOptimalSolutionFilePMPClusterBased!=NULL)free(outputOptimalSolutionFilePMPClusterBased);

    if(populationPMP1!=NULL)free(populationPMP1);
    if(populationPMP2!=NULL)free(populationPMP2);
    if(solutionDominantePMP!=NULL)free(solutionDominantePMP);
///=======================================================================================================
    fclose(outputFileVTC);
    if(outputFileVTC!=NULL)free(outputFileVTC);
    fclose(outputFilePopVTC);
    if(outputFilePopVTC!=NULL)free(outputFilePopVTC);
    fclose(outputOptimalSolutionFileVTC);
    if(outputOptimalSolutionFileVTC!=NULL)free(outputOptimalSolutionFileVTC);

    if(populationVTC1!=NULL)free(populationVTC1);
    if(populationVTC2!=NULL)free(populationVTC2);
    if(solutionDominanteVTC!=NULL)free(solutionDominanteVTC);
///=======================================================================================================
    fclose(outputFileCGE);
    if(outputFileCGE!=NULL)free(outputFileCGE);
    fclose(outputFilePopCGE);
    if(outputFilePopCGE!=NULL)free(outputFilePopCGE);
    fclose(outputOptimalSolutionFileCGE);
    if(outputOptimalSolutionFileCGE!=NULL)free(outputOptimalSolutionFileCGE);

    if(populationCGE1!=NULL)free(populationCGE1);
    if(populationCGE2!=NULL)free(populationCGE2);
    if(solutionDominanteCGE!=NULL)free(solutionDominanteCGE);

///=======================================================================================================

#endif // libre
    return 0;
}

