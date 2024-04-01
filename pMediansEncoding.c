#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "typeDeclaration.h"
#include "./communesFunctions.h"
#include "./pMediansEncoding.h"
#include "./compilationConditionnelle.h"
#define affectationArretes 0

extern int taillePopulation,tauxReproduction,sommeTotalFlux,nbrNoeuds,
       nbr_constraint,nbrArretes,nbrParties,nbrNoeudsInCluster, regeneration;
extern double tauxMutation ,averageTraffic;
extern int fluxMatrix[1000][1000],*fluxVector, maximumDistanceMatrix[1000][1000],shortestPath[1000][1000], maxFlow, sommeTotalFluxReal;
extern edge *edgeVector,cohabitationConstraintes[1000],nonCohabitationConstraintes[100];
extern int max_clusters , min_clusters,min_sizeCluster ,max_sizeCluster ;
extern int bestSolution, nbrApparition, nbrCohabitationConstraintes,nbrNonCohabitationConstraintes;
extern int nbrRun;


///***************************************************************************
///                     les corps des fonctions
///***************************************************************************

///****************************************************************************************************************************
///génération de la populationPMP initiale
void generatePopulationPMP(partitionPMP* populationPMP, int indiceFirstElt)
{

    int i,j,indice;
    /// génération des zero vertices
    for(i=indiceFirstElt; i< taillePopulation; i++)
    {
        (populationPMP+i)->id =i;
        /// affecter la valeur 0 à tous les noeuds de la solution
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationPMP+i)->genotype[j] = 0;
        }
        ///(populationPMP+i)->medians =2;
        (populationPMP+i)->medians =rnd(0,nbrNoeuds);
        ///(populationPMP+i)->medians =rnd(0,max_clusters); // laisser le nombre des medians libre donne des resumtats meilleurs que de le limité
        ///printf("le nombre de medains dans cette solutions : %d \n",(populationPMP+i)->medians);
        /// choisir aléatoirement les noeuds medians
        for(j=0; j<(populationPMP+i)->medians; j++) /// <= pour qu'on aura un élément en plus pour le nombre des parties
        {
            do
            {
                indice = rnd(0,nbrNoeuds);
            }
            while((populationPMP+i)->genotype[indice] == 1);
            (populationPMP+i)->genotype[indice] = 1;
        }
        (populationPMP+i)->nonMedians = nbrNoeuds - (populationPMP+i)->medians;
    }

}
///****************************************************************************************************************************
/// affichage des population
void affichePopulationPMP(partitionPMP* populationPMP)
{
    int i,j;
    for(i=0; i<taillePopulation; i++)
    {
        if((populationPMP+i)->contrainteViole == 0)
        {
            printf("\nid = %d \t", (populationPMP+i)->id);
            printf("la genotypes est  = \t");
            for(j=0; j<nbrNoeuds; j++)
            {
                printf("%d ",(populationPMP+i)->genotype[j]);
            }
            printf("\n");
            /// affichage des medians de cette solution
            printf("  les medians sont  = \t");
            for(j=0; j<(populationPMP+i)->medians; j++)
            {
                printf("%d ",(populationPMP+i)->mediansVertices[j]);
            }
            printf("\n");

            printf("  les NOT  medians sont  = \t");
            for(j=0; j<(populationPMP+i)->nonMedians; j++)
            {
                printf("%d ",(populationPMP+i)->nonMediansVertices[j]);
            }
            printf("\n");

            printf("\n la phenotypes est  = \t");
            for(j=0; j<nbrNoeuds; j++)
            {
                printf("%d ",(populationPMP+i)->phenotype[j]);
            }
            printf("\n");
            printf("la fitness de la %d solution est : %0.2f",i,(populationPMP+i)->fitness );
            printf("\n**************************************************************************\n");
        }
        else if((populationPMP+i)->contrainteViole == 1)
        {
            printf("la %d solution a viole la contrainte de capacite \n",i);
        }
    }

}
///****************************************************************************************************************************
void calculerMediansNonMediansVertices(partitionPMP* populationPMP)
{

    int i,j,testValueGenotype;
#if mouchard
    printf("calculerMediansNonMediansVertices : ...\n");
#endif // mouchard
    for(i=0; i<taillePopulation; i++)
    {
        (populationPMP+i)->medians = 0;
        (populationPMP+i)->nonMedians = 0;
        for(j=0; j<nbrNoeuds; j++)
        {
            testValueGenotype = (populationPMP+i)->genotype[j];
            switch (testValueGenotype)
            {
            case 0 :
                (populationPMP+i)->nonMediansVertices[(populationPMP+i)->nonMedians] = j;
                (populationPMP+i)->nonMedians++;
                break;
            case 1 :
                (populationPMP+i)->mediansVertices[(populationPMP+i)->medians] = j;
                (populationPMP+i)->medians++;
                break;
            default :
                printf("Cette valeur n\'est pas prise en consédiration \n");
                break;
            }
        }
    }
}
///****************************************************************************************************************************
void calculcoutCoupeeEtFitnessPMP(partitionPMP* populationPMP)
{
    int i,j,k,sommeTotalCoutCoupe = 0;

    for(i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (populationPMP+i)->coutCoupe = 0; /// réinitialisation de la valriable après chaque itération (individus de la population)
        for(j=0; j<nbrNoeuds-1; j++)
        {
            for(k=j+1; k<nbrNoeuds; k++)
            {
                if((populationPMP+i)->phenotype[j] == (populationPMP+i)->phenotype[k] && fluxMatrix[j][k] > 0)
                {
                    (populationPMP+i)->coutCoupe = (populationPMP+i)->coutCoupe +  fluxMatrix[j][k];
                }
            }
        }
        (populationPMP+i)->coutCoupeNormalise = (populationPMP+i)->coutCoupe + ((nbr_constraint - (populationPMP+i)->contrainteViole)*sommeTotalFlux);
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationPMP+i)->coutCoupeNormalise;
    }

#if scaling

    ///======================================================================================================
    /// Scaling Fintnss : sigma scaling Melanie mitchelle référence Goldberg | le 26/11/2015
    ///pour régler le problème des valeurs négatives des expected values ==> utiliser sigma truncation de goldberg
    /// g(f) = f + (moyenne - c * sigma )
    /// 1 <= c <= 5, généralement initilisé à 2
    ///======================================================================================================
    float varianceCoutDeCoupeNormalise = 0, moyenneCoutDeCoupeNormalise, sigma, sommeTotalOfExpectdValues = 0;
    int c= 2;

    moyenneCoutDeCoupeNormalise = sommeTotalCoutCoupe / taillePopulation;
    for(i=0; i<taillePopulation; i++)
    {
        varianceCoutDeCoupeNormalise =
            varianceCoutDeCoupeNormalise + pow(((populationPMP+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);
    /// tester la valeur ajouté
    for(i=0; i<taillePopulation; i++)
    {

        /// calcul de sigma truncation
        (populationPMP+i)->expectedValue = (populationPMP+i)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationPMP+i)->expectedValue;
    }


    /// calcule du fitness
    /// le 26/11/2015 : calcul du fitness des solution avec leur expectedValues
    for(i=0; i<taillePopulation ; i++)
    {
        (populationPMP+i)->fitness = (float)((populationPMP+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);
    }
#else

    FILE *file;
    file = fopen("D:/WorkingFolderAgProgramming/benchmark/pMedianEncoding/expectedValuePMP.txt","w");
    for(i=0; i<taillePopulation ; i++)
    {
        (populationPMP+i)->fitness = (float)((populationPMP+i)->coutCoupeNormalise)/(float)(sommeTotalCoutCoupe);
    }
#endif// scaling

}

///****************************************************************************************************************************
/// Cette fonction assure la séléction des individu pour la pprocedure d'appariement
/// elle émite la roue de lottrie biaisé
void naturalSelectionPMP(partitionPMP* populationPMP1,partitionPMP* populationPMP2)
{
    int i,j,maxFitness;
    int tmpVectIndice[taillePopulation];
    for(i=0;i<taillePopulation;i++) tmpVectIndice[i]=-1;
    float lotterie,sommeFitness;

    for(i=0; i<tauxReproduction; i++){
        maxFitness = 0;
        for(j=0;j<taillePopulation;j++){
            if((populationPMP1+maxFitness)->coutCoupeNormalise < (populationPMP1+j)->coutCoupeNormalise
               && tmpVectIndice[j] == -1 ){
                maxFitness = j;

            }
        }
        tmpVectIndice[maxFitness] = 0;
        *(populationPMP2+i) = *(populationPMP1+maxFitness);
        ///printf("(populationPMP2+%d)->coutCoupeNormalise = %d \n",i,(populationPMP2+i)->coutCoupeNormalise);
        ///mon_sleep(pause);
    }
    ///printf("\n\n\n");

    ///qsort (populationPMP2, taillePopulation, sizeof *populationPMP2, compareCroissantFitnessEE);
    ///sortingPopulationByFintness(populationPMP2,tmpPopulation);

    for(i=tauxReproduction; i<taillePopulation; i++)
    {
#if Windows
        lotterie= drand48ForWindows(0,1001);
#else
        lotterie= drand48();
#endif // Windows
        sommeFitness = 0;

        for(j=0; j<taillePopulation; j++)
        {
            sommeFitness = sommeFitness + (populationPMP1+j)->fitness;
            if (lotterie<=sommeFitness)
            {
                *(populationPMP2+i) = *(populationPMP1+j);
                break;
            }
            else if( j==taillePopulation-1)
            {
                *(populationPMP2+i) = *(populationPMP1+j);
            }
        }
    }
    ///free(tmpPopulation);

}
///****************************************************************************************************************************
void crossOverPMP(partitionPMP* populationPMP1, partitionPMP* populationPMP2)
{
    int i = 0,j;
    int choixInd1, choixInd2, choixLocus;
    ///***********************************************************************************
    for(i=0; i<tauxReproduction; i++)
    {
        *(populationPMP2+i) = *(populationPMP1+i);
        (populationPMP2+i)->id = i;
    }
    ///***********************************************************************************
#if NEW_CROSSOVER
    while (i < taillePopulation-regeneration)
    {
        ///choixInd1 = rnd(tauxReproduction,taillePopulation);
        choixInd1 = rnd(0,taillePopulation);
        do
        {
            /// ce traitement évite l'appariement d'un element avec lui même/
            ///choixInd2 = rnd(tauxReproduction,taillePopulation);
            choixInd2 = rnd(0,taillePopulation);
        }
        while(choixInd2 == choixInd1 );

        choixLocus = rnd(0,nbrNoeuds-1); /// -2 parce qu'on a commencer à partir du 0

        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationPMP2+i)->genotype[j] = (populationPMP1+choixInd1)->genotype[j];
            (populationPMP2+i+1)->genotype[j] = (populationPMP1+choixInd2)->genotype[j];

        }
        ///*************************************************************************************************
        for(j=choixLocus+1; j<nbrNoeuds; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationPMP2+i)->genotype[j] = (populationPMP1+choixInd2)->genotype[j];
            (populationPMP2+i+1)->genotype[j] = (populationPMP1+choixInd1)->genotype[j];
        }
        /// affectation de nouveau indice
        (populationPMP2+i)->id = i;
        (populationPMP2+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }
    generatePopulationPMP(populationPMP2, taillePopulation-regeneration );
#else
    while (i < taillePopulation)
    {
        ///choixInd1 = rnd(tauxReproduction,taillePopulation);
        choixInd1 = rnd(0,taillePopulation);
        do
        {
            /// ce traitement évite l'appariement d'un element avec lui même/
            /// choixInd2 = rnd(tauxReproduction,taillePopulation);
            choixInd2 = rnd(0,taillePopulation);
        }
        while(choixInd2 == choixInd1 );

        choixLocus = rnd(0,nbrNoeuds-1); /// -2 parce qu'on a commencer à partir du 0

        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationPMP2+i)->genotype[j] = (populationPMP1+choixInd1)->genotype[j];
            (populationPMP2+i+1)->genotype[j] = (populationPMP1+choixInd2)->genotype[j];

        }
        ///*************************************************************************************************
        for(j=choixLocus+1; j<nbrNoeuds; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationPMP2+i)->genotype[j] = (populationPMP1+choixInd2)->genotype[j];
            (populationPMP2+i+1)->genotype[j] = (populationPMP1+choixInd1)->genotype[j];
        }
        /// affectation de nouveau indice
        (populationPMP2+i)->id = i;
        (populationPMP2+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }

#endif
}
///****************************************************************************************************************************
int findTheBestSolutionPMP(partitionPMP *populationPMP)
{
    float maxFitness = populationPMP->fitness;
    int i,indice=0;
    for(i=0; i<taillePopulation ; i++)
    {
        if(maxFitness < (populationPMP+i)->fitness)
        {
            maxFitness = (populationPMP+i)->fitness;
            indice = i;
        }
    }
#if mouchard
    printf("\n maxFitness = %0.4f \t indice = %d \n",maxFitness, indice);
#endif // mouchard
    return indice;
}
///****************************************************************************************************************************
void mutationPMP(partitionPMP* populationPMP)
{

    int i,numeroGeneMute;
    double applicationOfMutation;
    for(i=tauxReproduction; i<taillePopulation; i++)
    {
#if Windows
        applicationOfMutation= drand48ForWindows(0,1001);
        ///printf("applicationOfMutation = %0.4f \n",applicationOfMutation);
        ///mon_sleep(pause);
#else
        applicationOfMutation= drand48();
#endif // Windows
        if(applicationOfMutation <= tauxMutation)
        {
            /// choisir un gène aléatoirement pour lui appliquer la mutation
            numeroGeneMute = rnd(0,nbrNoeuds);
            ///((populationPMP+i)->genotype[numeroGeneMute] == 1)? 0 : 1;
            if((populationPMP+i)->genotype[numeroGeneMute] == 1)
            {
                (populationPMP+i)->genotype[numeroGeneMute] =0;
            }
            else
            {
                (populationPMP+i)->genotype[numeroGeneMute] =1;
            }

        }
    }
}
///****************************************************************************************************************************
float testerLaSommeDesFitnessPMP(partitionPMP* populationPMP)
{
    printf("testerLaSommeDesFitness ...\n");
    int i;
    float sommeFitness = 0;
    for(i=0; i<taillePopulation ; i++)
    {
        sommeFitness = sommeFitness + (populationPMP+i)->fitness;
    }
    return sommeFitness;
}
///****************************************************************************************************************************
void displayTheBestSolutionPMP(partitionPMP* solutionDominante)
{
    int i;
    printf("id = %d | normalised cut size = %d | cut size intrA = %d cut size intEr = %d ",
           solutionDominante->id,solutionDominante->coutCoupeNormalise, solutionDominante->coutCoupe, (sommeTotalFlux-solutionDominante->coutCoupe));
    printf("le contraintes violees = %d \n", solutionDominante->contrainteViole);
    for(i=0; i<4; i++)
    {
        switch (i)
        {

        case 0:
            if(solutionDominante->constraintVector[i] ==0) printf("contrainte sur le nombre de clusters : OK \n");
            else printf("contrainte sur le nombre de clusters : NOT OK \n");
            break;
        case 1:
            if(solutionDominante->constraintVector[i] ==0) printf("contrainte sur la taille des clusters : OK \n");
            else printf("contrainte sur la taille des clusters : NOT OK \n");
            break;
        case 2:
            if(solutionDominante->constraintVector[i] ==0) printf("contrainte de cohabitation : OK \n");
            else printf("contrainte de cohabitation : NOT OK \n");
            break;
        case 3:
            if(solutionDominante->constraintVector[i] ==0) printf("contrainte Non Cohabitation : OK \n");
            else printf("contrainte Non Cohabitation : NOT OK \n");
            break;
        default:
            printf("Cette contrainte n'est pas prise en charge par le système \n");

        }

    }
}

///************************************************************************************************************
void writeSolutionInFilePMP(partitionPMP *populationPMP, FILE *outputFilePop,int iteration)
{

    int i,j;
    for(i=0; i<taillePopulation; i++)
    {

        ///écriture de la meilleur solution dans le fichier output
        fprintf(outputFilePop,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,(populationPMP+i)->id,
                (populationPMP+i)->coutCoupe,(populationPMP+i)->fitness,(populationPMP+i)->coutCoupeNormalise,
                (populationPMP+i)->contrainteViole,(populationPMP+i)->medians);
        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(outputFilePop,"%d ",(populationPMP+i)->genotype[j]);
        }
        fprintf(outputFilePop,"\t\t\t\t");
        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(outputFilePop,"%d ",(populationPMP+i)->phenotype[j]);
        }
        fprintf(outputFilePop,"\n");
    }
    /// fprintf(outputFilePop,"\n===============================================\n");
    fprintf(outputFilePop,"\n\n");

}
///************************************************************************************************************
void writeBestSolutionInFilePMP(partitionPMP *solutionDominante, FILE *outputFile,int iteration)
{
    int i;
    fprintf(outputFile,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,solutionDominante->id,
            solutionDominante->coutCoupe,solutionDominante->fitness,solutionDominante->coutCoupeNormalise,
            solutionDominante->contrainteViole,solutionDominante->medians);
    for(i=0; i<nbrNoeuds; i++)
    {
        fprintf(outputFile,"%d ",solutionDominante->genotype[i]);
    }
    fprintf(outputFile,"\t\t\t\t");
    for(i=0; i<nbrNoeuds; i++)
    {
        fprintf(outputFile,"%d ",solutionDominante->phenotype[i]);
    }
    fprintf(outputFile,"\n");

}
///**************************************************************************************
void checkContrainstAndFitnessPenalizationPMP(partitionPMP *populationPMP)
{

    int i,j,k;

    for(i=0; i<taillePopulation; i++)
    {
        ///(populationPMP+i)->nbrCluster=0; /// le nombre de cluster doit etre initialisé à 1 pour prendre en consédiration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes violées
        (populationPMP+i)->contrainteViole=0;
        for(j=0; j<4; j++)
        {
            (populationPMP+i)->constraintVector[j]=0;
        }
        ///**********************************************************************************


        ///trier le tableau résulant
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationPMP+i)->clustersSize[j] =0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if ((populationPMP+i)->phenotype[k] == j)
                {
                    (populationPMP+i)->clustersSize[j]++;
                }
            }
            /**if ((population+i)->clustersSize[j] !=0)
            {
                (populationPMP+i)->nbrCluster++;
            }*/
        }

        ///if((populationPMP+i)->nbrCluster > max_clusters || (populationPMP+i)->nbrCluster < min_clusters)
        if((populationPMP+i)->medians > max_clusters || (populationPMP+i)->medians < min_clusters)
        {
            (populationPMP+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationPMP+i)->contrainteViole++;
        }

        for(j=0; j<nbrNoeuds; j++)
        {
            if((populationPMP+i)->clustersSize[j]!=0)
            {
                if((populationPMP+i)->clustersSize[j]>max_sizeCluster || (populationPMP+i)->clustersSize[j]< min_sizeCluster)
                {
                    (populationPMP+i)->constraintVector[1]++;

                }
            }

        }

        if((populationPMP+i)->constraintVector[1] != 0)
        {
            (populationPMP+i)->contrainteViole++;
        }
    }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0)
    {
        checkCohabitationAndNonCohabitationConstraintePMP(populationPMP);
    }

}
///************************************************************************************************************
void checkCohabitationAndNonCohabitationConstraintePMP(partitionPMP *populationPMP)
{

    int i,j,noeud1,noeud2;
    for(i=0; i<taillePopulation; i++)
    {
	    /// Cette tache est déjà réalisé au niveau de la fonction checkCohabitationAndNonCohabitationConstrainteBE
		/**
        (populationPMP+i)->constraintVector[2]=0;
        (populationPMP+i)->constraintVector[3]=0;
		*/
        if(nbrCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrCohabitationConstraintes; j++)
            {
                noeud1=cohabitationConstraintes[j].nouedDepart;
                noeud2=cohabitationConstraintes[j].nouedArrive;

                if((populationPMP+i)->phenotype[noeud1]!= (populationPMP+i)->phenotype[noeud2])
                {
                    (populationPMP+i)->constraintVector[2]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationPMP+i)->constraintVector[2]!=0)
            {
                (populationPMP+i)->contrainteViole++;
            }
        }


        if(nbrNonCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrNonCohabitationConstraintes; j++)
            {
                noeud1=nonCohabitationConstraintes[j].nouedDepart;
                noeud2=nonCohabitationConstraintes[j].nouedArrive;

                if((populationPMP+i)->phenotype[noeud1]== (populationPMP+i)->phenotype[noeud2])
                {
                    (populationPMP+i)->constraintVector[3]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationPMP+i)->constraintVector[3]!=0)
            {
                (populationPMP+i)->contrainteViole++;
            }
        }


    }

}

///************************************************************************************************************
void pMedianEncodingAffectationParArretes(int nbrGeneration ,FILE *outputFilePMP,FILE *outputFilePopPMP,FILE *outputOptimalSolutionFilePMP,
                                          partitionPMP *populationPMP1,partitionPMP *populationPMP2, partitionPMP *solutionDominante,int iteration , int *bestSolutionIteration , int *nbrApparition)
{

        naturalSelectionPMP(populationPMP1,populationPMP2);
        crossOverPMP(populationPMP2, populationPMP1);
        mutationPMP(populationPMP1);
        calculerMediansNonMediansVertices(populationPMP1);
        affectationDesNonMediansParArretes(populationPMP1);
        checkContrainstAndFitnessPenalizationPMP(populationPMP1);
        calculcoutCoupeeEtFitnessPMP(populationPMP1);

#if writingPopulationInFile
        writeSolutionInFilePMP(populationPMP1,outputFilePopPMP,iteration);
#endif // writingPopulationInFile
///=====================================================================================================
        /// dans ce qui suit, on va juste localiser la meilleur solution
///=====================================================================================================

        bestSolution=findTheBestSolutionPMP(populationPMP1);
#if writingPopulationInFile
        writeBestSolutionInFilePMP((populationPMP1+bestSolution),outputFilePMP,iteration);
#endif // writingPopulationInFile
        if((populationPMP1+bestSolution)->coutCoupeNormalise > solutionDominante->coutCoupeNormalise)
        {
            *(solutionDominante) = *(populationPMP1+bestSolution);
            nbrApparition=1;
            bestSolutionIteration =iteration;
        }
        else
        {
            nbrApparition++;
        }

}


///************************************************************************************************************
void affectationDesNonMediansParArretes(partitionPMP* populationPMP){

    int i,j,l,k,maxDistance,sommeDistances,indiceNotSeed,indiceSeed;
    clock_t t1,t2;
    double temps;
    int distanceArray[nbrNoeuds],nbrNotSeed,tmpNotSeed[nbrNoeuds],nbrFullCluster;

    for(i=0; i<taillePopulation; i++){
//***********************************************************************************************************
        for(l=0; l<(populationPMP+i)->medians; l++){
            (populationPMP+i)->clustersSize[(populationPMP+i)->mediansVertices[l]] = 1;
        }
//***********************************************************************************************************
        for(l=0; l<(populationPMP+i)->nonMedians; l++){
            tmpNotSeed[l] = 1;
        }
//***********************************************************************************************************
		for(j=0;j<nbrNoeuds;j++){
			(populationPMP+i)->phenotype[j] = j;
		}
//***********************************************************************************************************
        nbrNotSeed = (populationPMP+i)->nonMedians;
        nbrFullCluster =0;

        do{
            for(j=0; j < (populationPMP+i)->nonMedians;j++){
                if(tmpNotSeed[j]>0){
                    indiceNotSeed = (populationPMP+i)->nonMediansVertices[j];
                    for(k=0; k<(populationPMP+i)->medians; k++){
                        indiceSeed = (populationPMP+i)->mediansVertices[k];
                        distanceArray[k] = 0;
                        if( (populationPMP+i)->clustersSize[indiceSeed] < max_sizeCluster){
                            for(l=0; l<nbrNoeuds; l++){
                                if((populationPMP+i)->phenotype[l] == indiceSeed){
                                    if(fluxMatrix[indiceNotSeed][l] >= 0 ){
                                        distanceArray[k] += fluxMatrix[indiceNotSeed][l];
                                    }
                                }
                            }
                        }
                        else distanceArray[k] = -1;
                    }
                    maxDistance = distanceArray[0];
                    indiceSeed=(populationPMP+i)->mediansVertices[0];
                    for(k=1; k<(populationPMP+i)->medians; k++){
                        if(maxDistance < distanceArray[k] ) { maxDistance = distanceArray[k]; indiceSeed = (populationPMP+i)->mediansVertices[k]; }
                    }
                    if(maxDistance != -1){
                        (populationPMP+i)->phenotype[indiceNotSeed] = indiceSeed;
                        (populationPMP+i)->clustersSize[indiceSeed]++;
                        if((populationPMP+i)->clustersSize[indiceSeed] >= max_sizeCluster)nbrFullCluster++;
                        tmpNotSeed[j] = -1;
                        nbrNotSeed--;
                    }
                }
            }
            if(nbrFullCluster ==(populationPMP+i)->medians ) nbrNotSeed = 0;
        }while(nbrNotSeed > 0);
 	}
}


///************************************************************************************************************
void pMedianEncodingBaseSurClusterSize(int nbrGeneration ,FILE *outputFilePMP,FILE *outputFilePopPMP,FILE *outputOptimalSolutionFilePMP,partitionPMP *populationPMP1,
                                       partitionPMP *populationPMP2, partitionPMP *solutionDominante,int iteration , int *bestSolutionIteration , int *nbrApparition)
{

        naturalSelectionPMP(populationPMP1,populationPMP2);
        crossOverPMP(populationPMP2, populationPMP1);
        mutationPMP(populationPMP1);
        calculerMediansNonMediansVertices(populationPMP1);
        affectationDesNonMediansBaseSurClusterSize(populationPMP1);
        checkContrainstAndFitnessPenalizationPMP(populationPMP1);
        calculcoutCoupeeEtFitnessPMP(populationPMP1);

#if writingPopulationInFile
        writeSolutionInFilePMP(populationPMP1,outputFilePopPMP,iteration);
#endif // writingPopulationInFile
///=====================================================================================================
        /// dans ce qui suit, on va juste localiser la meilleur solution
///=====================================================================================================

        bestSolution=findTheBestSolutionPMP(populationPMP1);
#if writingPopulationInFile
        writeBestSolutionInFilePMP((populationPMP1+bestSolution),outputFilePMP,iteration);
#endif // writingPopulationInFile
        if((populationPMP1+bestSolution)->coutCoupeNormalise > solutionDominante->coutCoupeNormalise)
        {
            *(solutionDominante) = *(populationPMP1+bestSolution);
            nbrApparition=1;
            bestSolutionIteration=iteration;

        }
        else
        {
            nbrApparition++;
        }
}

///****************************************************************************************************************************
/// cette fonction affecte les noeuds non medians au noeuds medians, permettant la création des
/// PHENOTYPE DES SOLUTION CE QUI PERMET PAR LA SUITE LE CALCULE DES FITNESS
///void affectationDesNonMediansBaseSurLesArretes(partition* population,int taillePopulation, int nbrNoeuds, int fluxMatrix[nbrNoeuds][nbrNoeuds],int max_sizeCluster){
void affectationDesNonMediansBaseSurClusterSize(partitionPMP* populationPMP)
{

    int i,j,l,k,m,maxDistance,indiceNotSeed,indiceSeed,nbrSommets,nonCohabit;
    for(i=0; i<taillePopulation; i++)
    {
///***************************************************************************************************
        for(l=0; l<nbrNoeuds; l++)
        {
            (populationPMP+i)->phenotype[l] = l;
        }

        for(l=0; l<(populationPMP+i)->medians; l++)
        {
            (populationPMP+i)->clustersSize[(populationPMP+i)->mediansVertices[l]] = 1;
        }
///***************************************************************************************************

        if((populationPMP+i)->nonMedians > 0 && (populationPMP+i)->medians > 0) // sometimes we have some solutions that contains no medians or not medians vertices
        {
            for(j=0; j<(populationPMP+i)->nonMedians;j++)
            {
                maxDistance = -1;
                indiceNotSeed= (populationPMP+i)->nonMediansVertices[j];

                for(k=0; k<(populationPMP+i)->medians; k++)
                {
                    indiceSeed = (populationPMP+i)->mediansVertices[k];
                    nonCohabit = 0;
                    for(l=0;l<nbrNonCohabitationConstraintes;l++){
                        if((nonCohabitationConstraintes[l].nouedDepart == indiceSeed && nonCohabitationConstraintes[l].nouedArrive == indiceNotSeed)||
                           (nonCohabitationConstraintes[l].nouedDepart == indiceNotSeed && nonCohabitationConstraintes[l].nouedArrive == indiceSeed))
                        {nonCohabit = 1; break;}

                    }

                    if( fluxMatrix[indiceNotSeed][indiceSeed] >= 0
					&&
					(populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceSeed]] < max_sizeCluster
                    &&
                    !nonCohabit ) /// nonCohabit = 0 ==> les deux peuvent cohabiter
                    {
                        if(maxDistance == -1)
                        {
                            maxDistance = fluxMatrix[indiceNotSeed][indiceSeed];
                            (populationPMP+i)->phenotype[indiceNotSeed] = (populationPMP+i)->phenotype[indiceSeed];

                        }
                        else  if(maxDistance < fluxMatrix[indiceNotSeed][indiceSeed])
                        {
                            maxDistance = fluxMatrix[indiceNotSeed][indiceSeed];
                            (populationPMP+i)->phenotype[indiceNotSeed] = (populationPMP+i)->phenotype[indiceSeed];
                        }
                    }
                }
                if(maxDistance >= 0)
                {
                    indiceSeed = (populationPMP+i)->phenotype[indiceNotSeed] ;
                    (populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceSeed]]++;
                }
                else
                {
                    for(k=0; k<(populationPMP+i)->medians; k++)
                    {
                        indiceSeed = (populationPMP+i)->mediansVertices[k];
                        if((populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceSeed]] < max_sizeCluster && !nonCohabit)
                        {
                            maxDistance = fluxMatrix[indiceNotSeed][indiceSeed];
                            (populationPMP+i)->phenotype[indiceNotSeed] = (populationPMP+i)->phenotype[indiceSeed];
                            (populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceNotSeed]]++;
                            break;
                        }
                    }
                    if(maxDistance == -1)
                    {
                        (populationPMP+i)->phenotype[indiceNotSeed] = (populationPMP+i)->phenotype[(populationPMP+i)->mediansVertices[(populationPMP+i)->medians-1]];
                        (populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceNotSeed]]++;
                    }
                }
            }
        }
    }

}


///**********************************************************************************************************
void writeOptimalSolutionInFilePMP(partitionPMP *solutionDominante,FILE* outputOptimalSolutionFilePMP,
                                   int nbrRun, int bestSolutionIteration, float runTime, int ES)
{
    int i;
    fprintf(outputOptimalSolutionFilePMP,"PMEB | RUN NUMBER = %d | ITERATION N = %d | SOLUTION IDENTIFIER = %d | CUT SIZE INTRA = %d | CUT SIZE INTER = %d | SOMME TOTAL DES FLUX = %d | RUN TIME = %f | ES = %d | NBR VIOLATED CONSTRAINTE = %d |",
            nbrRun, bestSolutionIteration, solutionDominante->id, solutionDominante->coutCoupe,(sommeTotalFlux - solutionDominante->coutCoupe),sommeTotalFlux ,runTime, ES, solutionDominante->contrainteViole);


    fprintf(outputOptimalSolutionFilePMP," max_clusters = %d | max_sizeCluster = %d  |  cohabitation = %d | non cohabitation = %d\n",
            solutionDominante->constraintVector[0],solutionDominante->constraintVector[1],solutionDominante->constraintVector[2],
            solutionDominante->constraintVector[3]);



}
///**********************************************************************************************************

int compareCroissantFitnessPMP(void const *a, void const *b)
{

    partitionPMP const *pa = a;
    partitionPMP const *pb = b;

    return  pb->coutCoupeNormalise - pa->coutCoupeNormalise;;
}
/**
///************************************************************************************************************
void writeOptimalSolutionInFilePMP(partitionPMP *solutionDominante,FILE* outputOptimalSolutionFilePMP){
    int i;
    fprintf(outputOptimalSolutionFilePMP,"id = %d | normalised cut size = %d | cut size intrA = %d cut size intEr = %d ",
           solutionDominante->id,solutionDominante->coutCoupeNormalise, solutionDominante->coutCoupe, (sommeTotalFlux-solutionDominante->coutCoupe));
    fprintf(outputOptimalSolutionFilePMP,"le contraintes violees = %d \n", solutionDominante->contrainteViole);
    for(i=0;i<4;i++){
            switch (i){

            case 0:
                if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFilePMP,"contrainte sur le nombre de clusters : OK \n");
                else fprintf(outputOptimalSolutionFilePMP,"contrainte sur le nombre de clusters : NOT OK \n");
                break;
            case 1:
                if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFilePMP,"contrainte sur la taille des clusters : OK \n");
                else fprintf(outputOptimalSolutionFilePMP,"contrainte sur la taille des clusters : NOT OK \n");
                break;
            case 2:
                if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFilePMP,"contrainte de cohabitation : OK \n");
                else fprintf(outputOptimalSolutionFilePMP,"contrainte de cohabitation : NOT OK \n");
                break;
            case 3:
                if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFilePMP,"contrainte Non Cohabitation : OK \n");
                else fprintf(outputOptimalSolutionFilePMP,"contrainte Non Cohabitation : NOT OK \n");
                break;
            default:
                fprintf(outputOptimalSolutionFilePMP,"Cette contrainte n'est pas prise en charge par le système \n");

            }

    }
///************************************************************************************************************
void checkContrainstAndFitnessPenalizationPMP(partitionPMP *populationPMP)
{
    int i,j;
    for(i=0; i<taillepopulationPMP; i++)
    {
        ///initialisation de vecteur des contraintes et de la variable de contraintes violées
        (populationPMP+i)->contrainteViole=0;
        for(j=0; j<nbr_constraint; j++)
        {
            (populationPMP+i)->constraintVector[j]=0;
        }

        for(j=0; j<(populationPMP+i)->medians; j++)
        {
            if((populationPMP+i)->clustersSize[j]>max_sizeCluster || (populationPMP+i)->clustersSize[j]< min_sizeCluster)
            {
                ///24.07.2015 ---> l'ajout de vecteur de contraintes = le nombre de cluster dont la taille a dépassé la borne
                (populationPMP+i)->constraintVector[1]++;

            }
        }



        if((populationPMP+i)->medians > max_clusters || (populationPMP+i)->medians < min_clusters)
        {
            (populationPMP+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationPMP+i)->contrainteViole++;
        }

        if((populationPMP+i)->constraintVector[1] !=0)
        {
            (populationPMP+i)->contrainteViole++;
        }
    }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0)
    {
        checkCohabitationAndNonCohabitationConstraintePMP(populationPMP);
    }

}
}*/


///****************************************************************************************************************************
/// cette fonction affecte les noeuds non medians au noeuds medians, permettant la création des
/// PHENOTYPE DES SOLUTION CE QUI PERMET PAR LA SUITE LE CALCULE DES FITNESS
///void affectationDesNonMediansBaseSurLesArretes(partition* population,int taillePopulation, int nbrNoeuds, int fluxMatrix[nbrNoeuds][nbrNoeuds],int max_sizeCluster){
/**
void affectationDesNonMediansParArretes(partitionPMP* populationPMP)
{

    int i,j,l,k,maxDistance,indiceNotSeed,indiceSeed,_medians, _mediansTmp,_nonMedians, _nonMediansTmp,nbrSommets;
    int _nonMediansVertices[1000]= {0},_mediansVertices[1000]= {0};

    for(i=0; i<taillePopulation; i++)
    {
        /// initialisation de vecteur des capacité avec 1 pour chaque médians
        /// initialisation des phénotype par la valeurs

        //(populationPMP+i)->medians = 3;(populationPMP+i)->mediansVertices[0] = 2; (populationPMP+i)->mediansVertices[1] = 4;(populationPMP+i)->mediansVertices[2] = 5;
        //(populationPMP+i)->nonMedians = 3; (populationPMP+i)->nonMediansVertices[0] = 0; (populationPMP+i)->nonMediansVertices[1] = 1; (populationPMP+i)->nonMediansVertices[2] = 3;
        //(populationPMP+i)->genotype[0]= 0; (populationPMP+i)->genotype[1]= 0; (populationPMP+i)->genotype[2]= 1;
        //(populationPMP+i)->genotype[3]= 0;(populationPMP+i)->genotype[4]= 1;(populationPMP+i)->genotype[5]= 1;


        _medians = _mediansTmp = (populationPMP+i)->medians;;
        _nonMedians =  _nonMediansTmp = (populationPMP+i)->nonMedians;;

        for(l=0; l<_nonMedians; l++)
        {
            _nonMediansVertices[l] = (populationPMP+i)->nonMediansVertices[l];
            (populationPMP+i)->phenotype[_nonMediansVertices[l]] = _nonMediansVertices[l];
            (populationPMP+i)->clustersSize[l]= 0;
        }

        for(l=0; l<_medians; l++)
        {
            _mediansVertices[l] = (populationPMP+i)->mediansVertices[l];
            (populationPMP+i)->phenotype[_mediansVertices[l]] = _mediansVertices[l];
            (populationPMP+i)->clustersSize[_mediansVertices[l]] = 1;

        }

        while(_nonMedians > 0 && _medians > 0)
        {
            j =0;
            _nonMediansTmp = _nonMedians ;
            _mediansTmp = _medians;
            maxDistance = -1;
            ///mon_sleep(1);
            while(j < _nonMediansTmp)
            {
                maxDistance = -1;
                indiceNotSeed= _nonMediansVertices[j];

                for(k=0; k<_mediansTmp; k++)
                {
                    indiceSeed = _mediansVertices[k];
                    if(fluxMatrix[indiceNotSeed][indiceSeed] >= 0 )
                    {
                        if(maxDistance == -1)
                        {
                            maxDistance = fluxMatrix[indiceNotSeed][indiceSeed];
                            (populationPMP+i)->phenotype[indiceNotSeed] = (populationPMP+i)->phenotype[indiceSeed];
                            /// il faut faire la mise à jours des taille des cluster à temps
                            (populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceSeed]]++;


                        }
                        else if((populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceSeed]] < max_sizeCluster)
                        {
                            ///*********************************
                            if((populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceNotSeed]] > max_sizeCluster)
                            {
                                maxDistance = fluxMatrix[indiceNotSeed][indiceSeed];
                                (populationPMP+i)->phenotype[indiceNotSeed] = (populationPMP+i)->phenotype[indiceSeed];
                                /// Mise à jour des tailles des cluster 10/11/2016
                                (populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceSeed]]++;
                                (populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceNotSeed]]--;
                            }
                            else if((populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceNotSeed]] < max_sizeCluster)
                            {
                                if(maxDistance < fluxMatrix[indiceNotSeed][indiceSeed])
                                {
                                    maxDistance = fluxMatrix[indiceNotSeed][indiceSeed];
                                    (populationPMP+i)->phenotype[indiceNotSeed] = (populationPMP+i)->phenotype[indiceSeed];
                                    /// Mise à jour des tailles des cluster 10/11/2016
                                    (populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceSeed]]++;
                                    (populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceNotSeed]]--;
                                }
                            }
                           ///*********************************
                        }
                    }
                }

                if(maxDistance >= 0)
                {

                    ///indiceSeed = (populationPMP+i)->phenotype[indiceNotSeed] ;
                    ///(populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceSeed]]++;
                    _mediansVertices[_medians] = indiceNotSeed;
                    _medians++;
                    ///_mediansTmp++;

                    for(l=j; l<_nonMediansTmp; l++)
                    {
                        _nonMediansVertices[l] = _nonMediansVertices[l+1];
                    }
                    _nonMedians--;
                    _nonMediansTmp--;
                }
                else
                {
                    j++;
                }

            }
        }

    }

}

*/
