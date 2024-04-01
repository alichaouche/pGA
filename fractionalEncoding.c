#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "typeDeclaration.h"
#include "./communesFunctions.h"
#include "./fractionalEncoding.h"
#include "./compilationConditionnelle.h"

extern int taillePopulation,tauxReproduction,sommeTotalFlux,nbrNoeuds,
       nbr_constraint,nbrArretes,nbrParties,nbrNoeudsInCluster,regeneration;
extern double tauxMutation ,averageTraffic;
extern int fluxMatrix[1000][1000],*fluxVector, maxFlow, sommeTotalFluxReal;
extern edge *edgeVector,cohabitationConstraintes[1000],nonCohabitationConstraintes[1000];
extern int max_clusters , min_clusters,min_sizeCluster ,max_sizeCluster ;
extern int bestSolution, nbrApparition, nbrCohabitationConstraintes,nbrNonCohabitationConstraintes;
extern int nbrRun;
///**************************************************************************************
///génération de la population initiale
void generatePopulationFVTC(partitionFVTC* populationFVTC, int indiceFirstElt )
{

    int i,j;
    for(i= indiceFirstElt; i<taillePopulation; i++)
    {
        ///srand(time(NULL));
        for(j=0; j<=nbrNoeuds; j++) /// <= pour qu'on aura un élément en plus pour le nombre des parties
        {
#if Windows
            (populationFVTC+i)->genotype[j] = drand48ForWindowsFVTC(0,101);
#else
            (populationFVTC+i)->genotype[j] = drand48();
#endif // Windows

        }
        (populationFVTC+i)->id =i;
    }
}
///**************************************************************************************
/// affichage des population
void affichePopulationFVTC(partitionFVTC* populationFVTC)
{
    int i,j;
    for(i=0; i<taillePopulation; i++)
    {
        printf("\nid = %d \t cout de coupe = %d \t fitness = %0.4f \n",
               (populationFVTC+i)->id,(populationFVTC+i)->coutCoupe, (populationFVTC+i)->fitness);
        printf("le cout de coupe normalise = %d \n",(populationFVTC+i)->coutCoupeNormalise);
        printf("le nombre des contraintes violees = %d \n",(populationFVTC+i)->contrainteViole);
        printf("le nombre des cluster de la partition = %d \n",(populationFVTC+i)->nbrCluster);
        printf("la genotype est  = \t");
        for(j=0; j<=nbrNoeuds; j++)
        {
            printf("%0.2f ",(populationFVTC+i)->genotype[j]);
        }
        printf("\n");
        printf("la phenotype est  = \t");
        for(j=0; j<=nbrNoeuds; j++)
        {
            printf("%d ",(populationFVTC+i)->phenotype[j]);
        }
        printf("\n_________________________________________________________________________\n");
    }
}
///************************************************************************************************************
/// calculer le cout de coupe de chaque partition ainsi que son fitness
///***************************************************************************************
/// on va appliquer la formule suivante pour le calcule des cout de coupe
/// T'(C) = B - T(C)            <=> la minimisation devient une maximisation
/// Z(C)=T'(C)+(u-v(C)) ×B      <=>  calculer la nouvelle fitness par rapport aux contrainte
/// T'(C) : cout de coupe (maximiser les flux INTRA), u ; nombre de contrainte totale
/// v(C): nombre de contrainte violé,  B: la somme totale des flux
///***************************************************************************************

void calculCoutCoupeEtFitnessFVTC(partitionFVTC* populationFVTC, int indiceFirst, int indiceLast)
{

    int i,j,k,sommeTotalCoutCoupe = 0;
    for(i=indiceFirst; i<indiceLast; i++) /// parcourir tous les individus de la population
    {
        (populationFVTC+i)->coutCoupe  =0; /// réinitialisation de la valriable après chaque itération (individus de la population)

        for(j=0; j<nbrNoeuds-1; j++)
        {
            for(k=j+1; k<nbrNoeuds; k++)
            {
                if((populationFVTC+i)->phenotype[j] == (populationFVTC+i)->phenotype[k] && fluxMatrix[j][k] > 0)
                {
                    (populationFVTC+i)->coutCoupe  = (populationFVTC+i)->coutCoupe  + fluxMatrix[j][k];
                }
            }
        }

        (populationFVTC+i)->coutCoupeNormalise = (populationFVTC+i)->coutCoupe + ((nbr_constraint - (populationFVTC+i)->contrainteViole)*sommeTotalFlux);

/// désormis on va utiliser le cout de coupe normaliser pour le calcule des fitness

        ///sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationFVTC+i)->coutCoupe;
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationFVTC+i)->coutCoupeNormalise;
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
            varianceCoutDeCoupeNormalise + pow(((populationFVTC+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);
    /// tester la valeur ajouté
    for(i=0; i<taillePopulation; i++)
    {

        /// calcul de sigma truncation
        (populationFVTC+i)->expectedValue = (populationFVTC+i)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationFVTC+i)->expectedValue;
    }

    /// calcule du fitness
    /// le 26/11/2015 : calcul du fitness des solution avec leur expectedValues
    for(i=0; i<taillePopulation ; i++)
    {
        (populationFVTC+i)->fitness = (float)((populationFVTC+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);
    }

#else
    for(i=0; i<taillePopulation ; i++)
    {
        (populationFVTC+i)->fitness = (float)((populationFVTC+i)->coutCoupeNormalise)/(float)(sommeTotalCoutCoupe);
    }
#endif// scaling
}
///**************************************************************************************
void naturalSelectionFVTC(partitionFVTC* populationFVTC1,partitionFVTC* populationFVTC2)
{
    int i,j,maxFitness;
    int tmpVectIndice[taillePopulation];
    for(i=0;i<taillePopulation;i++) tmpVectIndice[i]=-1;
    float lotterie,sommeFitness;

    for(i=0; i<tauxReproduction; i++){
        maxFitness = 0;
        for(j=0;j<taillePopulation;j++){
            if((populationFVTC1+maxFitness)->coutCoupeNormalise < (populationFVTC1+j)->coutCoupeNormalise
               && tmpVectIndice[j] == -1 ){
                maxFitness = j;

            }
        }
        tmpVectIndice[maxFitness] = 0;
        *(populationFVTC2+i) = *(populationFVTC1+maxFitness);
        ///printf("(populationFVTC2+%d)->coutCoupeNormalise = %d \n",i,(populationFVTC2+i)->coutCoupeNormalise);
        ///mon_sleep(pause);
    }
    ///printf("\n\n\n");

    ///qsort (populationFVTC2, taillePopulation, sizeof *populationFVTC2, compareCroissantFitnessEE);
    ///sortingPopulationByFintness(populationFVTC2,tmpPopulation);

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
            sommeFitness = sommeFitness + (populationFVTC1+j)->fitness;
            if (lotterie<=sommeFitness)
            {
                *(populationFVTC2+i) = *(populationFVTC1+j);
                break;
            }
            else if( j==taillePopulation-1)
            {
                *(populationFVTC2+i) = *(populationFVTC1+j);
            }
        }
    }
    ///free(tmpPopulation);
}
///**************************************************************************************
void crossOverFVTC(partitionFVTC* populationFVTC1, partitionFVTC* populationFVTC2)
{
    int i ,j;
    int choixInd1, choixInd2, choixLocus;
    ///***********************************************************************************
    for(i=0; i<tauxReproduction; i++)
    {
        *(populationFVTC2+i) = *(populationFVTC1+i);
        (populationFVTC2+i)->id = i;
    }
    ///***********************************************************************************
#if NEW_CROSSOVER
    while (i < taillePopulation-regeneration)
    {
        ///choixInd1 = rnd(tauxReproduction,taillePopulation); /** l'intervalle sera [0;taillePopulation-1] */
        choixInd1 = rnd(0,taillePopulation);
        do
        {
            /// ce traitement évite l'appariement d'un element avec lui même/
            ///choixInd2 = rnd(tauxReproduction,taillePopulation);
            choixInd2 = rnd(0,taillePopulation);
        }
        while(choixInd2 == choixInd1 );

        choixLocus = rnd(0,nbrNoeuds-1); /** l'intervalle sera [1;nbrNoeuds-2] */

        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationFVTC2+i)->genotype[j] = (populationFVTC1+choixInd1)->genotype[j];
            (populationFVTC2+i+1)->genotype[j] = (populationFVTC1+choixInd2)->genotype[j];
        }
        ///********************************************************************************
        for(j=choixLocus+1; j<=nbrNoeuds; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationFVTC2+i)->genotype[j] = (populationFVTC1+choixInd2)->genotype[j];
            (populationFVTC2+i+1)->genotype[j] = (populationFVTC1+choixInd1)->genotype[j];
        }

        /// affectation de nouveau indice
        (populationFVTC2+i)->id = i;
        (populationFVTC2+i+1)->id = i+1;
        i+=2;
    }
    generatePopulationFVTC(populationFVTC2, taillePopulation-regeneration);
#else

    while (i < taillePopulation)
    {
        ///choixInd1 = rnd(tauxReproduction,taillePopulation); /** l'intervalle sera [0;taillePopulation-1] */
       choixInd1 = rnd(0,taillePopulation); /** l'intervalle sera [0;taillePopulation-1] */
        do
        {
            /// ce traitement évite l'appariement d'un element avec lui même/
            ///choixInd2 = rnd(tauxReproduction,taillePopulation);
            choixInd2 = rnd(0,taillePopulation);
        }
        while(choixInd2 == choixInd1 );

        choixLocus = rnd(0,nbrNoeuds-1); /** l'intervalle sera [1;nbrNoeuds-2] */

        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationFVTC2+i)->genotype[j] = (populationFVTC1+choixInd1)->genotype[j];
            (populationFVTC2+i+1)->genotype[j] = (populationFVTC1+choixInd2)->genotype[j];
        }
        ///********************************************************************************
        for(j=choixLocus+1; j<=nbrNoeuds; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationFVTC2+i)->genotype[j] = (populationFVTC1+choixInd2)->genotype[j];
            (populationFVTC2+i+1)->genotype[j] = (populationFVTC1+choixInd1)->genotype[j];
        }

        /// affectation de nouveau indice
        (populationFVTC2+i)->id = i;
        (populationFVTC2+i+1)->id = i+1;
        i+=2;
    }

#endif

}
///**************************************************************************************
void calculPhenotypeFVTC(partitionFVTC* populationFVTC)
{

    /**<
    à partir d'une solution générée aléatoirement on déduit la partition
    1- multiplier la dernière valeur * max_clusters pour avoir le nombre de cluster en prennant la borne sup
    2- la valeur obtenu dans 1 est utilisée pour déterminer l'affectation des noeuds aux clusters
     */
    int i,j,nombreDesParties;
    for(i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        nombreDesParties =ceil((populationFVTC+i)->genotype[nbrNoeuds]*(max_clusters-1)); /// réinitialisation de la valriable après chaque itération (individus de la population)
        (populationFVTC+i)->phenotype[nbrNoeuds] = nombreDesParties;
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationFVTC+i)->phenotype[j] = floor((populationFVTC+i)->genotype[j]*nombreDesParties) +1;
        }
    }
}
///**************************************************************************************
int findTheBestSolutionFVTC(partitionFVTC *populationFVTC)
{
    float maxFitness = 0;
    int i,indice=0;
    for(i=0; i<taillePopulation ; i++)
    {
        if(maxFitness < (populationFVTC+i)->fitness)
        {
            maxFitness = (populationFVTC+i)->fitness;
            indice = i;
        }
    }
    return indice;
}
///**************************************************************************************
void mutationFVTC(partitionFVTC* populationFVTC)
{

    int i,numeroGeneMute=-1;
    double applicationOfMutation;
    for(i=tauxReproduction; i<taillePopulation; i++)
    {
#if Windows
        applicationOfMutation = drand48ForWindows(0,1001);
#else
        applicationOfMutation = drand48();
#endif // Windows
        if(applicationOfMutation <= tauxMutation)
        {
            /// choisir un gène aléatoirement pour lui appliquer la mutation
            numeroGeneMute = rnd(0,nbrNoeuds+1);
#if Windows
            (populationFVTC+i)->genotype[numeroGeneMute] = drand48ForWindows(0,101);
#else
            (populationFVTC+i)->genotype[numeroGeneMute] = drand48();
#endif // Windows

        }
    }

}
///**************************************************************************************
float testerLaSommeDesFitnessFVTC(partitionFVTC* populationFVTC)
{

#if mouchard
    printf("testerLaSommeDesFitness ...\n");
#endif // mouchard
    int i;
    float sommeFitness = 0;
    for(i=0; i<taillePopulation ; i++)
    {
        sommeFitness = sommeFitness + (populationFVTC+i)->fitness;
    }
    return sommeFitness;
}
///**************************************************************************************
void displayTheBestSolutionFVTC(partitionFVTC* solutionDominante)
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

///**************************************************************************************
void writeSolutionInFileFVTC(partitionFVTC *populationFVTC, FILE *outputFilePop,int iteration)
{

    int i,j;
    for(i=0; i<taillePopulation; i++)
    {

        ///écriture de la meilleur solution dans le fichier output
        fprintf(outputFilePop,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,(populationFVTC+i)->id,
                (populationFVTC+i)->coutCoupe,(populationFVTC+i)->fitness,(populationFVTC+i)->coutCoupeNormalise,
                (populationFVTC+i)->contrainteViole,(populationFVTC+i)->nbrCluster);
        for(j=0; j<=nbrNoeuds; j++)
        {
            fprintf(outputFilePop,"%0.2f ",(populationFVTC+i)->genotype[j]);
        }
        fprintf(outputFilePop,"\t\t\t\t");
        for(j=0; j<=nbrNoeuds; j++)
        {
            fprintf(outputFilePop,"%d ",(populationFVTC+i)->phenotype[j]);
        }
        fprintf(outputFilePop,"\n");
    }
    /// fprintf(outputFilePop,"\n===============================================\n");
    fprintf(outputFilePop,"\n\n");

}
///**************************************************************************************
void writeBestSolutionInFileFVTC(partitionFVTC *solutionDominante, FILE *outputFile,int iteration)
{

    int i;
    fprintf(outputFile,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,solutionDominante->id,
            solutionDominante->coutCoupe,solutionDominante->fitness,solutionDominante->coutCoupeNormalise,
            solutionDominante->contrainteViole,solutionDominante->nbrCluster);
    for(i=0; i<=nbrNoeuds; i++)
    {
        fprintf(outputFile,"%0.2f ",solutionDominante->genotype[i]);
    }
    fprintf(outputFile,"\t\t\t\t");
    for(i=0; i<=nbrNoeuds; i++)
    {
        fprintf(outputFile,"%d ",solutionDominante->phenotype[i]);
    }
    fprintf(outputFile,"\n");

}

///=======================================================================================
void fractionalEncoding(int nbrGeneration,FILE *outputFileFVTC,FILE *outputFilePopFVTC,FILE *outputOptimalSolutionFileFVTC,partitionFVTC *populationFVTC1,
                                        partitionFVTC *populationFVTC2,partitionFVTC *solutionDominante,int iteration , int *bestSolutionIteration , int *nbrApparition)
{

        naturalSelectionFVTC(populationFVTC1,populationFVTC2);
        crossOverFVTC(populationFVTC2,populationFVTC1); /// nbrNoeuds représente la taille des genotype
        mutationFVTC(populationFVTC1);
        calculPhenotypeFVTC(populationFVTC1);
        checkContrainstAndFitnessPenalizationFVTC(populationFVTC1);
        calculCoutCoupeEtFitnessFVTC(populationFVTC1,0,taillePopulation);


#if writingPopulationInFile
        writeSolutionInFileFVTC(populationFVTC1,outputFilePopFVTC,iteration);
#endif // writingPopulationInFile
///=====================================================================================================
        /// dans ce qui suit, on va juste localiser la meilleur solution
///=====================================================================================================

        bestSolution=findTheBestSolutionFVTC(populationFVTC1);
#if writingPopulationInFile
        writeBestSolutionInFileFVTC((populationFVTC1+bestSolution),outputFileFVTC,iteration);
#endif // writingPopulationInFile
        if((populationFVTC1+bestSolution)->coutCoupeNormalise > solutionDominante->coutCoupeNormalise)
        {
            *(solutionDominante) = *(populationFVTC1+bestSolution);
            nbrApparition=1;
            bestSolutionIteration=iteration;
        }
        else
        {
            nbrApparition++;
        }
        ///affichePopulation(populationFVTC1,taillePopulation,nbrNoeuds);


}

///**************************************************************************************
void checkContrainstAndFitnessPenalizationFVTC(partitionFVTC *populationFVTC)
{

    int i,j,k;

    for(i=0; i<taillePopulation; i++)
    {
        (populationFVTC+i)->nbrCluster=0; /// le nombre de cluster doit etre initialisé à 1 pour prendre en consédiration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes violées
        (populationFVTC+i)->contrainteViole=0;
        for(j=0; j<4; j++)
        {
            (populationFVTC+i)->constraintVector[j]=0;
        }
        ///**********************************************************************************



        for(j=0; j<nbrNoeuds; j++)
        {
            (populationFVTC+i)->clustersSize[j] =0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if ((populationFVTC+i)->phenotype[k] == j)
                {
                    (populationFVTC+i)->clustersSize[j]++;
                }
            }
            if ((populationFVTC+i)->clustersSize[j] !=0)
            {
                (populationFVTC+i)->nbrCluster++;
            }
        }

        if((populationFVTC+i)->nbrCluster > max_clusters || (populationFVTC+i)->nbrCluster < min_clusters)
        {
            (populationFVTC+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationFVTC+i)->contrainteViole++;
        }

        for(j=0; j<nbrNoeuds; j++)
        {
            if((populationFVTC+i)->clustersSize[j]!=0)
            {
                if((populationFVTC+i)->clustersSize[j]>max_sizeCluster || (populationFVTC+i)->clustersSize[j]< min_sizeCluster)
                {
                    (populationFVTC+i)->constraintVector[1]++;

                }
            }

        }

        if((populationFVTC+i)->constraintVector[1] != 0)
        {
            (populationFVTC+i)->contrainteViole++;
        }
    }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0)
    {
        checkCohabitationAndNonCohabitationConstrainteFVTC(populationFVTC);
    }

}

///**************************************************************************************
void checkCohabitationAndNonCohabitationConstrainteFVTC(partitionFVTC *populationFVTC)
{

    int i,j,noeud1,noeud2;
    for(i=0; i<taillePopulation; i++)
    {
        ///24/11/2015 normalement cette taches est déjà faite dans le fonction : checkContrainstAndFitnessPenalizationFVTC
        /**
        (populationFVTC+i)->constraintVector[2]=0;
        (populationFVTC+i)->constraintVector[3]=0;
        */
        if(nbrCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrCohabitationConstraintes; j++)
            {
                noeud1=cohabitationConstraintes[j].nouedDepart;
                noeud2=cohabitationConstraintes[j].nouedArrive;

                if((populationFVTC+i)->phenotype[noeud1]!= (populationFVTC+i)->phenotype[noeud2])
                {
                    (populationFVTC+i)->constraintVector[2]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationFVTC+i)->constraintVector[2]!=0)
            {
                (populationFVTC+i)->contrainteViole++;
            }
        }


        if(nbrNonCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrNonCohabitationConstraintes; j++)
            {
                noeud1=nonCohabitationConstraintes[j].nouedDepart;
                noeud2=nonCohabitationConstraintes[j].nouedArrive;

                if((populationFVTC+i)->phenotype[noeud1]== (populationFVTC+i)->phenotype[noeud2])
                {
                    (populationFVTC+i)->constraintVector[3]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationFVTC+i)->constraintVector[3]!=0)
            {
                (populationFVTC+i)->contrainteViole++;
            }
        }


    }

}

///********************************************************************************************
void writeOptimalSolutionInFileFVTC(partitionFVTC *solutionDominante,FILE* outputOptimalSolutionFileFVTC,
                                  int nbrRun, int bestSolutionIteration, float runTime, int ES)
{

    fprintf(outputOptimalSolutionFileFVTC,"FVTC | RUN NUMBER = %d | ITERATION N = %d | SOLUTION IDENTIFIER = %d | CUT SIZE INTRA = %d | CUT SIZE INTER = %d | SOMME TOTAL DES FLUX = %d | RUN TIME = %f | ES = %d | NBR VIOLATED CONSTRAINTE = %d | ",
            nbrRun, bestSolutionIteration, solutionDominante->id, solutionDominante->coutCoupe,(sommeTotalFlux - solutionDominante->coutCoupe),sommeTotalFlux ,runTime, ES, solutionDominante->contrainteViole);

    fprintf(outputOptimalSolutionFileFVTC," max_clusters = %d | max_sizeCluster = %d  |  cohabitation = %d | non cohabitation = %d\n",

            solutionDominante->constraintVector[0],
            solutionDominante->constraintVector[1],
            solutionDominante->constraintVector[2],
            solutionDominante->constraintVector[3]);


}
///******************************************************************************************
int compareCroissantFitnessFVTC (void const *a, void const *b)
{

    partitionFVTC const *pa = a;
    partitionFVTC const *pb = b;
    return  pb->coutCoupeNormalise - pa->coutCoupeNormalise;
}

double drand48ForWindowsFVTC(int v1, int v2)
{
    ///return (double)(rand()) / (double)(RAND_MAX);
    return (double)(rnd(v1,v2)) /100.0;
}


/**


void writeOptimalSolutionInFileFVTC(partitionFVTC *solutionDominante,FILE* outputOptimalSolutionFileFVTC)
{
    int i;
    fprintf(outputOptimalSolutionFileFVTC,"id = %d | normalised cut size = %d | cut size intrA = %d cut size intEr = %d ",
            solutionDominante->id,solutionDominante->coutCoupeNormalise, solutionDominante->coutCoupe, (sommeTotalFlux-solutionDominante->coutCoupe));
    fprintf(outputOptimalSolutionFileFVTC,"le contraintes violees = %d \n", solutionDominante->contrainteViole);
    for(i=0; i<4; i++)
    {
        switch (i)
        {

        case 0:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileFVTC,"contrainte sur le nombre de clusters : OK \n");
            else fprintf(outputOptimalSolutionFileFVTC,"contrainte sur le nombre de clusters : NOT OK \n");
            break;
        case 1:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileFVTC,"contrainte sur la taille des clusters : OK \n");
            else fprintf(outputOptimalSolutionFileFVTC,"contrainte sur la taille des clusters : NOT OK \n");
            break;
        case 2:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileFVTC,"contrainte de cohabitation : OK \n");
            else fprintf(outputOptimalSolutionFileFVTC,"contrainte de cohabitation : NOT OK \n");
            break;
        case 3:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileFVTC,"contrainte Non Cohabitation : OK \n");
            else fprintf(outputOptimalSolutionFileFVTC,"contrainte Non Cohabitation : NOT OK \n");
            break;
        default:
            fprintf(outputOptimalSolutionFileFVTC,"Cette contrainte n'est pas prise en charge par le système \n");

        }

    }

}
void checkContrainstAndFitnessPenalizationFVTC(partitionFVTC *populationFVTC)
{

    int i,j,k,m=nbrNoeuds,tabTmp[m],tmp;

    for(i=0; i<taillePopulation; i++)
    {
        (populationFVTC+i)->nbrCluster=1; /// le nombre de cluster doit etre initialisé à 1 pour prendre en consédiration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes violées, pour chaque solution
        (populationFVTC+i)->contrainteViole=0;
        for(j=0;j<nbr_constraint;j++){(populationFVTC+i)->constraintVector[j]=0;}


        m=nbrNoeuds;
        ///copier le phenotype dans un deuxième tableau
        for(j=0; j<nbrNoeuds; j++)
        {
            tabTmp[j]=(populationFVTC+i)->phenotype[j];
        }
        ///trier le tableau résulant pour savoir le nombre de clusters
        /// par exemple : 2 2 3 4 1 1 2 3 ==> 1 1 2 2 2 3 3 4
        /// l'étape suivante : 1 2 3 4 ==> 4 cluster, et combien de noeuds dans chaque cluster
        /// ceci pour vérifier le respect des contraintes.

        for(j=0; j<m-1; j++)
        {
            for(k=j+1; k<m; k++)
            {
                if (tabTmp[k] < tabTmp[j])
                {
                    tmp = tabTmp[j];
                    tabTmp[j] = tabTmp[k];
                    tabTmp[k] = tmp;
                }
            }
        }
        ///calculer le nombre des cluser
        for(j=1; j<nbrNoeuds; j++)
        {
            if(tabTmp[j-1] != tabTmp[j])
            {
                (populationFVTC+i)->nbrCluster++;
            }
        }

        if((populationFVTC+i)->nbrCluster > max_clusters || (populationFVTC+i)->nbrCluster < min_clusters)
        {
            (populationFVTC+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationFVTC+i)->contrainteViole++;
        }

        /// vérifier la taille des clusters

        for(j=0; j<m-1; j++)  /// je dois me limitter à m-1 pour que t arrive jusqu'à m
        {
            while(tabTmp[j] == tabTmp[j+1])
            {
                for(k=j; k<m; k++) /// la taille du vecteur change à chaque fois => m et non pas nbrNoeuds
                {
                    tabTmp[k] = tabTmp[k+1];
                }
                m--;
            }
        }
        for(j=0; j<m; j++)
        {
            (populationFVTC+i)->clustersSize[j]=0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if (tabTmp[j] == (populationFVTC+i)->phenotype[k])
                {
                    (populationFVTC+i)->clustersSize[j]=(populationFVTC+i)->clustersSize[j]+1;
                }
            }
            if((populationFVTC+i)->clustersSize[j]>max_sizeCluster || (populationFVTC+i)->clustersSize[j]< min_sizeCluster)
            {
                ///24.07.2015 ---> l'ajout de vecteur de contraintes = le nombre de cluster dont la taille a dépassé la borne
                (populationFVTC+i)->constraintVector[1]++;

            }
        }

        if((populationFVTC+i)->constraintVector[1] != 0)
        {
            (populationFVTC+i)->contrainteViole++;
        }
   }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0){
        checkCohabitationAndNonCohabitationConstrainteFVTC(populationFVTC);
    }

}

*/
