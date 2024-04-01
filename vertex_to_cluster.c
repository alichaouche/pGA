#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "typeDeclaration.h"
#include "./communesFunctions.h"
#include "./vertex_to_cluster.h"
#include "./compilationConditionnelle.h"


extern int taillePopulation,tauxReproduction,sommeTotalFlux,nbrNoeuds,nbr_constraint,
            nbrArretes,nbrParties,nbrNoeudsInCluster, regeneration;
extern double tauxMutation ,averageTraffic;
extern int fluxMatrix[1000][1000],*fluxVector, maxFlow, sommeTotalFluxReal;
extern edge *edgeVector,cohabitationConstraintes[1000],nonCohabitationConstraintes[1000];
extern int max_clusters , min_clusters,min_sizeCluster ,max_sizeCluster ;
extern int bestSolution, nbrApparition, nbrCohabitationConstraintes,nbrNonCohabitationConstraintes;
extern int nbrRun;

///**************************************************************************************
///génération de la population initiale
void generatePopulationVTC(partitionVTC* populationVTC, int indiceFirstElt)
{

    int i,j;
    for(i=indiceFirstElt; i<taillePopulation; i++)
    {
        for(j=0; j<nbrNoeuds; j++)
        {
            /// on va essayer de limiter le nombre de clusters pour voir ce que ça va donner
            ///(populationVTC+i)->genotype[j] = rnd(0,nbrParties);
            (populationVTC+i)->genotype[j] = rnd(0,nbrNoeuds); /// [0;nbrNoeuds-1]
        }
        (populationVTC+i)->id =i;
    }
}
///**************************************************************************************
/// affichage des population
void affichePopulationVTC(partitionVTC* populationVTC)
{
    int i,j;
    for(i=0; i<taillePopulation; i++)
    {
        printf("\n id = %d \t cout de coupe = %d \t fitness = %0.4f \t",
               (populationVTC+i)->id,(populationVTC+i)->coutCoupe, (populationVTC+i)->fitness);
        printf("le cout de coupe normalise = %d \n",(populationVTC+i)->coutCoupeNormalise);
        printf("la partition est  = \t");
        for(j=0; j<nbrNoeuds; j++)
        {
            printf("%d ",(populationVTC+i)->genotype[j]);
        }
        printf("\n");
        printf("le nombre des clusters est %d \n",(populationVTC+i)->nbrCluster);
        printf("le nombre des contraintes violees est  = %d \n",(populationVTC+i)->contrainteViole);
        printf("affichage des tailles des clusters \n");
        for(j=0; j<(populationVTC+i)->nbrCluster; j++)
        {
            printf("%d ",(populationVTC+i)->clustersSize[j]);
        }
        printf("\n============================================================\n");

    }
}
///************************************************************************************************************
/// calculer le cout de coupe de chaque partition ainsi que son fitness
///***************************************************************************************
/// on va appliquer la formule suivante pour le calcul des couts de coupe
/// T'(C) = B - T(C)            <=> la minimisation devient une maximisation
/// Z(C)=T'(C)+(u-v(C)) ×B      <=>  calculer la nouvelle fitness par rapport aux contrainte
/// T'(C) : cout de coupe (maximiser les flux INTRA), u ; nombre de contrainte totale
/// v(C): nombre de contrainte violé,  B: la somme totale des flux
///***************************************************************************************
void calculCoutCoupeEtFitnessVTC(partitionVTC* populationVTC)
{
    int i,j,k,sommeTotalCoutCoupe = 0;

    for(i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (populationVTC+i)->coutCoupe = 0; /// réinitialisation de la valriable après chaque itération (individus de la population)
        for(j=0; j<nbrNoeuds-1; j++)
        {
            for(k=j+1; k<nbrNoeuds; k++)
            {
                /// on peut avoir deux noeuds du même cluster mais qui sont par relier par une arrêtes
                if((populationVTC+i)->genotype[j] == (populationVTC+i)->genotype[k] && fluxMatrix[j][k] >= 0)
                {
                    (populationVTC+i)->coutCoupe = (populationVTC+i)->coutCoupe +  fluxMatrix[j][k];
                }
            }
        }
        (populationVTC+i)->coutCoupeNormalise = (populationVTC+i)->coutCoupe + ((nbr_constraint - (populationVTC+i)->contrainteViole)*sommeTotalFlux);
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationVTC+i)->coutCoupeNormalise;
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
            varianceCoutDeCoupeNormalise + pow(((populationVTC+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);
    /// tester la valeur ajouté
    for(i=0; i<taillePopulation; i++)
    {

        /// calcul de sigma truncation
        (populationVTC+i)->expectedValue = (populationVTC+i)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationVTC+i)->expectedValue;
    }

    /// calcule du fitness
    /// le 26/11/2015 : calcul du fitness des solution avec leur expectedValues
    for(i=0; i<taillePopulation ; i++)
    {
        (populationVTC+i)->fitness = (float)((populationVTC+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);
    }

#else
    FILE *file;
    file = fopen("D:/WorkingFolderAgProgramming/benchmark/vertexToClusterEncoding/expectedValueVTC.txt","w");
    for(i=0; i<taillePopulation ; i++)
    {
/**
        if(iteration == 1)
        {
            fprintf(file,"%d \t %d \t %0.2f \n",i,(populationVTC+i)->coutCoupeNormalise,(populationVTC+i)->expectedValue);
        }
*/
        (populationVTC+i)->fitness = (float)((populationVTC+i)->coutCoupeNormalise)/(float)(sommeTotalCoutCoupe);
    }
#endif// scaling
}
///**************************************************************************************
void naturalSelectionVTC(partitionVTC* populationVTC1,partitionVTC* populationVTC2)
{
    int i,j,maxFitness;
    int tmpVectIndice[taillePopulation];
    for(i=0;i<taillePopulation;i++) tmpVectIndice[i]=-1;
    float lotterie,sommeFitness;

    for(i=0; i<tauxReproduction; i++){
        maxFitness = 0;
        for(j=0;j<taillePopulation;j++){
            if((populationVTC1+maxFitness)->coutCoupeNormalise < (populationVTC1+j)->coutCoupeNormalise
               && tmpVectIndice[j] == -1 ){
                maxFitness = j;

            }
        }
        tmpVectIndice[maxFitness] = 0;
        *(populationVTC2+i) = *(populationVTC1+maxFitness);
        ///printf("(populationVTC2+%d)->coutCoupeNormalise = %d \n",i,(populationVTC2+i)->coutCoupeNormalise);
        ///mon_sleep(pause);
    }
    ///printf("\n\n\n");

    ///qsort (populationVTC2, taillePopulation, sizeof *populationVTC2, compareCroissantFitnessEE);
    ///sortingPopulationByFintness(populationVTC2,tmpPopulation);

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
            sommeFitness = sommeFitness + (populationVTC1+j)->fitness;
            if (lotterie<=sommeFitness)
            {
                *(populationVTC2+i) = *(populationVTC1+j);
                break;
            }
            else if( j==taillePopulation-1)
            {
                *(populationVTC2+i) = *(populationVTC1+j);
            }
        }
    }
    ///free(tmpPopulation);

}
///**************************************************************************************
void crossOverVTC(partitionVTC* populationVTC1, partitionVTC* populationVTC2)
{
    int i = 0,j =0;
    int choixInd1, choixInd2, choixLocus;
    ///***********************************************************************************
    for(i=0; i<tauxReproduction; i++)
    {
        *(populationVTC2+i) = *(populationVTC1+i);
        (populationVTC2+i)->id = i;
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
        choixLocus = rnd(0,nbrNoeuds-1); /// l'intervalle sera entre 1 et nbrNoeuds-2 (la borne superieur -1)
        for(j=0; j<=choixLocus; j++)
        {
            (populationVTC2+i)->genotype[j] = (populationVTC1+choixInd1)->genotype[j];
            (populationVTC2+i+1)->genotype[j] = (populationVTC1+choixInd2)->genotype[j];
        }
        ///**********************************************************************
        for(j=choixLocus+1; j<nbrNoeuds; j++)
        {
            (populationVTC2+i)->genotype[j] = (populationVTC1+choixInd2)->genotype[j];
            (populationVTC2+i+1)->genotype[j] = (populationVTC1+choixInd1)->genotype[j];
        }
        /// affectation de nouveau indice
        (populationVTC2+i)->id = i;
        (populationVTC2+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }
    generatePopulationVTC( populationVTC2, taillePopulation-regeneration);
#else
    while (i < taillePopulation)
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
        choixLocus = rnd(0,nbrNoeuds-1); /// l'intervalle sera entre 1 et nbrNoeuds-2 (la borne superieur -1)
        for(j=0; j<=choixLocus; j++)
        {
            (populationVTC2+i)->genotype[j] = (populationVTC1+choixInd1)->genotype[j];
            (populationVTC2+i+1)->genotype[j] = (populationVTC1+choixInd2)->genotype[j];
        }
        ///**********************************************************************
        for(j=choixLocus+1; j<nbrNoeuds; j++)
        {
            (populationVTC2+i)->genotype[j] = (populationVTC1+choixInd2)->genotype[j];
            (populationVTC2+i+1)->genotype[j] = (populationVTC1+choixInd1)->genotype[j];
        }
        /// affectation de nouveau indice
        (populationVTC2+i)->id = i;
        (populationVTC2+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }

#endif
}
///**************************************************************************************
int findTheBestSolutionVTC(partitionVTC *populationVTC)
{
    /// la fitness est calculé à partir des couts de coupe normalisés
    float maxFitness = populationVTC->fitness;
    int i,indice=0;
    for(i=0; i<taillePopulation ; i++)
    {
        if(maxFitness < (populationVTC+i)->fitness)
        {
            maxFitness = (populationVTC+i)->fitness;
            indice = i;
        }
    }
    return indice;
}
///**************************************************************************************
void mutationVTC(partitionVTC* populationVTC)
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
            numeroGeneMute = rnd(0,nbrParties);
            (populationVTC+i)->genotype[numeroGeneMute] = rnd(0,nbrParties);

        }
    }

}
///**************************************************************************************
float testerLaSommeDesFitnessVTC(partitionVTC* populationVTC)
{
#if mouchard
    printf("testerLaSommeDesFitness ...\n");
#endif // mouchard
    int i;
    float sommeFitness = 0;
    for(i=0; i<taillePopulation ; i++)
    {
        sommeFitness = sommeFitness + (populationVTC+i)->fitness;
    }
    return sommeFitness;
}
///************************************************************************************************************
void displayTheBestSolutionVTC(partitionVTC* solutionDominante)
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

///**********************************************************************************
void writeSolutionInFileVTC(partitionVTC *populationVTC, FILE *outputFilePop, int iteration)
{

    int i,j;
    for(i=0; i<taillePopulation; i++)
    {

        ///écriture de la meilleur solution dans le fichier output
        fprintf(outputFilePop,"%d\t%d\t%d\t%0.4f\t%d\t%d\t%d\t",iteration,(populationVTC+i)->id,
                (populationVTC+i)->coutCoupe,(populationVTC+i)->fitness,(populationVTC+i)->coutCoupeNormalise,
                (populationVTC+i)->contrainteViole,(populationVTC+i)->nbrCluster);
        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(outputFilePop,"%d ",(populationVTC+i)->genotype[j]);
        }
        fprintf(outputFilePop,"\t\t");

        for(j=0; j<nbr_constraint; j++)
        {
            fprintf(outputFilePop,"%d ",(populationVTC+i)->constraintVector[j]);
        }
        fprintf(outputFilePop,"\n");
    }
    fprintf(outputFilePop,"\n\n");

}
///**************************************************************************************
void writeBestSolutionInFileVTC(partitionVTC *solutionDominante, FILE *outputFile,int iteration)
{

    int i;
    fprintf(outputFile,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,solutionDominante->id,
            solutionDominante->coutCoupe,solutionDominante->fitness,solutionDominante->coutCoupeNormalise,
            solutionDominante->contrainteViole,solutionDominante->nbrCluster);
    for(i=0; i<nbrNoeuds; i++)
    {
        fprintf(outputFile,"%d ",solutionDominante->genotype[i]);
    }
    fprintf(outputFile,"\n");
}
///======================================================================================
void vertexToClusterEncoding(int nbrGeneration ,FILE *outputFileVTC,FILE *outputFilePopVTC,FILE *outputOptimalSolutionFileVTC,
                             partitionVTC *populationVTC1,partitionVTC *populationVTC2,partitionVTC *solutionDominante,
                             int iteration , int *bestSolutionIteration , int *nbrApparition)
{

///=====================================================================================================
        /// application des trois opérateurs de bases de génétiques
        naturalSelectionVTC(populationVTC1,populationVTC2);
        crossOverVTC(populationVTC2,populationVTC1); /// nbrNoeuds représente la taille des genotype
        mutationVTC(populationVTC1);
        checkContrainstAndFitnessPenalizationVTC(populationVTC1);
        calculCoutCoupeEtFitnessVTC(populationVTC1);
#if writingPopulationInFile
        writeSolutionInFileVTC(populationVTC1,outputFilePopVTC,iteration);
#endif // writingPopulationInFile
///=====================================================================================================
        /// dans ce qui suit, on va juste localiser la meilleur solution
///=====================================================================================================
        bestSolution=findTheBestSolutionVTC(populationVTC1);
#if writingPopulationInFile
        writeBestSolutionInFileVTC((populationVTC1+bestSolution),outputFileVTC,iteration);
#endif // writingPopulationInFile
        if((populationVTC1+bestSolution)->coutCoupeNormalise > solutionDominante->coutCoupeNormalise)
        {
            *(solutionDominante) = *(populationVTC1+bestSolution);
            *nbrApparition=1;
            *bestSolutionIteration = iteration;
        }
        else
        {
            *nbrApparition++;
        }
}


///**************************************************************************************
void checkContrainstAndFitnessPenalizationVTC(partitionVTC *populationVTC)
{

    int i,j,k;

    for(i=0; i<taillePopulation; i++)
    {
        (populationVTC+i)->nbrCluster=0; /// le nombre de cluster doit etre initialisé à 1 pour prendre en consédiration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes violées
        (populationVTC+i)->contrainteViole=0;
        for(j=0; j<4; j++)
        {
            (populationVTC+i)->constraintVector[j]=0;
        }
        ///**********************************************************************************


        ///trier le tableau résulant
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationVTC+i)->clustersSize[j] =0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if ((populationVTC+i)->genotype[k] == j)
                {
                    (populationVTC+i)->clustersSize[j]++;
                }
            }
            if ((populationVTC+i)->clustersSize[j] !=0)
            {
                (populationVTC+i)->nbrCluster++;
            }
        }

        if((populationVTC+i)->nbrCluster > max_clusters || (populationVTC+i)->nbrCluster < min_clusters)
        {
            (populationVTC+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationVTC+i)->contrainteViole++;
        }

        for(j=0; j<nbrNoeuds; j++)
        {
            if((populationVTC+i)->clustersSize[j]!=0)
            {
                if((populationVTC+i)->clustersSize[j]>max_sizeCluster || (populationVTC+i)->clustersSize[j]< min_sizeCluster)
                {
                    (populationVTC+i)->constraintVector[1]++;

                }
            }

        }

        if((populationVTC+i)->constraintVector[1] != 0)
        {
            (populationVTC+i)->contrainteViole++;
        }
    }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0)
    {
        checkCohabitationAndNonCohabitationConstrainteVTC(populationVTC);
    }

}
///**************************************************************************************
void checkCohabitationAndNonCohabitationConstrainteVTC(partitionVTC *populationVTC)
{

    int i,j,noeud1,noeud2;
    for(i=0; i<taillePopulation; i++)
    {
        /// Cette tache est déjà réalisé au niveau de la fonction checkCohabitationAndNonCohabitationConstrainteBE
		/**
        (populationVTC+i)->constraintVector[2]=0;
        (populationVTC+i)->constraintVector[3]=0;
		*/
        if(nbrCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrCohabitationConstraintes; j++)
            {
                noeud1=cohabitationConstraintes[j].nouedDepart;
                noeud2=cohabitationConstraintes[j].nouedArrive;

                if((populationVTC+i)->genotype[noeud1]!= (populationVTC+i)->genotype[noeud2])
                {
                    (populationVTC+i)->constraintVector[2]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationVTC+i)->constraintVector[2]!=0)
            {
                (populationVTC+i)->contrainteViole++;
            }
        }


        if(nbrNonCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrNonCohabitationConstraintes; j++)
            {
                noeud1=nonCohabitationConstraintes[j].nouedDepart;
                noeud2=nonCohabitationConstraintes[j].nouedArrive;
                if((populationVTC+i)->genotype[noeud1]== (populationVTC+i)->genotype[noeud2])
                {
                    (populationVTC+i)->constraintVector[3]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationVTC+i)->constraintVector[3]!=0)
            {
                (populationVTC+i)->contrainteViole++;
            }
        }


    }

}

///********************************************************************************************
void writeOptimalSolutionInFileVTC(partitionVTC *solutionDominante,FILE* outputOptimalSolutionFileVTC,
                                  int nbrRun, int bestSolutionIteration, float runTime, int ES){

    fprintf(outputOptimalSolutionFileVTC,"RUN NUMBER = %d | ITERATION N = %d | SOLUTION IDENTIFIER = %d | CUT SIZE INTRA = %d | CUT SIZE INTER = %d | SOMME TOTAL DES FLUX = %d | RUN TIME = %f | ES = %d | NBR VIOLATED CONSTRAINTE = %d | ",
            nbrRun, bestSolutionIteration, solutionDominante->id, solutionDominante->coutCoupe,(sommeTotalFlux - solutionDominante->coutCoupe),sommeTotalFlux ,runTime, ES, solutionDominante->contrainteViole);
fprintf(outputOptimalSolutionFileVTC," max_clusters = %d | max_sizeCluster = %d  |  cohabitation = %d | non cohabitation = %d\n",
            solutionDominante->constraintVector[0],solutionDominante->constraintVector[1],solutionDominante->constraintVector[2],
            solutionDominante->constraintVector[3]);

}
///******************************************************************************************
int compareCroissantFitnessVTC (void const *a, void const *b)
{

    partitionVTC const *pa = a;
    partitionVTC const *pb = b;
    return  pb->coutCoupeNormalise - pa->coutCoupeNormalise;
}
/**************************************************************************************
void writeOptimalSolutionInFileVTC(partitionVTC *solutionDominante,FILE* outputOptimalSolutionFileVTC)
{
    int i;
    fprintf(outputOptimalSolutionFileVTC,"id = %d | normalised cut size = %d | cut size intrA = %d cut size intEr = %d ",
            solutionDominante->id,solutionDominante->coutCoupeNormalise, solutionDominante->coutCoupe, (sommeTotalFlux-solutionDominante->coutCoupe));
    fprintf(outputOptimalSolutionFileVTC,"le contraintes violees = %d \n", solutionDominante->contrainteViole);
    for(i=0; i<4; i++)
    {
        switch (i)
        {

        case 0:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileVTC,"contrainte sur le nombre de clusters : OK \n");
            else fprintf(outputOptimalSolutionFileVTC,"contrainte sur le nombre de clusters : NOT OK \n");
            break;
        case 1:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileVTC,"contrainte sur la taille des clusters : OK \n");
            else fprintf(outputOptimalSolutionFileVTC,"contrainte sur la taille des clusters : NOT OK \n");
            break;
        case 2:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileVTC,"contrainte de cohabitation : OK \n");
            else fprintf(outputOptimalSolutionFileVTC,"contrainte de cohabitation : NOT OK \n");
            break;
        case 3:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileVTC,"contrainte Non Cohabitation : OK \n");
            else fprintf(outputOptimalSolutionFileVTC,"contrainte Non Cohabitation : NOT OK \n");
            break;
        default:
            fprintf(outputOptimalSolutionFileVTC,"Cette contrainte n'est pas prise en charge par le système \n");

        }

    }

}

void checkContrainstAndFitnessPenalizationVTC(partitionVTC *populationVTC)
{

    int i,j,k,m=nbrNoeuds,tabTmp[m],tmp;

    for(i=0; i<taillePopulation; i++)
    {
        (populationVTC+i)->nbrCluster=1; /// le nombre de cluster doit etre initialisé à 1 pour prendre en consédiration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes violées
        (populationVTC+i)->contrainteViole=0;
        for(j=0;j<nbr_constraint;j++){(populationVTC+i)->constraintVector[j]=0;}


        m=nbrNoeuds;
        ///copier le genotype dans un deuxième tableau
        for(j=0; j<nbrNoeuds; j++)
        {
            tabTmp[j]=(populationVTC+i)->genotype[j];
        }
        ///trier le tableau résulant
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
                (populationVTC+i)->nbrCluster++;
            }
        }
        if((populationVTC+i)->nbrCluster > max_clusters || (populationVTC+i)->nbrCluster < min_clusters)
        {
            (populationVTC+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationVTC+i)->contrainteViole++;
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
            (populationVTC+i)->clustersSize[j]=0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if (tabTmp[j] == (populationVTC+i)->genotype[k])
                {
                    (populationVTC+i)->clustersSize[j]=(populationVTC+i)->clustersSize[j]+1;
                }
            }
            if((populationVTC+i)->clustersSize[j]>max_sizeCluster || (populationVTC+i)->clustersSize[j]< min_sizeCluster)
            {
                ///24.07.2015 ---> l'ajout de vecteur de contraintes = le nombre de cluster dont la taille a dépassé la borne
                (populationVTC+i)->constraintVector[1]++;

            }
        }

        if((populationVTC+i)->constraintVector[1] != 0)
        {
            (populationVTC+i)->contrainteViole++;
        }
   }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0){
        checkCohabitationAndNonCohabitationConstrainteVTC(populationVTC);
    }

}
*/
