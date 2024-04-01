#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "typeDeclaration.h"
#include "./communesFunctions.h"
#include "./geneticClusterEncoding.h"
#include "./compilationConditionnelle.h"


extern int taillePopulation,tauxReproduction,sommeTotalFlux,nbrNoeuds,nbr_constraint,nbrArretes,nbrParties,nbrNoeudsInCluster;
extern double tauxMutation ,averageTraffic;
extern int fluxMatrix[1000][1000],*fluxVector, maxFlow, sommeTotalFluxReal;
extern edge *edgeVector,cohabitationConstraintes[1000],nonCohabitationConstraintes[1000];
extern int max_clusters , min_clusters,min_sizeCluster ,max_sizeCluster ;
extern int bestSolution, nbrApparition, nbrCohabitationConstraintes,nbrNonCohabitationConstraintes;
extern int nbrRun;

///**************************************************************************************
///génération de la population initiale
void generatePopulationGCE(partitionGCE* populationGCE)
{

    int i,j,assignmentVector[nbrNoeuds],allele,k,tmp;


    for(i=0; i<taillePopulation; i++)
    {
        for(j=0;j<nbrNoeuds;j++) assignmentVector[j]=-1;
        for(j=0; j<nbrNoeuds; j++)
        {
            do{
                allele = rnd(0,nbrNoeuds);
            }while(assignmentVector[allele] == 0);
            (populationGCE+i)->genotype[j] = allele; /// [0;nbrNoeuds-1]
             assignmentVector[allele] = 0;
        }
        /// affectation des sommets au cluster
        for(j=0;j<nbrNoeuds;j++) assignmentVector[j]=-1;
        (populationGCE+i)->nbrCluster = rnd(1,nbrNoeuds);

        for(j=0;j<(populationGCE+i)->nbrCluster-1; j++){
            do{
                allele = rnd(0,nbrNoeuds);
            }while(assignmentVector[allele] == 0);
            (populationGCE+i)->delimiterVector[j] = allele; /// [0;nbrNoeuds-1]
             assignmentVector[allele] = 0;
        }

        for(j=0;j< (populationGCE+i)->nbrCluster-1; j++){
                for(k=j+1;k< (populationGCE+i)->nbrCluster ; k++){
                    if((populationGCE+i)->delimiterVector[j] >  (populationGCE+i)->delimiterVector[k]){{
                        tmp = (populationGCE+i)->delimiterVector[j];
                        (populationGCE+i)->delimiterVector[j] = (populationGCE+i)->delimiterVector[k];
                        (populationGCE+i)->delimiterVector[k] =tmp;
                    }
                }
        }

        (populationGCE+i)->id =i;
    }
        getPhenotypeFromGenotype(populationGCE);
    }
}
///**************************************************************************************
void getPhenotypeFromGenotype(partitionGCE* populationGCE){

    int i,j,k,indiceVertex;
    for(i=0;i<taillePopulation;i++){
        indiceVertex = 0;
        for(j=0;j<(populationGCE+i)->nbrCluster-1;j++){
            for(k=0;k<(populationGCE+i)->delimiterVector[j];k++){
                (populationGCE+i)->phenotype[indiceVertex] = j;
                indiceVertex++;
            }
        }
    }
}
///**************************************************************************************
/// affichage des population
void affichePopulationGCE(partitionGCE* populationGCE)
{
    int i,j;
    for(i=0; i<taillePopulation; i++)
    {
        printf("\n id = %d \t cout de coupe = %d \t fitness = %0.4f \t",
               (populationGCE+i)->id,(populationGCE+i)->coutCoupe, (populationGCE+i)->fitness);
        printf("le cout de coupe normalise = %d \n",(populationGCE+i)->coutCoupeNormalise);
        printf("la partition est  = \t");
        for(j=0; j<nbrNoeuds; j++)
        {
            printf("%d ",(populationGCE+i)->genotype[j]);
        }
        printf("\n");
        printf("le nombre des clusters est %d \n",(populationGCE+i)->nbrCluster);
        printf("le nombre des contraintes violees est  = %d \n",(populationGCE+i)->contrainteViole);
        printf("affichage des tailles des clusters \n");
        for(j=0; j<(populationGCE+i)->nbrCluster; j++)
        {
            printf("%d ",(populationGCE+i)->clustersSize[j]);
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
void calculCoutCoupeEtFitnessGCE(partitionGCE* populationGCE)
{
    int i,j,k,sommeTotalCoutCoupe = 0;

    for(i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (populationGCE+i)->coutCoupe = 0; /// réinitialisation de la valriable après chaque itération (individus de la population)
        for(j=0; j<nbrNoeuds-1; j++)
        {
            for(k=j+1; k<nbrNoeuds; k++)
            {
                /// on peut avoir deux noeuds du même cluster mais qui sont par relier par une arrêtes
                if((populationGCE+i)->genotype[j] == (populationGCE+i)->genotype[k] && fluxMatrix[j][k] >= 0)
                {
                    (populationGCE+i)->coutCoupe = (populationGCE+i)->coutCoupe +  fluxMatrix[j][k];
                }
            }
        }
        (populationGCE+i)->coutCoupeNormalise = (populationGCE+i)->coutCoupe + ((nbr_constraint - (populationGCE+i)->contrainteViole)*sommeTotalFlux);
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationGCE+i)->coutCoupeNormalise;
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
            varianceCoutDeCoupeNormalise + pow(((populationGCE+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);
    /// tester la valeur ajouté
    for(i=0; i<taillePopulation; i++)
    {

        /// calcul de sigma truncation
        (populationGCE+i)->expectedValue = (populationGCE+i)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationGCE+i)->expectedValue;
    }

    /// calcule du fitness
    /// le 26/11/2015 : calcul du fitness des solution avec leur expectedValues
    for(i=0; i<taillePopulation ; i++)
    {
        (populationGCE+i)->fitness = (float)((populationGCE+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);
    }

#else
    FILE *file;
    file = fopen("D:/WorkingFolderAgProgramming/benchmark/vertexToClusterEncoding/expectedValueGCE.txt","w");
    for(i=0; i<taillePopulation ; i++)
    {
/**        if(iteration == 1)
        {
            fprintf(file,"%d \t %d \t %0.2f \n",i,(populationGCE+i)->coutCoupeNormalise,(populationGCE+i)->expectedValue);
        }
*/
        (populationGCE+i)->fitness = (float)((populationGCE+i)->coutCoupeNormalise)/(float)(sommeTotalCoutCoupe);
    }
#endif// scaling
}
///**************************************************************************************
void naturalSelectionGCE(partitionGCE* populationGCE1,partitionGCE* populationGCE2)
{
    int i,j,maxFitness;
    int tmpVectIndice[taillePopulation];
    for(i=0;i<taillePopulation;i++) tmpVectIndice[i]=-1;
    float lotterie,sommeFitness;

    for(i=0; i<tauxReproduction; i++){
        maxFitness = 0;
        for(j=0;j<taillePopulation;j++){
            if((populationGCE1+maxFitness)->coutCoupeNormalise < (populationGCE1+j)->coutCoupeNormalise
               && tmpVectIndice[j] == -1 ){
                maxFitness = j;

            }
        }
        tmpVectIndice[maxFitness] = 0;
        *(populationGCE2+i) = *(populationGCE1+maxFitness);
        ///printf("(populationEE2+%d)->coutCoupeNormalise = %d \n",i,(populationEE2+i)->coutCoupeNormalise);
        ///mon_sleep(pause);
    }
    ///printf("\n\n\n");

    ///qsort (populationEE2, taillePopulation, sizeof *populationEE2, compareCroissantFitnessEE);
    ///sortingPopulationByFintness(populationEE2,tmpPopulation);

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
            sommeFitness = sommeFitness + (populationGCE1+j)->fitness;
            if (lotterie<=sommeFitness)
            {
                *(populationGCE2+i) = *(populationGCE1+j);
                break;
            }
            else if( j==taillePopulation-1)
            {
                *(populationGCE2+i) = *(populationGCE1+j);
            }
        }
    }
    ///free(tmpPopulation);

}
///**************************************************************************************
void crossOverGCE(partitionGCE* populationGCE1, partitionGCE* populationGCE2)
{
    int i = 0,j =0;
    int choixInd1, choixInd2, choixLocus;
    ///***********************************************************************************
    for(i=0; i<tauxReproduction; i++)
    {
        *(populationGCE2+i) = *(populationGCE1+i);
        (populationGCE2+i)->id = i;
    }
    ///***********************************************************************************


    while (i < taillePopulation)
    {
        choixInd1 = rnd(tauxReproduction,taillePopulation);
        do
        {
            /// ce traitement évite l'appariement d'un element avec lui même/
            choixInd2 = rnd(tauxReproduction,taillePopulation);
        }
        while(choixInd2 == choixInd1 );
        choixLocus = rnd(0,nbrNoeuds-1); /// l'intervalle sera entre 1 et nbrNoeuds-2 (la borne superieur -1)
        for(j=0; j<=choixLocus; j++)
        {
            (populationGCE2+i)->genotype[j] = (populationGCE1+choixInd1)->genotype[j];
            (populationGCE2+i+1)->genotype[j] = (populationGCE1+choixInd2)->genotype[j];
        }
        ///**********************************************************************
        for(j=choixLocus+1; j<nbrNoeuds; j++)
        {
            (populationGCE2+i)->genotype[j] = (populationGCE1+choixInd2)->genotype[j];
            (populationGCE2+i+1)->genotype[j] = (populationGCE1+choixInd1)->genotype[j];
        }
        /// affectation de nouveau indice
        (populationGCE2+i)->id = i;
        (populationGCE2+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }
}
///**************************************************************************************
int findTheBestSolutionGCE(partitionGCE *populationGCE)
{
    /// la fitness est calculé à partir des couts de coupe normalisés
    float maxFitness = populationGCE->fitness;
    int i,indice=0;
    for(i=0; i<taillePopulation ; i++)
    {
        if(maxFitness < (populationGCE+i)->fitness)
        {
            maxFitness = (populationGCE+i)->fitness;
            indice = i;
        }
    }
    return indice;
}
///**************************************************************************************
void mutationGCE(partitionGCE* populationGCE)
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
            (populationGCE+i)->genotype[numeroGeneMute] = rnd(0,nbrParties);

        }
    }

}
///**************************************************************************************
float testerLaSommeDesFitnessGCE(partitionGCE* populationGCE)
{
#if mouchard
    printf("testerLaSommeDesFitness ...\n");
#endif // mouchard
    int i;
    float sommeFitness = 0;
    for(i=0; i<taillePopulation ; i++)
    {
        sommeFitness = sommeFitness + (populationGCE+i)->fitness;
    }
    return sommeFitness;
}
///************************************************************************************************************
void displayTheBestSolutionGCE(partitionGCE* solutionDominante)
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
void writeSolutionInFileGCE(partitionGCE *populationGCE, FILE *outputFilePop, int iteration)
{

    int i,j;
    for(i=0; i<taillePopulation; i++)
    {

        ///écriture de la meilleur solution dans le fichier output
        fprintf(outputFilePop,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,(populationGCE+i)->id,
                (populationGCE+i)->coutCoupe,(populationGCE+i)->fitness,(populationGCE+i)->coutCoupeNormalise,
                (populationGCE+i)->contrainteViole,(populationGCE+i)->nbrCluster);
        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(outputFilePop,"%d ",(populationGCE+i)->genotype[j]);
        }
        fprintf(outputFilePop,"\t\t");

        for(j=0; j<nbr_constraint; j++)
        {
            fprintf(outputFilePop,"%d ",(populationGCE+i)->constraintVector[j]);
        }
        fprintf(outputFilePop,"\n");
    }
    fprintf(outputFilePop,"\n\n");

}
///**************************************************************************************
void writeBestSolutionInFileGCE(partitionGCE *solutionDominante, FILE *outputFile,int iteration)
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
void geneticClusterEncoding(int nbrGeneration ,FILE *outputFileGCE,FILE *outputFilePopGCE,FILE *outputOptimalSolutionFileGCE, partitionGCE *populationGCE1,partitionGCE *populationGCE2,partitionGCE *solutionDominante)
{

    int bestSolution,iteration=1 ,bestSolutionIteration,ES = 0;
    clock_t t1,t2;
    double temps;

    t1=clock();


    generatePopulationGCE(populationGCE1);
    checkContrainstAndFitnessPenalizationGCE(populationGCE1);
    calculCoutCoupeEtFitnessGCE(populationGCE1);
#if writingPopulationInFile
    writeSolutionInFileGCE(populationGCE1,outputFilePopGCE,iteration);
#endif // writingPopulationInFile


    bestSolution=findTheBestSolutionGCE(populationGCE1);
    *solutionDominante=*(populationGCE1+bestSolution);
    bestSolutionIteration= 1;
#if writingPopulationInFile
    writeBestSolutionInFileGCE(solutionDominante,outputFileGCE,iteration);
#endif // writingPopulationInFile
    nbrApparition=1;
    for(iteration =2 ; iteration <= nbrGeneration; iteration++){
///=====================================================================================================
        /// application des trois opérateurs de bases de génétiques
        naturalSelectionGCE(populationGCE1,populationGCE2);
        crossOverGCE(populationGCE2,populationGCE1); /// nbrNoeuds représente la taille des genotype
        mutationGCE(populationGCE1);
        ///*************************************************************
        //getPhenotypeFromGenotype(populationGCE);
        ///*************************************************************
        checkContrainstAndFitnessPenalizationGCE(populationGCE1);
        calculCoutCoupeEtFitnessGCE(populationGCE1);
#if writingPopulationInFile
        writeSolutionInFileGCE(populationGCE1,outputFilePopGCE,iteration);
#endif // writingPopulationInFile
///=====================================================================================================
        /// dans ce qui suit, on va juste localiser la meilleur solution
///=====================================================================================================
        bestSolution=findTheBestSolutionGCE(populationGCE1);
#if writingPopulationInFile
        writeBestSolutionInFileGCE((populationGCE1+bestSolution),outputFileGCE,iteration);
#endif // writingPopulationInFile
        if((populationGCE1+bestSolution)->coutCoupeNormalise > solutionDominante->coutCoupeNormalise)
        {
            *(solutionDominante) = *(populationGCE1+bestSolution);
            nbrApparition=1;
            bestSolutionIteration = iteration;
        }
        else
        {
            nbrApparition++;
        }

    }

    t2=clock();

    ///displayTheBestSolutionGCE(solutionDominante);

    fprintf(outputFileGCE,"\n\n");
    writeBestSolutionInFileGCE(solutionDominante,outputFileGCE,bestSolutionIteration);
    fprintf(outputFileGCE,"\n\nLe nombre d'apparition de la meilleur solution est  = %d\n",nbrApparition);
    fprintf(outputFileGCE,"\n\n==============================================================================\n");
    fprintf(outputFilePopGCE,"\n\n==============================================================================\n");
    temps = (double)(t2-t1)/CLOCKS_PER_SEC;

    if(bestSolutionIteration>=2){
        ES = ((taillePopulation-tauxReproduction)*(bestSolutionIteration-2))+(taillePopulation+ solutionDominante->id -tauxReproduction+1);
    }
    else {
        ES = solutionDominante->id +1;
    }

    writeOptimalSolutionInFileGCE(solutionDominante,outputOptimalSolutionFileGCE,nbrRun,bestSolutionIteration,temps, ES);
    printf("le temps d execution est %lf \n", temps);



}

///**************************************************************************************
void checkContrainstAndFitnessPenalizationGCE(partitionGCE *populationGCE)
{

    int i,j,k;

    for(i=0; i<taillePopulation; i++)
    {
        (populationGCE+i)->nbrCluster=0; /// le nombre de cluster doit etre initialisé à 1 pour prendre en consédiration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes violées
        (populationGCE+i)->contrainteViole=0;
        for(j=0; j<4; j++)
        {
            (populationGCE+i)->constraintVector[j]=0;
        }
        ///**********************************************************************************


        ///trier le tableau résulant
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationGCE+i)->clustersSize[j] =0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if ((populationGCE+i)->genotype[k] == j)
                {
                    (populationGCE+i)->clustersSize[j]++;
                }
            }
            if ((populationGCE+i)->clustersSize[j] !=0)
            {
                (populationGCE+i)->nbrCluster++;
            }
        }

        if((populationGCE+i)->nbrCluster > max_clusters || (populationGCE+i)->nbrCluster < min_clusters)
        {
            (populationGCE+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationGCE+i)->contrainteViole++;
        }

        for(j=0; j<nbrNoeuds; j++)
        {
            if((populationGCE+i)->clustersSize[j]!=0)
            {
                if((populationGCE+i)->clustersSize[j]>max_sizeCluster || (populationGCE+i)->clustersSize[j]< min_sizeCluster)
                {
                    (populationGCE+i)->constraintVector[1]++;

                }
            }

        }

        if((populationGCE+i)->constraintVector[1] != 0)
        {
            (populationGCE+i)->contrainteViole++;
        }
    }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0)
    {
        checkCohabitationAndNonCohabitationConstrainteGCE(populationGCE);
    }

}
///**************************************************************************************
void checkCohabitationAndNonCohabitationConstrainteGCE(partitionGCE *populationGCE)
{

    int i,j,noeud1,noeud2;
    for(i=0; i<taillePopulation; i++)
    {
        /// Cette tache est déjà réalisé au niveau de la fonction checkCohabitationAndNonCohabitationConstrainteBE
		/**
        (populationGCE+i)->constraintVector[2]=0;
        (populationGCE+i)->constraintVector[3]=0;
		*/
        if(nbrCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrCohabitationConstraintes; j++)
            {
                noeud1=cohabitationConstraintes[j].nouedDepart;
                noeud2=cohabitationConstraintes[j].nouedArrive;

                if((populationGCE+i)->genotype[noeud1]!= (populationGCE+i)->genotype[noeud2])
                {
                    (populationGCE+i)->constraintVector[2]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationGCE+i)->constraintVector[2]!=0)
            {
                (populationGCE+i)->contrainteViole++;
            }
        }


        if(nbrNonCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrNonCohabitationConstraintes; j++)
            {
                noeud1=nonCohabitationConstraintes[j].nouedDepart;
                noeud2=nonCohabitationConstraintes[j].nouedArrive;
                if((populationGCE+i)->genotype[noeud1]== (populationGCE+i)->genotype[noeud2])
                {
                    (populationGCE+i)->constraintVector[3]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationGCE+i)->constraintVector[3]!=0)
            {
                (populationGCE+i)->contrainteViole++;
            }
        }


    }

}

///********************************************************************************************
void writeOptimalSolutionInFileGCE(partitionGCE *solutionDominante,FILE* outputOptimalSolutionFileGCE,
                                  int nbrRun, int bestSolutionIteration, float runTime, int ES){

    fprintf(outputOptimalSolutionFileGCE,"RUN NUMBER = %d | ITERATION N = %d | SOLUTION IDENTIFIER = %d | CUT SIZE INTRA = %d | CUT SIZE INTER = %d | SOMME TOTAL DES FLUX = %d | RUN TIME = %f | ES = %d | NBR VIOLATED CONSTRAINTE = %d \n",
            nbrRun, bestSolutionIteration, solutionDominante->id, solutionDominante->coutCoupe,(sommeTotalFlux - solutionDominante->coutCoupe),sommeTotalFlux ,runTime, ES, solutionDominante->contrainteViole);

}
///******************************************************************************************
int compareCroissantFitnessGCE (void const *a, void const *b)
{

    partitionGCE const *pa = a;
    partitionGCE const *pb = b;
    return  pb->coutCoupeNormalise - pa->coutCoupeNormalise;
}
