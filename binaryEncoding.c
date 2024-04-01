#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "typeDeclaration.h"
#include "./communesFunctions.h"
#include "./binaryEncoding.h"
#include "./compilationConditionnelle.h"

extern int taillePopulation,tauxReproduction,sommeTotalFlux,nbrNoeuds,nbr_constraint,nbrArretes,nbrParties,
nbrNoeudsInCluster,regeneration;
extern double tauxMutation;
extern int fluxMatrix[1000][1000],*fluxVector,maximumDistanceMatrix[1000][1000], maxFlow, sommeTotalFluxReal;
extern edge *edgeVector,cohabitationConstraintes[100],nonCohabitationConstraintes[1000],  *edgeVectorIntra;
extern int max_clusters , min_clusters,min_sizeCluster ,max_sizeCluster ;
extern int bestSolution,nbrCohabitationConstraintes,nbrNonCohabitationConstraintes;
extern int nbrRun;



///************************************************************************************************************
void generatePopulationEA(partitionEA* populationEA, int indiceFirstElt)
{

    /// la probabilité P1 correspond aux arrêtes intraClusters tandis que P2 correspond
    /// aux arrêtes interClusters tel que P1 = (m-1)/|E| et P2= 1-P1
    /// m : nombre de noeuds et |E|: cardinalité de l'ensemble E = nombre d'arrêtes dans le graphe
    float P1=0.3;
    int i,j,k,allele,nbrZeroArretes,nbrOneArretes;

    nbrZeroArretes = ceil(P1*nbrArretes);
    ///printf("Le nombre d'arretes total = %d le nombre d'arrêtes a zero est : %d\n",nbrArretes,nbrZeroArretes);
    for(i=indiceFirstElt; i<taillePopulation; i++)
    {
        (populationEA+i)->genotype= (char*)malloc(nbrArretes*sizeof(char));
        if((populationEA+i)->genotype == NULL){fprintf(stderr,"Memory allocation failled for genotypeEA\n");exit(EXIT_FAILURE);}
        for(j=0; j<nbrArretes; j++)
        {
            (populationEA+i)->genotype[j] = 1;
        }
        for(j=0; j<nbrZeroArretes; j++)
        {
            do
            {
                allele = rnd(0,nbrArretes);
            }
            while((populationEA+i)->genotype[allele]==0);
            (populationEA+i)->genotype[allele] = 0;
        }
        (populationEA+i)->id =i;

        for(j=0; j<nbrCohabitationConstraintes; j++)
            {
                for(k=0; k<nbrArretes; k++)
                {
                    if((cohabitationConstraintes[j].nouedDepart == edgeVector[k].nouedDepart &&
                            cohabitationConstraintes[j].nouedArrive == edgeVector[k].nouedArrive) ||
                            (cohabitationConstraintes[j].nouedDepart == edgeVector[k].nouedArrive &&
                            cohabitationConstraintes[j].nouedArrive == edgeVector[k].nouedDepart))
                    {
                        (populationEA+i)->genotype[k] = 0;
                    }
                }
            }
         for(j=0; j<nbrNonCohabitationConstraintes; j++)
            {
                for(k=0; k<nbrArretes; k++)
                {
                    if((nonCohabitationConstraintes[j].nouedDepart == edgeVector[k].nouedDepart &&
                            nonCohabitationConstraintes[j].nouedArrive == edgeVector[k].nouedArrive) ||
                       (nonCohabitationConstraintes[j].nouedDepart == edgeVector[k].nouedArrive &&
                            nonCohabitationConstraintes[j].nouedArrive == edgeVector[k].nouedDepart) )
                    {
                        (populationEA+i)->genotype[k] = 1;
                    }
                }
            }

    }
}
///***************************************************************************************************
void generatePopulationRandomlyEA(partitionEA* populationEA)
{
    int i,j;
    for(i=0; i<taillePopulation; i++)
    {
        for(j=0; j<nbrArretes; j++)
        {
            (populationEA+i)->genotype[j] = rnd(0,2);
        }
        (populationEA+i)->id =i;
    }
}
///************************************************************************************************************
void affichePopulationEA(partitionEA* populationEA)
{

    printf("affichePopulation ...\n");
    int i,j;
    for(i=0; i<taillePopulation; i++)
    {
        printf("affichage du %d ieme individu \n",i);
        for(j=0; j<nbrArretes; j++)
        {
            printf("%d\t",(populationEA+i)->genotype[j]);
        }
        printf("\n");
        ///*********************************************************************
        for(j=0; j<nbrArretes; j++)
        {
            printf("%d %d|",edgeVector[j].nouedDepart,edgeVector[j].nouedArrive);
        }

        printf("\nid = %d \t cout de coupe = %d \t fitness = %0.4f \n",
               (populationEA+i)->id,(populationEA+i)->coutCoupeNormalise, (populationEA+i)->fitness);
        printf("la partition est : \n");
        for(j=0; j<nbrNoeuds; j++)
        {
            printf("%d ",(populationEA+i)->phenotype[j]);
        }
        printf("\n le nombre de contrainte viole = %d \n",(populationEA+i)->contrainteViole);


        printf("\n\n_____________________________________________________________________________\n\n");
    }
}
///************************************************************************************************************
void naturalSelectionEA(partitionEA* populationEA1,partitionEA* populationEA2)
{
    int i,j,maxFitness;
    int tmpVectIndice[taillePopulation];
    for(i=0;i<taillePopulation;i++) tmpVectIndice[i]=-1;
    float lotterie,sommeFitness;

    for(i=0; i<tauxReproduction; i++){
        maxFitness = 0;
        for(j=0;j<taillePopulation;j++){
            if((populationEA1+maxFitness)->coutCoupeNormalise < (populationEA1+j)->coutCoupeNormalise
               && tmpVectIndice[j] == -1 ){
                maxFitness = j;

            }
        }
        tmpVectIndice[maxFitness] = 0;
        *(populationEA2+i) = *(populationEA1+maxFitness);
        ///printf("(populationEA2+%d)->coutCoupeNormalise = %d \n",i,(populationEA2+i)->coutCoupeNormalise);
        ///mon_sleep(pause);
    }
    ///printf("\n\n\n");
/*
    for(i=0; i<taillePopulation; i++)*(populationEA2+i) = *(populationEA1+i);
    qsort (populationEA2, taillePopulation, sizeof *populationEA2, compareCroissantFitnessEA);
*/
    ///sortingPopulationByFintness(populationEA2,tmpPopulation);

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
            sommeFitness = sommeFitness + (populationEA1+j)->fitness;
            if (lotterie<=sommeFitness)
            {
                *(populationEA2+i) = *(populationEA1+j);
                break;
            }
            else if( j==taillePopulation-1)
            {
                *(populationEA2+i) = *(populationEA1+j);
            }
        }
    }
}
///************************************************************************************************************
void crossOverEA(partitionEA* populationEA1, partitionEA* populationEA2)
{

    int i ,j;
    int choixInd1, choixInd2, choixLocus;

    ///***********************************************************************************
    for(i=0; i<tauxReproduction; i++)
    {
        *(populationEA1+i) = *(populationEA2+i);
        (populationEA1+i)->id = i;
    }
    ///***********************************************************************************

#if NEW_CROSSOVER

    while (i < (taillePopulation-regeneration))
    {

        choixInd1 = rnd(0,taillePopulation);
        ///choixInd1 = rnd(tauxReproduction,taillePopulation);
        do
        {
            /// ce traitement évite l'appariement d'un element avec lui même/
            choixInd2 = rnd(0,taillePopulation);
            ///choixInd2 = rnd(tauxReproduction,taillePopulation);
        }
        while(choixInd2 == choixInd1 );


        choixLocus = rnd(0,nbrArretes-2); /// -2 parce qu'on a commencer à partir du 0

        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationEA1+i)->genotype[j] = (populationEA2+choixInd1)->genotype[j];
            (populationEA1+i+1)->genotype[j] = (populationEA2+choixInd2)->genotype[j];
        }
        for(j=choixLocus+1; j<nbrArretes; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationEA1+i)->genotype[j] = (populationEA2+choixInd2)->genotype[j];
            (populationEA1+i+1)->genotype[j] = (populationEA2+choixInd1)->genotype[j];
        }

        /// affectation de nouveau indice
        (populationEA1+i)->id = i;
        (populationEA1+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }

     generatePopulationEA(populationEA1, taillePopulation-regeneration+1);

#else
    while (i < taillePopulation)
    {
        choixInd1 = rnd(0,taillePopulation);
        ///choixInd1 = rnd(tauxReproduction,taillePopulation);

        do
        {
            /// ce traitement évite l'appariement d'un element avec lui même/
            choixInd2 = rnd(0,taillePopulation);
            ///choixInd2 = rnd(tauxReproduction,taillePopulation);
        }
        while(choixInd2 == choixInd1 );

        choixLocus = rnd(0,nbrArretes-2); /// -2 parce qu'on a commencer à partir du 0

        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationEA1+i)->genotype[j] = (populationEA2+choixInd1)->genotype[j];
            (populationEA1+i+1)->genotype[j] = (populationEA2+choixInd2)->genotype[j];
        }
        for(j=choixLocus+1; j<nbrArretes; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationEA1+i)->genotype[j] = (populationEA2+choixInd2)->genotype[j];
            (populationEA1+i+1)->genotype[j] = (populationEA2+choixInd1)->genotype[j];
        }
        /// affectation de nouveau indice
        (populationEA1+i)->id = i;
        (populationEA1+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }
#endif// NEW_CROSSOVER
}
///************************************************************************************************************
int findTheBestSolutionEA(partitionEA *populationEA)
{

    float maxFitness = (populationEA+0)->fitness;
    int i,indice=0;
    for(i=0; i<taillePopulation ; i++)
    {
        if(maxFitness < (populationEA+i)->fitness)
        {
            maxFitness = (populationEA+i)->fitness;
            indice = i;
        }
    }
    return indice;

}
///************************************************************************************************************
void mutationEA(partitionEA* populationEA)
{
    int i,numeroGeneMute;
    double applicationOfMutation;
    for(i=tauxReproduction; i<taillePopulation; i++)
    {
#if Windows
        applicationOfMutation = drand48ForWindows(0,1001);
#else
        applicationOfMutation = drand48();
#endif
        if(applicationOfMutation <= tauxMutation)
        {
           /// printf("mutation appliquee, i=%d\n",i);
            /// choisir un gène aléatoirement pour lui appliquer la mutation
            numeroGeneMute = rnd(0,nbrArretes);
            if((populationEA+i)->genotype[numeroGeneMute] == 1)
            {
                (populationEA+i)->genotype[numeroGeneMute] =0;
            }
            else
            {
                (populationEA+i)->genotype[numeroGeneMute] =1;
            }
        }
    }
}
///************************************************************************************************************
float testerLaSommeDesFitnessEA(partitionEA* populationEA)
{
    int i;
    float sommeFitness = 0;
    for(i=0; i<taillePopulation ; i++)
    {
        sommeFitness = sommeFitness + (populationEA+i)->fitness;
    }
    return sommeFitness;
}
///************************************************************************************************************
void displayTheBestSolutionEA(partitionEA* solutionDominante)
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
void getPartitionFromSolutionEA(partitionEA *populationEA)
{

    int i,j,k,nd,na,numeroClustre,nbrIntra; /// nd = noued de départ , na = noued d'arrivé

    for(i=0; i<taillePopulation; i++)
    {
        /// initialisation des phénotype = partition par les numero des noeuds correspondant

        for(j=0; j<nbrNoeuds; j++)
        {
            (populationEA+i)->phenotype[j] = j;
        }
        /// 01/10/2017 : amélioration du temps d'exécution du codage binaire
        nbrIntra = 0;
        for(j=0; j<nbrArretes; j++){
            if((populationEA+i)->genotype[j]==0)
            {
               edgeVectorIntra[nbrIntra] = edgeVector[j];
               nbrIntra++;
            }

        }
        /// détermination de la parition
        for(j=0; j<nbrIntra; j++){

                nd = (edgeVectorIntra+j)->nouedDepart;
                na = (edgeVectorIntra+j)->nouedArrive;
               if((populationEA+i)->phenotype[na]!= na) /// il a été modifier auparavant
                {

                    if((populationEA+i)->phenotype[na] < (populationEA+i)->phenotype[nd])
                    {
                        ///dans le cas ou le noued d'arrivé est inférieur au noued depart ==> toutes les occurence du noued
                        /// de départ prendront ceux du noeud d'arrivé
                        numeroClustre = (populationEA+i)->phenotype[nd];
                        for(k=0; k<nbrNoeuds; k++)
                        {
                            if((populationEA+i)->phenotype[k] == numeroClustre)
                            {
                                (populationEA+i)->phenotype[k] = (populationEA+i)->phenotype[na];
                            }
                        }
                    }
                    else if((populationEA+i)->phenotype[na] > (populationEA+i)->phenotype[nd])
                    {
                        numeroClustre = (populationEA+i)->phenotype[na];
                        for(k=0; k<nbrNoeuds; k++)
                        {
                            if((populationEA+i)->phenotype[k] == numeroClustre)
                            {
                                (populationEA+i)->phenotype[k] = (populationEA+i)->phenotype[nd];
                            }
                        }
                    }
                }
                else
                {
                    (populationEA+i)->phenotype[na] = (populationEA+i)->phenotype[nd];
                }
            }
/**
/// L'ancienne version qui était basée sur tout les arêtes, la nouvelle utilise seulement les arêtes intra
/// pour éviter de parcourir tout le vecteur des arêtes
        for(j=0; j<nbrArretes; j++)
        {
            /// ajouter le 23/12/2015 pour déminuer l'impact de l'arrête fictive sur les solution générer
            ///if((population+i)->genotype[j]==0 && fluxVector[j]>0){
            if((populationEA+i)->genotype[j]==0)
            {
                nd = (edgeVector+j)->nouedDepart;
                na = (edgeVector+j)->nouedArrive;
               if((populationEA+i)->phenotype[na]!= na) /// il a été modifier auparavant
                {

                    if((populationEA+i)->phenotype[na] < (populationEA+i)->phenotype[nd])
                    {
                        ///dans le cas ou le noued d'arrivé est inférieur au noued depart ==> toutes les occurence du noued
                        /// de départ prendront ceux du noeud d'arrivé
                        numeroClustre = (populationEA+i)->phenotype[nd];
                        for(k=0; k<nbrNoeuds; k++)
                        {
                            if((populationEA+i)->phenotype[k] == numeroClustre)
                            {
                                (populationEA+i)->phenotype[k] = (populationEA+i)->phenotype[na];
                            }
                        }
                    }
                    else if((populationEA+i)->phenotype[na] > (populationEA+i)->phenotype[nd])
                    {
                        numeroClustre = (populationEA+i)->phenotype[na];
                        for(k=0; k<nbrNoeuds; k++)
                        {
                            if((populationEA+i)->phenotype[k] == numeroClustre)
                            {
                                (populationEA+i)->phenotype[k] = (populationEA+i)->phenotype[nd];
                            }
                        }
                    }
                }
                else
                {
                    (populationEA+i)->phenotype[na] = (populationEA+i)->phenotype[nd];
                }
            }
    }
}

*/
    }
}
///************************************************************************************************************
void writeSolutionInFileEA(partitionEA *populationEA, FILE *outputFilePop,int iteration)
{

    int i,j;
    for(i=0; i<taillePopulation; i++)
    {

        ///écriture de la meilleur solution dans le fichier output
        fprintf(outputFilePop,"%d\t%d\t%d\t%0.4f\t%d\t%d\t%d\t",iteration,(populationEA+i)->id,
                (populationEA+i)->coutCoupe,(populationEA+i)->fitness,(populationEA+i)->coutCoupeNormalise,
                (populationEA+i)->contrainteViole,(populationEA+i)->nbrCluster);

        for(j=0; j<nbrArretes; j++)
        {
            fprintf(outputFilePop,"%d ",(populationEA+i)->genotype[j]);
        }
        fprintf(outputFilePop,"\t\t\t\t");

        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(outputFilePop,"%d ",(populationEA+i)->phenotype[j]);
        }
        fprintf(outputFilePop,"\n");
    }
    fprintf(outputFilePop,"\n\n");

}
///************************************************************************************************************
void writeBestSolutionInFileEA(partitionEA *solutionDominante, FILE *outputFile,int iteration)
{

    int i;
    fprintf(outputFile,"%d\t%d\t%d\t%0.4f\t%d\t%d\t%d\t",iteration,solutionDominante->id,
            solutionDominante->coutCoupe,solutionDominante->fitness,solutionDominante->coutCoupeNormalise,
            solutionDominante->contrainteViole,solutionDominante->nbrCluster);
    for(i=0; i<nbrArretes; i++)
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

///************************************************************************************************************
/// calculer le cout de coupe de chaque partition ainsi que son fitness
///************************************************************************************************************
/// on va appliquer la formule suivante pour le calcule des cout de coupe
/// T'(C) = B - T(C)            <=> la minimisation devient une maximisation
/// Z(C)=T'(C)+(u-v(C)) ×B      <=>  calculer la nouvelle fitness par rapport aux contrainte
/// T'(C) : cout de coupe (maximiser les flux INTRA), u ; nombre de contrainte totale
/// v(C): nombre de contrainte violé,  B: la somme totale des flux
///************************************************************************************************************

void calculCoutCoupeEtFitnessEA(partitionEA* populationEA)
{
    int i,j,k,sommeTotalCoutCoupe = 0;

    for(i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (populationEA+i)->coutCoupe = 0; /// réinitialisation de la valriable après chaque itération (individus de la population)
        for(j=0; j<nbrNoeuds-1; j++)
        {
            for(k=j+1; k<nbrNoeuds; k++)
            {
                if((populationEA+i)->phenotype[j] == (populationEA+i)->phenotype[k] && fluxMatrix[j][k] > 0)
                {
                    (populationEA+i)->coutCoupe = (populationEA+i)->coutCoupe +  fluxMatrix[j][k];
                }
            }
        }

        (populationEA+i)->coutCoupeNormalise = (populationEA+i)->coutCoupe + ((nbr_constraint - (populationEA+i)->contrainteViole)*sommeTotalFlux);
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationEA+i)->coutCoupeNormalise;
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
            varianceCoutDeCoupeNormalise + pow(((populationEA+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);
    /// tester la valeur ajouté
    for(i=0; i<taillePopulation; i++)
    {
        /// calcul de sigma truncation
        (populationEA+i)->expectedValue = (populationEA+i)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationEA+i)->expectedValue;
    }

    /// calcule du fitness
    /// le 26/11/2015 : calcul du fitness des solution avec leur expectedValues
    for(i=0; i<taillePopulation ; i++)
    {
        (populationEA+i)->fitness = (float)((populationEA+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);
    }

#else
    for(i=0; i<taillePopulation ; i++)
    {
        (populationEA+i)->fitness = (float)((populationEA+i)->coutCoupeNormalise)/(float)(sommeTotalCoutCoupe);
    }
#endif// scaling
}
///************************************************************************************************************
void calculCoutCoupeEtFitnessWithFlowVectorEA(partitionEA* populationEA)
{
#if mouchard
    printf("calculCoutCoupeEtFitness ...\n");
#endif // mouchard
    int i,j,sommeTotalCoutCoupe = 0;

    for(i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (populationEA+i)->coutCoupe = 0;
        for(j=0; j<nbrArretes; j++) ///parcourir le genotype =  solution codée en binaire
        {
            if((populationEA+i)->genotype[j]==0)  /// on cherche à maximiser les intraCluster donc les arrêtes 0 != 1
            {
                (populationEA+i)->coutCoupe = (populationEA+i)->coutCoupe + fluxVector[j];
            }
        }

        (populationEA+i)->coutCoupeNormalise = (populationEA+i)->coutCoupe + ((nbr_constraint - (populationEA+i)->contrainteViole)*sommeTotalFlux);

        /// désormis on va utiliser le cout de coupe normaliser pour le calcule des fitness
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationEA+i)->coutCoupeNormalise;
    }
#if scaling

    FILE *file;
    file = fopen("D:/WorkingFolderAgProgramming/benchmark/binaryEncoding/expectedValueEA.txt","w");
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
            varianceCoutDeCoupeNormalise + pow(((populationEA+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);
    /// tester la valeur ajouté
    for(i=0; i<taillePopulation; i++)
    {
        /// calcul de sigma truncation
        (populationEA+i)->expectedValue = (populationEA+i)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);
        /*            if(iteration == 1){
                        fprintf(file,"%d \t %d\n",(populationEA+i)->contrainteViole,(populationEA+i)->coutCoupeNormalise);
                    }*/
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationEA+i)->expectedValue;
    }

    /// calcule du fitness
    /// le 26/11/2015 : calcul du fitness des solution avec leur expectedValues
    for(i=0; i<taillePopulation ; i++)
    {
        (populationEA+i)->fitness = (float)((populationEA+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);
    }

#else
    for(i=0; i<taillePopulation ; i++)
    {
        (populationEA+i)->fitness = (float)((populationEA+i)->coutCoupeNormalise)/(float)(sommeTotalCoutCoupe);
    }
#endif// scaling

}
///************************************************************************************************************

void binaryEncoding(int nbrGeneration ,FILE *outputFileEA,FILE *outputFilePopEA,FILE *outputOptimalSolutionFileEA,partitionEA *populationEA1,partitionEA *populationEA2
                    , partitionEA *solutionDominante,int iteration , int *bestSolutionIteration , int *nbrApparition){

        naturalSelectionEA(populationEA1,populationEA2);
        crossOverEA(populationEA2,populationEA1); /// nbrArretes représente la taille des genotype
        mutationEA(populationEA1);
        getPartitionFromSolutionEA(populationEA1);
        ///getPartitionFromSolutionWithoutRepetitionEA(populationEA1);
        checkContrainstAndFitnessPenalizationEA(populationEA1);
        calculCoutCoupeEtFitnessEA(populationEA1);
        ///calculCoutCoupeEtFitnessWithFlowVectorEA(populationEA1);

#if writingPopulationInFile
        writeSolutionInFileEA(populationEA1,outputFilePopEA,iteration);
#endif // writingPopulationInFile
///=====================================================================================================
        /// dans ce qui suit, on va juste localiser la meilleur solution
///=====================================================================================================

        bestSolution=findTheBestSolutionEA(populationEA1);
#if writingPopulationInFile
        writeBestSolutionInFileEA((populationEA1+bestSolution),outputFileEA,iteration);
#endif // writingPopulationInFile
        if((populationEA1+bestSolution)->coutCoupeNormalise > solutionDominante->coutCoupeNormalise )
                /** le 23/12/2015 : pour avoir la meilleur solution avec le moins de contrainte violées*/
                /// && (populationEA1+bestSolution)->contrainteViole   < solutionDominante->contrainteViole)
        {
            *(solutionDominante) = *(populationEA1+bestSolution);
            nbrApparition=1;
            bestSolutionIteration = iteration;
        }
        else
        {
            nbrApparition++;
        }

}
///************************************************************************************************************
void checkContrainstAndFitnessPenalizationEA(partitionEA *populationEA)
{

    int i,j,k;

    for(i=0; i<taillePopulation; i++)
    {
        (populationEA+i)->nbrCluster=0; /// le nombre de cluster doit etre initialisé à 1 pour prendre en consédiration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes violées
        (populationEA+i)->contrainteViole=0;
        for(j=0; j<4; j++)
        {
            (populationEA+i)->constraintVector[j]=0;
        }
        ///**********************************************************************************


        ///trier le tableau résulant
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationEA+i)->clustersSize[j] =0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if ((populationEA+i)->phenotype[k] == j)
                {
                    (populationEA+i)->clustersSize[j]++;
                }
            }
            if ((populationEA+i)->clustersSize[j] !=0)
            {
                (populationEA+i)->nbrCluster++;
            }
        }

        if((populationEA+i)->nbrCluster > max_clusters || (populationEA+i)->nbrCluster < min_clusters)
        {
            (populationEA+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationEA+i)->contrainteViole++;
        }

        for(j=0; j<nbrNoeuds; j++)
        {
            if((populationEA+i)->clustersSize[j]!=0)
            {
                if((populationEA+i)->clustersSize[j]>max_sizeCluster || (populationEA+i)->clustersSize[j]< min_sizeCluster)
                {
                    (populationEA+i)->constraintVector[1]++;

                }
            }

        }

        if((populationEA+i)->constraintVector[1] != 0)
        {
            (populationEA+i)->contrainteViole++;
        }
    }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0)
    {
        checkCohabitationAndNonCohabitationConstrainteEA(populationEA);
    }

}
///************************************************************************************************************
void checkCohabitationAndNonCohabitationConstrainteEA(partitionEA *populationEA)
{

    int i,j,noeud1,noeud2;
    for(i=0; i<taillePopulation; i++)
    {
	    ///24/11/2015 normalement cette taches est déjà faite dans le fonction : checkContrainstAndFitnessPenalizationFC
        /**
        (populationEA+i)->constraintVector[2]=0;
        (populationEA+i)->constraintVector[3]=0;
		*/
        if(nbrCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrCohabitationConstraintes; j++)
            {
                noeud1=cohabitationConstraintes[j].nouedDepart;
                noeud2=cohabitationConstraintes[j].nouedArrive;

                if((populationEA+i)->phenotype[noeud1]!= (populationEA+i)->phenotype[noeud2])
                {
                    (populationEA+i)->constraintVector[2]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationEA+i)->constraintVector[2]!=0)
            {
                (populationEA+i)->contrainteViole++;
            }
        }


        if(nbrNonCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrNonCohabitationConstraintes; j++)
            {
                noeud1=nonCohabitationConstraintes[j].nouedDepart;
                noeud2=nonCohabitationConstraintes[j].nouedArrive;

                if((populationEA+i)->phenotype[noeud1]== (populationEA+i)->phenotype[noeud2])
                {
                    (populationEA+i)->constraintVector[3]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationEA+i)->constraintVector[3]!=0)
            {
                (populationEA+i)->contrainteViole++;
            }
        }


    }

}

///********************************************************************************************
void writeOptimalSolutionInFileEA(partitionEA *solutionDominante,FILE* outputOptimalSolutionFileEA,
                                  int nbrRun, int bestSolutionIteration, float runTime,int ES)
{

    fprintf(outputOptimalSolutionFileEA,"EA | RUN NUMBER = %d | ITERATION N = %d | SOLUTION IDENTIFIER = %d | CUT SIZE INTRA = %d | CUT SIZE INTER = %d | SOMME TOTAL DES FLUX = %d | RUN TIME = %f | ES = %d | NBR VIOLATED CONSTRAINTE = %d |",
            nbrRun, bestSolutionIteration, solutionDominante->id, solutionDominante->coutCoupe,(sommeTotalFlux - solutionDominante->coutCoupe),sommeTotalFlux ,runTime, ES, solutionDominante->contrainteViole);
 fprintf(outputOptimalSolutionFileEA," max_clusters = %d | max_sizeCluster = %d  |  cohabitation = %d | non cohabitation = %d\n",
            solutionDominante->constraintVector[0],solutionDominante->constraintVector[1],solutionDominante->constraintVector[2],
            solutionDominante->constraintVector[3]);

}



///******************************************************************************************
int compareCroissantFitnessEA (void const *a, void const *b)
{

    partitionEA const *pa = a;
    partitionEA const *pb = b;

    return  pb->coutCoupeNormalise - pa->coutCoupeNormalise;;
}


///******************************************************************************************


/**

void writeOptimalSolutionInFileEA(partitionEA *solutionDominante,FILE* outputOptimalSolutionFileEA){
    int i;
    fprintf(outputOptimalSolutionFileEA,"id = %d | normalised cut size = %d | cut size intrA = %d cut size intEr = %d ",
            solutionDominante->id,solutionDominante->coutCoupeNormalise, solutionDominante->coutCoupe, (sommeTotalFlux-solutionDominante->coutCoupe));
    fprintf(outputOptimalSolutionFileEA,"le contraintes violees = %d \n", solutionDominante->contrainteViole);
    for(i=0; i<4; i++)
    {
        switch (i)
        {

        case 0:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileEA,"contrainte sur le nombre de clusters : OK \n");
            else fprintf(outputOptimalSolutionFileEA,"contrainte sur le nombre de clusters : NOT OK \n");
            break;
        case 1:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileEA,"contrainte sur la taille des clusters : OK \n");
            else fprintf(outputOptimalSolutionFileEA,"contrainte sur la taille des clusters : NOT OK \n");
            break;
        case 2:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileEA,"contrainte de cohabitation : OK \n");
            else fprintf(outputOptimalSolutionFileEA,"contrainte de cohabitation : NOT OK \n");
            break;
        case 3:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileEA,"contrainte Non Cohabitation : OK \n");
            else fprintf(outputOptimalSolutionFileEA,"contrainte Non Cohabitation : NOT OK \n");
            break;
        default:
            fprintf(outputOptimalSolutionFileEA,"Cette contrainte n'est pas prise en charge par le système \n");

        }

    }

}


void checkContrainstAndFitnessPenalizationEA(partitionEA *populationEA){

    int i,j,k,m=nbrNoeuds,tabTmp[1000],tmp;

    for(i=0; i<taillePopulation; i++)
    {
        (populationEA+i)->nbrCluster=1; /// le nombre de cluster doit etre initialisé à 1 pour prendre en consédiration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes violées
        (populationEA+i)->contrainteViole=0;
        for(j=0;j<nbr_constraint;j++){(populationEA+i)->constraintVector[j]=0;}


        m=nbrNoeuds;
        ///copier le phenotype dans un deuxième tableau
        for(j=0; j<nbrNoeuds; j++)
        {
            tabTmp[j]=(populationEA+i)->phenotype[j];
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
                (populationEA+i)->nbrCluster++;
            }
        }
        if((populationEA+i)->nbrCluster > max_clusters || (populationEA+i)->nbrCluster < min_clusters)
        {
            (populationEA+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationEA+i)->contrainteViole++;
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
            (populationEA+i)->clustersSize[j]=0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if (tabTmp[j] == (populationEA+i)->phenotype[k])
                {
                    (populationEA+i)->clustersSize[j]=(populationEA+i)->clustersSize[j]+1;
                }
            }
            if((populationEA+i)->clustersSize[j]>max_sizeCluster || (populationEA+i)->clustersSize[j]< min_sizeCluster)
            {
                ///24.07.2015 ---> l'ajout de vecteur de contraintes = le nombre de cluster dont la taille a dépassé la borne
                (populationEA+i)->constraintVector[1]++;

            }
        }

        if((populationEA+i)->constraintVector[1] != 0)
        {
            (populationEA+i)->contrainteViole++;
        }
   }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0){
        checkCohabitationAndNonCohabitationConstrainteEA(populationEA);
    }

}
*/
