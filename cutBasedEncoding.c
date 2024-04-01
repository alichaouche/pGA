/**================= CutBased Encoding Genetic Algorithm ====================================================
* 18/04/2015 :
* la population initiale est maintenant prête
* la tache du jour :
*       1- programmee dans un premier temps la fonction qui permet de calculer le cout de coupe et le fitness
*       2- programmer la fonction de séléction naturelle
*       3- Remarque : la taille de la population doit être toujours paire
*============================================================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "typeDeclaration.h"
#include "./communesFunctions.h"
#include "./cutBasedEncoding.h"
#include "./compilationConditionnelle.h"

extern char *cocyclesDeBase;
extern int taillePopulation,tauxReproduction,sommeTotalFlux,nbrNoeuds,
       nbr_constraint,nbrArretes,nbrParties,nbrNoeudsInCluster,regeneration;
extern double tauxMutation ,averageTraffic;
extern int fluxMatrix[1000][1000],*fluxVector, maximumDistanceMatrix[1000][1000], maxFlow, sommeTotalFluxReal;
///extern ul cocyclesDeBaseEntier[1000];
extern mpz_t gmpCocyclesDeBaseEntier[1000];
extern edge *edgeVector,cohabitationConstraintes[1000],nonCohabitationConstraintes[1000],*edgeVectorIntra;
extern int max_clusters , min_clusters,min_sizeCluster ,max_sizeCluster ;
extern int bestSolution, nbrApparition, nbrCohabitationConstraintes,nbrNonCohabitationConstraintes;
extern int nbrRun;

///**************************************************************************************
void generatePopulationDC(partitionDC* populationDC , int indiceFirstElt)
{
    int i,j,k,nbrCocycles=nbrParties-1,cocyclesSize=nbrCocycles*(nbrNoeuds-1);


   for(i=indiceFirstElt;i<taillePopulation;i++){
        (populationDC+i)->cocycles = (int*)malloc(cocyclesSize*sizeof(int));
        if((populationDC+i)->cocycles == NULL) {printf("Memory allocation failled for cocycles DC\n");exit(EXIT_FAILURE);}
    }

    for(i=indiceFirstElt;i<taillePopulation;i++){
        (populationDC+i)->genotype = (int*)malloc(nbrArretes*sizeof(int));
        if((populationDC+i)->genotype == NULL) {printf("Memory allocation failled for genotype DC\n");exit(EXIT_FAILURE);}
    }

    for(i=indiceFirstElt; i<taillePopulation; i++)
    {
        ///printf("i = %d (populationDC+i)->cocycles= %p \n",i,(populationDC+i)->cocycles);
        for(j=0; j<nbrCocycles; j++)
        {
            for(k=(j*(nbrNoeuds-1)); k<((j+1)*(nbrNoeuds-1)); k++)
            {
                (populationDC+i)->cocycles[k] = rnd(0,2); /// les valeurs générer seront entre 0 et 1
                ///if((populationDC+i)->cocycles[k] != 0 && (populationDC+i)->cocycles[k] !=1) printf("cannot access memory \n");
            }
        }
        (populationDC+i)->id =i;
        (populationDC+i)->coutCoupe =0;
        (populationDC+i)->fitness =0;

    }
        calculerGenotypeDC(populationDC, indiceFirstElt);
}
///**************************************************************************************
void affichePopulationDC(partitionDC* populationDC)
{
    printf("affichePopulation...\n");
    int i,j,k,l;
    for(i=0; i<taillePopulation; i++)
    {

        printf("\nid = %d | cout de coupe = %d | fitness = %0.4f",
               (populationDC+i)->id,(populationDC+i)->coutCoupe, (populationDC+i)->fitness);

        /// affichage de vecteur de la partition
        /**       printf("\nla partition est  = \n");
               for(l=0;l<nbrParties;l++){
                   printf("%d ",(populationDC+i)->genotypeEntier[l]);
               }
        */

        printf("\nla solution generer est  = \n");
        for(j=0; j<nbrParties-1; j++)
        {
            for(k=0; k<nbrNoeuds-1; k++)
            {
                printf("%d ",(populationDC+i)->cocycles[(j*(nbrNoeuds-1))+k]);
            }
            printf("\n");
        }
        printf("\n");
        /// affichage de vecteur de la partition
        printf("\nle genotype est  = \n");
        for(l=0; l<nbrArretes; l++)
        {
            printf("%d ",(populationDC+i)->genotype[l]);
        }
        printf("\n===============================================================================\n");

    }
}
///**************************************************************************************
void calculerGenotypeDC(partitionDC* populationDC, int indiceFirstElt)
{

    int i,j,k,numeroCocyclXor;
    mpz_t cocyclesXor[1000];  for(j=0;j<nbrParties-1;j++)mpz_init(cocyclesXor[j]);
    mpz_t genotypeTmp; mpz_init(genotypeTmp); /// il est de type mpz_t car il peut avoir des valeur allant jusqu'à [ pow(2, nbrNoueds-1)- 1]
    mpz_t binaryTmp; mpz_init(binaryTmp);

    for(i=indiceFirstElt; i<taillePopulation; i++)
    {
       /// affichePopulationDC(populationDC);
        for(j=0; j<nbrArretes; j++)
        {
            (populationDC+i)->genotype[j]=0;
        }
        ///***********************************************************************************************
        for(j=0; j<nbrParties-1; j++)
        {
            numeroCocyclXor=0;
            for(k=0; k<nbrNoeuds-1; k++){
            ///for(k=nbrNoeuds-2; k>=0; k--){
                if((populationDC+i)->cocycles[(j*(nbrNoeuds-1))+k]==1)
                {
                    if(numeroCocyclXor==0)
                    {
                        mpz_set(cocyclesXor[j] , gmpCocyclesDeBaseEntier[k]);
                        numeroCocyclXor++;

                    }
                    else
                    {
                        mpz_xor(cocyclesXor[j] , cocyclesXor[j] , gmpCocyclesDeBaseEntier[k]);
                    }
                }
            }
            if(j==0)
            {
                mpz_set(genotypeTmp , cocyclesXor[j]);
            }
            else
            {
                mpz_ior(genotypeTmp , genotypeTmp, cocyclesXor[j]);
            }
        }
        ///******************************************************************************
        j=nbrArretes-1;
        do
        {
            mpz_mod_ui(binaryTmp,genotypeTmp,2);
            (populationDC+i)->genotype[j]=mpz_get_ui(binaryTmp);
            mpz_fdiv_q_ui(genotypeTmp,genotypeTmp,2);
            j--;
        }
        while(j>=0);
    }

}
///**************************************************************************************
/// la séléction naturelle des individus pour la nouvelle génération
void naturalSelectionDC(partitionDC* populationDC1,partitionDC* populationDC2)
{
    int i,j,maxFitness;
    int tmpVectIndice[taillePopulation];
    for(i=0;i<taillePopulation;i++) tmpVectIndice[i]=-1;
    float lotterie,sommeFitness;

    for(i=0; i<tauxReproduction; i++){
        maxFitness = 0;
        for(j=0;j<taillePopulation;j++){
            if((populationDC1+maxFitness)->coutCoupeNormalise < (populationDC1+j)->coutCoupeNormalise
               && tmpVectIndice[j] == -1 ){
                maxFitness = j;

            }
        }
        tmpVectIndice[maxFitness] = 0;
        *(populationDC2+i) = *(populationDC1+maxFitness);
        ///printf("(populationDC2+%d)->coutCoupeNormalise = %d \n",i,(populationDC2+i)->coutCoupeNormalise);
        ///mon_sleep(pause);
    }
    ///printf("\n\n\n");

    ///qsort (populationDC2, taillePopulation, sizeof *populationDC2, compareCroissantFitnessEE);
    ///sortingPopulationByFintness(populationDC2,tmpPopulation);

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
            sommeFitness = sommeFitness + (populationDC1+j)->fitness;
            if (lotterie<=sommeFitness)
            {
                *(populationDC2+i) = *(populationDC1+j);
                break;
            }
            else if( j==taillePopulation-1)
            {
                *(populationDC2+i) = *(populationDC1+j);
            }
        }
    }
    ///free(tmpPopulation);

}
///**************************************************************************************
///le croisement des solution séléctionner
void crossOverDC(partitionDC* populationDC1, partitionDC* populationDC2)
{

    int i = 0,j;
    int choixInd1, choixInd2, choixLocus=0,tailleCocycles;
    tailleCocycles = (nbrParties-1)*(nbrNoeuds-1);

    ///***********************************************************************************
    for(i=0; i<tauxReproduction; i++)
    {
        *(populationDC1+i) = *(populationDC2+i);
        (populationDC1+i)->id = i;
    }
    ///***********************************************************************************
#if NEW_CROSSOVER
    while (i < taillePopulation-regeneration)
    {
        choixInd1 = rnd(0,taillePopulation);
        ///choixInd1 = rnd(tauxReproduction,taillePopulation);
        do
        {
            /// ce traitement évite l'appariement d'un element avec lui même/
            ///choixInd2 = rnd(tauxReproduction,taillePopulation);
            choixInd2 = rnd(0,taillePopulation);
        }
        while(choixInd2 == choixInd1 );

#if choixNbrPartiePourCrossoverDC
        choixLocus = rnd(0,nbrParties-1); /// -1 pour s'arrêter à moins 2
        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            for(k=0; k<nbrNoeuds-1; k++)
            {
                (populationDC1+i)->genotype[j*(nbrNoeuds-1)+k] = (populationDC2+choixInd1)->genotype[j*(nbrNoeuds-1)+k];
                (populationDC1+i+1)->genotype[j*(nbrNoeuds-1)+k] = (populationDC2+choixInd2)->genotype[j*(nbrNoeuds-1)+k];
            }
        }

        ///************************************************************************************

        for(j=choixLocus+1; j<nbrParties-1; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            for(k=0; k<nbrNoeuds-1; k++)
            {
                (populationDC1+i)->genotype[j*(nbrNoeuds-1)+k] = (populationDC2+choixInd2)->genotype[j*(nbrNoeuds-1)+k];
                (populationDC1+i+1)->genotype[j*(nbrNoeuds-1)+k] = (populationDC2+choixInd1)->genotype[j*(nbrNoeuds-1)+k];
            }
        }
#else
        choixLocus = rnd(0,tailleCocycles-1); /// -1 pour s'arrêter à moins 2
#endif
        ///****************************************************************************************************
        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationDC1+i)->cocycles[j] = (populationDC2+choixInd1)->cocycles[j];
            (populationDC1+i+1)->cocycles[j] = (populationDC2+choixInd2)->cocycles[j];
        }
        ///****************************************************************************************************
        for(j=choixLocus+1; j<tailleCocycles; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationDC1+i)->cocycles[j] = (populationDC2+choixInd2)->cocycles[j];
            (populationDC1+i+1)->cocycles[j] = (populationDC2+choixInd1)->cocycles[j];
        }
        /// affectation de nouveau indice
        (populationDC1+i)->id = i;
        (populationDC1+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }
    generatePopulationDC(populationDC1 , taillePopulation-regeneration);
#else

    while (i < taillePopulation)
    {
       /// choixInd1 = rnd(tauxReproduction,taillePopulation);
        choixInd1 = rnd(0,taillePopulation);
        do
        {
            /// ce traitement évite l'appariement d'un element avec lui même/
            ///choixInd2 = rnd(tauxReproduction,taillePopulation);
            choixInd2 = rnd(0,taillePopulation);
        }
        while(choixInd2 == choixInd1 );

#if choixNbrPartiePourCrossoverDC
        choixLocus = rnd(0,nbrParties-1); /// -1 pour s'arrêter à moins 2
        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            for(k=0; k<nbrNoeuds-1; k++)
            {
                (populationDC1+i)->genotype[j*(nbrNoeuds-1)+k] = (populationDC2+choixInd1)->genotype[j*(nbrNoeuds-1)+k];
                (populationDC1+i+1)->genotype[j*(nbrNoeuds-1)+k] = (populationDC2+choixInd2)->genotype[j*(nbrNoeuds-1)+k];
            }
        }

        ///************************************************************************************

        for(j=choixLocus+1; j<nbrParties-1; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            for(k=0; k<nbrNoeuds-1; k++)
            {
                (populationDC1+i)->genotype[j*(nbrNoeuds-1)+k] = (populationDC2+choixInd2)->genotype[j*(nbrNoeuds-1)+k];
                (populationDC1+i+1)->genotype[j*(nbrNoeuds-1)+k] = (populationDC2+choixInd1)->genotype[j*(nbrNoeuds-1)+k];
            }
        }
#else
        choixLocus = rnd(0,tailleCocycles-1); /// -2 pour s'arrêter à moins 2
#endif
        ///****************************************************************************************************
        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationDC1+i)->cocycles[j] = (populationDC2+choixInd1)->cocycles[j];
            (populationDC1+i+1)->cocycles[j] = (populationDC2+choixInd2)->cocycles[j];
        }
        ///****************************************************************************************************
        for(j=choixLocus+1; j<tailleCocycles; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationDC1+i)->cocycles[j] = (populationDC2+choixInd2)->cocycles[j];
            (populationDC1+i+1)->cocycles[j] = (populationDC2+choixInd1)->cocycles[j];
        }
        /// affectation de nouveau indice
        (populationDC1+i)->id = i;
        (populationDC1+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }

#endif
}
///**************************************************************************************
void mutationDC(partitionDC* populationDC)
{
    int i,numeroGeneMute,tailleCocycles;
    tailleCocycles = (nbrParties-1)*(nbrNoeuds-1);
    double applicationOfMutation;
    for(i=tauxReproduction; i<taillePopulation; i++)
    {
#if Windows
        applicationOfMutation= drand48ForWindows(0,1001);
#else
        applicationOfMutation= drand48();
#endif // Windows
        if(applicationOfMutation <= tauxMutation)
        {
            /// choisir un gène aléatoirement pour lui appliquer la mutation
            numeroGeneMute = rnd(0,tailleCocycles);
            ///((populationDC+i)->genotype[numeroGeneMute] == 1)? 0 : 1;
            if((populationDC+i)->cocycles[numeroGeneMute] == 1)
            {
                (populationDC+i)->cocycles[numeroGeneMute] =0;
            }
            else
            {
                (populationDC+i)->cocycles[numeroGeneMute] =1;
            }
        }
    }
        getPartitionFromSolutionDC(populationDC);
        checkContrainstAndFitnessPenalizationDC(populationDC);
        calculCoutCoupeEtFitnessDC(populationDC);

}
///**************************************************************************************
float testerLaSommeDesFitnessDC(partitionDC* populationDC)
{
    printf("testerLaSommeDesFitness ...\n");
    int i;
    float sommeFitness = 0;
    for(i=0; i<taillePopulation ; i++)
    {
        sommeFitness = sommeFitness + (populationDC+i)->fitness;
    }
    return sommeFitness;
}
///**************************************************************************************
int findTheBestSolutionDC(partitionDC *populationDC)
{

    float maxFitness = (populationDC+0)->fitness;
    int i,indice=0;
    for(i=0; i<taillePopulation ; i++)
    {
        if(maxFitness < (populationDC+i)->fitness)
        {
            maxFitness = (populationDC+i)->fitness;
            indice = i;
        }
    }
    return indice;

}
///**************************************************************************************
void displayTheBestSolutionDC(partitionDC* solutionDominante)
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
void writeSolutionInFileDC(partitionDC *populationDC, FILE *outputFilePop,int iteration)
{
    int i,j;
    for(i=0; i<taillePopulation; i++)
    {

        ///écriture de la meilleur solution dans le fichier output
        fprintf(outputFilePop,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,(populationDC+i)->id,
                (populationDC+i)->coutCoupe,(populationDC+i)->fitness,(populationDC+i)->coutCoupeNormalise,
                (populationDC+i)->contrainteViole,(populationDC+i)->nbrCluster);
        for(j=0; j<nbrArretes; j++)
        {
            fprintf(outputFilePop,"%d ",(populationDC+i)->genotype[j]);
        }
        fprintf(outputFilePop,"\t\t\t\t");
        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(outputFilePop,"%d ",(populationDC+i)->phenotype[j]);
        }
        fprintf(outputFilePop,"\n");
    }
    /// fprintf(outputFilePop,"\n===============================================\n");
    fprintf(outputFilePop,"\n\n");

}
///************************************************************************************************************
void writeBestSolutionInFileDC(partitionDC *solutionDominante, FILE *outputFile,int iteration)
{

    int i;
    fprintf(outputFile,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,solutionDominante->id,
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
void getPartitionFromSolutionDC(partitionDC *populationDC)
{
        int i,j,k,nd,na,numeroClustre,nbrIntra; /// nd = noued de départ , na = noued d'arrivé

    for(i=0; i<taillePopulation; i++)
    {
        /// initialisation des phénotype = partition par les numero des noeuds correspondant

        for(j=0; j<nbrNoeuds; j++)
        {
            (populationDC+i)->phenotype[j] = j;
        }
        /// 01/10/2017 : amélioration du temps d'exécution du codage binaire
        nbrIntra = 0;
        for(j=0; j<nbrArretes; j++){
            if((populationDC+i)->genotype[j]==0)
            {
               edgeVectorIntra[nbrIntra] = edgeVector[j];
               nbrIntra++;
            }

        }
        /// détermination de la parition
        for(j=0; j<nbrIntra; j++){

                nd = (edgeVectorIntra+j)->nouedDepart;
                na = (edgeVectorIntra+j)->nouedArrive;
               if((populationDC+i)->phenotype[na]!= na) /// il a été modifier auparavant
                {

                    if((populationDC+i)->phenotype[na] < (populationDC+i)->phenotype[nd])
                    {
                        ///dans le cas ou le noued d'arrivé est inférieur au noued depart ==> toutes les occurence du noued
                        /// de départ prendront ceux du noeud d'arrivé
                        numeroClustre = (populationDC+i)->phenotype[nd];
                        for(k=0; k<nbrNoeuds; k++)
                        {
                            if((populationDC+i)->phenotype[k] == numeroClustre)
                            {
                                (populationDC+i)->phenotype[k] = (populationDC+i)->phenotype[na];
                            }
                        }
                    }
                    else if((populationDC+i)->phenotype[na] > (populationDC+i)->phenotype[nd])
                    {
                        numeroClustre = (populationDC+i)->phenotype[na];
                        for(k=0; k<nbrNoeuds; k++)
                        {
                            if((populationDC+i)->phenotype[k] == numeroClustre)
                            {
                                (populationDC+i)->phenotype[k] = (populationDC+i)->phenotype[nd];
                            }
                        }
                    }
                }
                else
                {
                    (populationDC+i)->phenotype[na] = (populationDC+i)->phenotype[nd];
                }
            }
    }

/**   19/11/2017
    int i,j,k,nd,na,numeroClustre; /// nd = noued de départ , na = noued d'arrivé
    for(i=0; i<taillePopulation; i++)
    {
        /// initialisation des phénotype = partition par les numero des noeuds correspondant
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationDC+i)->phenotype[j] = j;
        }

        /// détermination de la parition
        ///for(j=0; j<nbrArretes; j++){
        for(j=nbrArretes-1; j>=0; j--)
        {
            if((populationDC+i)->genotype[j]==0)
            {
                na = edgeVector[j].nouedArrive;
                nd = edgeVector[j].nouedDepart;

                if((populationDC+i)->phenotype[na]!= na) /// il a été modifier auparavant
                {

                    if((populationDC+i)->phenotype[na] < (populationDC+i)->phenotype[nd])
                    {
                        ///dans le cas ou le noued d'arrivé est inférieur au noued d'arrivé ==> toutes les occurence du noued
                        /// de départ prendront ceux du noeud d'arrivé
                        numeroClustre = (populationDC+i)->phenotype[nd];
                        for(k=0; k<nbrNoeuds; k++)
                        {
                            if((populationDC+i)->phenotype[k] == numeroClustre)
                            {
                                (populationDC+i)->phenotype[k] = (populationDC+i)->phenotype[na];
                            }
                        }
                    }
                    else if((populationDC+i)->phenotype[na] > (populationDC+i)->phenotype[nd])
                    {
                        numeroClustre = (populationDC+i)->phenotype[na];
                        for(k=0; k<nbrNoeuds; k++)
                        {
                            if((populationDC+i)->phenotype[k] == numeroClustre)
                            {
                                (populationDC+i)->phenotype[k] = (populationDC+i)->phenotype[nd];
                            }
                        }
                    }
                }
                else
                {
                    (populationDC+i)->phenotype[na] = (populationDC+i)->phenotype[nd];
                }
            }
        }
    }
    */
}
///**************************************************************************************
/// calculer le cout de coupe de chaque partition ainsi que son fitness
///**************************************************************************************
/// on va appliquer la formule suivante pour le calcule des cout de coupe
/// T'(C) = B - T(C)            <=> la minimisation devient une maximisation
/// Z(C)=T'(C)+(u-v(C)) ×B      <=>  calculer la nouvelle fitness par rapport aux contrainte
/// T'(C) : cout de coupe (maximiser les flux INTRA), u ; nombre de contrainte totale
/// v(C): nombre de contrainte violé,  B: la somme totale des flux
///**************************************************************************************

void calculCoutCoupeEtFitnessDC(partitionDC* populationDC)
{
    int i,j,k,sommeTotalCoutCoupe = 0;

    for(i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (populationDC+i)->coutCoupe = 0; /// réinitialisation de la valriable après chaque itération (individus de la population)
        for(j=0; j<nbrNoeuds-1; j++)
        {
            for(k=j+1; k<nbrNoeuds; k++)
            {
                if((populationDC+i)->phenotype[j] == (populationDC+i)->phenotype[k] && fluxMatrix[j][k] > 0)
                {
                    (populationDC+i)->coutCoupe = (populationDC+i)->coutCoupe +  fluxMatrix[j][k];
                }
            }
        }
        (populationDC+i)->coutCoupeNormalise = (populationDC+i)->coutCoupe + ((nbr_constraint - (populationDC+i)->contrainteViole)*sommeTotalFlux);
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationDC+i)->coutCoupeNormalise;
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
            varianceCoutDeCoupeNormalise + pow(((populationDC+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);
    /// tester la valeur ajouté
    for(i=0; i<taillePopulation; i++)
    {

        /// calcul de sigma truncation
        (populationDC+i)->expectedValue = (populationDC+i)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationDC+i)->expectedValue;
    }

    /// calcule du fitness
    /// le 26/11/2015 : calcul du fitness des solution avec leur expectedValues
    for(i=0; i<taillePopulation ; i++)
    {
        (populationDC+i)->fitness = (float)((populationDC+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);
    }

#else
    FILE *file;
    file = fopen("D:/WorkingFolderAgProgramming/benchmark/cutBasedEncoding/expectedValueDC.txt","w");

    for(i=0; i<taillePopulation ; i++)
    {
        (populationDC+i)->fitness = (float)((populationDC+i)->coutCoupeNormalise)/(float)(sommeTotalCoutCoupe);
    }
#endif// scaling

}

///**************************************************************************************
/// calculer le cout de coupe de chacune des solution dans la population
void calculCoutCoupeEtFitnessWithFlowVectorDC(partitionDC* populationDC)
{
    int i,j,sommeTotalCoutCoupe = 0;

    for(i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (populationDC+i)->coutCoupe = 0;

        for(j=0; j<nbrArretes; j++) ///parcourir le vecteurSolution =  solution codée en binaire
        {
            if((populationDC+i)->genotype[j]==0)  /// on cherche à maximiser les intraCluster donc les arrêtes 0 != 1
            {
                (populationDC+i)->coutCoupe = (populationDC+i)->coutCoupe + fluxVector[j];
            }
        }
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationDC+i)->coutCoupe;

    }
    /// calcule du fitness
    for(i=0; i<taillePopulation ; i++)
    {
        (populationDC+i)->fitness = (float)((populationDC+i)->coutCoupe)/(float)(sommeTotalCoutCoupe);
    }
}

///======================================================================================
void cutBasedEncoding(int nbrGeneration ,FILE *outputFileDC,FILE *outputFilePopDC,FILE *outputOptimalSolutionFileDC,partitionDC *populationDC1,
                      partitionDC *populationDC2,partitionDC *solutionDominante,int iteration , int *bestSolutionIteration , int *nbrApparition)
{


        naturalSelectionDC(populationDC1,populationDC2);
        crossOverDC(populationDC2, populationDC1);
        mutationDC(populationDC1);
        ///calculerGenotypeDC(populationDC1);
        getPartitionFromSolutionDC(populationDC1);
        checkContrainstAndFitnessPenalizationDC(populationDC1);
        calculCoutCoupeEtFitnessDC(populationDC1);

#if writingPopulationInFile
        writeSolutionInFileDC(populationDC1,outputFilePopDC,iteration);
#endif // writingPopulationInFile
        bestSolution=findTheBestSolutionDC(populationDC1);
#if writingPopulationInFile
        writeBestSolutionInFileDC((populationDC1+bestSolution),outputFileDC,iteration);
#endif // writingPopulationInFile
        if((populationDC1+bestSolution)->coutCoupeNormalise > solutionDominante->coutCoupeNormalise)
        {
            *(solutionDominante) = *(populationDC1+bestSolution);
            nbrApparition=1;
            bestSolutionIteration = iteration;
        }
        else
        {
            nbrApparition++;
        }

}
///**************************************************************************************
void checkContrainstAndFitnessPenalizationDC(partitionDC *populationDC)
{

    int i,j,k;

    for(i=0; i<taillePopulation; i++)
    {
        (populationDC+i)->nbrCluster=0; /// le nombre de cluster doit etre initialisé à 1 pour prendre en consédiration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes violées
        (populationDC+i)->contrainteViole=0;
        for(j=0; j<4; j++)
        {
            (populationDC+i)->constraintVector[j]=0;
        }
        ///**********************************************************************************


        ///trier le tableau résulant
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationDC+i)->clustersSize[j] =0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if ((populationDC+i)->phenotype[k] == j)
                {
                    (populationDC+i)->clustersSize[j]++;
                }
            }
            if ((populationDC+i)->clustersSize[j] !=0)
            {
                (populationDC+i)->nbrCluster++;
            }
        }

        if((populationDC+i)->nbrCluster > max_clusters || (populationDC+i)->nbrCluster < min_clusters)
        {
            (populationDC+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationDC+i)->contrainteViole++;
        }

        for(j=0; j<nbrNoeuds; j++)
        {
            if((populationDC+i)->clustersSize[j]!=0)
            {
                if((populationDC+i)->clustersSize[j]>max_sizeCluster || (populationDC+i)->clustersSize[j]< min_sizeCluster)
                {
                    (populationDC+i)->constraintVector[1]++;

                }
            }

        }

        if((populationDC+i)->constraintVector[1] != 0)
        {
            (populationDC+i)->contrainteViole++;
        }
    }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0)
    {
        checkCohabitationAndNonCohabitationConstrainteDC(populationDC);
    }

}
///**************************************************************************************
void checkCohabitationAndNonCohabitationConstrainteDC(partitionDC *populationDC)
{

    int i,j,noeud1,noeud2;
    for(i=0; i<taillePopulation; i++)
    {
        /// Cette tache est déjà réalisé au niveau de la fonction checkCohabitationAndNonCohabitationConstrainteBE
		/**
		(populationDC+i)->constraintVector[2]=0;
        (populationDC+i)->constraintVector[3]=0;
		*/
        if(nbrCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrCohabitationConstraintes; j++)
            {
                noeud1=cohabitationConstraintes[j].nouedDepart;
                noeud2=cohabitationConstraintes[j].nouedArrive;

                if((populationDC+i)->phenotype[noeud1]!= (populationDC+i)->phenotype[noeud2])
                {
                    (populationDC+i)->constraintVector[2]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationDC+i)->constraintVector[2]!=0)
            {
                (populationDC+i)->contrainteViole++;
            }
        }


        if(nbrNonCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrNonCohabitationConstraintes; j++)
            {
                noeud1=nonCohabitationConstraintes[j].nouedDepart;
                noeud2=nonCohabitationConstraintes[j].nouedArrive;

                if((populationDC+i)->phenotype[noeud1]== (populationDC+i)->phenotype[noeud2])
                {
                    (populationDC+i)->constraintVector[3]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationDC+i)->constraintVector[3]!=0)
            {
                (populationDC+i)->contrainteViole++;
            }
        }


    }

}

///********************************************************************************************
void writeOptimalSolutionInFileDC(partitionDC *solutionDominante,FILE* outputOptimalSolutionFileDC,
                                   int nbrRun, int bestSolutionIteration, float runTime,int ES)
{

    fprintf(outputOptimalSolutionFileDC,"DC | RUN NUMBER = %d | ITERATION N = %d | SOLUTION IDENTIFIER = %d | CUT SIZE INTRA = %d | CUT SIZE INTER = %d | SOMME TOTAL DES FLUX = %d | RUN TIME = %f | ES = %d | NBR VIOLATED CONSTRAINTE = %d | ",
            nbrRun, bestSolutionIteration, solutionDominante->id, solutionDominante->coutCoupe,(sommeTotalFlux - solutionDominante->coutCoupe),sommeTotalFlux ,runTime, ES, solutionDominante->contrainteViole);
 fprintf(outputOptimalSolutionFileDC," max_clusters = %d | max_sizeCluster = %d  |  cohabitation = %d | non cohabitation = %d\n",
            solutionDominante->constraintVector[0],solutionDominante->constraintVector[1],solutionDominante->constraintVector[2],
            solutionDominante->constraintVector[3]);


}

///********************************************************************************************
int compareCroissantFitnessDC (void const *a, void const *b)
{

    partitionDC const *pa = a;
    partitionDC const *pb = b;
    return  pb->coutCoupeNormalise - pa->coutCoupeNormalise;
}
/**********************************************************************************************




void writeOptimalSolutionInFileDC(partitionDC *solutionDominante,FILE* outputOptimalSolutionFileDC){
    int i;
    fprintf(outputOptimalSolutionFileDC,"id = %d | normalised cut size = %d | cut size intrA = %d cut size intEr = %d ",
            solutionDominante->id,solutionDominante->coutCoupeNormalise, solutionDominante->coutCoupe, (sommeTotalFlux-solutionDominante->coutCoupe));
    fprintf(outputOptimalSolutionFileDC,"le contraintes violees = %d \n", solutionDominante->contrainteViole);
    for(i=0; i<4; i++)
    {
        switch (i)
        {

        case 0:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileDC,"contrainte sur le nombre de clusters : OK \n");
            else fprintf(outputOptimalSolutionFileDC,"contrainte sur le nombre de clusters : NOT OK \n");
            break;
        case 1:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileDC,"contrainte sur la taille des clusters : OK \n");
            else fprintf(outputOptimalSolutionFileDC,"contrainte sur la taille des clusters : NOT OK \n");
            break;
        case 2:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileDC,"contrainte de cohabitation : OK \n");
            else fprintf(outputOptimalSolutionFileDC,"contrainte de cohabitation : NOT OK \n");
            break;
        case 3:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileDC,"contrainte Non Cohabitation : OK \n");
            else fprintf(outputOptimalSolutionFileDC,"contrainte Non Cohabitation : NOT OK \n");
            break;
        default:
            fprintf(outputOptimalSolutionFileDC,"Cette contrainte n'est pas prise en charge par le système \n");

        }

    }

}

void checkContrainstAndFitnessPenalizationDC(partitionDC *populationDC)
{

    int i,j,k,m=nbrNoeuds,tabTmp[m],tmp;

    for(i=0; i<taillePopulation; i++)
    {
        (populationDC+i)->nbrCluster=1; /// le nombre de cluster doit etre initialisé à 1 pour prendre en consédiration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes violées
        (populationDC+i)->contrainteViole=0;
        for(j=0;j<nbr_constraint;j++){(populationDC+i)->constraintVector[j]=0;}


        m=nbrNoeuds;
        ///copier le phenotype dans un deuxième tableau
        for(j=0; j<nbrNoeuds; j++)
        {
            tabTmp[j]=(populationDC+i)->phenotype[j];
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
                (populationDC+i)->nbrCluster++;
            }
        }
        if((populationDC+i)->nbrCluster > max_clusters || (populationDC+i)->nbrCluster < min_clusters)
        {
            (populationDC+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationDC+i)->contrainteViole++;
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
            (populationDC+i)->clustersSize[j]=0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if (tabTmp[j] == (populationDC+i)->phenotype[k])
                {
                    (populationDC+i)->clustersSize[j]=(populationDC+i)->clustersSize[j]+1;
                }
            }
            if((populationDC+i)->clustersSize[j]>max_sizeCluster || (populationDC+i)->clustersSize[j]< min_sizeCluster)
            {
                ///24.07.2015 ---> l'ajout de vecteur de contraintes = le nombre de cluster dont la taille a dépassé la borne
                (populationDC+i)->constraintVector[1]++;

            }
        }

        if((populationDC+i)->constraintVector[1] != 0)
        {
            (populationDC+i)->contrainteViole++;
        }
   }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0){
        checkCohabitationAndNonCohabitationConstrainteDC(populationDC);
    }

}


*/
