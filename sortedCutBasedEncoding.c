/**================= CutBased Encoding Genetic Algorithm ====================================================
* 18/04/2015 :
* la population initiale est maintenant prête
* la tache du jour :
*       1- programmee dans un premier temps la fonction qui permet de calculer le cout de coupe et le fitness
*       2- programmer la fonction de séléction naturelle
*       3- Remarque : la taille de la population doit être toujours paire
*
*       - Si : taillePopulation > pow(2, nbrPartie-1), avec une population initial sans redondance
*       on ne pourra pas attiendre la nombre de solutions nécessaires.
*============================================================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "typeDeclaration.h"
#include "./communesFunctions.h"
#include "./sortedCutBasedEncoding.h"
#include "./compilationConditionnelle.h"

extern char *cocyclesDeBase;
extern int taillePopulation,tauxReproduction,sommeTotalFlux,nbrNoeuds,sommeTotalFluxReal,
    nbr_constraint,nbrArretes,nbrParties,nbrNoeudsInCluster,regeneration;
///extern ul cocyclesDeBaseEntier[1000];
extern mpz_t gmpCocyclesDeBaseEntier[1000];
extern double tauxMutation ,averageTraffic;
extern int fluxMatrix[1000][1000],*fluxVector, maximumDistanceMatrix[1000][1000], maxFlow;
extern edge *edgeVector,cohabitationConstraintes[1000],nonCohabitationConstraintes[1000],*edgeVectorIntra;
extern int max_clusters , min_clusters,min_sizeCluster ,max_sizeCluster ;
extern int bestSolution, nbrApparition, nbrCohabitationConstraintes,nbrNonCohabitationConstraintes;
extern int nbrRun;



///**************************************************************************************
void affichePopulationIC(partitionIC* populationIC)
{
    int i,j,k,l;
    for(i=0; i<taillePopulation; i++)
    {

        printf("\nid = %d | cout de coupe = %d | fitness = %0.4f",
               (populationIC+i)->id,(populationIC+i)->coutCoupe, (populationIC+i)->fitness);

        /// affichage de vecteur de la partition
        /**       printf("\nla partition est  = \n");
               for(l=0;l<nbrParties;l++){
                   printf("%d ",(Ba+i)->genotypeEntier[l]);
               }
        */

        printf("\nla solution generer est  = \n");
        for(j=0; j<nbrParties-1; j++)
        {
            for(k=0; k<nbrNoeuds-1; k++)
            {
                printf("%d ",(populationIC+i)->cocycles[k+(j*-1)]);
            }
            printf("\n");
        }
        printf("\n");
        /// affichage de vecteur de la partition
        printf("\nla partition est  = \n");
        for(l=0; l<nbrArretes; l++)
        {
            printf("%d ",(populationIC+i)->genotype[l]);
        }
        printf("\n===============================================================================\n");

    }
}
///**************************************************************************************
/**===============================================================
Les cocycle de base : leur taille est [nbrArretes * (nbrNoeuds-1)]
Les cocycles : ce sont les individus de la population leur taille
               est donnée par [(nbrParties-1) * (nbrNoeuds-1)]
Les cocyclesXor : c'est le résultat de la somme XOR entre les
                 les cocycles de base correspondant au numero
                 1 dans les cocycles générés leur taille est
                 donnée par [(nbrParties-1)*nbrArretes]
 les genotype :c'est le résultat des traitement précédents
              donné par la somme OR entre les cocyclesXor
              la taille est [nbrArretes]
==================================================================*/
void calculerGenotypeIC(partitionIC * populationIC, int indiceFirstElt)
{

    int i,j,k,l,numeroCocyclXor;
    mpz_t cocyclesXor[1000];
    mpz_t genotypeTmp; mpz_init2(genotypeTmp,1000); ///mpz_init2(genotypeTmp,99999999999999999999999999999);
    mpz_t binaryTmp; mpz_init2(binaryTmp,1000); ///mpz_init2(binaryTmp,99999999999999999999999999999);
    for(j=0;j<nbrParties-1;j++) mpz_init2(cocyclesXor[j],1000); ///mpz_init2(cocyclesXor[j],99999999999999999999999999999);
    for(i=indiceFirstElt; i<taillePopulation; i++)
    {
        for(j=0; j<nbrArretes; j++)
        {
            (populationIC+i)->genotype[j]=0;
        }
        ///***********************************************************************************************
        for(j=0; j<nbrParties-1; j++)
        {
            numeroCocyclXor=0;
            for(k=0; k<nbrNoeuds-1; k++)
            ///for(k=nbrNoeuds-2; k>=0; k--)
            {
                if((populationIC+i)->cocycles[(j*(nbrNoeuds-1))+k]==1)
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
                mpz_ior(genotypeTmp , genotypeTmp , cocyclesXor[j]);
            }
        }
        j=nbrArretes-1;
        do
        {
            mpz_mod_ui(binaryTmp,genotypeTmp,2);
            (populationIC+i)->genotype[j]=mpz_get_ui(binaryTmp);
            mpz_fdiv_q_ui(genotypeTmp,genotypeTmp,2);
            j--;
        }
        while(j>=0);
    }


}
///**************************************************************************************
/// la séléction naturelle des individus pour la nouvelle génération
void naturalSelectionIC(partitionIC* populationIC1,partitionIC* populationIC2)
{
   int i,j,maxFitness;
    int tmpVectIndice[taillePopulation];
    for(i=0;i<taillePopulation;i++) tmpVectIndice[i]=-1;
    float lotterie,sommeFitness;

    for(i=0; i<tauxReproduction; i++){
        maxFitness = 0;
        for(j=0;j<taillePopulation;j++){
            if((populationIC1+maxFitness)->coutCoupeNormalise < (populationIC1+j)->coutCoupeNormalise
               && tmpVectIndice[j] == -1 ){
                maxFitness = j;

            }
        }
        tmpVectIndice[maxFitness] = 0;
        *(populationIC2+i) = *(populationIC1+maxFitness);
        ///printf("(populationIC2+%d)->coutCoupeNormalise = %d \n",i,(populationIC2+i)->coutCoupeNormalise);
        ///mon_sleep(pause);
    }
    ///printf("\n\n\n");

    ///qsort (populationIC2, taillePopulation, sizeof *populationIC2, compareCroissantFitnessEE);
    ///sortingPopulationByFintness(populationIC2,tmpPopulation);

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
            sommeFitness = sommeFitness + (populationIC1+j)->fitness;
            if (lotterie<=sommeFitness)
            {
                *(populationIC2+i) = *(populationIC1+j);
                break;
            }
            else if( j==taillePopulation-1)
            {
                *(populationIC2+i) = *(populationIC1+j);
            }
        }
    }
    ///free(tmpPopulation);

}
///**************************************************************************************
///le croisement des solution séléctionner
void crossOverIC(partitionIC* populationIC1, partitionIC* populationIC2)
{

    int i = 0,j;
    int choixInd1, choixInd2, choixLocus=0,nbrCocyles;
    ///tailleCocycles = (nbrParties-1)*(nbrNoeuds-1);
    nbrCocyles = nbrParties-1;
    ///***********************************************************************************
    for(i=0; i<tauxReproduction; i++)
    {
        *(populationIC2+i) = *(populationIC1+i);
        (populationIC2+i)->id = i;
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
        /// c'est nbrPartie = 2 ce qui représente la valeur minimale que cette variable peut prendre
        /// donc nbrCocyles = 1 ==> rnd(0, nbrCocyles-1) <=> rnd(0, 0) ce qui va générer erreur
        if(nbrCocyles != 1)
        {
            choixLocus = rnd(0,nbrCocyles-1);
        } /// -1 pour s'arrêter à moins 2}

        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            mpz_set((populationIC2+i)->genotypeEntier[j] , (populationIC1+choixInd1)->genotypeEntier[j]);
            mpz_set((populationIC2+i+1)->genotypeEntier[j] , (populationIC1+choixInd2)->genotypeEntier[j]);
        }
        ///*****************************************************************************************
        for(j=choixLocus+1; j<nbrCocyles; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            mpz_set((populationIC2+i)->genotypeEntier[j] , (populationIC1+choixInd2)->genotypeEntier[j]);
            mpz_set((populationIC2+i+1)->genotypeEntier[j] , (populationIC1+choixInd1)->genotypeEntier[j]);
        }
        /// affectation de nouveau indice
        (populationIC2+i)->id = i;
        (populationIC2+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }
    generatePopulationWithoutRedandancyIC(populationIC2, taillePopulation-regeneration);
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
        /// c'est nbrPartie = 2 ce qui représente la valeur minimale que cette variable peut prendre
        /// donc nbrCocyles = 1 ==> rnd(0, nbrCocyles-1) <=> rnd(0, 0) ce qui va générer erreur
        if(nbrCocyles != 1)
        {
            choixLocus = rnd(0,nbrCocyles-1);
        } /// -1 pour s'arrêter à moins 2}

        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            mpz_set((populationIC2+i)->genotypeEntier[j] , (populationIC1+choixInd1)->genotypeEntier[j]);
            mpz_set((populationIC2+i+1)->genotypeEntier[j] , (populationIC1+choixInd2)->genotypeEntier[j]);
        }
        ///*****************************************************************************************
        for(j=choixLocus+1; j<nbrCocyles; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            mpz_set((populationIC2+i)->genotypeEntier[j] , (populationIC1+choixInd2)->genotypeEntier[j]);
            mpz_set((populationIC2+i+1)->genotypeEntier[j] , (populationIC1+choixInd1)->genotypeEntier[j]);
        }
        /// affectation de nouveau indice
        (populationIC2+i)->id = i;
        (populationIC2+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }

#endif
}
///**************************************************************************************
void mutationIC(partitionIC* populationIC)
{
    int i,numeroGeneMute,nbeCocyles=nbrParties-1;
    double applicationOfMutation;
    mpz_t borneSup; mpz_init2(borneSup,1000); ///mpz_init2(borneSup,99999999999999999999999999999);
    mpz_ui_pow_ui(borneSup,2,nbrNoeuds-1);

    gmp_randstate_t gmpRandState; /* Random generator state object */
    gmp_randinit_default(gmpRandState);
    gmp_randseed_ui(gmpRandState, time(NULL));

    for(i=tauxReproduction; i<taillePopulation; i++)
    {
#if Windows
        applicationOfMutation = drand48ForWindows(0,1001);
#else
        applicationOfMutation = drand48();
#endif
        if(applicationOfMutation <= tauxMutation)
        {
            /// choisir un gène aléatoirement pour lui appliquer la mutation
            numeroGeneMute = rnd(0,nbeCocyles);
            mpz_urandomm((populationIC+i)->genotypeEntier[numeroGeneMute], gmpRandState, borneSup);
        }
    }
}
///**************************************************************************************
float testerLaSommeDesFitnessIC(partitionIC* populationIC)
{
    int i;
    float sommeFitness = 0;
    for(i=0; i<taillePopulation ; i++)
    {
        sommeFitness = sommeFitness + (populationIC+i)->fitness;
    }
    return sommeFitness;
}
///**************************************************************************************
int findTheBestSolutionIC(partitionIC *populationIC)
{
    float maxFitness = (populationIC+0)->fitness;
    int i,indice=0;
    for(i=0; i<taillePopulation ; i++)
    {
        if(maxFitness < (populationIC+i)->fitness)
        {
            maxFitness = (populationIC+i)->fitness;
            indice = i;
        }
    }
    return indice;

}
///**************************************************************************************
void displayTheBestSolutionIC(partitionIC* solutionDominante)
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
void writeSolutionInFileIC(partitionIC *populationIC, FILE *outputFilePop ,int iteration)
{

    int i,j,k;
    for(i=0; i<taillePopulation; i++)
    {

        ///écriture de la meilleur solution dans le fichier output
        fprintf(outputFilePop,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,(populationIC+i)->id,
                (populationIC+i)->coutCoupe,(populationIC+i)->fitness,(populationIC+i)->coutCoupeNormalise,
                (populationIC+i)->contrainteViole,(populationIC+i)->nbrCluster);
        for(j=0; j<nbrArretes; j++)
        {
            fprintf(outputFilePop,"%d ",(populationIC+i)->genotype[j]);
        }
        fprintf(outputFilePop,"\t\t\t\t");
        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(outputFilePop,"%d ",(populationIC+i)->phenotype[j]);
        }
        fprintf(outputFilePop,"\t\t\t\t");
        for(j=0; j<nbrParties-1; j++)
            for(k=0;k<nbrNoeuds-1;k++)
                fprintf(outputFilePop,"%d ",(populationIC+i)->cocycles[j*(nbrNoeuds-1)+k]);
        fprintf(outputFilePop,"\n");
    }
    fprintf(outputFilePop,"\n\n");
}
///**************************************************************************************
void writeBestSolutionInFileIC(partitionIC *solutionDominante, FILE *outputFile, int iteration)
{

    int i,j;
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
    fprintf(outputFile,"\t\t\t\t");
    for(i=0; i<nbrParties-1; i++)
        for(j=0;j<nbrNoeuds-1;j++)
            fprintf(outputFile,"%d ",solutionDominante->cocycles[i*(nbrNoeuds-1)+j]);
    fprintf(outputFile,"\n");

}

///**************************************************************************************
void getPartitionFromSolutionIC(partitionIC *populationIC)
{
    int i,j,k,nd,na,numeroClustre,nbrIntra; /// nd = noued de départ , na = noued d'arrivé

    for(i=0; i<taillePopulation; i++)
    {
        /// initialisation des phénotype = partition par les numero des noeuds correspondant

        for(j=0; j<nbrNoeuds; j++)
        {
            (populationIC+i)->phenotype[j] = j;
        }
        /// 01/10/2017 : amélioration du temps d'exécution du codage binaire
        nbrIntra = 0;
        for(j=0; j<nbrArretes; j++){
            if((populationIC+i)->genotype[j]==0)
            {
               edgeVectorIntra[nbrIntra] = edgeVector[j];
               nbrIntra++;
            }

        }
        /// détermination de la parition
        for(j=0; j<nbrIntra; j++){

                nd = (edgeVectorIntra+j)->nouedDepart;
                na = (edgeVectorIntra+j)->nouedArrive;
               if((populationIC+i)->phenotype[na]!= na) /// il a été modifier auparavant
                {

                    if((populationIC+i)->phenotype[na] < (populationIC+i)->phenotype[nd])
                    {
                        ///dans le cas ou le noued d'arrivé est inférieur au noued depart ==> toutes les occurence du noued
                        /// de départ prendront ceux du noeud d'arrivé
                        numeroClustre = (populationIC+i)->phenotype[nd];
                        for(k=0; k<nbrNoeuds; k++)
                        {
                            if((populationIC+i)->phenotype[k] == numeroClustre)
                            {
                                (populationIC+i)->phenotype[k] = (populationIC+i)->phenotype[na];
                            }
                        }
                    }
                    else if((populationIC+i)->phenotype[na] > (populationIC+i)->phenotype[nd])
                    {
                        numeroClustre = (populationIC+i)->phenotype[na];
                        for(k=0; k<nbrNoeuds; k++)
                        {
                            if((populationIC+i)->phenotype[k] == numeroClustre)
                            {
                                (populationIC+i)->phenotype[k] = (populationIC+i)->phenotype[nd];
                            }
                        }
                    }
                }
                else
                {
                    (populationIC+i)->phenotype[na] = (populationIC+i)->phenotype[nd];
                }
            }
    }

/**  19/11/2017
    int i,j,k,nd,na,numeroClustre; /// nd = noued de départ , na = noued d'arrivé

    for(i=0; i<taillePopulation; i++)
    {
        /// initialisation des phénotype = partition par les numero des noeuds correspondant
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationIC+i)->phenotype[j] = j;
        }

        /// détermination de la parition
        ///for(j=0; j<nbrArretes; j++)
        for(j=nbrArretes-1; j>=0; j--)
        {
            if((populationIC+i)->genotype[j]==0)
            {
                na = edgeVector[j].nouedArrive;
                nd = edgeVector[j].nouedDepart;

                if((populationIC+i)->phenotype[na]!= na) /// il a été modifier auparavant
                {

                    if((populationIC+i)->phenotype[na] < (populationIC+i)->phenotype[nd])
                    {
                        ///dans le cas ou le noued d'arrivé est inférieur au noued d'arrivé ==> toutes les occurence du noued
                        /// de départ prendront ceux du noeud d'arrivé
                        numeroClustre = (populationIC+i)->phenotype[nd];
                        for(k=0; k<nbrNoeuds; k++)
                        {
                            if((populationIC+i)->phenotype[k] == numeroClustre)
                            {
                                (populationIC+i)->phenotype[k] = (populationIC+i)->phenotype[na];
                            }
                        }
                    }
                    else if((populationIC+i)->phenotype[na] > (populationIC+i)->phenotype[nd])
                    {
                        numeroClustre = (populationIC+i)->phenotype[na];
                        for(k=0; k<nbrNoeuds; k++)
                        {
                            if((populationIC+i)->phenotype[k] == numeroClustre)
                            {
                                (populationIC+i)->phenotype[k] = (populationIC+i)->phenotype[nd];
                            }
                        }
                    }
                }
                else
                {
                    (populationIC+i)->phenotype[na] = (populationIC+i)->phenotype[nd];
                }
            }
        }
    }
    */
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

void calculCoutCoupeEtFitnessIC(partitionIC* populationIC)
{
    int i,j,k,sommeTotalCoutCoupe = 0;

    for(i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (populationIC+i)->coutCoupe = 0; /// réinitialisation de la valriable après chaque itération (individus de la population)
        for(j=0; j<nbrNoeuds-1; j++)
        {
            for(k=j+1; k<nbrNoeuds; k++)
            {
                if((populationIC+i)->phenotype[j] == (populationIC+i)->phenotype[k] && fluxMatrix[j][k] > 0)
                {
                    (populationIC+i)->coutCoupe = (populationIC+i)->coutCoupe +  fluxMatrix[j][k];
                }
            }
        }
        (populationIC+i)->coutCoupeNormalise = (populationIC+i)->coutCoupe + ((nbr_constraint - (populationIC+i)->contrainteViole)*sommeTotalFlux);
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationIC+i)->coutCoupeNormalise;
    }

#if scaling

    FILE *file;
    file = fopen("D:/WorkingFolderAgProgramming/benchmark/sortedCutBasedEncoding/expectedValueIC.txt","w");
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
            varianceCoutDeCoupeNormalise + pow(((populationIC+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);
    /// tester la valeur ajouté
    for(i=0; i<taillePopulation; i++)
    {

        /// calcul de sigma truncation
        (populationIC+i)->expectedValue = (populationIC+i)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);

        /**
                    if(iteration == 1){
                        fprintf(file,"%d \t %d \t %0.2f \n",i,(populationIC+i)->coutCoupeNormalise,(populationIC+i)->expectedValue);
                    }
        */

        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationIC+i)->expectedValue;
    }

    /// calcule du fitness
    /// le 26/11/2015 : calcul du fitness des solution avec leur expectedValues
    for(i=0; i<taillePopulation ; i++)
    {
        (populationIC+i)->fitness = (float)((populationIC+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);
    }

#else
    for(i=0; i<taillePopulation ; i++)
    {
        (populationIC+i)->fitness = (float)((populationIC+i)->coutCoupeNormalise)/(float)(sommeTotalCoutCoupe);
    }
#endif// scaling


}

///**************************************************************************************
/// calculer le cout de coupe de chacune des solution dans la population
void calculCoutCoupeEtFitnessWithFlowVectorIC(partitionIC* populationIC)
{
#if mouchard
    printf("calculCoutCoupeEtFitness ...\n");
#endif // mouchard
    int i,j,sommeTotalCoutCoupe = 0;
    for(i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (populationIC+i)->coutCoupe = 0;

        for(j=0; j<nbrArretes; j++) ///parcourir le vecteurSolution =  solution codée en binaire
        {
            if((populationIC+i)->genotype[j]==0)  /// on cherche à maximiser les intraCluster donc les arrêtes 0 != 1
            {
                (populationIC+i)->coutCoupe = (populationIC+i)->coutCoupe + fluxVector[j];
            }
        }
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationIC+i)->coutCoupe;
    }
    /// calcule du fitness
    for(i=0; i<taillePopulation ; i++)
    {
        ///(populationIC+i)->fitness = (float)((populationIC+i)->coutCoupe)/sommeTotalCoutCoupe;
        (populationIC+i)->fitness = (float)((populationIC+i)->coutCoupe)/(float)(sommeTotalCoutCoupe);
    }
}


///=====================================================================================
/// Cette partie concerne la génération des partie sans répétition, le choix de la méthode
/// de génération se fait au niveau de la fonction main
///=====================================================================================
///**************************************************************************************
void generatePopulationWithoutRedandancyIC(partitionIC* populationIC, int indiceFirstElt)
{
    genererSolutionEntiere(populationIC, indiceFirstElt);
    genererLesSolutionBinaireIC(populationIC, indiceFirstElt);
    calculerGenotypeIC(populationIC, indiceFirstElt);
}
/**
void generatePopulationWithoutRedandancyIC(partitionIC* populationIC)
{
    int i;
    for(i=0; i<taillePopulation; i++)
    {
        do
        {
            genererSolutionEntiere((populationIC+i)->genotypeEntier);
            trieGenotypeEntier((populationIC+i)->genotypeEntier);
        }
        while(existanceDeSolution(i, populationIC));

        (populationIC+i)->id =i;
        (populationIC+i)->coutCoupe =0;
        (populationIC+i)->fitness =0;


    }
    genererLesSolutionBinaireIC(populationIC);
}*/
///**************************************************************************************
void genererSolutionEntiere(partitionIC* populationIC, int indiceFirstElt)
{
    int i=0,j,k;
    mpz_t tmp; mpz_init2(tmp,1000);///mpz_init2(tmp,99999999999999999999999999999);
    mpz_t borneSup; mpz_init2(borneSup,1000); ///mpz_init2(borneSup,99999999999999999999999999999);
    mpz_ui_pow_ui(borneSup,2,nbrNoeuds-1);
    int nbrCocycles = nbrParties-1;
    gmp_randstate_t gmpRandState; /* Random generator state object */
    ///gmp_randinit_default(gmpRandState);
    gmp_randinit_mt(gmpRandState);
    gmp_randseed_ui(gmpRandState, time(NULL));


    for(i=indiceFirstElt;i<taillePopulation;i++){
        (populationIC+i)->id =i;
        (populationIC+i)->coutCoupe =0;
        (populationIC+i)->fitness =0;
        j=0;
        while(j<nbrCocycles){
            mpz_urandomm((populationIC+i)->genotypeEntier[j], gmpRandState, borneSup);
            mpz_urandomm(tmp, gmpRandState, borneSup);
            for(k=0;k<j;k++){
                if(mpz_cmp((populationIC+i)->genotypeEntier[j], tmp) == 0){break;}
            }
            if(k==j){
                    mpz_set((populationIC+i)->genotypeEntier[j],tmp);
                    j++;
            }
        }
    }

    ///mpz_clears(tmp,borneSup);
    ///trieGenotypeEntier(A);
}
/**
void genererSolutionEntiere(ul A[])
{
    int i;
    ul borneSup = pow(2,nbrNoeuds-1)-1;

    for(i=0; i<nbrParties-1; i++)
    {
        A[i] = rndUl(0,borneSup);
    }
}
*/
///**************************************************************************************
int existanceDeSolution(int indexNewSolution, partitionIC* populationIC )
{

    int i=0,j,Exist = 0,theSameElement;

    while(i<indexNewSolution && Exist == 0)
    {
        theSameElement = 0;
        for(j=0; j<nbrParties-1; j++)
        {
            if(mpz_cmp((populationIC+i)->genotypeEntier[j] ,(populationIC+indexNewSolution)->genotypeEntier[j]) ==0)
            {
                theSameElement++;
            }
        }
        if(theSameElement==nbrParties-1)
        {
            Exist = 1;
            break;
        }
        i++;
    }
    return Exist;
}
///**************************************************************************************
void genererLesSolutionBinaireIC(partitionIC *populationIC, int indiceFirstElt)
{
    int i,j,k,nbrCocyle=nbrParties-1, cocyclesSize = nbrCocyle*(nbrNoeuds-1);
    mpz_t tmp; mpz_init2(tmp,1000); ///mpz_init2(tmp,99999999999999999999999999999);
    mpz_t binaryTmp; mpz_init2(binaryTmp,1000); ///mpz_init2(binaryTmp,99999999999999999999999999999);

    for(i=indiceFirstElt;i<taillePopulation;i++){
        if((populationIC+i)->genotype == NULL){
            (populationIC+i)->genotype = (int*)malloc(nbrArretes*sizeof(int));
            if((populationIC+i)->genotype == NULL){printf("Memory allocation for genotypeIC failled \n"); exit(EXIT_FAILURE);}
        }
    }

    for(i=indiceFirstElt;i<taillePopulation;i++){
        if((populationIC+i)->cocycles == NULL){
            (populationIC+i)->cocycles = (int*)malloc(cocyclesSize*sizeof(int));
            if((populationIC+i)->cocycles == NULL){printf("Memory allocation for cocyclesIC failled\n");exit(EXIT_FAILURE);}
        }
    }



    for(i=indiceFirstElt; i<taillePopulation; i++)
    {
        for(j=0; j<nbrCocyle; j++)
        {
            mpz_set(tmp,(populationIC+i)->genotypeEntier[j]);
            for(k=0; k<nbrNoeuds-1; k++)
            {
                /**
                    nbrNoeuds-2-k : on la taille de chaque partie du cocyle = nbrNoeuds-1
                    ça veut dire de 0-->nbrNoeuds-2, et pour commencer à partir de la droite
                    nbrNoeuds-2-k ou k représente l'indice de la boucle.
                */
                mpz_mod_ui(binaryTmp,tmp,2);
                (populationIC+i)->cocycles[(j*(nbrNoeuds-1))+k] = mpz_get_ui(binaryTmp);
                mpz_fdiv_q_ui(tmp,tmp,2);
            }
        }
    }
}

///======================================================================================
/// the main programm
void cutBasedEncodingWithoutRedandancyIC(int nbrGeneration,FILE *outputFileIC,FILE *outputFilePopIC,FILE *outputOptimalSolutionFileIC,partitionIC *populationIC1,partitionIC *populationIC2, partitionIC *solutionDominante,
                                         int iteration , int *bestSolutionIteration , int *nbrApparition)
{

        naturalSelectionIC(populationIC1,populationIC2);
        crossOverIC(populationIC2, populationIC1);
        mutationIC(populationIC1);
        ///************************************************
        genererLesSolutionBinaireIC(populationIC1,0);
        ///************************************************
        calculerGenotypeIC(populationIC1,0);
        getPartitionFromSolutionIC(populationIC1);
        checkContrainstAndFitnessPenalizationIC(populationIC1);
        calculCoutCoupeEtFitnessIC(populationIC1);

#if writingPopulationInFile
        writeSolutionInFileIC(populationIC1,outputFilePopIC,iteration);
#endif // writingPopulationInFile

        bestSolution=findTheBestSolutionIC(populationIC1);
#if writingPopulationInFile
        writeBestSolutionInFileIC((populationIC1+bestSolution),outputFileIC,iteration);
#endif // writingPopulationInFile
        if((populationIC1+bestSolution)->coutCoupeNormalise > solutionDominante->coutCoupeNormalise)
        {
            *(solutionDominante) = *(populationIC1+bestSolution);
            nbrApparition=1;
            bestSolutionIteration = iteration;
        }
        else
        {
            nbrApparition++;
        }

}

///**************************************************************************************
void checkContrainstAndFitnessPenalizationIC(partitionIC *populationIC)
{

    int i,j,k;

    for(i=0; i<taillePopulation; i++)
    {
        (populationIC+i)->nbrCluster=0; /// le nombre de cluster doit etre initialisé à 1 pour prendre en consédiration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes violées
        (populationIC+i)->contrainteViole=0;
        for(j=0; j<4; j++)
        {
            (populationIC+i)->constraintVector[j]=0;
        }
        ///**********************************************************************************


        ///trier le tableau résulant
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationIC+i)->clustersSize[j] =0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if ((populationIC+i)->phenotype[k] == j)
                {
                    (populationIC+i)->clustersSize[j]++;
                }
            }
            if ((populationIC+i)->clustersSize[j] !=0)
            {
                (populationIC+i)->nbrCluster++;
            }
        }

        if((populationIC+i)->nbrCluster > max_clusters || (populationIC+i)->nbrCluster < min_clusters)
        {
            (populationIC+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationIC+i)->contrainteViole++;
        }

        for(j=0; j<nbrNoeuds; j++)
        {
            if((populationIC+i)->clustersSize[j]!=0)
            {
                if((populationIC+i)->clustersSize[j]>max_sizeCluster || (populationIC+i)->clustersSize[j]< min_sizeCluster)
                {
                    (populationIC+i)->constraintVector[1]++;

                }
            }

        }

        if((populationIC+i)->constraintVector[1] != 0)
        {
            (populationIC+i)->contrainteViole++;
        }
    }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0)
    {
        checkCohabitationAndNonCohabitationConstrainteIC(populationIC);
    }

}

///**************************************************************************************
void checkCohabitationAndNonCohabitationConstrainteIC(partitionIC *populationIC)
{

    int i,j,noeud1,noeud2;
    for(i=0; i<taillePopulation; i++)
    {
        /// Cette tache est déjà réalisé au niveau de la fonction checkCohabitationAndNonCohabitationConstrainteBE
		/**
		(populationIC+i)->constraintVector[2]=0;
        (populationIC+i)->constraintVector[3]=0;
		*/
        if(nbrCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrCohabitationConstraintes; j++)
            {
                noeud1=cohabitationConstraintes[j].nouedDepart;
                noeud2=cohabitationConstraintes[j].nouedArrive;

                if((populationIC+i)->phenotype[noeud1]!= (populationIC+i)->phenotype[noeud2])
                {
                    (populationIC+i)->constraintVector[2]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationIC+i)->constraintVector[2]!=0)
            {
                (populationIC+i)->contrainteViole++;
            }
        }


        if(nbrNonCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrNonCohabitationConstraintes; j++)
            {
                noeud1=nonCohabitationConstraintes[j].nouedDepart;
                noeud2=nonCohabitationConstraintes[j].nouedArrive;
                if((populationIC+i)->phenotype[noeud1]== (populationIC+i)->phenotype[noeud2])
                {
                    (populationIC+i)->constraintVector[3]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationIC+i)->constraintVector[3]!=0)
            {
                (populationIC+i)->contrainteViole++;
            }
        }


    }

}

///********************************************************************************************
void writeOptimalSolutionInFileIC(partitionIC *solutionDominante,FILE* outputOptimalSolutionFileIC,
                                    int nbrRun, int bestSolutionIteration, float runTime, int ES)
{
    int i;
    fprintf(outputOptimalSolutionFileIC,"RUN NUMBER = %d | ITERATION N = %d | SOLUTION IDENTIFIER = %d | CUT SIZE INTRA = %d | CUT SIZE INTER = %d | SOMME TOTAL DES FLUX = %d | RUN TIME = %f | ES = %d | NBR VIOLATED CONSTRAINTE = %d | ",
            nbrRun, bestSolutionIteration, solutionDominante->id, solutionDominante->coutCoupe,(sommeTotalFlux - solutionDominante->coutCoupe),sommeTotalFlux ,runTime, ES, solutionDominante->contrainteViole);
fprintf(outputOptimalSolutionFileIC," max_clusters = %d | max_sizeCluster = %d  |  cohabitation = %d | non cohabitation = %d\n",
            solutionDominante->constraintVector[0],solutionDominante->constraintVector[1],solutionDominante->constraintVector[2],
            solutionDominante->constraintVector[3]);

}

///******************************************************************************************
int compareCroissantFitnessIC (void const *a, void const *b)
{

    partitionIC const *pa = a;
    partitionIC const *pb = b;
    return  pb->coutCoupeNormalise - pa->coutCoupeNormalise;
}
/**************************************************************************************

///******************************************************************************************
void writeOptimalSolutionInFileIC(partitionIC *solutionDominante,FILE* outputOptimalSolutionFileIC){
    int i;
    fprintf(outputOptimalSolutionFileIC,"id = %d | normalised cut size = %d | cut size intrA = %d cut size intEr = %d ",
            solutionDominante->id,solutionDominante->coutCoupeNormalise, solutionDominante->coutCoupe, (sommeTotalFlux-solutionDominante->coutCoupe));
    fprintf(outputOptimalSolutionFileIC,"le contraintes violees = %d \n", solutionDominante->contrainteViole);
    for(i=0; i<4; i++)
    {
        switch (i)
        {

        case 0:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileIC,"contrainte sur le nombre de clusters : OK \n");
            else fprintf(outputOptimalSolutionFileIC,"contrainte sur le nombre de clusters : NOT OK \n");
            break;
        case 1:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileIC,"contrainte sur la taille des clusters : OK \n");
            else fprintf(outputOptimalSolutionFileIC,"contrainte sur la taille des clusters : NOT OK \n");
            break;
        case 2:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileIC,"contrainte de cohabitation : OK \n");
            else fprintf(outputOptimalSolutionFileIC,"contrainte de cohabitation : NOT OK \n");
            break;
        case 3:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileIC,"contrainte Non Cohabitation : OK \n");
            else fprintf(outputOptimalSolutionFileIC,"contrainte Non Cohabitation : NOT OK \n");
            break;
        default:
            fprintf(outputOptimalSolutionFileIC,"Cette contrainte n'est pas prise en charge par le système \n");

        }

    }

}


void checkContrainstAndFitnessPenalizationIC(partitionIC *populationIC)
{

    int i,j,k,m=nbrNoeuds,tabTmp[m],tmp;

    for(i=0; i<taillePopulation; i++)
    {
        (populationIC+i)->nbrCluster=1; /// le nombre de cluster doit etre initialisé à 1 pour prendre en consédiration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes violées
        (populationIC+i)->contrainteViole=0;
        for(j=0;j<nbr_constraint;j++){(populationIC+i)->constraintVector[j]=0;}


        m=nbrNoeuds;
        ///copier le phenotype dans un deuxième tableau
        for(j=0; j<nbrNoeuds; j++)
        {
            tabTmp[j]=(populationIC+i)->phenotype[j];
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
                (populationIC+i)->nbrCluster++;
            }
        }
        if((populationIC+i)->nbrCluster > max_clusters || (populationIC+i)->nbrCluster < min_clusters)
        {
            (populationIC+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationIC+i)->contrainteViole++;
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
            (populationIC+i)->clustersSize[j]=0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if (tabTmp[j] == (populationIC+i)->phenotype[k])
                {
                    (populationIC+i)->clustersSize[j]=(populationIC+i)->clustersSize[j]+1;
                }
            }
            if((populationIC+i)->clustersSize[j]>max_sizeCluster || (populationIC+i)->clustersSize[j]< min_sizeCluster)
            {
                ///24.07.2015 ---> l'ajout de vecteur de contraintes = le nombre de cluster dont la taille a dépassé la borne
                (populationIC+i)->constraintVector[1]++;

            }
        }

        if((populationIC+i)->constraintVector[1] != 0)
        {
            (populationIC+i)->contrainteViole++;
        }
   }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0){
        checkCohabitationAndNonCohabitationConstrainteIC(populationIC);
    }

}

*/

