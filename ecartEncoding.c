/**================= CutBased Encoding Genetic Algorithm ====================================================
* 18/04/2015 :
* la population initiale est maintenant prête
* la tache du jour :
*       1- programmCGE dans un premier temps la fonction qui permet de calculer le cout de coupe et le fitness
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
#include "./ecartEncoding.h"
#include "./compilationConditionnelle.h"
#define tousLeRestAZero 1


extern char *cocyclesDeBase;
extern int taillePopulation,tauxReproduction,sommeTotalFlux,nbrNoeuds,
       nbr_constraint,nbrArretes,nbrParties,nbrNoeudsInCluster,regeneration;
extern double tauxMutation ,averageTraffic;
extern int fluxMatrix[1000][1000],*fluxVector,maximumDistanceMatrix[1000][1000], maxFlow, sommeTotalFluxReal;
///extern ul cocyclesDeBasCGEntier[1000];
extern mpz_t gmpCocyclesDeBaseEntier[1000];
extern edge *edgeVector,cohabitationConstraintes[1000],nonCohabitationConstraintes[1000],*edgeVectorIntra;
extern int max_clusters , min_clusters,min_sizeCluster ,max_sizeCluster ;
extern int bestSolution, nbrApparition, nbrCohabitationConstraintes,nbrNonCohabitationConstraintes;
extern int nbrRun;



///**************************************************************************************
void affichePopulationCGE(partitionCGE* populationCGE){
    printf("affichePopulation...\n");
    int i,j,k,l;
    for(i=0; i<taillePopulation; i++)
    {

        printf("\nid = %d | cout de coupe = %d | fitness = %0.4f",
               (populationCGE+i)->id,(populationCGE+i)->coutCoupe, (populationCGE+i)->fitness);

        /// affichage de vecteur de la partition
        /**       printf("\nla partition est  = \n");
               for(l=0;l<nbrParties;l++){
                   printf("%d ",(populationCGE+i)->genotypeEcart[l]);
               }
        */

        printf("\nla solution generer est  = \n");
        for(j=0; j<nbrParties; j++)
        {
            for(k=0; k<nbrNoeuds-1; k++)
            {
                printf("%d ",(populationCGE+i)->cocycles[k+(j*nbrParties)]);
            }
            printf("\n");
        }
        printf("\n");
        /// affichage de vecteur de la partition
        printf("\nla partition est  = \n");
        for(l=0; l<nbrArretes; l++)
        {
            printf("%d ",(populationCGE+i)->genotype[l]);
        }
        printf("\n===============================================================================\n");

    }
}
///**************************************************************************************
void calculerGenotypeCGE(partitionCGE * populationCGE){

    int i,j,k,numeroCocyclXor,nbrCocycles=nbrParties-1;
    ///ul cocyclesXor[1000]={0},genotypeTmp;
    static mpz_t cocyclesXor[1000];
    static mpz_t genotypeTmp;mpz_init(genotypeTmp);
    static mpz_t binaryTmp;mpz_init(binaryTmp);

    for(j=0;j<nbrParties-1;j++)mpz_init(cocyclesXor[j]);
    for(i=0; i<taillePopulation; i++){

        for(j=0; j<nbrArretes; j++){
            (populationCGE+i)->genotype[j]=0;
        }
        ///***********************************************************************************************
        for(j=0; j<nbrCocycles; j++){
            numeroCocyclXor=0;
            for(k=nbrNoeuds-2; k>=0; k--){
                if((populationCGE+i)->cocycles[(j*(nbrNoeuds-1))+k]==1){
                    if(numeroCocyclXor==0){
                        mpz_set(cocyclesXor[j],gmpCocyclesDeBaseEntier[k]);
                        numeroCocyclXor++;
                    }
                    else{
                        mpz_xor(cocyclesXor[j] , cocyclesXor[j] , gmpCocyclesDeBaseEntier[k] );
                    }
                }
            }

            if(j==0){
                    mpz_set(genotypeTmp ,cocyclesXor[j]);
            }
            else{
                    mpz_ior(genotypeTmp , genotypeTmp , cocyclesXor[j]);
            }

        }

        j=nbrArretes-1;
        do
        {
            mpz_mod_ui(binaryTmp,genotypeTmp,2);
            (populationCGE+i)->genotype[j]=mpz_get_ui(binaryTmp);
            mpz_fdiv_q_ui(genotypeTmp,genotypeTmp,2);
            j--;
        }
        while(j>=0);
    }

}
///**************************************************************************************
/// la séléction naturelle des individus pour la nouvelle génération
void naturalSelectionCGE(partitionCGE* populationCGE1,partitionCGE* populationCGE2){
    int i,j,maxFitness;
    int tmpVectIndice[taillePopulation];
    for(i=0;i<taillePopulation;i++) tmpVectIndice[i]=-1;
    float lotterie,sommeFitness;

    for(i=0; i<tauxReproduction; i++){
        maxFitness = 0;
        for(j=0;j<taillePopulation;j++){
            if((populationCGE1+maxFitness)->coutCoupeNormalise < (populationCGE1+j)->coutCoupeNormalise
               && tmpVectIndice[j] == -1 ){
                maxFitness = j;

            }
        }
        tmpVectIndice[maxFitness] = 0;
        *(populationCGE2+i) = *(populationCGE1+maxFitness);
        ///printf("(populationCGE2+%d)->coutCoupeNormalise = %d \n",i,(populationCGE2+i)->coutCoupeNormalise);
        ///mon_slCGEp(pause);
    }
    ///printf("\n\n\n");

    ///qsort (populationCGE2, taillePopulation, sizeof *populationCGE2, compareCroissantFitnessCGE);
    ///sortingPopulationByFintness(populationCGE2,tmpPopulation);

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
            sommeFitness = sommeFitness + (populationCGE1+j)->fitness;
            if (lotterie<=sommeFitness)
            {
                *(populationCGE2+i) = *(populationCGE1+j);
                break;
            }
            else if( j==taillePopulation-1)
            {
                *(populationCGE2+i) = *(populationCGE1+j);
            }
        }
    }
    ///frCGE(tmpPopulation);

}
///**************************************************************************************
///le croisement des solution séléctionner
void crossOverCGE(partitionCGE* populationCGE1, partitionCGE* populationCGE2){

    int i = 0,j;
    int choixInd1, choixInd2, choixLocus=0,nbrCocyles;
    ///tailleCocycles = (nbrParties-1)*(nbrNoeuds-1);
    nbrCocyles = nbrParties-1;
    ///***********************************************************************************
    for(i=0; i<tauxReproduction; i++)
    {
        *(populationCGE2+i) = *(populationCGE1+i);
        (populationCGE2+i)->id = i;
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

        if(nbrCocyles != 1){choixLocus = rnd(0,nbrCocyles-1);} /// -1 pour s'arrêter à moins 2

        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            mpz_set((populationCGE2+i)->genotypeEcart[j], (populationCGE1+choixInd1)->genotypeEcart[j]);
            mpz_set((populationCGE2+i+1)->genotypeEcart[j] , (populationCGE1+choixInd2)->genotypeEcart[j]);
        }
        ///*****************************************************************
        for(j=choixLocus+1; j<nbrCocyles; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            mpz_set((populationCGE2+i)->genotypeEcart[j] , (populationCGE1+choixInd2)->genotypeEcart[j]);
            mpz_set((populationCGE2+i+1)->genotypeEcart[j] , (populationCGE1+choixInd1)->genotypeEcart[j]);
        }
        /// affectation de nouveau indice
        (populationCGE2+i)->id = i;
        (populationCGE2+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }
    generatePopulationWithoutRedandancyCGE(populationCGE2, taillePopulation-regeneration);
#else
    while (i < taillePopulation)
    {
        //choixInd1 = rnd(tauxReproduction,taillePopulation);
        choixInd1 = rnd(0,taillePopulation);
        do
        {
            /// ce traitement évite l'appariement d'un element avec lui même/
            ///choixInd2 = rnd(tauxReproduction,taillePopulation);
            choixInd2 = rnd(0,taillePopulation);
        }
        while(choixInd2 == choixInd1 );

        if(nbrCocyles != 1){choixLocus = rnd(0,nbrCocyles-1);} /// -1 pour s'arrêter à moins 2

        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            mpz_set((populationCGE2+i)->genotypeEcart[j], (populationCGE1+choixInd1)->genotypeEcart[j]);
            mpz_set((populationCGE2+i+1)->genotypeEcart[j] , (populationCGE1+choixInd2)->genotypeEcart[j]);
        }
        ///*****************************************************************
        for(j=choixLocus+1; j<nbrCocyles; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            mpz_set((populationCGE2+i)->genotypeEcart[j] , (populationCGE1+choixInd2)->genotypeEcart[j]);
            mpz_set((populationCGE2+i+1)->genotypeEcart[j] , (populationCGE1+choixInd1)->genotypeEcart[j]);
        }
        /// affectation de nouveau indice
        (populationCGE2+i)->id = i;
        (populationCGE2+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }

#endif
}
///**************************************************************************************
void mutationCGE(partitionCGE* populationCGE){

    int i,numeroGeneMute,nbeCocyles=nbrParties-1;
    double applicationOfMutation;

    static mpz_t borneSup; mpz_init(borneSup);
    mpz_ui_pow_ui(borneSup,2,nbrNoeuds-1);

    gmp_randstate_t gmpRandState; /* Random generator state object */
    gmp_randinit_default(gmpRandState);
    gmp_randseed_ui(gmpRandState, time(NULL));

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
            numeroGeneMute = rnd(0,nbeCocyles);
            mpz_urandomm((populationCGE+i)->genotypeEcart[numeroGeneMute], gmpRandState, borneSup);
        }
    }
}
///**************************************************************************************
float testerLaSommeDesFitnessCGE(partitionCGE* populationCGE){
    int i;
    float sommeFitness = 0;
    for(i=0; i<taillePopulation ; i++)
    {
        sommeFitness = sommeFitness + (populationCGE+i)->fitness;
    }
    return sommeFitness;
}
///**************************************************************************************
int findTheBestSolutionCGE(partitionCGE *populationCGE){
    float maxFitness = (populationCGE+0)->fitness;
    int i,indice=0;
    for(i=0; i<taillePopulation ; i++)
    {
        if(maxFitness < (populationCGE+i)->fitness)
        {
            maxFitness = (populationCGE+i)->fitness;
            indice = i;
        }
    }
    return indice;

}
///**************************************************************************************
void displayTheBestSolutionCGE(partitionCGE* solutionDominante){
    int i;
    printf("id = %d | normalised cut size = %d | cut size intrA = %d cut size intEr = %d ",
           solutionDominante->id,solutionDominante->coutCoupeNormalise, solutionDominante->coutCoupe, (sommeTotalFlux-solutionDominante->coutCoupe));
    printf("le contraintes violCGEs = %d \n", solutionDominante->contrainteViole);
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
void writeSolutionInFileCGE(partitionCGE *populationCGE, FILE *outputFilePop,int iteration){

    int i,j;
    for(i=0; i<taillePopulation; i++)
    {

        ///écriture de la meilleur solution dans le fichier output
        fprintf(outputFilePop,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,(populationCGE+i)->id,
                (populationCGE+i)->coutCoupe,(populationCGE+i)->fitness,(populationCGE+i)->coutCoupeNormalise,
                (populationCGE+i)->contrainteViole,(populationCGE+i)->nbrCluster);
        for(j=0; j<nbrArretes; j++)
        {
            fprintf(outputFilePop,"%d ",(populationCGE+i)->genotype[j]);
        }
        fprintf(outputFilePop,"\t\t\t\t");
        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(outputFilePop,"%d ",(populationCGE+i)->phenotype[j]);
        }
        fprintf(outputFilePop,"\t\t\t\t");
        for(j=0; j<nbrParties-1; j++)
        {
            gmp_fprintf(outputFilePop,"%Zd ",(populationCGE+i)->genotypeEcart[j]);
        }
        fprintf(outputFilePop,"\t\t\t\t");

        for(j=0; j<nbrParties-1; j++)
        {
            gmp_fprintf(outputFilePop,"%Zd ",(populationCGE+i)->genotypeEntier[j]);
        }
        fprintf(outputFilePop,"\t\t\t\t");
        for(j=0; j<(nbrParties-1)*(nbrNoeuds-1); j++)
        {
            fprintf(outputFilePop,"%d ",(populationCGE+i)->cocycles[j]);
        }
        fprintf(outputFilePop,"\n");

        for(j=0; j<nbr_constraint; j++)
        {
            fprintf(outputFilePop,"%d ",(populationCGE+i)->constraintVector[j]);
        }

    }
    fprintf(outputFilePop,"\n\n");
}
///**************************************************************************************
void writeBestSolutionInFileCGE(partitionCGE *solutionDominante, FILE *outputFile,int iteration){

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

///**************************************************************************************
void getPartitionFromSolutionCGE(partitionCGE *populationCGE){

    int i,j,k,nd,na,numeroClustre,nbrIntra; /// nd = noued de départ , na = noued d'arrivé

    for(i=0; i<taillePopulation; i++)
    {
        /// initialisation des phénotype = partition par les numero des noeuds correspondant

        for(j=0; j<nbrNoeuds; j++)
        {
            (populationCGE+i)->phenotype[j] = j;
        }
        /// 01/10/2017 : amélioration du temps d'exécution du codage binaire
        nbrIntra = 0;
        for(j=0; j<nbrArretes; j++){
            if((populationCGE+i)->genotype[j]==0)
            {
               edgeVectorIntra[nbrIntra] = edgeVector[j];
               nbrIntra++;
            }

        }
        /// détermination de la parition
        for(j=0; j<nbrIntra; j++){

                nd = (edgeVectorIntra+j)->nouedDepart;
                na = (edgeVectorIntra+j)->nouedArrive;
               if((populationCGE+i)->phenotype[na]!= na) /// il a été modifier auparavant
                {

                    if((populationCGE+i)->phenotype[na] < (populationCGE+i)->phenotype[nd])
                    {
                        ///dans le cas ou le noued d'arrivé est inférieur au noued depart ==> toutes les occurence du noued
                        /// de départ prendront ceux du noeud d'arrivé
                        numeroClustre = (populationCGE+i)->phenotype[nd];
                        for(k=0; k<nbrNoeuds; k++)
                        {
                            if((populationCGE+i)->phenotype[k] == numeroClustre)
                            {
                                (populationCGE+i)->phenotype[k] = (populationCGE+i)->phenotype[na];
                            }
                        }
                    }
                    else if((populationCGE+i)->phenotype[na] > (populationCGE+i)->phenotype[nd])
                    {
                        numeroClustre = (populationCGE+i)->phenotype[na];
                        for(k=0; k<nbrNoeuds; k++)
                        {
                            if((populationCGE+i)->phenotype[k] == numeroClustre)
                            {
                                (populationCGE+i)->phenotype[k] = (populationCGE+i)->phenotype[nd];
                            }
                        }
                    }
                }
                else
                {
                    (populationCGE+i)->phenotype[na] = (populationCGE+i)->phenotype[nd];
                }
            }
        }

/** 19/11/2017
    int i,j,k,nd,na,numeroClustre; /// nd = noued de départ , na = noued d'arrivé

    for(i=0; i<taillePopulation; i++)
    {
        /// initialisation des phénotype = partition par les numero des noeuds correspondant
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationCGE+i)->phenotype[j] = j;
        }

        /// détermination de la parition
        ///for(j=0; j<nbrArretes; j++)
        for(j=nbrArretes-1; j>=0; j--){
            if((populationCGE+i)->genotype[j]==0){
                na = edgeVector[j].nouedArrive;
                nd = edgeVector[j].nouedDepart;

                if((populationCGE+i)->phenotype[na]!= na) /// il a été modifier auparavant
                {

                    if((populationCGE+i)->phenotype[na] < (populationCGE+i)->phenotype[nd])
                    {
                        ///dans le cas ou le noued d'arrivé est inférieur au noued d'arrivé ==> toutes les occurence du noued
                        /// de départ prendront ceux du noeud d'arrivé
                        numeroClustre = (populationCGE+i)->phenotype[nd];
                        for(k=0; k<nbrNoeuds; k++)
                        {
                            if((populationCGE+i)->phenotype[k] == numeroClustre)
                            {
                                (populationCGE+i)->phenotype[k] = (populationCGE+i)->phenotype[na];
                            }
                        }
                    }
                    else if((populationCGE+i)->phenotype[na] > (populationCGE+i)->phenotype[nd])
                    {
                        numeroClustre = (populationCGE+i)->phenotype[na];
                        for(k=0; k<nbrNoeuds; k++)
                        {
                            if((populationCGE+i)->phenotype[k] == numeroClustre)
                            {
                                (populationCGE+i)->phenotype[k] = (populationCGE+i)->phenotype[nd];
                            }
                        }
                    }
                }
                else
                {
                    (populationCGE+i)->phenotype[na] = (populationCGE+i)->phenotype[nd];
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

void calculCoutCoupeEtFitnessCGE(partitionCGE* populationCGE){
    int i,j,k,sommeTotalCoutCoupe = 0;

    for(i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (populationCGE+i)->coutCoupe = 0; /// réinitialisation de la valriable après chaque itération (individus de la population)
        for(j=0; j<nbrNoeuds-1; j++)
        {
            for(k=j+1; k<nbrNoeuds; k++)
            {
                if((populationCGE+i)->phenotype[j] == (populationCGE+i)->phenotype[k] && fluxMatrix[j][k] > 0)
                {
                    (populationCGE+i)->coutCoupe = (populationCGE+i)->coutCoupe +  fluxMatrix[j][k];
                }
            }
        }
        (populationCGE+i)->coutCoupeNormalise = (populationCGE+i)->coutCoupe + ((nbr_constraint - (populationCGE+i)->contrainteViole)*sommeTotalFlux);
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationCGE+i)->coutCoupeNormalise;
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
            varianceCoutDeCoupeNormalise + pow(((populationCGE+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);
    /// tester la valeur ajouté
    for(i=0; i<taillePopulation; i++)
    {

        /// calcul de sigma truncation
        (populationCGE+i)->expectedValue = (populationCGE+i)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationCGE+i)->expectedValue;
    }

    /// calcule du fitness
    /// le 26/11/2015 : calcul du fitness des solution avec leur expectedValues
    for(i=0; i<taillePopulation ; i++)
    {
        (populationCGE+i)->fitness = (float)((populationCGE+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);
    }

#else
    FILE *file;
    file = fopen("D:/WorkingFolderAgProgramming/benchmark/ecartEncoding/expectedValueCGE.txt","w");

    for(i=0; i<taillePopulation ; i++)
    {
        (populationCGE+i)->fitness = (float)((populationCGE+i)->coutCoupeNormalise)/(float)(sommeTotalCoutCoupe);
    }
#endif// scaling

}

///**************************************************************************************
/// calculer le cout de coupe de chacune des solution dans la population
void calculCoutCoupCGEtFitnessWithFlowVectorCGE(partitionCGE* populationCGE){

    int i,j,sommeTotalCoutCoupe = 0;
    for(i=0; i<taillePopulation; i++){
        (populationCGE+i)->coutCoupe = 0;
        for(j=0; j<nbrArretes; j++) ///parcourir le vecteurSolution =  solution codée en binaire
        {
            if((populationCGE+i)->genotype[j]==0)  /// on cherche à maximiser les intraCluster donc les arrêtes 0 != 1
            {
                (populationCGE+i)->coutCoupe = (populationCGE+i)->coutCoupe + fluxVector[j];
            }
        }
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationCGE+i)->coutCoupe;
    }
    /// calcule du fitness
    for(i=0; i<taillePopulation ; i++)
    {
        ///(populationCGE+i)->fitness = (float)((populationCGE+i)->coutCoupe)/sommeTotalCoutCoupe;
        (populationCGE+i)->fitness = (float)((populationCGE+i)->coutCoupe)/(float)(sommeTotalCoutCoupe);
    }
}

///=====================================================================================
/// Cette partie concerne la génération des partie sans répétition, le choix de la méthode
/// de génération se fait au niveau de la fonction main
///=====================================================================================
///**************************************************************************************
void generatePopulationWithoutRedandancyCGE(partitionCGE* populationCGE, int indiceFirstElt ){

    genererLesEcartCGE(populationCGE, indiceFirstElt);
    genererSolutionEntiereCGE(populationCGE, indiceFirstElt);
    genererLesSolutionBinaireCGE(populationCGE , indiceFirstElt);

}
///**************************************************************************************

void genererLesEcartCGE(partitionCGE* populationCGE, int indiceFirstElt){

    int i,j,nbrCocycles = nbrParties-1;
    static mpz_t borneSup;
    mpz_init2(borneSup,4000);
    mpz_ui_pow_ui(borneSup,2,nbrNoeuds-1);

    gmp_randstate_t gmpRandState; /* Random genegenotypeEcartrator state object */
    gmp_randinit_default(gmpRandState);
    gmp_randseed_ui(gmpRandState, time(NULL));

    for(i=indiceFirstElt;i<taillePopulation;i++){
        mpz_init2((populationCGE+i)->genotypeEcart[0],4000);
        mpz_urandomm((populationCGE+i)->genotypeEcart[0], gmpRandState, borneSup);
        for(j=1;j<nbrCocycles;j++){
            mpz_init2((populationCGE+i)->genotypeEcart[j],4000);
            mpz_urandomm((populationCGE+i)->genotypeEcart[j], gmpRandState, borneSup);
        }
        (populationCGE+i)->id =i;
        (populationCGE+i)->coutCoupe =0;
        (populationCGE+i)->fitness =0;
    }
}

///**************************************************************************************
void genererSolutionEntiereCGE(partitionCGE* populationCGE,int indiceFirstElt){

    int i,j,k,nbrCocycles=nbrParties-1;
    mpz_t borneSup,genotypeEntier;
    mpz_init2(borneSup,4000);
    mpz_ui_pow_ui(borneSup,2,nbrNoeuds-1);
    ///printf("la taille de la variable borneSup = %d \n",mpz_sizeinbase (borneSup, 10));


    mpz_init2(genotypeEntier,4000);///mpz_set(genotypeEntier,borneSup);

    ///printf("la taille de la variable genotypeEntier = %d \n",mpz_sizeinbase (genotypeEntier, 10));

    for(i=indiceFirstElt;i<taillePopulation;i++){
        mpz_init2((populationCGE+i)->genotypeEcart[0],4000);
        mpz_set((populationCGE+i)->genotypeEntier[0],(populationCGE+i)->genotypeEcart[0]);
        for(j=1; j<nbrCocycles; j++){
            mpz_add(genotypeEntier,(populationCGE+i)->genotypeEcart[j-1],(populationCGE+i)->genotypeEcart[j]);
            if(mpz_cmp(genotypeEntier,borneSup) <= 0 ){
                 mpz_init2((populationCGE+i)->genotypeEntier[j],4000);
                 mpz_set((populationCGE+i)->genotypeEntier[j],genotypeEntier);
            }
            else {
                for(k=j;k<nbrCocycles;k++){
                        mpz_init2((populationCGE+i)->genotypeEcart[k],4000);
                        mpz_set_ui((populationCGE+i)->genotypeEcart[k],0);
                }
                break;
            }
        }
    }
}
///**************************************************************************************
void genererLesSolutionBinaireCGE(partitionCGE *populationCGE, int indiceFirstElt){

    int i,j,k,nbrCocyles = nbrParties-1;
    int cocyclesSize = nbrParties*nbrNoeuds;
    int bit;
    static mpz_t tmp;
    mpz_init(tmp);
    static mpz_t binaryTmp; mpz_init(binaryTmp);

    for(i=indiceFirstElt; i<taillePopulation; i++){
        for(j=0; j<nbrCocyles; j++){
            mpz_set(tmp,(populationCGE+i)->genotypeEcart[j]);
            k=nbrNoeuds-2;
            do{
                mpz_mod_ui(binaryTmp,tmp,2);
                (populationCGE+i)->cocycles[(j*(nbrNoeuds-1))+k] = mpz_get_ui(binaryTmp);
                mpz_fdiv_q_ui(tmp,tmp,2);
                k--;
            }
            while(k>=0);
        }
    }
}

///**************************************************************************************

///======================================================================================
/// the main programm
void ecartEncoding(int nbrGeneration,FILE *outputFileCGE,FILE *outputFilePopCGE,FILE *outputOptimalSolutionFileCGE,
                   partitionCGE *populationCGE1,partitionCGE *populationCGE2,partitionCGE *solutionDominante,int iteration , int *bestSolutionIteration , int *nbrApparition){

        naturalSelectionCGE(populationCGE1,populationCGE2);
        crossOverCGE(populationCGE2, populationCGE1);
        mutationCGE(populationCGE1);
        ///==================================================
        genererSolutionEntiereCGE(populationCGE1,0);

        genererLesSolutionBinaireCGE(populationCGE1,0);
        ///==================================================
        calculerGenotypeCGE(populationCGE1);
        getPartitionFromSolutionCGE(populationCGE1);
        checkContrainstAndFitnessPenalizationCGE(populationCGE1);
        calculCoutCoupeEtFitnessCGE(populationCGE1);


#if writingPopulationInFile
        writeSolutionInFileCGE(populationCGE1,outputFilePopCGE,iteration);
#endif // writingPopulationInFile

        bestSolution=findTheBestSolutionCGE(populationCGE1);
#if writingPopulationInFile
        writeBestSolutionInFileCGE((populationCGE1+bestSolution),outputFileCGE,iteration);
#endif // writingPopulationInFile
        if((populationCGE1+bestSolution)->coutCoupeNormalise > solutionDominante->coutCoupeNormalise)
        {
            ///*(solutionDominante) = *(populationCGE1+bestSolution);
            nbrApparition=1;
            bestSolutionIteration = iteration;
        }
        else
        {
            nbrApparition++;
        }

}


///**************************************************************************************
void checkContrainstAndFitnessPenalizationCGE(partitionCGE *populationCGE){

    int i,j,k;

    for(i=0; i<taillePopulation; i++)
    {
        (populationCGE+i)->nbrCluster=0; /// le nombre de cluster doit etre initialisé à 1 pour prendre en consédiration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes violées
        (populationCGE+i)->contrainteViole=0;
        for(j=0; j<4; j++)
        {
            (populationCGE+i)->constraintVector[j]=0;
        }
        ///**********************************************************************************


        ///trier le tableau résulant
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationCGE+i)->clustersSize[j] =0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if ((populationCGE+i)->phenotype[k] == j)
                {
                    (populationCGE+i)->clustersSize[j]++;
                }
            }
            if ((populationCGE+i)->clustersSize[j] !=0)
            {
                (populationCGE+i)->nbrCluster++;
            }
        }

        if((populationCGE+i)->nbrCluster > max_clusters || (populationCGE+i)->nbrCluster < min_clusters)
        {
            (populationCGE+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationCGE+i)->contrainteViole++;
        }

        for(j=0; j<nbrNoeuds; j++)
        {
            if((populationCGE+i)->clustersSize[j]!=0)
            {
                if((populationCGE+i)->clustersSize[j]>max_sizeCluster || (populationCGE+i)->clustersSize[j]< min_sizeCluster)
                {
                    (populationCGE+i)->constraintVector[1]++;

                }
            }

        }

        if((populationCGE+i)->constraintVector[1] != 0)
        {
            (populationCGE+i)->contrainteViole++;
        }
    }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0)
    {
        checkCohabitationAndNonCohabitationConstraintCGE(populationCGE);
    }

}

///**************************************************************************************
void checkCohabitationAndNonCohabitationConstraintCGE(partitionCGE *populationCGE){

    int i,j,noeud1,noeud2;
    for(i=0; i<taillePopulation; i++)
    {
        /// Cette tache est déjà réalisé au niveau de la fonction checkCohabitationAndNonCohabitationConstrainteBE
		/**
		(populationCGE+i)->constraintVector[2]=0;
        (populationCGE+i)->constraintVector[3]=0;
		*/
        if(nbrCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrCohabitationConstraintes; j++)
            {
                noeud1=cohabitationConstraintes[j].nouedDepart;
                noeud2=cohabitationConstraintes[j].nouedArrive;
                if((populationCGE+i)->phenotype[noeud1]!= (populationCGE+i)->phenotype[noeud2])
                {
                    (populationCGE+i)->constraintVector[2]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationCGE+i)->constraintVector[2]!=0)
            {
                (populationCGE+i)->contrainteViole++;
            }
        }


        if(nbrNonCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrNonCohabitationConstraintes; j++)
            {
                noeud1=nonCohabitationConstraintes[j].nouedDepart;
                noeud2=nonCohabitationConstraintes[j].nouedArrive;

                if((populationCGE+i)->phenotype[noeud1]== (populationCGE+i)->phenotype[noeud2])
                {
                    (populationCGE+i)->constraintVector[3]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationCGE+i)->constraintVector[3]!=0)
            {
                (populationCGE+i)->contrainteViole++;
            }
        }


    }

}


///********************************************************************************************
void writeOptimalSolutionInFileCGE(partitionCGE *solutionDominante,FILE* outputOptimalSolutionFilCGE,
                                  int nbrRun, int bestSolutionIteration, float runTime,int ES){

    fprintf(outputOptimalSolutionFilCGE,"RUN NUMBER = %d | ITERATION N = %d | SOLUTION IDENTIFIER = %d | CUT SIZE INTRA = %d | CUT SIZE INTER = %d | SOMME TOTAL DES FLUX = %d | RUN TIME = %f | ES = %d | NBR VIOLATED CONSTRAINTE = %d | ",
            nbrRun, bestSolutionIteration, solutionDominante->id, solutionDominante->coutCoupe,(sommeTotalFlux - solutionDominante->coutCoupe),sommeTotalFlux ,runTime, ES, solutionDominante->contrainteViole);
     fprintf(outputOptimalSolutionFilCGE," max_clusters = %d | max_sizeCluster = %d  |  cohabitation = %d | non cohabitation = %d\n",
            solutionDominante->constraintVector[0],solutionDominante->constraintVector[1],solutionDominante->constraintVector[2],
            solutionDominante->constraintVector[3]);


}
///******************************************************************************************
int compareCroissantFitnessCGE (void const *a, void const *b){
    partitionCGE const *pa = a;
    partitionCGE const *pb = b;
    return  pb->coutCoupeNormalise - pa->coutCoupeNormalise;
}

void sortingPopulationByFintness(partitionCGE* populationCGE2, partitionCGE* tmpPopulation){


    if(tmpPopulation != NULL){
        int i,j;
        ///for(i=0;i<taillePopulation-1;i++){
        for(i=0;i<tauxReproduction;i++){
            for(j=i+1;j<taillePopulation;j++){
                if((populationCGE2+i)->coutCoupeNormalise < (populationCGE2+j)->coutCoupeNormalise){
                        *tmpPopulation = *(populationCGE2+i);
                        *(populationCGE2+i)= *(populationCGE2+j);
                        *(populationCGE2+j) = *tmpPopulation;

                }
            }
            ///printf("*(populationCGE2+%d)->coutCoupeNormalise = %d\n",i,(populationCGE2+i)->coutCoupeNormalise);
            ///mon_slCGEp(pause);
        }
        ///printf("\n\n\n");
    }
    else{
        printf("tmpPopulation = %p \n",tmpPopulation);
    }
}

///**************************************************************************************
/***
void writeOptimalSolutionInFilCGE(partitionCGE *solutionDominante,FILE* outputOptimalSolutionFilCGE){
    int i;
    fprintf(outputOptimalSolutionFilCGE,"id = %d | normalised cut size = %d | cut size intrA = %d cut size intEr = %d ",
            solutionDominante->id,solutionDominante->coutCoupeNormalise, solutionDominante->coutCoupe, (sommeTotalFlux-solutionDominante->coutCoupe));
    fprintf(outputOptimalSolutionFilCGE,"le contraintes violCGEs = %d \n", solutionDominante->contrainteViole);


    for(i=0; i<4; i++)
    {
        switch (i)
        {

        case 0:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFilCGE,"contrainte sur le nombre de clusters : OK \n");
            else fprintf(outputOptimalSolutionFilCGE,"contrainte sur le nombre de clusters : NOT OK \n");
            break;
        case 1:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFilCGE,"contrainte sur la taille des clusters : OK \n");
            else fprintf(outputOptimalSolutionFilCGE,"contrainte sur la taille des clusters : NOT OK \n");
            break;
        case 2:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFilCGE,"contrainte de cohabitation : OK \n");
            else fprintf(outputOptimalSolutionFilCGE,"contrainte de cohabitation : NOT OK \n");
            break;
        case 3:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFilCGE,"contrainte Non Cohabitation : OK \n");
            else fprintf(outputOptimalSolutionFilCGE,"contrainte Non Cohabitation : NOT OK \n");
            break;
        default:
            fprintf(outputOptimalSolutionFilCGE,"Cette contrainte n'est pas prise en charge par le système \n");

        }

    }

}

void checkContrainstAndFitnessPenalizationCGE(partitionCGE *populationCGE)
{

    int i,j,k,m=nbrNoeuds,tabTmp[m],tmp;

    for(i=0; i<taillepopulationCGE; i++)
    {
        (populationCGE+i)->nbrCluster=1; /// le nombre de cluster doit etre initialisé à 1 pour prendre en consédiration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes violées
        (populationCGE+i)->contrainteViole=0;
        for(j=0;j<nbr_constraint;j++){(populationCGE+i)->constraintVector[j]=0;}


        m=nbrNoeuds;
        ///copier le phenotype dans un deuxième tableau
        for(j=0; j<nbrNoeuds; j++)
        {
            tabTmp[j]=(populationCGE+i)->phenotype[j];
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
                (populationCGE+i)->nbrCluster++;
            }
        }
        if((populationCGE+i)->nbrCluster > max_clusters || (populationCGE+i)->nbrCluster < min_clusters)
        {
            (populationCGE+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationCGE+i)->contrainteViole++;
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
            (populationCGE+i)->clustersSize[j]=0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if (tabTmp[j] == (populationCGE+i)->phenotype[k])
                {
                    (populationCGE+i)->clustersSize[j]=(populationCGE+i)->clustersSize[j]+1;
                }
            }
            if((populationCGE+i)->clustersSize[j]>max_sizeCluster || (populationCGE+i)->clustersSize[j]< min_sizeCluster)
            {
                ///24.07.2015 ---> l'ajout de vecteur de contraintes = le nombre de cluster dont la taille a dépassé la borne
                (populationCGE+i)->constraintVector[1]++;

            }
        }

        if((populationCGE+i)->constraintVector[1] != 0)
        {
            (populationCGE+i)->contrainteViole++;
        }
   }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0){
        checkCohabitationAndNonCohabitationConstraintCGE(populationCGE);
    }

}

///**************************************************************************************
void genererSolutionEntierCGE(partitionCGE* populationCGE){
    int i,j,k,nbrCocycles=nbrParties-1;
    ///ul borneSup = pow(2,nbrNoeuds-1)-1,genotypeEntier;
    static mpz_t borneSup,genotypeEntier;
    mpz_init2(borneSup,500);mpz_ui_pow_ui(borneSup,2,nbrNoeuds-1);
    mpz_init2(genotypeEntier,500);///mpz_set(genotypeEntier,borneSup);

    gmp_printf ("%s is an mpz %Zd\n", "borneSup", borneSup);
    ///gmp_printf ("%s is an mpz %Zd\n", "genotypeEntier", genotypeEntier);
    for(i=0;i<taillePopulation;i++){
        mpz_init2((populationCGE+i)->genotypeEcart[i],500);
    }

    for(i=0;i<taillePopulation;i++){
        printf("nbrRun = %d  individu = %d\n",nbrRun,i);
        gmp_printf ("%s is an mpz %Zd\n", "(populationCGE+i)->genotypeEcart[i] ", (populationCGE+i)->genotypeEcart[i]);
        ///mpz_init2((populationCGE+i)->genotypeEcart[i],500);
        mpz_set((populationCGE+i)->genotypeEcart[0],(populationCGE+i)->genotypeEcart[0]);
        ///gmp_printf ("%s is an mpz %Zd\n", "(populationCGE+i)->genotypeEcart[i] ", (populationCGE+i)->genotypeEcart[i]);
        for(j=1; j<nbrCocycles; j++){
            mpz_add(genotypeEntier,(populationCGE+i)->genotypeEcart[j-1],(populationCGE+i)->genotypeEcart[j]);
            gmp_printf("%s is an mpz %Zd\n", "genotypeEntier ", genotypeEntier);
            if(mpz_cmp(genotypeEntier,borneSup) < 0 )
                    mpz_set((populationCGE+i)->genotypeEcart[j],genotypeEntier);
            gmp_printf ("%s is an mpz %Zd\n", "(populationCGE+i)->genotypeEcart[i] ", (populationCGE+i)->genotypeEcart[i]);


            else {
                    for(k=j;k<nbrCocycles;k++){
                        ///mpz_init2((populationCGE+i)->genotypeEcart[k],500);
                        mpz_set((populationCGE+i)->genotypeEcart[k] ,0);
                    }
            break;

            }
        }
    }
}
*/
