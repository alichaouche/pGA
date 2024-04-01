
/**================= Sorted Binary Group Number Encoding ====================================================
* 04/12/2015 :
* la population initiale est maintenant prête
*============================================================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "typeDeclaration.h"
#include "./communesFunctions.h"
#include "./sortedBinaryGroupNumberEncoding.h"
#include "./compilationConditionnelle.h"

extern int taillePopulation,tauxReproduction,sommeTotalFlux,nbrNoeuds,
       nbr_constraint,nbrArretes,nbrParties,nbrNoeudsInCluster, regeneration;
extern double tauxMutation ,averageTraffic;
extern int fluxMatrix[1000][1000],*fluxVector, maxFlow, sommeTotalFluxReal;
extern edge *edgeVector,cohabitationConstraintes[1000],nonCohabitationConstraintes[1000];
extern int max_clusters , min_clusters,min_sizeCluster ,max_sizeCluster ;
extern int bestSolution, nbrApparition, nbrCohabitationConstraintes,nbrNonCohabitationConstraintes;
extern int nbrRun;





///======================================================================================
/// the main programm
void sortedBinaryGroupNumberEncoding(int nbrGeneration, FILE* outputFileVAE, FILE* outputFilePopVAE,FILE *outputOptimalSolutionFileVAE,
                                     partitionVAE *populationVAE1,partitionVAE *populationVAE2, partitionVAE *solutionDominante,int iteration , int *bestSolutionIteration , int *nbrApparition )
{

        naturalSelectionVAE(populationVAE1,populationVAE2);
        crossOverVAE(populationVAE2, populationVAE1);
        mutationVAE(populationVAE1);
        getPartitionFromSolutionVAE(populationVAE1);
        checkContrainstAndFitnessPenalizationVAE(populationVAE1);
        calculCoutCoupeEtFitnessVAE(populationVAE1);

#if writingPopulationInFile
        writeSolutionInFileVAE(populationVAE1,outputFilePopVAE,iteration);
#endif // writingPopulationInFile

        bestSolution=findTheBestSolutionVAE(populationVAE1);
#if writingPopulationInFile
        writeBestSolutionInFileVAE((populationVAE1+bestSolution),outputFileVAE,iteration);
#endif // writingPopulationInFile
        if((populationVAE1+bestSolution)->coutCoupeNormalise > solutionDominante->coutCoupeNormalise)
        {
            *(solutionDominante) = *(populationVAE1+bestSolution);
            nbrApparition=1;
            bestSolutionIteration =iteration;

        }
        else
        {
            nbrApparition++;
        }
}
///**************************************************************************************
void genererPopulationInitialeVAE(partitionVAE *populationVAE, int indiceFirstElt)
{

    int i,nbrSolutionRepeter;

    for(i=indiceFirstElt; i<taillePopulation; i++)
    {
        (populationVAE+i)->id =i;
        generateBinarySolutionVAE(populationVAE, i);
    }
}
/**
void genererPopulationInitialeVAE(partitionVAE *populationVAE)
{

    int i,nbrSolutionRepeter;
    ///printf("the genererPopulationInitialeVAE function ...\n");
    for(i=0; i<taillePopulation; i++)
    {
        (populationVAE+i)->id =i;
        ///printf("(populationVAE+i)->id = %d",(populationVAE+i)->id);
        /// affecter la valeur 0 à tous les noeuds de la solution

        nbrSolutionRepeter =0;

        do
        {
            generateBinarySolutionVAE(populationVAE, i);
            conversionSolutionBinaireToEntiereVAE(populationVAE,i);
            trieGenotypeEntieVAE(populationVAE, i);
            nbrSolutionRepeter++;
        }
        while(existanceDeSolutionVAE(populationVAE, i));
        ///printf("création du %d individus après %d tentatives...!\n",i,nbrSolutionRepeter);

    }
}*/
///**************************************************************************************
void genererPopulationInitialeRandomlyVAE(partitionVAE *populationVAE, int indiceFirstElt)
{

    int i;
    ///printf("the genererPopulationInitialeVAE function ...\n");
    for(i=indiceFirstElt; i<taillePopulation; i++)
    {
        (populationVAE+i)->id =i;
        generateBinarySolutionVAE(populationVAE, i);
    }
}
/// la génération des solution binaires
void generateBinarySolutionVAE(partitionVAE* populationVAE, int indiceIndividu)
{

    int numeroPartie,i,indiceNoeudAffecte;
    ///(populationVAE+indiceIndividu)->genotype = (char*)malloc((max_clusters * nbrNoeuds)*sizeof(char));
    (populationVAE+indiceIndividu)->genotype = (char*)malloc((max_clusters * nbrNoeuds)*sizeof(char));
    if((populationVAE+indiceIndividu)->genotype ==NULL){printf(stderr,"Memory allocation for genotype VAE failled\n");exit(EXIT_FAILURE);}
    ///for(i=0; i<(max_clusters * nbrNoeuds); i++)
    for(i=0; i<(max_clusters * nbrNoeuds); i++)
    {
        (populationVAE+indiceIndividu)->genotype[i] = 0;
    }

    for(i=0; i<nbrNoeuds; i++)
    {
        ///numeroPartie=rnd(0,max_clusters);
        numeroPartie=rnd(0,max_clusters);
        indiceNoeudAffecte = nbrNoeuds*numeroPartie + i;
        (populationVAE+indiceIndividu)->genotype[indiceNoeudAffecte] = 1;
    }

}
///********************************************************************************************
/**
    après la génération de la solution binaire, il va falloir procéder à son homogination par le fait
    que deux noeuds n'ont pas le droit d'appartenir à deux cluster au même temps, voici la procédure :
    - si un noeud n'est guère affecté ==> il sera affecté au dernier cluster
    - si un noued est affecté à plusieurs cluster ==> il sera affecter au cluster du plus petit rang

    Important : Cette technique sera utilisée lors de croiement
*/

void homoginiserLesSolutionBinaire(partitionVAE* populationVAE, int indiceIndividu)
{

    int i,j,k,choixCluster,noeudAffecte;

    for(k=0; k<nbrNoeuds; k++)
    {
        noeudAffecte = 0;
        ///for(j=0; j<max_clusters; j++)
        for(j=0; j<max_clusters; j++)
        {

            if((populationVAE+indiceIndividu)->genotype[(j*nbrNoeuds) + k]  == 1)
            {
                ///for(i=j+1; i<max_clusters; i++)
                for(i=j+1; i<max_clusters; i++)
                {
                    (populationVAE+indiceIndividu)->genotype[(i*nbrNoeuds) + k] =0;
                }
                noeudAffecte = 1;
                break;
            }
        }
        if(noeudAffecte == 0)
        {
            ///choixCluster = rnd(0,max_clusters);
            choixCluster = rnd(0,max_clusters);
            (populationVAE+indiceIndividu)->genotype[(choixCluster*nbrNoeuds) + k] =1;
        }
    }
}



/// conversion des solution binaires générées au solutions entiere pour illiminé la rendondance
///********************************************************************************************
void conversionSolutionBinaireToEntiereVAE(partitionVAE* populationVAE,int indiceIndividu)
{

    int i,j;
    ///for(i=0; i<max_clusters; i++)
    for(i=0; i<max_clusters; i++)
    {
        (populationVAE+indiceIndividu)->genotypeEntier[i] = 0;
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationVAE+indiceIndividu)->genotypeEntier[i] =
                (populationVAE+indiceIndividu)->genotypeEntier[i] +
                (pow(2,(nbrNoeuds-1-j))*(populationVAE+indiceIndividu)->genotype[(i*nbrNoeuds)+j]);
        }
    }
}

///**************************************************************************************
void trieGenotypeEntieVAE(partitionVAE* populationVAE, int indiceIndividu)
{

    int i,j,tmp;
    ///for(i=0; i<max_clusters-1; i++)
    for(i=0; i<max_clusters-1; i++)
    {
        ///for(j=i+1; j<max_clusters; j++)
        for(j=i+1; j<max_clusters; j++)
        {
            if((populationVAE+indiceIndividu)->genotypeEntier[i] > (populationVAE+indiceIndividu)->genotypeEntier[j])
            {
                tmp = (populationVAE+indiceIndividu)->genotypeEntier[i];
                (populationVAE+indiceIndividu)->genotypeEntier[i] = (populationVAE+indiceIndividu)->genotypeEntier[j];
                (populationVAE+indiceIndividu)->genotypeEntier[j] = tmp;
            }
        }
    }
}

///**************************************************************************************
int existanceDeSolutionVAE(partitionVAE* populationVAE, int indexNewSolution)
{

    int i=0,j,Exist = 0,theSameElement;

    while(i<indexNewSolution && Exist == 0)
    {
        theSameElement = 0;
        ///for(j=0; j<max_clusters; j++)
        for(j=0; j<max_clusters; j++)
        {
            if((populationVAE+i)->genotypeEntier[j] == (populationVAE+indexNewSolution)->genotypeEntier[j])
            {
                theSameElement++;
            }
        }
        ///if(theSameElement==max_clusters)
        if(theSameElement==max_clusters)
        {
            Exist = 1;
            break;
        }
        i++;
    }
    return Exist;
}


///**************************************************************************************
void getPartitionFromSolutionVAE(partitionVAE *populationVAE)
{

    int i,j,k;

    for(i=0; i<taillePopulation; i++)
    {

        /// détermination de la parition
        ///for(j=0; j<max_clusters; j++)
        for(j=0; j<max_clusters; j++)
        {
            for(k=0; k<nbrNoeuds; k++)
            {
                if((populationVAE+i)->genotype[(j*nbrNoeuds)+k] == 1)
                {
                    (populationVAE+i)->phenotype[k] = j;
                }
            }
        }
    }
}


///***************************************************************************************
/// calculer le cout de coupe de chaque partition ainsi que son fitness
///***************************************************************************************
/// on va appliquer la formule suivante pour le calcule des cout de coupe
/// T'(C) = B - T(C)            <=> la minimisation devient une maximisation
/// Z(C)=T'(C)+(u-v(C)) ×B      <=>  calculer la nouvelle fitness par rapport aux contrainte
/// T'(C) : cout de coupe (maximiser les flux INTRA), u ; nombre de contrainte totale
/// v(C): nombre de contrainte violé,  B: la somme totale des flux
///***************************************************************************************

void calculCoutCoupeEtFitnessVAE(partitionVAE* populationVAE)
{
    int i,j,k,sommeTotalCoutCoupe = 0;

    for(i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (populationVAE+i)->coutCoupe = 0; /// réinitialisation de la valriable après chaque itération (individus de la population)
        for(j=0; j<nbrNoeuds-1; j++)
        {
            for(k=j+1; k<nbrNoeuds; k++)
            {
                if((populationVAE+i)->phenotype[j] == (populationVAE+i)->phenotype[k] && fluxMatrix[j][k] > 0)
                {
                    (populationVAE+i)->coutCoupe = (populationVAE+i)->coutCoupe +  fluxMatrix[j][k];
                }
            }
        }
        (populationVAE+i)->coutCoupeNormalise = (populationVAE+i)->coutCoupe + ((nbr_constraint - (populationVAE+i)->contrainteViole)*sommeTotalFlux);
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationVAE+i)->coutCoupeNormalise;
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
            varianceCoutDeCoupeNormalise + pow(((populationVAE+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);
    /// tester la valeur ajouté
    for(i=0; i<taillePopulation; i++)
    {

        /// calcul de sigma truncation
        (populationVAE+i)->expectedValue = (populationVAE+i)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationVAE+i)->expectedValue;
    }

    for(i=0; i<taillePopulation ; i++)
    {
        (populationVAE+i)->fitness = (float)((populationVAE+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);
    }
#else
    FILE *file;
    file = fopen("D:/WorkingFolderAgProgramming/benchmark/sortedBinaryGroupNumberEncoding/expectedValueVAE.txt","w");

    for(i=0; i<taillePopulation ; i++)
    {
        (populationVAE+i)->fitness = (float)((populationVAE+i)->coutCoupeNormalise)/(float)(sommeTotalCoutCoupe);
    }
#endif// scaling



}
///**************************************************************************************
void checkContrainstAndFitnessPenalizationVAE(partitionVAE *populationVAE)
{

    int i,j,k;

    for(i=0; i<taillePopulation; i++)
    {
        (populationVAE+i)->nbrCluster=0; /// le nombre de cluster doit etre initialisé à 1 pour prendre en consédiration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes violées
        (populationVAE+i)->contrainteViole=0;
        for(j=0; j<4; j++)
        {
            (populationVAE+i)->constraintVector[j]=0;
        }
        ///**********************************************************************************


        ///trier le tableau résulant
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationVAE+i)->clustersSize[j] =0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if ((populationVAE+i)->phenotype[k] == j)
                {
                    (populationVAE+i)->clustersSize[j]++;
                }
            }
            if ((populationVAE+i)->clustersSize[j] !=0)
            {
                (populationVAE+i)->nbrCluster++;
            }
        }

        if((populationVAE+i)->nbrCluster > max_clusters || (populationVAE+i)->nbrCluster < min_clusters)
        {
            (populationVAE+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationVAE+i)->contrainteViole++;
        }

        for(j=0; j<nbrNoeuds; j++)
        {
            if((populationVAE+i)->clustersSize[j]!=0)
            {
                if((populationVAE+i)->clustersSize[j]>max_sizeCluster || (populationVAE+i)->clustersSize[j]< min_sizeCluster)
                {
                    (populationVAE+i)->constraintVector[1]++;

                }
            }

        }

        if((populationVAE+i)->constraintVector[1] != 0)
        {
            (populationVAE+i)->contrainteViole++;
        }
    }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0)
    {
        checkCohabitationAndNonCohabitationConstrainteVAE(populationVAE);
    }

}

///**************************************************************************************
void checkCohabitationAndNonCohabitationConstrainteVAE(partitionVAE *populationVAE)
{

    int i,j,noeud1,noeud2;
    for(i=0; i<taillePopulation; i++)
    {
	    /// Cette tache est déjà réalisé au niveau de la fonction checkCohabitationAndNonCohabitationConstrainteBE
		/**
        (populationVAE+i)->constraintVector[2]=0;
        (populationVAE+i)->constraintVector[3]=0;
		*/
        if(nbrCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrCohabitationConstraintes; j++)
            {
                noeud1=cohabitationConstraintes[j].nouedDepart;
                noeud2=cohabitationConstraintes[j].nouedArrive;

                if((populationVAE+i)->phenotype[noeud1]!= (populationVAE+i)->phenotype[noeud2])
                {
                    (populationVAE+i)->constraintVector[2]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationVAE+i)->constraintVector[2]!=0)
            {
                (populationVAE+i)->contrainteViole++;
            }
        }


        if(nbrNonCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrNonCohabitationConstraintes; j++)
            {
                noeud1=nonCohabitationConstraintes[j].nouedDepart;
                noeud2=nonCohabitationConstraintes[j].nouedArrive;

                if((populationVAE+i)->phenotype[noeud1]== (populationVAE+i)->phenotype[noeud2])
                {
                    (populationVAE+i)->constraintVector[3]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationVAE+i)->constraintVector[3]!=0)
            {
                (populationVAE+i)->contrainteViole++;
            }
        }


    }

}

///**************************************************************************************
/// la séléction naturelle des individus pour la nouvelle génération
void naturalSelectionVAE(partitionVAE* populationVAE1,partitionVAE* populationVAE2)
{
   int i,j,maxFitness;
    int tmpVectIndice[taillePopulation];
    for(i=0;i<taillePopulation;i++) tmpVectIndice[i]=-1;
    float lotterie,sommeFitness;

    for(i=0; i<tauxReproduction; i++){
        maxFitness = 0;
        for(j=0;j<taillePopulation;j++){
            if((populationVAE1+maxFitness)->coutCoupeNormalise < (populationVAE1+j)->coutCoupeNormalise
               && tmpVectIndice[j] == -1 ){
                maxFitness = j;

            }
        }
        tmpVectIndice[maxFitness] = 0;
        *(populationVAE2+i) = *(populationVAE1+maxFitness);
        ///printf("(populationVAE2+%d)->coutCoupeNormalise = %d \n",i,(populationVAE2+i)->coutCoupeNormalise);
        ///mon_sleep(pause);
    }
    ///printf("\n\n\n");

    ///qsort (populationVAE2, taillePopulation, sizeof *populationVAE2, compareCroissantFitnessEE);
    ///sortingPopulationByFintness(populationVAE2,tmpPopulation);

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
            sommeFitness = sommeFitness + (populationVAE1+j)->fitness;
            if (lotterie<=sommeFitness)
            {
                *(populationVAE2+i) = *(populationVAE1+j);
                break;
            }
            else if( j==taillePopulation-1)
            {
                *(populationVAE2+i) = *(populationVAE1+j);
            }
        }
    }
}
///**************************************************************************************
///le croisement des solution séléctionner
void crossOverVAE(partitionVAE* populationVAE1, partitionVAE* populationVAE2)
{

/// le croisement se fait au niveau des points d'articulation du génotype

    int i = 0,j,k;
    int choixInd1, choixInd2, choixLocus;
    ///***********************************************************************************
    for(i=0; i<tauxReproduction; i++)
    {
        *(populationVAE2+i) = *(populationVAE1+i);
        (populationVAE2+i)->id = i;
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
#if choixNbrPartiePourCrossoverVAE
        ///choixLocus = rnd(0,max_clusters-1); /// -1 pour s'arrêter à moins 2
        choixLocus = rnd(0,max_clusters-1); /// -1 pour s'arrêter à moins 2
        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            for(k=0; k<nbrNoeuds; k++)
            {
                (populationVAE2+i)->genotype[(j*nbrNoeuds)+k] = (populationVAE1+choixInd1)->genotype[(j*nbrNoeuds)+k];
                (populationVAE2+i+1)->genotype[(j*nbrNoeuds)+k] = (populationVAE1+choixInd2)->genotype[(j*nbrNoeuds)+k];
            }
        }

        ///************************************************************************************

        ///for(j=choixLocus+1; j<max_clusters; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        for(j=choixLocus+1; j<max_clusters; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            for(k=0; k<nbrNoeuds; k++)
            {
                (populationVAE2+i)->genotype[(j*nbrNoeuds)+k] = (populationVAE1+choixInd2)->genotype[(j*nbrNoeuds)+k];
                (populationVAE2+i+1)->genotype[(j*nbrNoeuds)+k] = (populationVAE1+choixInd1)->genotype[(j*nbrNoeuds)+k];
            }
        }
#else
        ///choixLocus = rnd(0,nbrNoeuds*max_clusters-1); /// -1 pour s'arrêter à moins 2
        choixLocus = rnd(0,nbrNoeuds*max_clusters-1); /// -1 pour s'arrêter à moins 2

        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationVAE2+i)->genotype[j] = (populationVAE1+choixInd1)->genotype[j];
            (populationVAE2+i+1)->genotype[j] = (populationVAE1+choixInd2)->genotype[j];
        }

        ///************************************************************************************

        ///for(j=choixLocus+1; j<max_clusters; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        for(j=choixLocus+1; j<max_clusters; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationVAE2+i)->genotype[j] = (populationVAE1+choixInd2)->genotype[j];
            (populationVAE2+i+1)->genotype[j] = (populationVAE1+choixInd1)->genotype[j];
        }


#endif

        /// affectation de nouveau indice
        (populationVAE2+i)->id = i;
        (populationVAE2+i+1)->id = i+1;

        homoginiserLesSolutionBinaire(populationVAE2, i);
        homoginiserLesSolutionBinaire(populationVAE2, i+1);
        ///*****************************************************************
        i+=2;
    }
    generateBinarySolutionVAE(populationVAE2, taillePopulation-regeneration);
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
        #if choixNbrPartiePourCrossoverVAE
                ///choixLocus = rnd(0,max_clusters-1); /// -1 pour s'arrêter à moins 2
                choixLocus = rnd(0,max_clusters-1); /// -1 pour s'arrêter à moins 2
                for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
                {
                    for(k=0; k<nbrNoeuds; k++)
                    {
                        (populationVAE2+i)->genotype[(j*nbrNoeuds)+k] = (populationVAE1+choixInd1)->genotype[(j*nbrNoeuds)+k];
                        (populationVAE2+i+1)->genotype[(j*nbrNoeuds)+k] = (populationVAE1+choixInd2)->genotype[(j*nbrNoeuds)+k];
                    }
                }

                ///************************************************************************************

                ///for(j=choixLocus+1; j<max_clusters; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
                for(j=choixLocus+1; j<max_clusters; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
                {
                    for(k=0; k<nbrNoeuds; k++)
                    {
                        (populationVAE2+i)->genotype[(j*nbrNoeuds)+k] = (populationVAE1+choixInd2)->genotype[(j*nbrNoeuds)+k];
                        (populationVAE2+i+1)->genotype[(j*nbrNoeuds)+k] = (populationVAE1+choixInd1)->genotype[(j*nbrNoeuds)+k];
                    }
                }
        #else
                ///choixLocus = rnd(0,nbrNoeuds*max_clusters-1); /// -1 pour s'arrêter à moins 2
                choixLocus = rnd(0,nbrNoeuds*max_clusters-1); /// -1 pour s'arrêter à moins 2

                for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
                {
                    (populationVAE2+i)->genotype[j] = (populationVAE1+choixInd1)->genotype[j];
                    (populationVAE2+i+1)->genotype[j] = (populationVAE1+choixInd2)->genotype[j];
                }

                ///************************************************************************************

                ///for(j=choixLocus+1; j<max_clusters; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
                for(j=choixLocus+1; j<max_clusters; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
                {
                    (populationVAE2+i)->genotype[j] = (populationVAE1+choixInd2)->genotype[j];
                    (populationVAE2+i+1)->genotype[j] = (populationVAE1+choixInd1)->genotype[j];
                }


        #endif

        /// affectation de nouveau indice
        (populationVAE2+i)->id = i;
        (populationVAE2+i+1)->id = i+1;

        homoginiserLesSolutionBinaire(populationVAE2, i);
        homoginiserLesSolutionBinaire(populationVAE2, i+1);
        ///*****************************************************************
        i+=2;
    }

#endif
}
///**************************************************************************************
int findTheBestSolutionVAE(partitionVAE *populationVAE)
{
    float maxFitness = (populationVAE+0)->fitness;
    int i,indice=0;
    for(i=0; i<taillePopulation ; i++)
    {
        if(maxFitness < (populationVAE+i)->fitness)
        {
            maxFitness = (populationVAE+i)->fitness;
            indice = i;
        }
    }
    return indice;

}
///**************************************************************************************
void mutationVAE(partitionVAE* populationVAE)
{

    int i,j,numeroGeneMute, choixCluster, choixNoeud;
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
            ///choixCluster = rnd(0,max_clusters);
            choixCluster = rnd(0,max_clusters);
            choixNoeud = rnd(0,nbrNoeuds);
            numeroGeneMute = (choixCluster*nbrNoeuds)+choixNoeud;
            if((populationVAE+i)->genotype[numeroGeneMute] == 1)
            {
                (populationVAE+i)->genotype[numeroGeneMute] =0;
            }
            else
            {
                (populationVAE+i)->genotype[numeroGeneMute] =1;
            }
            ///après l'application de la mutation, il va falloire remettre le genotype en ordre
            ///for(j=0; j<max_clusters; j++)
            for(j=0; j<max_clusters; j++)
            {
                if(j*nbrNoeuds+choixNoeud != numeroGeneMute && (populationVAE+i)->genotype[numeroGeneMute] == 1)
                {
                    (populationVAE+i)->genotype[j*nbrNoeuds+choixNoeud]=0;
                }
                else if(j*nbrNoeuds+choixNoeud != numeroGeneMute && (populationVAE+i)->genotype[numeroGeneMute] == 0)
                {
                    (populationVAE+i)->genotype[j*nbrNoeuds+choixNoeud]=1;
                }

            }

        }
    }
}
///**************************************************************************************
float testerLaSommeDesFitnessVAE(partitionVAE* populationVAE)
{
    int i;
    float sommeFitness = 0;
    for(i=0; i<taillePopulation ; i++)
    {
        sommeFitness = sommeFitness + (populationVAE+i)->fitness;
    }
    return sommeFitness;
}
///**************************************************************************************
void displayTheBestSolutionVAE(partitionVAE* solutionDominante)
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
void writeSolutionInFileVAE(partitionVAE *populationVAE, FILE *outputFilePop,int iteration)
{

    int i,j,k;
    for(i=0; i<taillePopulation; i++)
    {

        ///écriture de la meilleur solution dans le fichier output
        fprintf(outputFilePop,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,(populationVAE+i)->id,
                (populationVAE+i)->coutCoupe,(populationVAE+i)->fitness,(populationVAE+i)->coutCoupeNormalise,
                (populationVAE+i)->contrainteViole,(populationVAE+i)->nbrCluster);
        ///for(j=0; j<max_clusters; j++)
        for(j=0; j<max_clusters; j++)
        {
            for(k=0; k<nbrNoeuds; k++)
            {
                fprintf(outputFilePop,"%d ",(populationVAE+i)->genotype[j*nbrNoeuds+k]);
            }
        }
        fprintf(outputFilePop,"\t\t\t\t");
        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(outputFilePop,"%d ",(populationVAE+i)->phenotype[j]);
        }
        fprintf(outputFilePop,"\n");
    }
    fprintf(outputFilePop,"\n\n");
}
///**************************************************************************************
void writeBestSolutionInFileVAE(partitionVAE *solutionDominante, FILE *outputFile,int iteration)
{

    int i;
    fprintf(outputFile,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,solutionDominante->id,
            solutionDominante->coutCoupe,solutionDominante->fitness,solutionDominante->coutCoupeNormalise,
            solutionDominante->contrainteViole,solutionDominante->nbrCluster);
    ///for(i=0; i<max_clusters*nbrNoeuds; i++)
    for(i=0; i<max_clusters*nbrNoeuds; i++)
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
void affichePopulationVAE(partitionVAE* populationVAE)
{
    printf("affichePopulation...\n");
    int i,j,k,l;
    for(i=0; i<taillePopulation; i++)
    {

        printf("\nid = %d | cout de coupe = %d | fitness = %0.4f",
               (populationVAE+i)->id,(populationVAE+i)->coutCoupe, (populationVAE+i)->fitness);

        /// affichage de vecteur de la partition
        /**       printf("\nla partition est  = \n");
               for(l=0;l<max_clusters;l++){
                   printf("%d ",(populationVAE+i)->genotypeEntier[l]);
               }
        */

        printf("\nla solution generer est  = \n");
        ///for(j=0; j<max_clusters; j++)
        for(j=0; j<max_clusters; j++)
        {
            for(k=0; k<nbrNoeuds; k++)
            {
                ///printf("%d ",(populationVAE+i)->genotype[k+(j*max_clusters)]);
                printf("%d ",(populationVAE+i)->genotype[k+(j*max_clusters)]);
            }
            printf("\n");
        }
        printf("\n");
        /// affichage de vecteur de la partition
        printf("\nla partition est  = \n");
        for(l=0; l<nbrNoeuds; l++)
        {
            printf("%d ",(populationVAE+i)->phenotype[l]);
        }
        printf("\n===============================================================================\n");

    }
}

///********************************************************************************************
void writeOptimalSolutionInFileVAE(partitionVAE *solutionDominante,FILE* outputOptimalSolutionFileVAE,
                                     int nbrRun, int bestSolutionIteration, float runTime, int ES)
{
    int i;
    fprintf(outputOptimalSolutionFileVAE,"RUN NUMBER = %d | ITERATION N = %d | SOLUTION IDENTIFIER = %d | CUT SIZE INTRA = %d | CUT SIZE INTER = %d | SOMME TOTAL DES FLUX = %d | RUN TIME = %f | ES = %d | NBR VIOLATED CONSTRAINTE = %d | ",
            nbrRun, bestSolutionIteration, solutionDominante->id, solutionDominante->coutCoupe,(sommeTotalFlux - solutionDominante->coutCoupe),sommeTotalFlux ,runTime, ES, solutionDominante->contrainteViole);
 fprintf(outputOptimalSolutionFileVAE," max_clusters = %d | max_sizeCluster = %d  |  cohabitation = %d | non cohabitation = %d\n",
            solutionDominante->constraintVector[0],solutionDominante->constraintVector[1],solutionDominante->constraintVector[2],
            solutionDominante->constraintVector[3]);


}
///******************************************************************************************
int compareCroissantFitnessVAE (void const *a, void const *b)
{

    partitionVAE const *pa = a;
    partitionVAE const *pb = b;
    return  pb->coutCoupeNormalise - pa->coutCoupeNormalise;
}
/**************************************************************************************

///**************************************************************************************
void writeOptimalSolutionInFileVAE(partitionVAE *solutionDominante,FILE* outputOptimalSolutionFileVAE)
{
    int i;
    fprintf(outputOptimalSolutionFileVAE,"id = %d | normalised cut size = %d | cut size intrA = %d cut size intEr = %d ",
            solutionDominante->id,solutionDominante->coutCoupeNormalise, solutionDominante->coutCoupe, (sommeTotalFlux-solutionDominante->coutCoupe));
    fprintf(outputOptimalSolutionFileVAE,"le contraintes violees = %d \n", solutionDominante->contrainteViole);
    for(i=0; i<4; i++)
    {
        switch (i)
        {

        case 0:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileVAE,"contrainte sur le nombre de clusters : OK \n");
            else fprintf(outputOptimalSolutionFileVAE,"contrainte sur le nombre de clusters : NOT OK \n");
            break;
        case 1:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileVAE,"contrainte sur la taille des clusters : OK \n");
            else fprintf(outputOptimalSolutionFileVAE,"contrainte sur la taille des clusters : NOT OK \n");
            break;
        case 2:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileVAE,"contrainte de cohabitation : OK \n");
            else fprintf(outputOptimalSolutionFileVAE,"contrainte de cohabitation : NOT OK \n");
            break;
        case 3:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileVAE,"contrainte Non Cohabitation : OK \n");
            else fprintf(outputOptimalSolutionFileVAE,"contrainte Non Cohabitation : NOT OK \n");
            break;
        default:
            fprintf(outputOptimalSolutionFileVAE,"Cette contrainte n'est pas prise en charge par le système \n");

        }

    }

}

void checkContrainstAndFitnessPenalizationVAE(partitionVAE *populationVAE)
{

    int i,j,k,m=nbrNoeuds,tabTmp[m],tmp;

    for(i=0; i<taillePopulation; i++)
    {
        (populationVAE+i)->nbrCluster=1; /// le nombre de cluster doit etre initialisé à 1 pour prendre en consédiration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes violées
        (populationVAE+i)->contrainteViole=0;
        for(j=0;j<nbr_constraint;j++){(populationVAE+i)->constraintVector[j]=0;}


        m=nbrNoeuds;
        ///copier le phenotype dans un deuxième tableau
        for(j=0; j<nbrNoeuds; j++)
        {
            tabTmp[j]=(populationVAE+i)->phenotype[j];
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
                (populationVAE+i)->nbrCluster++;
            }
        }
        if((populationVAE+i)->nbrCluster > max_clusters || (populationVAE+i)->nbrCluster < min_clusters)
        {
            (populationVAE+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationVAE+i)->contrainteViole++;
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
            (populationVAE+i)->clustersSize[j]=0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if (tabTmp[j] == (populationVAE+i)->phenotype[k])
                {
                    (populationVAE+i)->clustersSize[j]=(populationVAE+i)->clustersSize[j]+1;
                }
            }
            if((populationVAE+i)->clustersSize[j]>max_sizeCluster || (populationVAE+i)->clustersSize[j]< min_sizeCluster)
            {
                ///24.07.2015 ---> l'ajout de vecteur de contraintes = le nombre de cluster dont la taille a dépassé la borne
                (populationVAE+i)->constraintVector[1]++;

            }
        }

        if((populationVAE+i)->constraintVector[1] != 0)
        {
            (populationVAE+i)->contrainteViole++;
        }
   }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0){
        checkCohabitationAndNonCohabitationConstrainteVAE(populationVAE);
    }

}

*/
