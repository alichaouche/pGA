#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "./communesFunctions.h"
#include "./pMediansProblemBinaryEncoding.h"
#include "./compilationConditionnelle.h"
#define affectationArretes 0

extern int taillePopulation,tauxReproduction,sommeTotalFlux,nbrNoeuds,
       nbr_constraint,nbrArretes,nbrParties,nbrNoeudsInCluster;
extern double tauxMutation ,averageTraffic;
extern int fluxMatrix[1000][1000],fluxVector[1000],cocyclesDeBase[1000],cocyclesDeBaseEntier[1000], maximumDistanceMatrix[1000][1000],shortestPath[1000][1000];
extern edge edgeVector[1000],cohabitationConstraintes[100],nonCohabitationConstraintes[100];
extern int max_clusters , min_clusters,min_sizeCluster ,max_sizeCluster ;
extern int bestSolution, nbrApparition, nbrCohabitationConstraintes,nbrNonCohabitationConstraintes;


///***************************************************************************
///                     les corps des fonctions
///***************************************************************************

///****************************************************************************************************************************
///génération de la population initiale
void generatePopulationPMP(partitionPMP* population)
{

    int i,j,indice;
    /// génération des zero vertices
    for(i=0; i< taillePopulation; i++)
    {
        (population+i)->id =i;
        /// affecter la valeur 0 à tous les noeuds de la solution
        for(j=0; j<nbrNoeuds; j++)
        {
            (population+i)->genotype[j] = 0;
        }
        (population+i)->medians =rnd(2,nbrNoeuds-2);
        ///printf("le nombre de medains dans cette solutions : %d \n",(population+i)->medians);
        /// choisir aléatoirement les noeuds medians
        for(j=0; j<(population+i)->medians; j++) /// <= pour qu'on aura un élément en plus pour le nombre des parties
        {
            do
            {
                indice = rnd(0,nbrNoeuds);
            }
            while((population+i)->genotype[indice] == 1);
            (population+i)->genotype[indice] = 1;
        }
        (population+i)->nonMedians = nbrNoeuds - (population+i)->medians;
    }
}
///****************************************************************************************************************************
/// affichage des population
void affichePopulationPMP(partitionPMP* population)
{
    int i,j;
    for(i=0; i<taillePopulation; i++)
    {
        if((population+i)->contrainteViole == 0)
        {
            printf("\nid = %d \t", (population+i)->id);
            printf("la genotypes est  = \t");
            for(j=0; j<nbrNoeuds; j++)
            {
                printf("%d ",(population+i)->genotype[j]);
            }
            printf("\n");
            /// affichage des medians de cette solution
            printf("  les medians sont  = \t");
            for(j=0; j<(population+i)->medians; j++)
            {
                printf("%d ",(population+i)->mediansVertices[j]);
            }
            printf("\n");

            printf("  les NOT  medians sont  = \t");
            for(j=0; j<(population+i)->nonMedians; j++)
            {
                printf("%d ",(population+i)->nonMediansVertices[j]);
            }
            printf("\n");

            printf("\n la phenotypes est  = \t");
            for(j=0; j<nbrNoeuds; j++)
            {
                printf("%d ",(population+i)->phenotype[j]);
            }
            printf("\n");
            printf("la fitness de la %d solution est : %0.2f",i,(population+i)->fitness );
            printf("\n**************************************************************************\n");
        }
        else if((population+i)->contrainteViole == 1)
        {
            printf("la %d solution a viole la contrainte de capacite \n",i);
        }
    }

}

void calculerMediansNonMediansVertices(partitionPMP* population)
{

    int i,j,testValueGenotype;
#if mouchard
    printf("calculerMediansNonMediansVertices : ...\n");
#endif // mouchard
    for(i=0; i<taillePopulation; i++)
    {
        (population+i)->medians = 0;
        (population+i)->nonMedians = 0;
        for(j=0; j<nbrNoeuds; j++)
        {
            testValueGenotype = (population+i)->genotype[j];
            switch (testValueGenotype)
            {
            case 0 :
                (population+i)->nonMediansVertices[(population+i)->nonMedians] = j;
                (population+i)->nonMedians++;
                break;
            case 1 :
                (population+i)->mediansVertices[(population+i)->medians] = j;
                (population+i)->medians++;
                break;
            default :
                printf("Cette valeur n\'est pas prise en consédiration \n");
                break;
            }
        }
    }
}
///****************************************************************************************************************************
/// cette fonction affecte les noeuds non medians au noeuds medians, permettant la création des
/// PHENOTYPE DES SOLUTION CE QUI PERMET PAR LA SUITE LE CALCULE DES FITNESS
///void affectationDesNonMediansBaseSurLesArretes(partition* population,int taillePopulation, int nbrNoeuds, int fluxMatrix[nbrNoeuds][nbrNoeuds],int max_sizeCluster){
void affectationDesNonMedians(partitionPMP* population)
{

    int i,j,l,k,maxDistance,indiceNotSeed,indiceSeed,_medians, _mediansTmp,_nonMedians, _nonMediansTmp,nbrSommets,nbrG=0;
    int _nonMediansVertices[100]= {0},_mediansVertices[100]= {0};
#if affectationArretes
    for(i=0; i<taillePopulation; i++)
    {
        /// initialisation de vecteur des capacité avec 1 pour chaque médians
        /// initialisation des phénotype par la valeurs
        _medians = _mediansTmp = (population+i)->medians;;
        _nonMedians =  _nonMediansTmp = (population+i)->nonMedians;;

        for(l=0; l<_nonMedians; l++)
        {
            _nonMediansVertices[l] = (population+i)->nonMediansVertices[l];
            (population+i)->phenotype[_nonMediansVertices[l]] = _nonMediansVertices[l];
            (population+i)->clustersSize[l]= 0;
        }

        for(l=0; l<_medians; l++)
        {
            _mediansVertices[l] = (population+i)->mediansVertices[l];
            (population+i)->phenotype[_mediansVertices[l]] = _mediansVertices[l];
            (population+i)->clustersSize[_mediansVertices[l]] = 1;

        }

        while(_nonMedians > 0 && _medians > 0)
        {
            j =0;
            _nonMediansTmp = _nonMedians ;
            _mediansTmp = _medians;
            maxDistance = 0;
            //mon_sleep(1);
            while(j < _nonMediansTmp)
            {
                maxDistance = 0;
                indiceNotSeed= _nonMediansVertices[j];

                for(k=0; k<_mediansTmp; k++)
                {
                    indiceSeed = _mediansVertices[k];
                    if(fluxMatrix[indiceNotSeed][indiceSeed] > 0 )
                    {
                        if(maxDistance == 0)
                        {
                            maxDistance = fluxMatrix[indiceNotSeed][indiceSeed];
                            (population+i)->phenotype[indiceNotSeed] = (population+i)->phenotype[indiceSeed];

                        }
                        else if((population+i)->clustersSize[(population+i)->phenotype[indiceSeed]] < max_sizeCluster)
                        {
                            if((population+i)->clustersSize[(population+i)->phenotype[indiceNotSeed]] < max_sizeCluster)
                            {
                                if(maxDistance < fluxMatrix[indiceNotSeed][indiceSeed])
                                {
                                    maxDistance = fluxMatrix[indiceNotSeed][indiceSeed];
                                    (population+i)->phenotype[indiceNotSeed] = (population+i)->phenotype[indiceSeed];
                                }
                            }
                            else if((population+i)->clustersSize[(population+i)->phenotype[indiceNotSeed]] > max_sizeCluster)
                            {
                                maxDistance = fluxMatrix[indiceNotSeed][indiceSeed];
                                (population+i)->phenotype[indiceNotSeed] = (population+i)->phenotype[indiceSeed];
                            }
                        }
                    }
                }

                if(maxDistance != 0)
                {

                    indiceSeed = (population+i)->phenotype[indiceNotSeed] ;
                    (population+i)->clustersSize[(population+i)->phenotype[indiceSeed]]++;
                    _mediansVertices[_medians] = indiceNotSeed;
                    _medians++;

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
        nbrSommets = nbrNoeuds;
        l=0;
        while(l<nbrSommets)
        {
            if((population+i)->clustersSize[l]==0)
            {
                for(k=l; k<nbrSommets; k++)
                {
                    (population+i)->clustersSize[k] = (population+i)->clustersSize[k+1];
                }
                nbrSommets--;
            }
            else
            {
                l++;
            }
        }

    }

#else
    for(i=0; i<taillePopulation; i++)
    {

        for(l=0; l<(population+i)->nonMedians; l++)
        {
            (population+i)->phenotype[(population+i)->nonMediansVertices[l]] = (population+i)->nonMediansVertices[l];
            (population+i)->clustersSize[l]= 0;
        }

        for(l=0; l<(population+i)->medians; l++)
        {
            (population+i)->phenotype[(population+i)->mediansVertices[l]] = (population+i)->mediansVertices[l];
            (population+i)->clustersSize[(population+i)->mediansVertices[l]] = 1;

        }

        if((population+i)->nonMedians > 0 && (population+i)->medians > 0)
        {
            j =0;
            while(j < (population+i)->nonMedians)
            {
                maxDistance = 0;
                indiceNotSeed= (population+i)->nonMediansVertices[j];

                for(k=0; k<(population+i)->medians; k++)
                {
                    indiceSeed = (population+i)->mediansVertices[k];
                    if(fluxMatrix[indiceNotSeed][indiceSeed] > 0 && (population+i)->clustersSize[(population+i)->phenotype[indiceSeed]] < max_sizeCluster)
                    {
                        if(maxDistance == 0)
                        {
                            maxDistance = fluxMatrix[indiceNotSeed][indiceSeed];
                            (population+i)->phenotype[indiceNotSeed] = (population+i)->phenotype[indiceSeed];

                        }
                        else  if(maxDistance < fluxMatrix[indiceNotSeed][indiceSeed])
                        {
                            maxDistance = fluxMatrix[indiceNotSeed][indiceSeed];
                            (population+i)->phenotype[indiceNotSeed] = (population+i)->phenotype[indiceSeed];
                        }
                    }
                }
                if(maxDistance != 0)
                {
                    indiceSeed = (population+i)->phenotype[indiceNotSeed] ;
                    (population+i)->clustersSize[(population+i)->phenotype[indiceSeed]]++;
                }
                else
                {
                    for(k=0; k<(population+i)->medians; k++)
                    {
                        indiceSeed = (population+i)->mediansVertices[k];
                        if((population+i)->clustersSize[(population+i)->phenotype[indiceSeed]] < max_sizeCluster)
                        {
                            maxDistance = fluxMatrix[indiceNotSeed][indiceSeed];
                            (population+i)->phenotype[indiceNotSeed] = (population+i)->phenotype[indiceSeed];
                            (population+i)->clustersSize[(population+i)->phenotype[indiceNotSeed]]++;
                            break;
                        }
                    }
                    if(maxDistance == 0)
                    {
                        (population+i)->phenotype[indiceNotSeed] = (population+i)->phenotype[(population+i)->mediansVertices[0]];
                        (population+i)->clustersSize[(population+i)->phenotype[indiceNotSeed]]++;
                    }
                }
                nbrSommets = nbrNoeuds;
                l=0;
                while(l<nbrSommets)
                {
                    if((population+i)->clustersSize[l]==0)
                    {
                        for(k=l; k<nbrSommets; k++)
                        {
                            (population+i)->clustersSize[k] = (population+i)->clustersSize[k+1];
                        }
                        nbrSommets--;
                    }
                    else
                    {
                        l++;
                    }
                }
                j++;
            }
        }
    }

#endif

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

void calculcoutCoupeeEtFitnessPMP(partitionPMP* population)
{
    int i,j,k,sommeTotalCoutCoupe = 0;

    for(i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (population+i)->coutCoupe = 0; /// réinitialisation de la valriable après chaque itération (individus de la population)
        for(j=0; j<nbrNoeuds-1; j++)
        {
            for(k=j+1; k<nbrNoeuds; k++)
            {
                if((population+i)->phenotype[j] == (population+i)->phenotype[k] && fluxMatrix[j][k] > 0)
                {
                    (population+i)->coutCoupe = (population+i)->coutCoupe +  fluxMatrix[j][k];
                }
            }
        }
        (population+i)->coutCoupeNormalise = (population+i)->coutCoupe + ((nbr_constraint - (population+i)->contrainteViole)*sommeTotalFlux);
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (population+i)->coutCoupeNormalise;
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
            varianceCoutDeCoupeNormalise + pow(((population+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);
    /// tester la valeur ajouté
    for(i=0; i<taillePopulation; i++)
    {

        /// calcul de sigma truncation
        (population+i)->expectedValue = (population+i)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (population+i)->expectedValue;
    }


    /// calcule du fitness
    /// le 26/11/2015 : calcul du fitness des solution avec leur expectedValues
    for(i=0; i<taillePopulation ; i++)
    {
        (population+i)->fitness = (float)((population+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);
    }
#else

    FILE *file;
    file = fopen("D:/WorkingFolderAgProgramming/benchmark/pMedianEncoding/expectedValuePMP.txt","w");
    for(i=0; i<taillePopulation ; i++)
    {
        (population+i)->fitness = (float)((population+i)->coutCoupeNormalise)/(float)(sommeTotalCoutCoupe);
    }
#endif// scaling

}

///****************************************************************************************************************************
/// Cette fonction assure la séléction des individu pour la pprocedure d'appariement
/// elle émite la roue de lottrie biaisé
void naturalSelectionPMP(partitionPMP* population1,partitionPMP* population2)
{
    int i,j;
    double lotterie,sommeFitness;
    /** je devrais récuperer les meilleurs solution est les injecter dans la nouvelle génération
        il faut que la solution soit générique pour n'importe quel nombre d'individus a préserver.
    */
    /** cette partie concerne la reproduction*/
    partitionPMP *tmpPopulation,*tmpInverserSolution;
    tmpPopulation = (partitionPMP*)malloc(taillePopulation*sizeof(partitionPMP));
    tmpInverserSolution = (partitionPMP*)malloc(sizeof(partitionPMP));

    /// affecter la population1 à tmpPopulation
    for(i=0; i<taillePopulation; i++)
    {
        *(tmpPopulation+i) = *(population1+i);
    }
    /// trier la nouvelle population pour en choisir les meilleur solution
    for(i=0; i<taillePopulation; i++)
    {
        for(j=i+1; j<taillePopulation; j++)
        {
            if((tmpPopulation+i)->fitness < (tmpPopulation+j)->fitness)
            {
                *tmpInverserSolution = *(tmpPopulation+i);
                *(tmpPopulation+i) = *(tmpPopulation+j);
                *(tmpPopulation+j) = *tmpInverserSolution;
            }
        }
    }
    qsort (tmpPopulation, taillePopulation, sizeof *tmpPopulation, compareCroissantFitnessPMP);
    /** injecter les k premières solutions dans population2*/

    for(i=0; i<tauxReproduction; i++)
    {
        *(population2+i) = *(tmpPopulation+i);
    }

    for(i=tauxReproduction; i<taillePopulation; i++)
    {
#if Windows
        lotterie= drand48ForWindows(0,101);
#else
        lotterie= drand48();
#endif // Windows

        sommeFitness = 0;

        for(j=0; j<taillePopulation; j++)
        {
            sommeFitness = sommeFitness + (population1+j)->fitness;
            if (lotterie<=sommeFitness)
            {
                *(population2+i) = *(population1+j);
                break;
            }
            else if( j==taillePopulation-1)
            {
                *(population2+i) = *(population1+j);
            }
        }
    }

    free(tmpPopulation);

}
///****************************************************************************************************************************
void crossOverPMP(partitionPMP* population1, partitionPMP* population2)
{
    int i = 0,j;
    int choixInd1, choixInd2, choixLocus;

    while (i < taillePopulation)
    {
        choixInd1 = rnd(0,taillePopulation);
        do
        {
            /// ce traitement évite l'appariement d'un element avec lui même/
            choixInd2 = rnd(0,taillePopulation);
        }
        while(choixInd2 == choixInd1 );

        choixLocus = rnd(0,nbrNoeuds-1); /// -2 parce qu'on a commencer à partir du 0

        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (population2+i)->genotype[j] = (population1+choixInd1)->genotype[j];
            (population2+i+1)->genotype[j] = (population1+choixInd2)->genotype[j];

        }
        ///*************************************************************************************************
        for(j=choixLocus+1; j<nbrNoeuds; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (population2+i)->genotype[j] = (population1+choixInd2)->genotype[j];
            (population2+i+1)->genotype[j] = (population1+choixInd1)->genotype[j];
        }
        /// affectation de nouveau indice
        (population2+i)->id = i;
        (population2+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }

}
///****************************************************************************************************************************
int findTheBestSolutionPMP(partitionPMP *population)
{
    float maxFitness = population->fitness;
    int i,indice=0;
    for(i=0; i<taillePopulation ; i++)
    {
        if(maxFitness < (population+i)->fitness)
        {
            maxFitness = (population+i)->fitness;
            indice = i;
        }
    }
#if mouchard
    printf("\n maxFitness = %0.4f \t indice = %d \n",maxFitness, indice);
#endif // mouchard
    return indice;
}
///****************************************************************************************************************************
void mutationPMP(partitionPMP* population)
{

    int i,numeroGeneMute;
    double applicationOfMutation;
    for(i=0; i<taillePopulation; i++)
    {
#if Windows
        applicationOfMutation = drand48ForWindows(0,101);
#else
        applicationOfMutation = drand48();
#endif // Windows
        if(applicationOfMutation <= tauxMutation)
        {
            /// choisir un gène aléatoirement pour lui appliquer la mutation
            numeroGeneMute = rnd(0,nbrNoeuds);
            ///((population+i)->genotype[numeroGeneMute] == 1)? 0 : 1;
            if((population+i)->genotype[numeroGeneMute] == 1)
            {
                (population+i)->genotype[numeroGeneMute] =0;
            }
            else
            {
                (population+i)->genotype[numeroGeneMute] =1;
            }

        }
    }
}
///****************************************************************************************************************************
float testerLaSommeDesFitnessPMP(partitionPMP* population)
{
    printf("testerLaSommeDesFitness ...\n");
    int i;
    float sommeFitness = 0;
    for(i=0; i<taillePopulation ; i++)
    {
        sommeFitness = sommeFitness + (population+i)->fitness;
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
void writeSolutionInFilePMP(partitionPMP *population, FILE *outputFilePop,int iteration)
{

    int i,j;
    for(i=0; i<taillePopulation; i++)
    {

        ///écriture de la meilleur solution dans le fichier output
        fprintf(outputFilePop,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,(population+i)->id,
                (population+i)->coutCoupe,(population+i)->fitness,(population+i)->coutCoupeNormalise,
                (population+i)->contrainteViole,(population+i)->medians);
        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(outputFilePop,"%d ",(population+i)->genotype[j]);
        }
        fprintf(outputFilePop,"\t\t\t\t");
        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(outputFilePop,"%d ",(population+i)->phenotype[j]);
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


///************************************************************************************************************
void pMedianEncoding(int nbrGeneration ,FILE *outputFilePMP,FILE *outputFilePopPMP,FILE *outputOptimalSolutionFilePMP)
{

    int bestSolution,iteration=0;
    clock_t t1,t2;
    double temps,averageTrafficPMP=0;


    partitionPMP *population1,*population2, *solutionDominante;
    population1 = (partitionPMP*)malloc(taillePopulation*sizeof (partitionPMP));
    population2 = (partitionPMP*)malloc(taillePopulation*sizeof(partitionPMP));
    solutionDominante = (partitionPMP*)malloc(sizeof(partitionPMP));

    t1=clock();

    generatePopulationPMP(population1);
    calculerMediansNonMediansVertices(population1);
    affectationDesNonMedians(population1);
    checkContrainstAndFitnessPenalizationPMP(population1);
    calculcoutCoupeeEtFitnessPMP(population1);
#if writingPopulationInFile
    writeSolutionInFilePMP(population1,outputFilePopPMP,iteration);
#endif // writingPopulationInFile


    ///initialisation de la solution dominante
    bestSolution=findTheBestSolutionPMP(population1);
    *solutionDominante=*(population1+bestSolution);
    averageTrafficPMP = averageTrafficPMP + (population1+bestSolution)->coutCoupeNormalise;
#if writingPopulationInFile
    writeBestSolutionInFilePMP(solutionDominante,outputFilePMP,iteration);
#endif // writingPopulationInFile
    nbrApparition=1;


    while(iteration<nbrGeneration && nbrApparition <=300 )
    {

        naturalSelectionPMP(population1,population2);
        crossOverPMP(population2, population1);
        mutationPMP(population1);
        calculerMediansNonMediansVertices(population1);
        affectationDesNonMedians(population1);
        checkContrainstAndFitnessPenalizationPMP(population1);
        calculcoutCoupeeEtFitnessPMP(population1);

#if writingPopulationInFile
        writeSolutionInFilePMP(population1,outputFilePopPMP,iteration);
#endif // writingPopulationInFile
///=====================================================================================================
        /// dans ce qui suit, on va juste localiser la meilleur solution
///=====================================================================================================

        bestSolution=findTheBestSolutionPMP(population1);
        averageTrafficPMP = averageTrafficPMP + (population1+bestSolution)->coutCoupeNormalise;
#if writingPopulationInFile
        writeBestSolutionInFilePMP((population1+bestSolution),outputFilePMP,iteration);
#endif // writingPopulationInFile
        if((population1+bestSolution)->coutCoupeNormalise > solutionDominante->coutCoupeNormalise)
        {
            *(solutionDominante) = *(population1+bestSolution);
            nbrApparition=1;
        }
        else
        {
            nbrApparition++;
        }


        iteration++;
    }

    t2=clock();
    displayTheBestSolutionPMP(solutionDominante);
    writeOptimalSolutionInFilePMP(solutionDominante,outputOptimalSolutionFilePMP);
    temps = (double)(t2-t1)/CLOCKS_PER_SEC;
    printf("le temps d execution est %lf \n", temps);

    averageTrafficPMP = averageTrafficPMP / nbrGeneration;
    fprintf(outputFilePMP,"\n\n");
    writeBestSolutionInFilePMP(solutionDominante,outputFilePMP,iteration);
    fprintf(outputFilePMP,"\n\nthe average traffic = %lf\n",averageTrafficPMP);
    fprintf(outputFilePMP,"\n\nLe nombre d'apparition de la meilleur solution est  = %d\n",nbrApparition);
    fprintf(outputFilePMP,"\n\n==============================================================================\n");
    fprintf(outputFilePopPMP,"\n\n==============================================================================\n");

    free(population1);
    free(population2);
    free(solutionDominante);


}


///**************************************************************************************
void checkCohabitationAndNonCohabitationConstraintePMP(partitionPMP *population)
{

    int i,j,noeud1,noeud2;
    for(i=0; i<taillePopulation; i++)
    {
        (population+i)->constraintVector[2]=0;
        (population+i)->constraintVector[3]=0;

        if(nbrCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrCohabitationConstraintes; j++)
            {
                noeud1=cohabitationConstraintes[j].nouedDepart;
                noeud2=cohabitationConstraintes[j].nouedArrive;

                if((population+i)->phenotype[noeud1]!= (population+i)->phenotype[noeud2])
                {
                    (population+i)->constraintVector[2]++;
                }
            }
        }
        /// cette instruction permet d'avoir le nombre de contraintes violées
        if((population+i)->constraintVector[2]!=0)
        {
            (population+i)->contrainteViole++;
        }

        if(nbrNonCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrNonCohabitationConstraintes; j++)
            {
                noeud1=nonCohabitationConstraintes[j].nouedDepart;
                noeud2=nonCohabitationConstraintes[j].nouedArrive;

                if((population+i)->phenotype[noeud1]== (population+i)->phenotype[noeud2])
                {
                    (population+i)->constraintVector[3]++;
                }
            }
        }
        /// cette instruction permet d'avoir le nombre de contraintes violées
        if((population+i)->constraintVector[3]!=0)
        {
            (population+i)->contrainteViole++;
        }

    }

}

///************************************************************************************************************
void checkContrainstAndFitnessPenalizationPMP(partitionPMP *population)
{
    int i,j;
    for(i=0; i<taillePopulation; i++)
    {
        ///initialisation de vecteur des contraintes et de la variable de contraintes violées
        (population+i)->contrainteViole=0;
        for(j=0; j<nbr_constraint; j++)
        {
            (population+i)->constraintVector[j]=0;
        }

        for(j=0; j<(population+i)->medians; j++)
        {
            if((population+i)->clustersSize[j]>max_sizeCluster || (population+i)->clustersSize[j]< min_sizeCluster)
            {
                ///24.07.2015 ---> l'ajout de vecteur de contraintes = le nombre de cluster dont la taille a dépassé la borne
                (population+i)->constraintVector[1]++;

            }
        }



        if((population+i)->medians > max_clusters || (population+i)->medians < min_clusters)
        {
            (population+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (population+i)->contrainteViole++;
        }

        if((population+i)->constraintVector[1] !=0)
        {
            (population+i)->contrainteViole++;
        }
    }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0)
    {
        checkCohabitationAndNonCohabitationConstraintePMP(population);
    }

}
///************************************************************************************************************
void writeOptimalSolutionInFilePMP(partitionPMP *solutionDominante,FILE* outputOptimalSolutionFilePMP)
{
    int i;
    fprintf(outputOptimalSolutionFilePMP,"id = %d | normalised cut size = %d | cut size intrA = %d cut size intEr = %d ",
            solutionDominante->id,solutionDominante->coutCoupeNormalise, solutionDominante->coutCoupe, (sommeTotalFlux-solutionDominante->coutCoupe));
    fprintf(outputOptimalSolutionFilePMP,"le contraintes violees = %d \n", solutionDominante->contrainteViole);
    for(i=0; i<4; i++)
    {
        switch (i)
        {

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

}

///******************************************************************************************
int compareCroissantFitnessPMP (void const *a, void const *b)
{
    partitionPMP const *pa = a;
    partitionPMP const *pb = b;
    return  pb->coutCoupeNormalise - pa->coutCoupeNormalise;
}
/*** La fonction d'affectationDesNonMedians() basé sur le plus court chemin et la matrices des distances max
void affectationDesNonMedians(partitionPMP* population)
{

    int i,j,l,k,maxDistance,indiceNotSeed,indiceSeed,nbrClusterSature;
    int affecter,plusCourtChemin; /// pour tester si le noeuds est affecter à un cluster ou pas


    for(i=0; i<taillePopulation; i++)
    {
        nbrClusterSature = 0;
        /// initialisation de vecteur des capacité avec 1 pour chaque médians
        /// initialisation des phénotype par la valeurs
        for(l=0; l<(population+i)->medians; l++)
        {
            (population+i)->clustersSize[l] = 1;
            (population+i)->phenotype[(population+i)->mediansVertices[l]] = (population+i)->mediansVertices[l];

        }
        /// for each individual in the population
        for(j=0; j<(population+i)->nonMedians; j++)
        {
            affecter = -1; /// il faut la réinitialiser pour chaque itération pour chaque !Seed
            indiceNotSeed= (population+i)->nonMediansVertices[j];

            for(l=0; l<(population+i)->medians; l++)
            {
                indiceSeed = (population+i)->mediansVertices[l];
                if((population+i)->clustersSize[l] < max_sizeCluster && fluxMatrix[indiceNotSeed][indiceSeed] > 0)
                {
                    (population+i)->phenotype[indiceNotSeed] = indiceSeed; /// initialiser l'affectation de !Seed au premier Seed
                    (population+i)->clustersSize[l]++;
                    ///03.07.2015_19.26 : quel est le numéro du cluster auquel on a affecté le noeud !Seed
                    maxDistance = fluxMatrix[indiceNotSeed][indiceSeed];
                    affecter = l;
                    break;
                }
            }

            if(affecter!=-1){
                /// je commence les teste à partir du median qui vient just après celui qu'on lui affecter le !Seed
                for(l=affecter+1; l<(population+i)->medians; l++)
                {
                    indiceSeed = (population+i)->mediansVertices[l];
                    /// les liens entre medians et non medians sont des liens intraCluster ==> il faut charcher les distance max
                    if(maxDistance < fluxMatrix[indiceNotSeed][indiceSeed] && (population+i)->clustersSize[l]<max_sizeCluster)  /// rajouter la contrainte de capacitée des medians
                    {
                        maxDistance = fluxMatrix[indiceNotSeed][indiceSeed];
                        (population+i)->phenotype[indiceNotSeed] = indiceSeed; /// indiceNotSeed = non medians
                        affecter = l;/// pour récuperer l'indice de nouveau cluster qui acceuillir le !Seed pour pouvoir l'utiliser la prochaine fois
                    }
                }
            }else {
                /// si le !Seed n'a pas été affecter : soit tous les cluster sont saturé, soit il n'est relié à aucun median
                if(nbrClusterSature < (population+i)->medians){
                    for(k=0; k<(population+i)->medians; k++){
                        if((population+i)->clustersSize[k] >= max_sizeCluster){nbrClusterSature++;}
                    }
                }
                if(nbrClusterSature >= (population+i)->medians){
                    for(k=0; k<(population+i)->medians; k++){
                        indiceSeed = (population+i)->mediansVertices[k];
                        ///affecter le sommet au dernier Median auquel il est lié
                        if(fluxMatrix[indiceNotSeed][indiceSeed] > 0){
                            (population+i)->phenotype[indiceNotSeed] = indiceSeed;
                            affecter = k;
                        }
                    }
                } else {

                        for(k=0; k<(population+i)->medians; k++){
                            indiceSeed = (population+i)->mediansVertices[k];
                            if((population+i)->clustersSize[k] < max_sizeCluster ){
                                (population+i)->phenotype[indiceNotSeed] = indiceSeed;
                                maxDistance = maximumDistanceMatrix[indiceNotSeed][indiceSeed];
                                plusCourtChemin = shortestPath[indiceNotSeed][indiceSeed];
                                affecter = k;
                                break;
                            }
                        }

                        for(k=affecter+1; k<(population+i)->medians; k++){
                            indiceSeed = (population+i)->mediansVertices[k];

                            if(plusCourtChemin == shortestPath[indiceNotSeed][indiceSeed]){

                                if( maxDistance < maximumDistanceMatrix[indiceNotSeed][indiceSeed] && (population+i)->clustersSize[k] < max_sizeCluster){
                                    (population+i)->phenotype[indiceNotSeed] = indiceSeed;
                                    maxDistance = maximumDistanceMatrix[indiceNotSeed][indiceSeed];
                                    affecter = k;
                                }
                            }
                            else if(plusCourtChemin > shortestPath[indiceNotSeed][indiceSeed]){
                                    (population+i)->phenotype[indiceNotSeed] = indiceSeed;
                                    maxDistance = maximumDistanceMatrix[indiceNotSeed][indiceSeed];
                                    affecter = k;
                            }
                        }
                }
            }
            (population+i)->clustersSize[affecter]++;
        }
    }
}
*/
