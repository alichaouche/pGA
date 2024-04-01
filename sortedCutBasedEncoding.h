#ifndef SORTEDCUTBASEDENCODING_H_INCLUDED
#define SORTEDCUTBASEDENCODING_H_INCLUDED



void affichePopulationIC(partitionIC* populationIC);
void calculerGenotypeIC(partitionIC * populationIC, int indiceFirstElt);
void calculCoutCoupeEtFitnessIC(partitionIC* populationIC);
void naturalSelectionIC(partitionIC* populationIC1,partitionIC* populationIC2);
void crossOverIC(partitionIC* populationIC1, partitionIC* populationIC2);
void mutationIC(partitionIC* populationIC);

float testerLaSommeDesFitnessIC(partitionIC* populationIC);
int findTheBestSolutionIC(partitionIC *populationIC);
void displayTheBestSolutionIC(partitionIC* solutionDominante);
void writeBestSolutionInFileIC(partitionIC *solutionDominante, FILE *outputFile,int iteration);
void writeSolutionInFileIC(partitionIC *populationIC, FILE *outputFilePop,int iteration);
void getPartitionFromSolutionIC(partitionIC *populationIC);
void calculCoutCoupeEtFitnessWithFlowVectorIC(partitionIC* populationIC);
void checkContrainstAndFitnessPenalizationIC(partitionIC *populationIC);

///========================================================================================================================
void generatePopulationWithoutRedandancyIC(partitionIC* populationIC, int indiceFirstElt);
///void genererSolutionEntiere(ul A[]);
void genererSolutionEntiere(partitionIC* populationIC, int indiceFirstElt);
int  existanceDeSolution(int indexNewSolution, partitionIC* populationIC);
void genererLesSolutionBinaireIC(partitionIC *populationIC, int indiceFirstElt);
///========================================================================================================================
///void cutBasedEncodingWithoutRedandancyIC(int nbrGeneration,FILE *outputFileIC,FILE *outputFilePopIC,FILE *outputOptimalSolutionFileIC,partitionIC *populationIC1,partitionIC *populationIC2,partitionIC *solutionDominante);
void cutBasedEncodingWithoutRedandancyIC(int nbrGeneration,FILE *outputFileIC,FILE *outputFilePopIC,FILE *outputOptimalSolutionFileIC,partitionIC *populationIC1,partitionIC *populationIC2, partitionIC *solutionDominante,
                                         int iteration , int *bestSolutionIteration , int *nbrApparition);


void checkCohabitationAndNonCohabitationConstrainteIC(partitionIC *populationIC);
///void writeOptimalSolutionInFileIC(partitionIC *solutionDominante,FILE *outputOptimalSolutionFileIC);
int compareCroissantFitnessIC (void const *a, void const *b);
void writeOptimalSolutionInFileIC(partitionIC *solutionDominante,FILE* outputOptimalSolutionFileIC,
                                    int nbrRun, int bestSolutionIteration, float runTime, int ES);

#endif // SORTEDCUTBASEDENCODING_H_INCLUDED
