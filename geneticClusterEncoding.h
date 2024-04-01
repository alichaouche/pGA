#ifndef GENETICCLUSTERENCODING_H_INCLUDED
#define GENETICCLUSTERENCODING_H_INCLUDED

void generatePopulationGCE(partitionGCE* populationGCE);

void affichePopulationGCE(partitionGCE* populationGCE);
void calculCoutCoupeEtFitnessGCE(partitionGCE* populationGCE);
void crossOverGCE(partitionGCE* populationGCE1, partitionGCE* populationGCE2);
void naturalSelectionGCE(partitionGCE* populationGCE1,partitionGCE* populationGCE2);
int findTheBestSolutionGCE(partitionGCE *populationGCE);
void mutationGCE(partitionGCE* populationGCE);
float testerLaSommeDesFitnessGCE(partitionGCE* populationGCE);
void displayTheBestSolutionGCE(partitionGCE* solutionDominante);
void checkContrainstAndFitnessPenalizationGCE(partitionGCE *populationGCE);


void writeSolutionInFileGCE(partitionGCE *populationGCE,FILE *outputFilePop,int iteration);
void writeBestSolutionInFileGCE(partitionGCE *bestSolution, FILE *outputFile,int iteration);
void vertexToClusterEncoding(int nbrGeneration ,FILE *outputFileGCE,FILE *outputFilePopGCE,FILE *outputOptimalSolutionFileGCE, partitionGCE *populationGCE1,partitionGCE *populationGCE2,partitionGCE *solutionDominante);
void checkCohabitationAndNonCohabitationConstrainteGCE(partitionGCE *populationGCE);
///void writeOptimalSolutionInFileGCE(partitionGCE *solutionDominante,FILE *outputOptimalSolutionFileGCE);
int compareCroissantFitnessGCE (void const *a, void const *b);
void writeOptimalSolutionInFileGCE(partitionGCE *solutionDominante,FILE* outputOptimalSolutionFileGCE,
                                  int nbrRun, int bestSolutionIteration, float runTime, int ES);
void getPhenotypeFromGenotype(partitionGCE* populationGCE);

#endif // VERTEX_TO_CLUSTER_H_INCLUDED


