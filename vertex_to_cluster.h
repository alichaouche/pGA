#ifndef VERTEX_TO_CLUSTER_H_INCLUDED
#define VERTEX_TO_CLUSTER_H_INCLUDED





void generatePopulationVTC(partitionVTC* populationVTC, int indiceFirstElt);

void affichePopulationVTC(partitionVTC* populationVTC);
void calculCoutCoupeEtFitnessVTC(partitionVTC* populationVTC);
void crossOverVTC(partitionVTC* populationVTC1, partitionVTC* populationVTC2);
void naturalSelectionVTC(partitionVTC* populationVTC1,partitionVTC* populationVTC2);
int findTheBestSolutionVTC(partitionVTC *populationVTC);
void mutationVTC(partitionVTC* populationVTC);
float testerLaSommeDesFitnessVTC(partitionVTC* populationVTC);
void displayTheBestSolutionVTC(partitionVTC* solutionDominante);
void checkContrainstAndFitnessPenalizationVTC(partitionVTC *populationVTC);


void writeSolutionInFileVTC(partitionVTC *populationVTC,FILE *outputFilePop,int iteration);
void writeBestSolutionInFileVTC(partitionVTC *bestSolution, FILE *outputFile,int iteration);
///void vertexToClusterEncoding(int nbrGeneration ,FILE *outputFileVTC,FILE *outputFilePopVTC,FILE *outputOptimalSolutionFileVTC, partitionVTC *populationVTC1,partitionVTC *populationVTC2,partitionVTC *solutionDominante);
void vertexToClusterEncoding(int nbrGeneration ,FILE *outputFileVTC,FILE *outputFilePopVTC,FILE *outputOptimalSolutionFileVTC,
                             partitionVTC *populationVTC1,partitionVTC *populationVTC2,partitionVTC *solutionDominante,
                             int iteration , int *bestSolutionIteration , int *nbrApparition);

void checkCohabitationAndNonCohabitationConstrainteVTC(partitionVTC *populationVTC);
///void writeOptimalSolutionInFileVTC(partitionVTC *solutionDominante,FILE *outputOptimalSolutionFileVTC);
int compareCroissantFitnessVTC (void const *a, void const *b);
void writeOptimalSolutionInFileVTC(partitionVTC *solutionDominante,FILE* outputOptimalSolutionFileVTC,
                                  int nbrRun, int bestSolutionIteration, float runTime, int ES);

#endif // VERTEX_TO_CLUSTER_H_INCLUDED
