#ifndef BINARYENCODING_H_INCLUDED
#define BINARYENCODING_H_INCLUDED


///*******************************************************************
///*************** Les entêtes des fonctions ************************///
///*******************************************************************


void generatePopulationEA(partitionEA* populationEA, int indiceFirstElt);
void generatePopulationRandomlyEA(partitionEA* populationEA);
void affichePopulationEA(partitionEA* populationEA);
void calculCoutCoupeEtFitnessWithFlowVectorEA(partitionEA* populationEA);

void naturalSelectionEA(partitionEA* populationEA1,partitionEA* populationEA2);
void crossOverEA(partitionEA* populationEA1, partitionEA* populationEA2);

int findTheBestSolutionEA(partitionEA *populationEA);
void mutationEA(partitionEA* populationEA);
float testerLaSommeDesFitnessEA(partitionEA* populationEA);
void displayTheBestSolutionEA(partitionEA* solutionDominante);
void getPartitionFromSolutionEA(partitionEA *populationEA);
void checkContrainstAndFitnessPenalizationEA(partitionEA *populationEA);
void writeSolutionInFileEA(partitionEA *populationEA, FILE *outputFilePop,int iteration);
void writeBestSolutionInFileEA(partitionEA *solutionDominante, FILE *outputFile,int iteration);
//void binaryEncoding(int nbrGeneration ,FILE *outputFileEA,FILE *outputFilePopEA,FILE *outputOptimalSolutionFileEA,partitionEA *populationEA1,partitionEA *populationEA2,partitionEA *solutionDominante);
void binaryEncoding(int nbrGeneration ,FILE *outputFileEA,FILE *outputFilePopEA,FILE *outputOptimalSolutionFileEA,partitionEA *populationEA1,partitionEA *populationEA2
                    , partitionEA *solutionDominante,int iteration , int *bestSolutionIteration , int *nbrApparition);

void checkCohabitationAndNonCohabitationConstrainteEA(partitionEA *populationEA);

void getPartitionFromSolutionWithoutRepetitionEA(partitionEA *populationEA);
///void writeOptimalSolutionInFileEA(partitionEA *solutionDominante,FILE *outputOptimalSolutionFileEA);
void writeOptimalSolutionInFileEA(partitionEA *solutionDominante,FILE* outputOptimalSolutionFileEA, int nbrRun, int BestSolutionIteration, float runTime,int ES);
int compareCroissantFitnessEA (void const *a, void const *b);

#endif // BINARYENCODING_H_INCLUDED
