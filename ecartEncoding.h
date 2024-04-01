#ifndef SORTEDCUTBASEDENCODING_H_INCLUDED
#define SORTEDCUTBASEDENCODING_H_INCLUDED




void affichePopulationCGE(partitionCGE* populationCGE);
void calculerGenotypeCGE(partitionCGE * populationCGE);
void calculCoutCoupeEtFitnessCGE(partitionCGE* populationCGE);
void naturalSelectionCGE(partitionCGE* populationCGE1,partitionCGE* populationCGE2);
void crossOverCGE(partitionCGE* populationCGE1, partitionCGE* populationCGE2);
void mutationCGE(partitionCGE* populationCGE);

float testerLaSommeDesFitnessCGE(partitionCGE* populationCGE);
int findTheBestSolutionCGE(partitionCGE *populationCGE);
void displayTheBestSolutionCGE(partitionCGE* solutionDominante);
void writeBestSolutionInFileCGE(partitionCGE *solutionDominante, FILE *outputFile,int iteration);
void writeSolutionInFileCGE(partitionCGE *populationCGE, FILE *outputFilePop,int iteration);
void getPartitionFromSolutionCGE(partitionCGE *populationCGE);
void calculCoutCoupeEtFitnessWithFlowVectorCGE(partitionCGE* populationCGE);
void checkContrainstAndFitnessPenalizationCGE(partitionCGE *populationCGE);

///========================================================================================================================
void generatePopulationWithoutRedandancyCGE(partitionCGE* populationCGE, int indiceFirstElt);
void genererLesEcartCGE(partitionCGE* populationCGE , int indiceFirstElt);
void genererSolutionEntiereCGE(partitionCGE* populationCGE , int indiceFirstElt);
///void trieGenotypeEntier(ul A[]);
int  existanceDeSolutionCGE(int indexNewSolution, partitionCGE* populationCGE );
void genererLesSolutionBinaireCGE(partitionCGE *populationCGE, int indiceFirstElt);
///========================================================================================================================
///void ecartEncoding(int nbrGeneration,FILE *outputFileCGE,FILE *outputFilePopCGE,FILE *outputOptimalSolutionFileCGE,partitionCGE *populationCGE1, partitionCGE *populationCGE2,partitionCGE *solutionDominante);
void ecartEncoding(int nbrGeneration,FILE *outputFileCGE,FILE *outputFilePopCGE,FILE *outputOptimalSolutionFileCGE,
                   partitionCGE *populationCGE1,partitionCGE *populationCGE2,partitionCGE *solutionDominante,int iteration , int *bestSolutionIteration , int *nbrApparition);

void checkCohabitationAndNonCohabitationConstrainteCGE(partitionCGE *populationCGE);
///void writeOptimalSolutionInFileCGE(partitionCGE *solutionDominante,FILE* outputOptimalSolutionFileCGE);
int compareCroissantFitnessCGE (void const *a, void const *b);
void writeOptimalSolutionInFileCGE(partitionCGE *solutionDominante,FILE* outputOptimalSolutionFileCGE,
                                  int nbrRun, int bestSolutionIteration, float runTime,int ES);
void sortingPopulationByFintness(partitionCGE* populationCGE2, partitionCGE* tmpPopulation);


#endif // SORTEDCUTBASEDENCODING_H_INCLUDED

