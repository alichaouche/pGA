#ifndef CUTBASEDENCODING_H_INCLUDED
#define CUTBASEDENCODING_H_INCLUDED





int calculerNombreArreteEtFluxVectorEtEdgeVectorDC();
void generatePopulationDC(partitionDC* populationDC , int indiceFirstElt);
void affichePopulationDC(partitionDC* populationDC);
void calculerGenotypeDC(partitionDC * populationDC  , int indiceFirstElt);
void calculCoutCoupeEtFitnessDC(partitionDC* populationDC);
void naturalSelectionDC(partitionDC* populationDC1,partitionDC* populationDC2);
void crossOverDC(partitionDC* populationDC1, partitionDC* populationDC2);
void mutationDC(partitionDC* populationDC);

float testerLaSommeDesFitnessDC(partitionDC* populationDC);
int findTheBestSolutionDC(partitionDC *populationDC);
void displayTheBestSolutionDC(partitionDC* solutionDominante);
void writeBestSolutionInFileDC(partitionDC *solutionDominante, FILE *outputFile,int iteration);
void writeSolutionInFileDC(partitionDC *populationDC, FILE *outputFilePop,int iteration);
void writeCocyclesDeBaseInFileDC(FILE *outputFileCocyclesDeBase);
void getPartitionFromSolutionDC(partitionDC *populationDC);
void calculCoutCoupeEtFitnessWithFlowVectorDC(partitionDC* populationDC);
void checkContrainstAndFitnessPenalizationDC(partitionDC *populationDC);

//void cutBasedEncoding(int nbrGeneration ,FILE *outputFileDC,FILE *outputFilePopDC,FILE *outputOptimalSolutionFileDC,partitionDC *populationDC1,partitionDC *populationDC2,partitionDC *solutionDominante);

void cutBasedEncoding(int nbrGeneration ,FILE *outputFileDC,FILE *outputFilePopDC,FILE *outputOptimalSolutionFileDC,partitionDC *populationDC1,
                      partitionDC *populationDC2,partitionDC *solutionDominante,int iteration , int *bestSolutionIteration , int *nbrApparition);

void checkCohabitationAndNonCohabitationConstrainteDC(partitionDC *populationDC);
///void writeOptimalSolutionInFileDC(partitionDC *solutionDominante,FILE *outputOptimalSolutionFileDC);
int compareCroissantFitnessDC (void const *a, void const *b);
void writeOptimalSolutionInFileDC(partitionDC *solutionDominante,FILE* outputOptimalSolutionFileDC,
                                   int nbrRun, int bestSolutionIteration, float runTime,int ES);


#endif // CUTBASEDENCODING_H_INCLUDED
