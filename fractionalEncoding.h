#ifndef FRACTIONALENCODING_H_INCLUDED
#define FRACTIONALENCODING_H_INCLUDED


void generatePopulationFVTC(partitionFVTC* populationFVTC,int indiceFirstElt);
void affichePopulationFVTC(partitionFVTC* populationFVTC);
void calculCoutCoupeEtFitnessFVTC(partitionFVTC* populationFVTC, int indiceFirst, int indiceLast);
void naturalSelectionFVTC(partitionFVTC* populationFVTC1,partitionFVTC* populationFVTC2);
void crossOverFVTC(partitionFVTC* populationFVTC1, partitionFVTC* populationFVTC2);
void calculPhenotypeFVTC(partitionFVTC* populationFVTC);
int findTheBestSolutionFVTC(partitionFVTC *populationFVTC);
void mutationFVTC(partitionFVTC* populationFVTC);
float testerLaSommeDesFitnessFVTC(partitionFVTC* populationFVTC);
void displayTheBestSolutionFVTC(partitionFVTC* solutionDominante);
void checkContrainstAndFitnessPenalizationFVTC(partitionFVTC *populationFVTC);
void writeSolutionInFileFVTC(partitionFVTC *populationFVTC,FILE *outputFilePop,int iteration);
void writeBestSolutionInFileFVTC(partitionFVTC *bestSolution, FILE *outputFile,int iteration);
///void fractionalEncoding(int nbrGeneration,FILE *outputFileFVTC,FILE *outputFilePopFVTC,FILE *outputOptimalSolutionFileFVTC,partitionFVTC *populationFVTC1,partitionFVTC *populationFVTC2,partitionFVTC *solutionDominante );
void fractionalEncoding(int nbrGeneration,FILE *outputFileFVTC,FILE *outputFilePopFVTC,FILE *outputOptimalSolutionFileFVTC,partitionFVTC *populationFVTC1,
                                        partitionFVTC *populationFVTC2,partitionFVTC *solutionDominante,int iteration , int *bestSolutionIteration , int *nbrApparition);

void checkCohabitationAndNonCohabitationConstrainteFVTC(partitionFVTC *populationFVTC);
///void writeOptimalSolutionInFileFVTC(partitionFVTC *solutionDominante,FILE *outputOptimalSolutionFileFVTC);
void writeOptimalSolutionInFileFVTC(partitionFVTC *solutionDominante,FILE* outputOptimalSolutionFileFVTC,
                                  int nbrRun, int bestSolutionIteration, float runTime,int ES);
int compareCroissantFitnessFVTC (void const *a, void const *b);
double drand48ForWindowsFVTC(int v1,int v2);



#endif // FRACTIONALENCODING_H_INCLUDED
