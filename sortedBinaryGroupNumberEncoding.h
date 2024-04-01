#ifndef SORTEDCUTBASEDENCODING_H_INCLUDED


///=========================== Définitions des entêtes des fonction =======================
///void sortedBinaryGroupNumberEncoding(int nbrGeneration, FILE* outputFileVAE, FILE* outputFilePopVAE,FILE *outputOptimalSolutionFileVAE,partitionVAE *populationVAE1,partitionVAE *populationVAE2, partitionVAE *solutionDominante);
void sortedBinaryGroupNumberEncoding(int nbrGeneration, FILE* outputFileVAE, FILE* outputFilePopVAE,FILE *outputOptimalSolutionFileVAE,
                                     partitionVAE *populationVAE1,partitionVAE *populationVAE2, partitionVAE *solutionDominante,int iteration , int *bestSolutionIteration , int *nbrApparition );

void genererPopulationInitialeVAE(partitionVAE *populationVAE , int indiceFirstElt);
void generateBinarySolutionVAE(partitionVAE* populationVAE, int indiceIndividu);
void homoginiserLesSolutionBinaire(partitionVAE* populaiton, int indiceIndividu);
void conversionSolutionBinaireToEntiereVAE(partitionVAE* populationVAE,int indiceIndividu);
void trieGenotypeEntieVAE(partitionVAE* populationVAE, int indiceIndividu);
int existanceDeSolutionVAE(partitionVAE* populationVAE, int indexNewSolution);
void getPartitionFromSolutionVAE(partitionVAE *populationVAE);
void calculCoutCoupeEtFitnessVAE(partitionVAE* populationVAE);
void checkContrainstAndFitnessPenalizationVAE(partitionVAE *populationVAE);
void checkCohabitationAndNonCohabitationConstrainteVAE(partitionVAE *populationVAE);
void naturalSelectionVAE(partitionVAE* populationVAE1,partitionVAE* populationVAE2);
void crossOverVAE(partitionVAE* populationVAE1, partitionVAE* populationVAE2);
int findTheBestSolutionVAE(partitionVAE *populationVAE);
void mutationVAE(partitionVAE* populationVAE);
float testerLaSommeDesFitnessVAE(partitionVAE* populationVAE);
void displayTheBestSolutionVAE(partitionVAE* solutionDominante);
void writeSolutionInFileVAE(partitionVAE *populationVAE, FILE *outputFilePop,int iteration);
void writeBestSolutionInFileVAE(partitionVAE *solutionDominante, FILE *outputFile,int iteration);
void affichePopulationVAE(partitionVAE* populationVAE);


void genererPopulationInitialeRandomlyVAE(partitionVAE *populationVAE, int indiceFirstElt);
///void writeOptimalSolutionInFileSBGBE(partitionVAE *solutionDominante,FILE *outputOptimalSolutionFileVAE);
void writeOptimalSolutionInFileVAE(partitionVAE *solutionDominante,FILE* outputOptimalSolutionFileVAE,
                                     int nbrRun, int bestSolutionIteration, float runTime, int ES);

int compareCroissantFitnessVAE (void const *a, void const *b);

#endif // SORTEDCUTBASEDENCODING_H_INCLUDED
