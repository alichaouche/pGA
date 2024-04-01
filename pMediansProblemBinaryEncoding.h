#ifndef PMEDIANSPROBLEMBINARYENCODING_H_INCLUDED
#define PMEDIANSPROBLEMBINARYENCODING_H_INCLUDED



typedef struct
{
    int id;
    int genotype[1000];
    int phenotype[1000];
    int mediansVertices[1000];
    int medians;
    int nonMediansVertices[1000];
    int nonMedians;
///************************************
    int contrainteViole;
    int clustersSize[1000];
    int coutCoupeNormalise;
    int constraintVector[10];
///************************************
    int coutCoupe;
    float fitness;
///************************************
    float expectedValue;
} partitionPMP;



void generatePopulationPMP(partitionPMP* population);
void affichePopulationPMP(partitionPMP* population);
void calculerMediansNonMediansVertices(partitionPMP* population);
void affectationDesNonMedians(partitionPMP* population);
void calculeDeCutSizeEtFitnessPMP(partitionPMP* population);
void naturalSelectionPMP(partitionPMP* population1,partitionPMP* population2);
void crossOverPMP(partitionPMP* population1, partitionPMP* population2);
int findTheBestSolutionPMP(partitionPMP *population);
void mutationPMP(partitionPMP* population);
float testerLaSommeDesFitnessPMP(partitionPMP* population);
void displayTheBestSolutionPMP(partitionPMP* solutionDominante);
void checkContrainstAndFitnessPenalizationPMP(partitionPMP *population);
void writeSolutionInFilePMP(partitionPMP *population, FILE *outputFilePop,int iteration);
void writeBestSolutionInFilePMP(partitionPMP *solutionDominante, FILE *outputFile,int iteration);
void pMedianEncoding(int nbrGeneration,FILE *outputFilePMP,FILE *outputFilePopPMP,FILE *outputOptimalSolutionFilePMP);
void checkContrainstAndFitnessPenalizationPMP(partitionPMP *population);
void writeOptimalSolutionInFilePMP(partitionPMP *solutionDominante,FILE *outputOptimalSolutionFilePMP);
int compareCroissantFitnessPMP (void const *a, void const *b);


#endif // PMEDIANSPROBLEMBINARYENCODING_H_INCLUDED
