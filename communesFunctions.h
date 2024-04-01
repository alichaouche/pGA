#ifndef COMMUNESFUNCTIONS_H_INCLUDED
#define COMMUNESFUNCTIONS_H_INCLUDED

typedef unsigned long ul;

typedef struct
{
    int numeroEdge;
    int nouedDepart;
    int nouedArrive;
} edge;

double drand48ForWindows(int v1, int v2);
int rnd(int v1, int v2);
void afficherMatriceFlux();
void remplissageMatriceFlux();
void connecterGraphe();
void afficherCocyclesDeBase();
void creeCocycleDeBase();
void writeCocyclesDeBaseInFile(FILE *outputFileCocyclesDeBase);
void agencementDeGene(int t[10000]);
void mon_sleep( int nbr_seconds );
int calculSommeTotalFlux();
int calculerNombreArreteEtFluxVectorEtEdgeVector();
void themaximumDistanceMatrix();
int compareCroissant (void const *a, void const *b);
int compareTrieDecroissant(void const *a, void const *b);
void trieGenotypeEntier(mpz_t A[]);
void writeDetailsProblemInSolutionsFile(FILE* solutionsFile);
void openingFile(int numeroGraphe);
int isConnex();
void voisinEtNonVoisinArray(int voisin[],int nonVoisin[],int *v,int *nv);
void graphConnexite(int voisin[],int nonVoisin[],int v,int nv);
void cohabitationFictivesEdges();
int  maxFlux();
int SVTC_Representation(int solution[]);
void conversionVersFVTC(int solution[], partitionFVTC* populationFVTC,int indiceFVTC);
void conversionVersEA(int solution[], partitionEA* populationEA, int indiceEA);
void conversionVersPMEB(int solution[], partitionPMP* populationPMP, int indicePMEB);
void conversionVersDC(int solution[], partitionDC* populationDC, int indiceDC);
void writeBestSolutionInFilePGA(partitionFVTC *solutionDominanteFVTC, partitionEA *solutionDominanteEA,
                                partitionDC *solutionDominanteDC, partitionPMP *solutionDominantePMEB , partitionPMP *solutionDominantePMCA , FILE *outputFilePGA );


int choixDuSommetMedian(int cluster[], int tailleCluster);
#endif // COMMUNESFUNCTIONS_H_INCLUDED
