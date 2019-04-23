#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "mpi.h"

char DATA_FILE_PATH[200] = "/home/nada/cpp/parallelElectionProcessManagement/";
int numberOfCandidates;
int numberOfVoters;
int rank;
int p;

int* initDynamicArray(int size){
    int* arr = malloc(size* sizeof(int));
    int i = 0;
    for(i = 0; i < size; ++i){
        arr[i] = 0;
    }
    return arr;
}

FILE* openFile(char* filePath, char* mode){
    FILE* fptr = fopen(filePath, mode);
    if (fptr == NULL)
    {
        printf("Error!! opening file");
        exit(1);
    }
    return fptr;
}

void generateData() {
    int numOfCandidates, numOfVoters;

    FILE *fptr = openFile(DATA_FILE_PATH, "w+");

    printf("Enter the number of candidates: ");
    scanf("%d", &numOfCandidates);
    fprintf(fptr, "%d\n", numOfCandidates);

    printf("Enter the number of voters: ");
    scanf("%d", &numOfVoters);
    fprintf(fptr, "%d\n", numOfVoters);

    int preference[numOfVoters][numOfCandidates];

    int i = 0;
    int j = 0;
    for (i = 0; i < numOfVoters; ++i) {
        int *selectedCandidates = initDynamicArray(numOfCandidates + 1);
        for (j = 0; j < numOfCandidates; ++j) {
            while (1) {
                int candidate = rand() % numOfCandidates + 1;
                if (selectedCandidates[candidate] == 0) {
                    selectedCandidates[candidate] = 1;
                    preference[i][j] = candidate;
                    break;
                }
            }
            fprintf(fptr, "%d ", preference[i][j]);
        }
        free(selectedCandidates);
        fprintf(fptr, "\n");
    }
    fclose(fptr);
}

void skipLines(FILE *stream, int numberOfLines) {
    char line[256];
    int i = 0;
    while (fgets(line, sizeof(line), stream)) {
        i++;
        if (i == numberOfLines) {
            break;
        }
    }
}

int* getTopTwo(const int* candidateVotes, int numberOfCandidates, int *maxVotes){
    int* winners = initDynamicArray(2);
    int * winnersVotes = initDynamicArray(2);

    int i = 0;
    for(i = 1; i < numberOfCandidates+1; ++i){
        if(candidateVotes[i] > winnersVotes[0]){
            winnersVotes[1] = winnersVotes[0];
            winners[1] = winners[0];
            winnersVotes[0] = candidateVotes[i];
            winners[0] = i;
        }else if(candidateVotes[i] > winnersVotes[1]){
            winnersVotes[1] = candidateVotes[i];
            winners[1] = i;
        }
    }
    *maxVotes = winnersVotes[0];
    free(winnersVotes);
    return winners;
}

int* runRound1ForPart(int partNum, int toSkipPortion, int toReadPortion) {
    FILE *f = openFile(DATA_FILE_PATH, "r");

    int firstLine = partNum * toSkipPortion + 2;
    skipLines(f, firstLine);

    int *candidatesVotes = initDynamicArray(numberOfCandidates + 1);

    int i = 0;
    for (i = 0; i < toReadPortion; ++i) {
        int candidate;
        fscanf(f, "%d", &candidate);
        candidatesVotes[candidate] += 1;
        skipLines(f, 1);
    }
    fclose(f);
    return candidatesVotes;

}

int* runRound2ForPart(int partNum, int toSkipPortion, int toReadPortion, const int *winners) {
    FILE *f = openFile(DATA_FILE_PATH, "r");

    int firstLine = partNum * toSkipPortion + 2;
    skipLines(f, firstLine);

    int *candidatesVotes = initDynamicArray(numberOfCandidates + 1);

    int i = 0;
    int j = 0;
    for (i = 0; i < toReadPortion; ++i) {
        for (j = 0; j < numberOfCandidates; ++j) {
            int candidate;
            fscanf(f, "%d", &candidate);
            if (candidate == winners[0] || candidate == winners[1]) {
                candidatesVotes[candidate] += 1;
                skipLines(f, 1);
                break;
            }
        }
    }
    fclose(f);
    return candidatesVotes;
}

void runRemainder(int* globalCandidateVotes, int votersPortion, const int* winners) {
    int sentVoters = votersPortion * p;
    int remVoters = numberOfVoters - sentVoters;

    int *localCandidateVotes;
    if (winners == NULL) {
        localCandidateVotes = runRound1ForPart(p, votersPortion, remVoters);
    } else {
        localCandidateVotes = runRound2ForPart(p, votersPortion, remVoters, winners);
    }

    int i = 0;
    for (i = 0; i < numberOfCandidates + 1; ++i) {
        globalCandidateVotes[i] += localCandidateVotes[i];
    }
    free(localCandidateVotes);
}

int* runRound1() {
    int votersPortion;
    if (rank == 0) {
        FILE *f = openFile(DATA_FILE_PATH, "r");
        fscanf(f, "%d", &numberOfCandidates);
        fscanf(f, "%d", &numberOfVoters);
        votersPortion = numberOfVoters / p;
        fclose(f);
    }

    MPI_Bcast(&numberOfCandidates, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&votersPortion, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int *localCandidateVotes = runRound1ForPart(rank, votersPortion, votersPortion);

    int *globalCandidateVotes = malloc((numberOfCandidates + 1) * sizeof(int));
    MPI_Reduce(localCandidateVotes, globalCandidateVotes, numberOfCandidates + 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    free(localCandidateVotes);

    if (rank == 0) {
        runRemainder(globalCandidateVotes, votersPortion, NULL);
    }
    return globalCandidateVotes;
}

int* runRound2(int rank, const int* winners) {
    int votersPortion;
    if (rank == 0) {
        votersPortion = numberOfVoters / p;
    }

    MPI_Bcast(&votersPortion, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int *localCandidateVotes = runRound2ForPart(rank, votersPortion, votersPortion, winners);

    int *globalCandidateVotes = malloc((numberOfCandidates + 1) * sizeof(int));
    MPI_Reduce(localCandidateVotes, globalCandidateVotes, numberOfCandidates + 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    free(localCandidateVotes);

    if (rank == 0) {
        runRemainder(globalCandidateVotes, votersPortion, winners);
    }
    return globalCandidateVotes;
}

void calcResult() {
    int *candidateVotes = runRound1();
    int i;
    int *winners = malloc(2 * sizeof(int));

    if (rank == 0) {
        for (i = 1; i < numberOfCandidates + 1; ++i) {
            int percentage = (int) ((candidateVotes[i] / (float) numberOfVoters) * 100);
            printf("Candidate [%d] got %d/%d which is %d%%\n", i, candidateVotes[i], numberOfVoters, percentage);
        }

        int maxVotes;
        winners = getTopTwo(candidateVotes, numberOfCandidates, &maxVotes);
        free(candidateVotes);

        if ((maxVotes / (float) numberOfVoters) >= 0.5) {
            printf("Candidate %d wins in round 1", winners[0]);
            exit(1);
        } else {
            printf("----------------------------------------------\n");
            printf("Second round will take place between candidates %d and %d\n", winners[0], winners[1]);
            printf("----------------------------------------------\n");
        }
    }
    MPI_Bcast(winners, 2, MPI_INT, 0, MPI_COMM_WORLD);

    int *finalVotes = runRound2(rank, winners);

    if (rank == 0) {
        for (i = 0; i < 2; ++i) {
            int percentage = (int) ((finalVotes[winners[i]] / (float) numberOfVoters) * 100);
            printf("Candidate [%d] got %d/%d which is %d%%\n", winners[i], finalVotes[winners[i]], numberOfVoters,
                   percentage);
        }
        if (finalVotes[winners[0]] > finalVotes[winners[1]]) {
            printf("Candidate %d wins in round 2", winners[0]);
        } else {
            printf("Candidate %d wins in round 2", winners[1]);
        }
    }
    free(finalVotes);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    if(rank == 0){
        printf("1- Create Data File\n2- Calculate Result\nEnter Your Choice: ");
        int n;
        scanf("%d", &n);
        if(n==1){
            strcat(DATA_FILE_PATH, "data.txt");
            generateData();
            printf("Data generated Successfully! \n");
            exit(0);
        }else{
            char fileName[50];
            printf("Enter Data file name: ");
            scanf("%s", fileName);
            strcat(DATA_FILE_PATH, fileName);
        }
    }
    MPI_Bcast(&DATA_FILE_PATH, 200, MPI_CHAR, 0, MPI_COMM_WORLD);
    calcResult();
    MPI_Finalize();

    return 0;
}