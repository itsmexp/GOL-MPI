#include<stdio.h>
#include<stdlib.h>
#include <mpi.h>  
#include <pthread.h>
#include <allegro.h>
#include <time.h>

#define v(r,c) (r)*(nCols+2)+(c)

clock_t start, end;

bool allegroManualStep = false; // true: step manuale, false: step automatico
int delayAllegro = 80;			// delay per la grafica allegro
int outputFormat = 1;  	  	 // 0: terminale, 1: allegro, 2: file
bool onlyEndOutput = false;  // true: stampa solo l'output finale, false: stampa tutti gli step
int graphicCellDim = 10;	 // dimensione di una cella nella grafica allegro

int size;			// numero di processi
int rank;			// id del processo
int nRowsTot;		// numero di righe totali
int nColsTot;		// numero di colonne totali
int nPartX;			// numero di partizioni orizzontali (asse x)
int nPartY;			// numero di partizioni verticali (asse y)
int nThreads;		// numero di thread
int nStep;			// numero di step

int nRows;			// numero di righe del processo
int nCols;			// numero di colonne del processo
int * nRowsPart;	// numero di righe per partizione  [array di nPartY elementi]
int * nColsPart;	// numero di colonne per partizione [array di nPartX elementi]
int rankX;			// rank del processo lungo l'asse x
int rankY;			// rank del processo lungo l'asse y
int rankUp;			// rank del processo superiore
int rankDown;		// rank del processo inferiore
int rankLeft;		// rank del processo di sinistra
int rankRight;		// rank del processo di destra
int nPartXproc;		// numero di partizioni orizzontali per processo
int nPartYproc;		// numero di partizioni verticali per processo

int beginX;			// indice di partenza delle colonne del processo
int beginY;			// indice di partenza delle righe del processo

int* arrBeginX;		// array di indici di partenza delle colonne dei processi
int* arrBeginY; 	// array di indici di partenza delle righe dei processi
int* arrNRows;		// array di numero di righe dei processi
int* arrNCols;    	// array di numero di colonne dei processi

int * readM;		// matrice di lettura
int * writeM;		// matrice di scrittura
int * mat;           // matrice di lavoro: viene usata per l'input iniziale e dal rank 0 per la scrittura finale



pthread_t * threads;	            // array di thread
int contBarrierExe, contBarrierCom;	// contatori per la barriera
pthread_mutex_t mutexBarrierExe;    // mutex per la barriera di esecuzione
pthread_mutex_t mutexBarrierCom;    // mutex per la barriera di comunicazione
pthread_cond_t condBarrierExe;      // condition per la barriera di esecuzione
pthread_cond_t condBarrierCom;      // condition per la barriera di comunicazione

FILE * configurationFile;    // file di configurazione
FILE * inputFile;            // file di input
FILE * outputFile;           // file di output

MPI_Datatype tipoColonna;					// tipo colonna
MPI_Datatype tipoMatriceSenzaHaloBorders;	// tipo matrice senza halo cells

// inizializzo i datatypes
void initDatatypes(){
	MPI_Type_vector(nRows+2, 1, nCols+2, MPI_INT, &tipoColonna);
	// va chiamato a partire da nCols + 3
	MPI_Type_vector(nRows, nCols, nCols+2, MPI_INT, &tipoMatriceSenzaHaloBorders);
	MPI_Type_commit(&tipoColonna);
	MPI_Type_commit(&tipoMatriceSenzaHaloBorders);
}

// libero i datatypes
void finalizeDatatypes(){
	MPI_Type_free(&tipoColonna);
	MPI_Type_free(&tipoMatriceSenzaHaloBorders);
}

void initGraphics(){
    allegro_init();
	install_keyboard();
    set_color_depth(16);
    set_gfx_mode(GFX_AUTODETECT_WINDOWED, nColsTot * graphicCellDim, nRowsTot * graphicCellDim, 0, 0);
	set_window_title("Game of Life");
}


void updateGraphics(int step){
	for(int i=0;i<nRowsTot;i++){
		for(int j=0;j<nColsTot;j++){
			int x = j * graphicCellDim;
			int y = i * graphicCellDim;
			if (mat[i*nColsTot+j] == 1)
				rectfill(screen, x, y, x + graphicCellDim, y + graphicCellDim, makecol(10, 10, 10));
			else
				rectfill(screen, x, y, x + graphicCellDim, y + graphicCellDim, makecol(255, 255, 255));
		}
	}
	for(int i = 0; i < nRowsTot; i++)
		line(screen, 0, i * graphicCellDim, nColsTot * graphicCellDim, i * graphicCellDim, makecol(240, 240, 240));
	for(int i = 0; i < nColsTot; i++)
		line(screen, i * graphicCellDim, 0, i * graphicCellDim, nRowsTot * graphicCellDim, makecol(240, 240, 240));

	int offsetPartX = 0;
	int offsetPartY = 0;

	for(int i = 0; i < nPartX; i++){
		if(i%nPartXproc == 0)
			line(screen, offsetPartX, 0, offsetPartX, nRowsTot * graphicCellDim, makecol(255, 0, 0)); 
		else
			line(screen, offsetPartX, 0, offsetPartX, nRowsTot * graphicCellDim, makecol(0, 255, 0));
		offsetPartX += nColsPart[i] * graphicCellDim;
	}

	for(int i = 0; i < nPartY; i++){
		if(i%nPartYproc == 0)
			line(screen, 0, offsetPartY, nColsTot * graphicCellDim, offsetPartY, makecol(255, 0, 0));
		else
			line(screen, 0, offsetPartY, nColsTot * graphicCellDim, offsetPartY, makecol(0, 255, 0));
		offsetPartY += nRowsPart[i] * graphicCellDim;
	}
	
	//print step number on allegro screen
	char str[10];
	sprintf(str, "%d", step);
	textout_ex(screen, font, str, 5, 5, makecol(255, 0, 0), -1);
	if(allegroManualStep)
		readkey();
	else
		usleep(delayAllegro * 1000);
}
void finalizeGraphics(){
	allegro_exit();
}

void init(){
	// calcolo prima come il processo è posizionato
	int i = 1;
	nRowsPart = new int[nPartY];
	nColsPart = new int[nPartX];
	bool cont = false;
	while(i <= nThreads && !cont){
		if(nThreads % i == 0){
			if(nPartX%(nThreads/i) == 0){
				cont = true;
				nPartXproc = nThreads/i;
				nPartYproc = i;
			} else if(nPartY%(nThreads/i) == 0){
				cont = true;
				nPartXproc = i;
				nPartYproc = nThreads/i;
			}
		}
		i++;
	}
	int sizeX = nPartX/nPartXproc; // numero di processi lungo l'asse x
	int sizeY = nPartY/nPartYproc; // numero di processi lungo l'asse y
	// calcolo il numero di righe e colonne per partizione
	for(int i=0;i<nPartY;i++){
		nRowsPart[i] = nRowsTot/nPartY;
		if(i<nRowsTot%nPartY && nRowsTot>nPartY)
			nRowsPart[i]++;
	}
	for(int i=0;i<nPartX;i++){
		nColsPart[i] = nColsTot/nPartX;
		if(i<nColsTot%nPartX && nColsTot>nPartX)
			nColsPart[i]++;
	}
	// calcolo il rank del processo lungo l'asse x e y
	rankX = rank % sizeX;
	rankY = rank / sizeX;
	// calcolo il rank del processo superiore e inferiore
	rankUp = ((((rankY - 1 + sizeY) % sizeY) * sizeX) + rankX);
	rankDown = ((((rankY + 1) % sizeY) * sizeX) + rankX);
	// calcolo il rank del processo di sinistra e destra
	rankLeft = (((rankX - 1 + sizeX) % sizeX) + (rankY * sizeX));
	rankRight = (((rankX + 1) % sizeX) + (rankY * sizeX));

	// calcolo il numero di righe e colonne per processo
	nRows = 0;
	nCols = 0;
	for(int i = 0;i<nPartYproc;i++){
		nRows = nRows + nRowsPart[rankY*nPartYproc+i];
	}
	for(int i = 0;i<nPartXproc;i++){
		nCols = nCols + nColsPart[rankX*nPartXproc+i];
	}
	readM = new int[(nRows+2)*(nCols+2)];
	writeM = new int[(nRows+2)*(nCols+2)];
}

// funzione per scambiare i bordi tra i processi
void exchangeBorders(){
	// così gli scambi funzionano solo se la send non è bloccante (usa il buffer)
	// per matrici grandi conviene usare una Bsend per forzare l'uso del buffer
	if(rankUp == rank && rankDown == rank){
		if(rankLeft == rank && rankRight == rank){
			// nessun partizionamento multiprocesso
			// non devo scambiare nulla
			// copio righe e colonne, anche se sono a memoria condivisa, perché così non ho problemi
			for(int i=0;i<nRows+2;i++){
				readM[v(i,0)] = readM[v(i, nCols)];
				readM[v(i,nCols+1)] = readM[v(i, 1)];
			}
			for(int i=0;i<nCols+2;i++){
				readM[v(0, i)] = readM[v(nRows, i)];
				readM[v(nRows+1, i)] = readM[v(1, i)];
			}
		} else {
			// devo scambiare solo le colonne
			// (con un elemento in più rispetto a quando devo scambiarle insieme alle righe)
			// copio le righe
			for(int i=1;i<nRows+1;i++){
				readM[v(0, i)] = readM[v(nRows, i)];
				readM[v(nRows+1, i)] = readM[v(1, i)];
			}
			MPI_Send(&readM[v(0, 1)],1,tipoColonna,rankLeft,0,MPI_COMM_WORLD);
			MPI_Send(&readM[v(0, nCols)],1,tipoColonna,rankRight,1,MPI_COMM_WORLD);
			MPI_Recv(&readM[v(0, 0)],1,tipoColonna, rankLeft,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Recv(&readM[v(0, nCols+1)],1,tipoColonna,rankRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
	} else {
		if(rankLeft == rank && rankRight == rank){
			// devo scambiare solo le righe
			// (con un elemento in più sotto o due in più sopra
			// rispetto a quando devo scambiarle insieme alle colonne)
			// copio le colonne
			for(int i=1;i<nRows+1;i++){
				readM[v(i, 0)] = readM[v(i, nCols)];
				readM[v(i, nCols+1)] = readM[v(i, 1)];
			}
			MPI_Send(&readM[v(1,0)],nCols+2,MPI_INT, rankUp,0,MPI_COMM_WORLD);
			MPI_Recv(&readM[v(nRows+1, 0)],nCols+2,MPI_INT,rankDown,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Send(&readM[v(nRows, 0)],nCols+2,MPI_INT,rankDown,1,MPI_COMM_WORLD);
			MPI_Recv(&readM[v(0, 0)],nCols+2,MPI_INT,rankUp,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		} else {
			// devo scambiare sia le righe che le colonne
			// scambio prima la riga superiore (ricevo quella inferiore)
			MPI_Send(&readM[v(1,1)],nCols,MPI_INT,rankUp,1,MPI_COMM_WORLD);
			MPI_Recv(&readM[v(nRows+1, 1)],nCols,MPI_INT,rankDown,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			// scambio la riga inferiore (ricevo quella superiore)
			MPI_Send(&readM[v(nRows, 1)], nCols, MPI_INT, rankDown, 0, MPI_COMM_WORLD);
			MPI_Recv(&readM[v(0,1)], nCols, MPI_INT, rankUp, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			// scambio le colonne a sinistra e a destra
			MPI_Send(&readM[v(0, 1)],1,tipoColonna,rankLeft,2,MPI_COMM_WORLD);
			MPI_Send(&readM[v(0, nCols)],1,tipoColonna,rankRight,3, MPI_COMM_WORLD);
			MPI_Recv(&readM[v(0, 0)],1,tipoColonna,rankLeft,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Recv(&readM[v(0, nCols+1)],1,tipoColonna,rankRight,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
	}
}


// legge il file configuration.txt
void readConfig(){
    configurationFile = fopen("configuration.txt", "r");
    fscanf(configurationFile, "%d\n", &nPartX);
    fscanf(configurationFile, "%d\n", &nPartY);
    fscanf(configurationFile, "%d\n", &nThreads);
    fscanf(configurationFile, "%d\n", &nStep);
    fclose(configurationFile);
}

//legge la dimensione della matrice di input
void readInputDim(){
    inputFile = fopen("input.txt", "r");
	char c;
	bool end = false;
	nRowsTot = 0;
	nColsTot = 0;
	while(!end){
		c = fgetc(inputFile);
		if(c == EOF){
			end = true;
		} else {
			if(c == '\n'){
				nRowsTot++;
			} else {
				nColsTot++;
			}
		}
	}
	nColsTot = nColsTot/nRowsTot;
	fclose(inputFile);
	return;
}

// legge la sottomatrice di input che deve essere calcolata dal processo
void readInput(){
	inputFile = fopen("input.txt", "r");
	mat = new int[nRowsTot*nColsTot];
	char c;
	int cont = 0;
	bool end = false;
	while(!end){
		c = fgetc(inputFile);
		if(c == EOF)
			end = true;
		else if(c != '\n'){
				mat[cont] = c - '0';
				cont++;
		}
	}
	
	// trovo da dove deve iniziare la sottomatrice
	int nPartXAnt = rankX * nPartXproc;
	int nPartYAnt = rankY * nPartYproc;
	beginY = 0;
	beginX = 0;
	for(int i = 0; i < nPartYAnt; i++)
		beginY += nRowsPart[i];
	for(int i = 0; i < nPartXAnt; i++)
		beginX += nColsPart[i];

	// copio la sottomatrice
	for(int i = 0; i < nRows + 2; i++)
		for(int j = 0; j < nCols + 2; j++)
			if(i == 0 || i == nRows + 1 || j == 0 || j == nCols + 1)
				readM[v(i,j)] = 0;
			else
				readM[v(i,j)] = mat[(i-1 + beginY)*nColsTot+(j-1+beginX)];

	//inizializzo writeM
	for(int i = 0; i < nRows + 2; i++)
		for(int j = 0; j < nCols + 2; j++)
			writeM[v(i,j)] = 0;
	fclose(inputFile);
	if(rank != 0)
		delete [] mat;
	return;
}
// scrive la sottomatrice di input calcolata dal processo
void printReadM(){
	for(int i=0;i<nRows+2;i++){
		for(int j=0;j<nCols+2;j++){
			printf("%d ",readM[v(i,j)]);
		}
		printf("\n");
	}
	printf("\n");
}

// scrive la sottomatrice di input calcolata dal processo ma senza le halo cells
void printReadM_NoHalo(){
	for(int i=1;i<nRows+1;i++){
		for(int j=1;j<nCols+1;j++){
			printf("%d ",readM[v(i,j)]);
		}
		printf("\n");
	}
	printf("\n");
}

// scrive la sottomatrice di output calcolata dal processo
void printWriteM(){
	for(int i=0;i<nRows+2;i++){
		for(int j=0;j<nCols+2;j++){
			printf("%d ",writeM[v(i,j)]);
		}
		printf("\n");
	}
	printf("\n");
}

// scrive la sottomatrice di output calcolata dal processo ma senza le halo cells
void printWriteM_NoHalo(){
	for(int i=1;i<nRows+1;i++){
		for(int j=1;j<nCols+1;j++){
			printf("%d ",writeM[v(i,j)]);
		}
		printf("\n");
	}
	printf("\n");
}

void transitionFunction(int x, int y){
	// implemento il gioco della vita:
	// se una cella è viva e ha meno di 2 o più di 3 vicini muore
	// se una cella è viva e ha 2 o 3 vicini sopravvive
	// se una cella è morta e ha 3 vicini nasce
	int nNeighbours = 0;
	for(int i = -1; i < 2; i++)
		for(int j = -1; j < 2; j++)
			nNeighbours += readM[v(x+i,y+j)];
	nNeighbours -= readM[v(x,y)];
	if(readM[v(x,y)] == 1){
		if(nNeighbours < 2 || nNeighbours > 3)
			writeM[v(x,y)] = 0;
		else
			writeM[v(x,y)] = 1;
	} else {
		if(nNeighbours == 3)
			writeM[v(x,y)] = 1;
		else
			writeM[v(x,y)] = 0;
	}
	return;
}

void receiveMatTot(){
	for(int i = 1; i < nRows+2; i++)
		for(int j = 1; j < nCols+2; j++)
			mat[(i-1)*nColsTot+(j-1)] = writeM[v(i,j)];
	for(int r = 1; r < size; r++)
	{
		int temp[(arrNRows[r]*arrNCols[r])];
		MPI_Recv(temp, arrNRows[r]*arrNCols[r], MPI_INT, r, 25, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for(int i = 0; i < arrNRows[r]; i++)
			for(int j = 0; j < arrNCols[r]; j++)
				mat[(i+arrBeginY[r])*nColsTot+(j+arrBeginX[r])] = temp[i*arrNCols[r]+j];
	}
	return;
}

void printMatTot(int step){
	printf("Step %d\n", step);
	for(int i = 0; i < nRowsTot; i++){
		for(int j = 0; j < nColsTot; j++){
			printf("%d ", mat[i*nColsTot+j]);
		}
		printf("\n");
	}
	return;
}

void openOutputFile(){
	outputFile = fopen("output.txt", "w");
	fprintf(outputFile, "Game of Life - MPI + POSIX Threads\n");
	return;
}

void writeMatTot(int step){
	fprintf(outputFile, "\nStep %d\n", step);
	for(int i = 0; i < nRowsTot; i++){
		for(int j = 0; j < nColsTot; j++){
			fprintf(outputFile, "%d", mat[i*nColsTot+j]);
		}
		fprintf(outputFile, "\n");
	}
	return;
}

void closeOutputFile(){
	double timeElapsed = (double)(end-start)/CLOCKS_PER_SEC;
	fprintf(outputFile, "\nTEMPO DI ESECUZIONE: %fs \nProgetto APSD 2023 - Università della Calabria\n© Galardo Emanuele - Davide Pirrò - Salvatore Biamonte", timeElapsed);
	fclose(outputFile);
	return;
}



void sendSubmat(){
	MPI_Send(&writeM[v(1,1)], 1, tipoMatriceSenzaHaloBorders, 0, 25, MPI_COMM_WORLD);
	return;
}

// Riceve le informazioni necessarie per la stampa grafica
void recvCommunicationInformation(){
	arrBeginX = new int[size];
	arrBeginY = new int[size];
	arrNRows = new int[size];
	arrNCols = new int[size];
	int info[4];
	arrBeginX[0] = beginX;
	arrBeginY[0] = beginY;
	arrNRows[0] = nRows;
	arrNCols[0] = nCols;
	for(int i = 1; i < size; i++){
		MPI_Recv(info, 4, MPI_INT, i, 50, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		arrNRows[i] = info[0];
		arrNCols[i] = info[1];
		arrBeginY[i] = info[2];	
		arrBeginX[i] = info[3];
	}
}

// Invia le informazioni necessarie per la stampa grafica
void sendCommunicationInformation(){
	int info[4] = {nRows, nCols, beginY, beginX};
	MPI_Send(info, 4, MPI_INT, 0, 50, MPI_COMM_WORLD);
}

void swap(){
	int* tmp = readM;
	readM = writeM;
	writeM = tmp;
	return;
}

void initBarrier(){
    contBarrierCom = 0;
    contBarrierExe = 0;
    pthread_mutex_init(&mutexBarrierCom, NULL);
    pthread_cond_init(&condBarrierCom, NULL);
    pthread_mutex_init(&mutexBarrierExe, NULL);
    pthread_cond_init(&condBarrierExe, NULL);
}

void barrierCom(){
    pthread_mutex_lock(&mutexBarrierCom);
    contBarrierCom++;
    if(contBarrierCom == nThreads){
        contBarrierCom = 0;
        pthread_cond_broadcast(&condBarrierCom);
    } else {
        pthread_cond_wait(&condBarrierCom, &mutexBarrierCom);
    }
    pthread_mutex_unlock(&mutexBarrierCom);
}

void barrierExe(){
    pthread_mutex_lock(&mutexBarrierExe);
    contBarrierExe++;
    if(contBarrierExe == nThreads){
        contBarrierExe = 0;
        pthread_cond_broadcast(&condBarrierExe);
    } else {
        pthread_cond_wait(&condBarrierExe, &mutexBarrierExe);
    }
    pthread_mutex_unlock(&mutexBarrierExe);
}

void finalizeBarrier(){
	pthread_mutex_destroy(&mutexBarrierCom);
	pthread_cond_destroy(&condBarrierCom);
	pthread_mutex_destroy(&mutexBarrierExe);
	pthread_cond_destroy(&condBarrierExe);
}

void compute(int beginPartX, int beginPartY, int sizePartX, int sizePartY, int id){
	for(int i = beginPartX; i < beginPartX + sizePartX; i++){
		for(int j = beginPartY; j < beginPartY + sizePartY; j++){
			transitionFunction(j, i);
		}
	}
}

void * threadFunction(void * arg){
	int id = *(int *)arg;
	int xId = id % nPartXproc;
	int yId = id / nPartXproc; 
	int nPartXAnt = rankX * nPartXproc + xId;
	int nPartYAnt = rankY * nPartYproc + yId;
	int sizePartX = nColsPart[nPartXAnt];
	int sizePartY = nRowsPart[nPartYAnt];
	int beginPartX = 1;
	int beginPartY = 1;
	for(int i = rankX * nPartXproc; i < nPartXAnt; i++)
		beginPartX += nColsPart[i];
	for(int i = rankY * nPartYproc; i < nPartYAnt; i++)
		beginPartY += nRowsPart[i];
    for(int s = 0; s < nStep; s++){
        barrierCom();
		compute(beginPartX, beginPartY, sizePartX, sizePartY, id);
        barrierExe();
    }
    return NULL;
}

void initThreads(){
    threads = new pthread_t[nThreads];
    for(int i = 1; i < nThreads; i++){
        int * id = new int(i);
        pthread_create(&threads[i], NULL, threadFunction, (void *) id);
    }
}

int main(int argc, char** argv){
	start = clock();
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	readConfig();
    if(size != nPartX * nPartY / nThreads)
        exit(-1);
	readInputDim();
	init();
	readInput();
	initDatatypes();
    initBarrier();
	if(rank == 0){
		switch (outputFormat){
		case 0:
			printMatTot(0);
			break;
		case 1:
			initGraphics();
			updateGraphics(0);
			break;
		case 2:
			openOutputFile();
			writeMatTot(0);
			break;
		}
	}
	if(rank == 0){
		recvCommunicationInformation();
	}
	else{
		sendCommunicationInformation();
	}
    initThreads();
	for (int i = 0; i < nStep; i++){
		exchangeBorders();
        barrierCom();
		compute(1, 1, nColsPart[rankX * nPartXproc], nRowsPart[rankY * nPartYproc], 0);
        barrierExe();
		if(!onlyEndOutput || i == nStep - 1)
			if(rank == 0){
				receiveMatTot();
				switch(outputFormat){
					case 0:
						printMatTot(i + 1);
						break;
					case 1:
						updateGraphics(i + 1);
						break;
					case 2:
						writeMatTot(i + 1);
						break;
				}
			}
			else
				sendSubmat();
		swap();
	}
	end = clock();
	if(rank == 0 && outputFormat == 1)
		finalizeGraphics();
	if(rank == 0 && outputFormat == 2)
		closeOutputFile();
	finalizeBarrier();
	finalizeDatatypes();
	MPI_Finalize();
	return 0;
}
