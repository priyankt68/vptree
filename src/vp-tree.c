/*
 * vptree.c
 *
 * VP-tree
 *
 * - 1) Dataset reading
 * - 2) Building VP-Tree
 * - 3) Read/Write VP-Tree from/to a file
 * - 4) K-NNSearch
 * - 5) Releasing Resource
 *
 * By Philip Fu
 *
 * Fri Jan 12 22:58:44 EST 2001
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <assert.h>

#include "vptree.h"
#include "qsort.h"
#include "standard.h"




////////////////////////////////////////////////////////////////
// Function Definitions
////////////////////////////////////////////////////////////////


//#define _DEBUG




/////////////////////////////////////////////////////////
// (1) Reading the dataset
/////////////////////////////////////////////////////////


DataSet *
readDataSet(int n, int dim, char *dataFile)
{
    FILE   *fp;
    double *ptr,f;

    DataSet *dataSet;
    char     buf[255];
    int      i;



    /////////////////////////////////////////////
    // (1) Open the dataFile

    fp = fopen(dataFile,"r");

    if (fp == NULL) {
	printf("Err : Cannot open file [%s] for reading!\n",dataFile);
	return NULL;
    }


    /////////////////////////////////////////////
    // (2) Initialize dataSet

    dataSet = (DataSet *) malloc(sizeof(DataSet));

    if (dataSet == NULL) {
	printf("Err : not enough memory for data allocation!\n\n");
	return NULL;
    }


    /////////////////////////////////////////////
    // (3) Read the header

    // read number of points
    // fgets(buf,255,fp);
    // sscanf(buf,"%d",&i);
    dataSet->nPoints = n;

    // read dimension
    // fgets(buf,255,fp);
    // sscanf(buf,"%d",&i);
    dataSet->dimension = dim;

    //printf("number of points : %d\n",dataSet->nPoints);
    //printf("dimension        : %d\n",dataSet->dimension);


    /////////////////////////////////////////////
    // (4) Read the points

    dataSet->points = (double *) 
	malloc(sizeof(double)*(dataSet->nPoints)*(dataSet->dimension));

    if (dataSet->points == NULL) {
	printf("Err : not enough memory for data allocation!\n\n");
	return NULL;
    }

    ptr = dataSet->points;

    for (i=0; i<(dataSet->nPoints)*(dataSet->dimension); i++) {
	fscanf(fp,"%lf",&f);
	*ptr++ = f;
    }


    /////////////////////////////////////////////
    // (5) Done

    fclose(fp);

    return dataSet;
}




/////////////////////////////////////////////////////////
// (2) Building the VP Tree
/////////////////////////////////////////////////////////


double
sqrDist(double *p1, double *p2, int dim)
{
    double dist;
    int    i;

    dist = pow( ((*p1++) - (*p2++)), 2.0 );
    for (i=1; i<dim; i++)
	dist += pow( ((*p1++) - (*p2++)), 2.0 );
    return dist;
}



// Building a VP Node

VPNode *
buildVPNode(double *points, int *dataID, int nPoints, 
            int numBranch, int dimension)
{
    VPNode *vpNode;
    int     i,j;

    double *distances;

    int  nSubset;
    int *subset;

    double stddev,bestStddev,dist;

    int    vp,bestVP;
    double *ptrVP, *ptrPt, *ptrDist;

    int    *sortList;
    double *tmpPt;
    int    *tmpID, *ptrID;

    int lastIndex,childIndex;


    #ifdef _DEBUG
    printf("\nEnter buildVPNode (nPoints : %d) and (address %d)\n",
	nPoints,points);
    printf("dataID : %d\n",dataID[0]);
    #endif


    /////////////////////////////////////////////
    // (1) Who am I?

    vpNode = (VPNode *) malloc(sizeof(VPNode));

    if (vpNode == NULL)
	errexit("Err : Not enough memory for allocation!\n\n");



    /////////////////////////////////////////////
    // (2) Fill in information

    // avoid leakage when free()
    vpNode->child   = NULL;
    vpNode->medians = NULL;
    vpNode->points  = NULL;
    vpNode->dataID  = NULL;

    // should we stop?
    if (nPoints <= _MAX_POINTS_NODE) {

	vpNode->isLeaf = _TRUE;

	// alloc. data points in leaf node
	vpNode->nPoints = nPoints;
	vpNode->points  = (double *) malloc(sizeof(double)*dimension*nPoints);
	vpNode->dataID  = (int *)    malloc(sizeof(int)*nPoints);

	if (!vpNode->points || !vpNode->dataID)
	    errexit("Err : Not enough memory for allocation!\n\n");

	memcpy(vpNode->points,
	       points,
	       sizeof(double)*dimension*nPoints);

	memcpy(vpNode->dataID,
	       dataID,
	       sizeof(int)*nPoints);

	return vpNode;

    } else {

	vpNode->isLeaf = _FALSE;

	// alloc. vp
	vpNode->nPoints = 1;
	vpNode->points  = (double *) malloc(sizeof(double)*dimension);
	vpNode->dataID  = (int *)    malloc(sizeof(int));

	if (!vpNode->points || !vpNode->dataID)
	    errexit("Err : Not enough memory for allocation!\n\n");

	// * vp to be found later
    }


    // Initialize distances
    distances = (double *) malloc(sizeof(double)*(nPoints-1));

    if (!distances)
	errexit("Err : Not enough memory for allocation!\n\n");



    /////////////////////////////////////////////
    // (3) Find the vantage point


    // first find a random "subset" in the dataset

    if (nPoints <= _NUM_PT_RANDOM_SUBSET) {

	// no need to find a random subset

	nSubset = nPoints;
	subset  = (int *) malloc(sizeof(int)*nSubset);

	if (subset == NULL)
	    errexit("Err : Not enough memory for allocation!\n\n");

	for (i=0; i < nSubset; i++)
	    subset[i] = i;

    } else {

	// find a random subset

	nSubset = _NUM_PT_RANDOM_SUBSET;
	subset  = (int *) malloc(sizeof(int)*nSubset);

	if (subset == NULL)
	    errexit("Err : Not enough memory for allocation!\n\n");

	srand(time(NULL));

	for (i=0; i < nSubset; i++)
	    subset[i] = rand() % nPoints;
    }


    // Randomize some candidate vps and the one 
    // with the largest standard dev. will be the VP

    srand(time(NULL));

    // initialize
    bestStddev = 0.0;
    bestVP     = -1;

    for (j=0; j  < _NUM_CAND_VP; j++) {

	// randomize a candidate vp
	vp    = rand() % nPoints;
	ptrVP = points + vp*dimension;

	// measure the standard dev. against this vp
	dist = 0.0;
	for (i=0; i  < nSubset; i++)
	    dist += sqrDist( ptrVP, 
	                     points + subset[i]*dimension, 
	                     dimension );
	stddev = sqrt(dist/nSubset);

	// largest?
	if (bestStddev  < stddev) {
	    bestStddev = stddev;
	    bestVP     = vp;
	}
    }

    free(subset);

    #ifdef _DEBUG
    printf("bestVP : %d at (%f %f %f) with stdev %f\n",
	bestVP,ptrVP[0],ptrVP[1],ptrVP[2],bestStddev);
    #endif



    /////////////////////////////////////////////
    // (4) Make the VP at the 1st
    //     position by swapping


    // copy from (bestVP)th to vpNode
    memcpy(vpNode->points,
           points + bestVP*dimension,
           sizeof(double)*dimension);
    vpNode->dataID[0] = dataID[bestVP];

    // copy from 0th element to (bestVP)th
    memcpy(points + bestVP*dimension,
           points,
           sizeof(double)*dimension);
    dataID[bestVP] = dataID[0];

    // copy from vpNode to 0th element
    memcpy(points,
           vpNode->points,
           sizeof(double)*dimension);
    dataID[0] = vpNode->dataID[0];



    /////////////////////////////////////////////
    // (5) Find the distance of points to the vp

    ptrPt   = points + dimension;	// skip the 1st vp
    ptrDist = distances;

    for (i=0; i< nPoints-1; i++) {

	*ptrDist++  = sqrt( sqrDist(points, ptrPt, dimension) );
	 ptrPt     += dimension;

	#ifdef _DEBUG
	printf("cal dist. %f (%f %f %f)\n",*(ptrDist-1),
		*(ptrPt-3),*(ptrPt-2),*(ptrPt-1));
	#endif
    }



    /////////////////////////////////////////////
    // (6) Sort the points in the array in the
    //     ascending order of their corr. distances


    // initialize

    sortList = (int *) malloc(sizeof(int)*(nPoints-1));

    if (sortList == NULL)
	errexit("Err : Not enough memory for allocation!\n\n");

    for (i=0; i<nPoints-1; i++)
	sortList[i] = i;

    // sorting (find the sort index)
    //bubbleSortByIndex(&sortList, distances, nPoints-1);
    q_sort(sortList,distances,0,nPoints-2,nPoints-1);

    // rearrange the points -> tmpPt and dataID -> tmpID
    tmpPt = (double *) malloc(sizeof(double)*(nPoints-1)*dimension);
    tmpID = (int *)    malloc(sizeof(int)*(nPoints-1));

    if (!tmpPt || !tmpID)
	errexit("Err : Not enough memory for allocation!\n\n");

    ptrPt = tmpPt;
    ptrID = tmpID;

    for (i=0; i<nPoints-1; i++) {

	memcpy( ptrPt,
	        points + (sortList[i]+1)*dimension,	// +1 to skip 1st vp
	        sizeof(double)*dimension );

	*ptrID++ = dataID[sortList[i]+1];

	ptrPt += dimension;
    }

    // update points and dataID
    memcpy(points+dimension,
           tmpPt,
           sizeof(double)*(nPoints-1)*dimension);
    memcpy(dataID+1,
           tmpID,
           sizeof(int)*(nPoints-1));

    free(tmpPt);
    free(tmpID);

    #ifdef _DEBUG
    for (i=0; i<nPoints-1; i++)
	printf("dist[%02d] = %f\n",i,distances[i]);
    for (i=0; i<nPoints-1; i++)
	printf("sort[%02d] = %d\n",i,sortList[i]);
    #endif



    /////////////////////////////////////////////
    // (7) Fill in medians


    // initialize medians

    vpNode->medians = (double *) malloc(sizeof(double)*(numBranch-1));

    if (!vpNode->medians)
	errexit("Err : Not enough memory for allocation!\n\n");

    childIndex = 1;

    for (i=0; i<numBranch-1; i++) {

	// last point in this branch
	lastIndex = 1 + (nPoints-1)*(i+1) / numBranch - 1;
	if (i == numBranch-1)
	    lastIndex = nPoints - 1;

	// set median
	vpNode->medians[i] = distances[lastIndex-1];

	// next branch
	childIndex = lastIndex+1;

	//printf("vpNode->medians[%d] = %f\n",i,vpNode->medians[i]);
    }

    free(distances);
    free(sortList);



    /////////////////////////////////////////////
    // (8) Build VPNode of each child


    // initialize childs

    vpNode->child = 
	(VPNode **) malloc(sizeof(VPNode *)*numBranch);

    if (! vpNode->child)
	errexit("Err : Not enough memory for allocation!\n\n");

    childIndex = 1;

    for (i=0; i<numBranch; i++) {

	// last point in this branch
	lastIndex = 1 + (nPoints-1)*(i+1) / numBranch - 1;
	if (i == numBranch-1)
	    lastIndex = nPoints - 1;

        #ifdef _DEBUG
	printf("from %d to %d\n",childIndex,lastIndex);
	#endif

	// set child
	vpNode->child[i] = buildVPNode(points+childIndex*dimension,
	                               dataID+childIndex,
	                               lastIndex-childIndex+1,
	                               numBranch,
	                               dimension);

	#ifdef _DEBUG
	if (i != numBranch-1)
	    printf("median : %f\n",vpNode->medians[i]);
	#endif

	// next branch
	childIndex = lastIndex+1;
    }


    return vpNode;
}



// External Function :
//
// Building the VP Tree
// - if numBranch    <= 1, use _DEFAULT_BRANCH
// - note that ordering of points in dataSet is changed after calling
//   (return the root node)

VPTree *
buildVPTree(DataSet *dataSet, int numBranch)
{
    VPTree *vpTree;
    int    *dataID, i;


    ///////////////////////////////////////////
    // Allocate memory

    vpTree = (VPTree *) malloc(sizeof(VPTree));
    dataID = (int *)    malloc(sizeof(int)*dataSet->nPoints);

    if (!vpTree || !dataID)
	errexit("Err : Not enough memory for allocation!\n\n");


    ///////////////////////////////////////////
    // (2) Initialization

    if (numBranch <= 1)
	numBranch  = _DEFAULT_BRANCH;

    vpTree->numBranch = numBranch;
    vpTree->dimension = dataSet->dimension;

    // dataIDs' keep track of where the points goes in the VP tree
    for (i=0; i<dataSet->nPoints; i++)
	dataID[i] = i;


    ///////////////////////////////////////////
    // (3) Build the VP-Tree recursively

    vpTree->root = buildVPNode(dataSet->points,dataID,dataSet->nPoints,
                               numBranch,dataSet->dimension);


    return vpTree;
}




/////////////////////////////////////////////////////////
// (3) Read/Write VPTree file
/////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////
// (3a) read the VP-Tree from a file


// Internal function : read a VP-node

VPNode *
readVPNode(FILE *fp, int numBranch, int dimension, int *nPoints)
{
    VPNode *vpNode;
    char    buf[255],dummy[255];
    int     i,j,k;

    double *ptrPt, f;
    int    *ptrID;


    /////////////////////////////////
    // (1) Initialization

    vpNode = (VPNode *) malloc(sizeof(VPNode));

    if (!vpNode)
	errexit("Err : Not enough memory for allocation!\n\n");

    // for safety
    vpNode->child   = NULL;
    vpNode->medians = NULL;
    vpNode->points  = NULL;
    vpNode->dataID  = NULL;


    /////////////////////////////////
    // (2) What kind of node?

    fgets(buf,255,fp);

    if (!strncmp(buf,_VPTREE_LEAF,strlen(_VPTREE_LEAF))) {

	// Leaf Node
	vpNode->isLeaf = _TRUE;

	// How many data points?
	sscanf(buf,"%s %d",dummy,nPoints);
	vpNode->nPoints = *nPoints;

	// memory allocation
	vpNode->points = (double *) malloc(sizeof(double)*(*nPoints)*dimension);
	vpNode->dataID = (int *)    malloc(sizeof(int)*(*nPoints));

	if (!vpNode->points || !vpNode->dataID)
	    errexit("Err : Not enough memory for allocation!\n\n");

	// for fast reference
	ptrPt = vpNode->points;
	ptrID = vpNode->dataID;

	// Read the points and dataID
	for (i=0; i<*nPoints; i++) {
	    for (j=0; j<dimension; j++) {
		fscanf(fp,"%lf",&f);
		*ptrPt++ = f;
	    }
	    fscanf(fp,"%d",ptrID++);

	    //printf("%f %f %f (%d)\n",
		//*(ptrPt-3),*(ptrPt-2),*(ptrPt-1),*(ptrID-1));

	    // clear this line in fp
	    fgets(buf,255,fp);
	}

    } else
    if (!strncmp(buf,_VPTREE_INTERNAL,strlen(_VPTREE_INTERNAL))) {

	// Internal Node
	vpNode->isLeaf = _FALSE;

	// Only one Vantage point
	vpNode->nPoints = *nPoints = 1;

	// memory allocation
	vpNode->points  = (double *)  malloc(sizeof(double)*dimension);
	vpNode->dataID  = (int *)     malloc(sizeof(int));
	vpNode->medians = (double *)  malloc(sizeof(double)*(numBranch-1));
	vpNode->child   = (VPNode **) malloc(sizeof(VPNode *)*numBranch);

	if ( !vpNode->points  || !vpNode->dataID
	  || !vpNode->medians || !vpNode->child )
	    errexit("Err : Not enough memory for allocation!\n\n");

	// for fast reference
	ptrPt = vpNode->points;

	// Read the vp and its dataID
	for (j=0; j<dimension; j++) {
	    fscanf(fp,"%lf",&f);
	    *ptrPt++ = f;
	}
	fscanf(fp,"%d",vpNode->dataID);

	// clear this line in fp
	fgets(buf,255,fp);

	//printf("%f %f %f (%d)\n",
		//*(ptrPt-3),*(ptrPt-2),*(ptrPt-1),vpNode->dataID[0]);

	// Read the children and medians
	for (i=0; i<numBranch; i++) {

	    // read child
	    vpNode->child[i] = readVPNode(fp,numBranch,dimension,&k);
	    *nPoints += k;

	    // read medians
	    if (i != numBranch-1) {
		fgets(buf,255,fp);
		sscanf(buf,"%s %s %lf",dummy,dummy,&f);
		vpNode->medians[i] = f;
	    }
	}

    } else {

	printf("UNKNOWN NODE TYPE [%s]\n",buf);

	return NULL;
    }

    return vpNode;
}


// on error, return NULL

VPTree *
readVPTree(char *treeFile, int *nPoints)
{
    VPTree *vpTree;
    FILE   *fp;
    char    buf[255],dummy[255];



    /////////////////////////////////
    // (1) open file

    fp = fopen(treeFile,"r");

    if (!fp) {
	fprintf(stderr,"ERR : cannot open [%s] for reading!\n\n",treeFile);
	return NULL;
    }

    vpTree = (VPTree *) malloc(sizeof(VPTree));

    if (!vpTree)
	errexit("Err : Not enough memory for allocation!\n\n");


    /////////////////////////////////
    // (2) write header

    fgets(buf,255,fp);
    sscanf(buf,"%s",dummy);

    if (strcmp(dummy,_VPTREE_ID)) {
	fprintf(stderr,"Err : [%s] is not a vptree file\n\n");
	return NULL;
    }

    fgets(buf,255,fp);
    sscanf(buf,"%s %s %d",dummy,dummy,&(vpTree->numBranch));
    
    fgets(buf,255,fp);
    sscanf(buf,"%s %s %d",dummy,dummy,&(vpTree->dimension));

    //printf("Branching : %d\n",vpTree->numBranch);
    //printf("Dimension : %d\n",vpTree->dimension);


    /////////////////////////////////
    // (3) read vpNode recursively

    vpTree->root =
	readVPNode(fp,vpTree->numBranch,vpTree->dimension,nPoints);


    /////////////////////////////////
    // (4) Done

    fclose(fp);


    return vpTree;
}



/////////////////////////////////////////////////////////
// (3b) write the VP-Tree to a file


// Internal function

int
writeVPNode(FILE *fp, VPNode *vpNode, int numBranch, int dimension)
{
    int     i,j,nPoints,*ptrID;
    double *ptrPt;


    if (vpNode->isLeaf) {

	////////////////////////
	// leaf node
	////////////////////////

	// initialization
	nPoints = vpNode->nPoints;
	ptrPt   = vpNode->points;
	ptrID   = vpNode->dataID;

	// identifier
	fprintf(fp,"%s %d\n",_VPTREE_LEAF,nPoints);

	// data points
	for (i=0; i<nPoints; i++) {
	    for (j=0; j<dimension; j++)
		fprintf(fp,"%f ",*ptrPt++);
	    fprintf(fp,"%d\n",*ptrID++);
	}

    } else {

	////////////////////////
	// Internal node
	////////////////////////

	// initialization
	nPoints = 1;
	ptrPt   = vpNode->points;

	// identifier
	fprintf(fp,"%s\n",_VPTREE_INTERNAL);

	// vantage point(s)
	for (i=0; i<dimension; i++)
	    fprintf(fp,"%f ",*ptrPt++);
	fprintf(fp,"%d\n",vpNode->dataID[0]);

	for (i=0; i<numBranch; i++) {
	    nPoints += writeVPNode(fp,vpNode->child[i],numBranch,dimension);
	    if (i != numBranch-1)
		fprintf(fp,"median : %f\n",vpNode->medians[i]);
	}
    }

    return nPoints;
}


// External function
// on error, return -1, else number of data points written

int
writeVPTree(char *treeFile, VPTree *vpTree)
{
    int   nPoints;
    FILE *fp;


    /////////////////////////////////
    // (1) open file

    fp = fopen(treeFile,"w");

    if (!fp) {
	fprintf(stderr,"ERR : cannot open [%s] for writing!\n\n",treeFile);
	return -1;
    }


    /////////////////////////////////
    // (2) write header

    fprintf(fp,"%s\n",_VPTREE_ID);
    fprintf(fp,"Branching : %d\n",vpTree->numBranch);
    fprintf(fp,"Dimension : %d\n",vpTree->dimension);


    /////////////////////////////////
    // (3) write vpNode recursively

    nPoints = writeVPNode( fp, vpTree->root,
                           vpTree->numBranch, vpTree->dimension );


    /////////////////////////////////
    // (4) Done

    fclose(fp);


    return nPoints;
}




/////////////////////////////////////////////////////////
// (4) K-Nearest Neighbor Search
/////////////////////////////////////////////////////////

int
searchKNNVPNode(VPNode *vpNode, double *queryPt,
                double *kMinDist, double *resultPt, int *resultID,
                int numBranch, int dimension, int k)
{
    double dist, *ptrPt;
    int    i, x, y;
    int    index, nodeCount;

    // visited me
    nodeCount = 1;

    // what kind of node?
    if (vpNode->isLeaf) {

	ptrPt = vpNode->points;

	for (i=0; i<vpNode->nPoints; i++) {

	    // compute distance
	    dist = sqrt( sqrDist(queryPt,ptrPt,dimension) );

	    // shorter?
	    if (kMinDist[k-1] > dist) {
	     	for (x=0; x<k-1; x++)
			if (kMinDist[x] > dist)
				break;
	     	index = x;
	     	for (x=k-1; x>index; x--) {
			kMinDist[x] = kMinDist[x-1];	
			resultID[x] = resultID[x-1];
			for (y=0; y<dimension; y++)
				resultPt[x*dimension+y] = resultPt[(x-1)*dimension+y];
	     	}

	     	kMinDist[index] = dist;
		memcpy(&(resultPt[index*dimension]),ptrPt,sizeof(double)*dimension);
	     	resultID[index] = vpNode->dataID[i];
	    }

	    // next
	    ptrPt += dimension;
	}

    } else {

	// compute distance
	dist = sqrt( sqrDist(queryPt,vpNode->points,dimension) );

	// shorter?
	if (kMinDist[k-1] > dist) {
	     for (i=0; i<k-1; i++)
		if (kMinDist[i] > dist)
			break;
	     index = i;
	     for (i=k-1; i>index; i--) {
		kMinDist[i] = kMinDist[i-1];	
		resultID[i] = resultID[i-1];
		for (x=0; x<dimension; x++)
			resultPt[i*dimension+x] = resultPt[(i-1)*dimension+x];
	     }

	     kMinDist[index] = dist;
	     memcpy(&(resultPt[index*dimension]),vpNode->points,sizeof(double)*dimension);
	     resultID[index] = vpNode->dataID[0];
	}

	// Prune the near part
	for (i=0; i<numBranch-1; i++)
	    if (kMinDist[k-1] + vpNode->medians[i] >= dist)
		break;

	// Prune the far part
	for (; i<numBranch; i++) {

	    nodeCount += 
		searchKNNVPNode( vpNode->child[i], queryPt,
		                 kMinDist, resultPt, resultID,
		                 numBranch, dimension, k );

	    if ( i != numBranch-1
	      && kMinDist[k-1] + dist < vpNode->medians[i] )
		break;
	}
    }

    return nodeCount;
}


// Perform K-Nearest Neighbor Search
// Input :
// - vpTree
// - queryPt
// - k
// Return :
// - number of nodes visited (int)
// - list of points -> resultPt is returned in argument
// - list of dataID -> resultID is returned in argument
// Note : - memory for resultID and resultPt should have been 
//        - allocated before calling.

int
knnsearch(VPTree *vpTree, double *queryPt, int k, 
          double *resultPt, int *resultID)
{
    double *kMinDist;
    int i,nodeCount;

    // Assertion
    assert(resultID);
    assert(resultPt);

    // initialize distances and nodeCount
    kMinDist = (double *) malloc(sizeof(double) * k);

    if (kMinDist == NULL) {
	printf("Err : not enough memory!\n\n");
	return;
    }

    for (i=0; i < k; i++)
	kMinDist[i] = _LARGE_VALUE;


    // Start searching
    nodeCount = searchKNNVPNode( vpTree->root, queryPt,
                                 kMinDist, resultPt, resultID,
                                 vpTree->numBranch, vpTree->dimension, k );

    return nodeCount;
}




/////////////////////////////////////////////////////////
// (5) Free the resource
/////////////////////////////////////////////////////////


void
freeDataSet(DataSet *dataSet)
{
    if (dataSet) {
	if (dataSet->points)
	    free(dataSet->points);
	free(dataSet);
    }
}


void
freeVPNode(VPNode *vpNode, int numBranch)
{
    int i;

    if (vpNode) {

	// free the child
	if (!vpNode->isLeaf && vpNode->child) {
	    for (i=0; i<numBranch; i++)
		if (vpNode->child[i])
		    freeVPNode(vpNode->child[i],numBranch);
	}

	// free myself
	if (vpNode->child)
	    free(vpNode->child);
	if (vpNode->medians)
	    free(vpNode->medians);
	if (vpNode->points)
	    free(vpNode->points);
	if (vpNode->dataID)
	    free(vpNode->dataID);

	free(vpNode);
    }
}


void
freeVPTree(VPTree *vpTree)
{
    if (vpTree) {

	// free the root
	if (vpTree->root)
	    freeVPNode(vpTree->root,vpTree->numBranch);

	// free myself
	free(vpTree);
    }
}



// Print point

void
printPt(FILE *out, double *pt, int dim)
{
    int i;

    fprintf(out,"%f",pt[0]);
    for (i=1; i<dim; i++)
	fprintf(out," %f",pt[i]);
}
