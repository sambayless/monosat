
#ifndef _IBFS_H__
#define _IBFS_H__

//Source is based on that from http://www.cs.tau.ac.il/~sagihed/ibfs/, from this paper:e "A.V. Goldberg, S.Hed, H. Kaplan, R.E. Tarjan, and R.F. Werneck, Maximum Flows By Incremental Breadth-First Search"

//#define STATS


#pragma warning(disable:4786)
#include <time.h>
#include <sys/timeb.h>


template <typename captype, typename tcaptype, typename flowtype> class IBFSGraph
{
public:
	typedef enum
	{
		SOURCE	= 0,
		SINK	= 1
	} termtype;
	typedef int node_id;

	IBFSGraph(int numNodes, int numEdges);
	~IBFSGraph();
	int add_node(int numNodes);
	void add_edge(int nodeIndexFrom, int nodeIndexTo, captype capacity, captype reverseCapacity);
	void add_tweights(int nodeIndex, tcaptype capacityFromSource, tcaptype CapacityToSink);

	// to separate the graph creation and maximum flow for measurements,
	// call prepareGraph and then call maxflowClean
	// prepareGraph is only required because of the limited API for building the graph
	// (specifically - the degree of nodes is not given)
	void prepareGraph();
	flowtype maxflowClean();

	flowtype maxflow()
	{
		prepareGraph();
		return maxflowClean();
	}

	termtype what_segment(int nodeIndex, termtype default_segm = SOURCE);


private:

	struct node;
	struct arc;

	struct arc
	{
		node*		head;
		arc*		sister;
		int			sister_rCap :1;
		captype		rCap :31;
	};

	struct node
	{
		arc			*firstArc;
		arc			*parent;
		node		*nextActive;
		node		*firstSon;
		int			nextSibling;
		int			label;		// distance to the source or the sink
								// label > 0: distance from src
								// label < 0: -distance from sink
		union
		{
			tcaptype	srcSinkCap;		// srcSinkCap > 0: capacity from the source
										// srcSinkCap < 0: -capacity to the sink
			node		*nextOrphan;
		};
	};

	struct AugmentationInfo
	{
		captype remainingDeficit;
		captype remainingExcess;
		captype flowDeficit;
		captype flowExcess;
	};

	node		*nodes, *nodeLast;
	arc			*arcs, *arcLast;
	flowtype	flow;

	void augment(arc *bridge, AugmentationInfo* augInfo);
	void adoptionSrc();
	void adoptionSink();

	node* orphanFirst;
	node* orphanLast;

	int activeLevel;
	node* activeFirst0;
	node* activeFirst1;
	node* activeLast1;


public:
	int nNodes;

#ifdef STATS
	double numAugs;
	double grownSinkTree;
	double grownSourceTree;
	double numOrphans;

	double growthArcs;
	double numPushes;
	double orphanArcs1;
	double orphanArcs2;
	double orphanArcs3;
	
	double numOrphans0;
	double numOrphans1;
	double numOrphans2;
	double augLenMin;
	double augLenMax;
#endif

};





template <typename captype, typename tcaptype, typename flowtype> inline void IBFSGraph<captype, tcaptype, flowtype>::add_tweights(int nodeIndex, tcaptype capacitySource, tcaptype capacitySink)
{
	flowtype f = nodes[nodeIndex].srcSinkCap;
	if (f > 0)
	{
		capacitySource += f;
	}
	else
	{
		capacitySink -= f;
	}
	if (capacitySource < capacitySink)
	{
		flow += capacitySource;
	}
	else
	{
		flow += capacitySink;
	}
	nodes[nodeIndex].srcSinkCap = capacitySource - capacitySink;
}

template <typename captype, typename tcaptype, typename flowtype> inline void IBFSGraph<captype, tcaptype, flowtype>::add_edge(int nodeIndexFrom, int nodeIndexTo, captype capacity, captype reverseCapacity)
{
	arc *aFwd = arcLast;
	arcLast++;
	arc *aRev = arcLast;
	arcLast++;

	node* x = nodes + nodeIndexFrom;
	x->label++;
	node* y = nodes + nodeIndexTo;
	y->label++;

	aRev->sister = aFwd;
	aFwd->sister = aRev;
	aFwd->rCap = capacity;
	aRev->rCap = reverseCapacity;
	aFwd->head = y;
	aRev->head = x;
}



template <typename captype, typename tcaptype, typename flowtype> inline typename IBFSGraph<captype, tcaptype, flowtype>::termtype IBFSGraph<captype, tcaptype, flowtype>::what_segment(int nodeIndex, termtype default_segm)
{
	if (nodes[nodeIndex].parent != NULL)
	{
		if (nodes[nodeIndex].label > 0)
		{
			return SOURCE;
		}
		else
		{
			return SINK;
		}
	}
	return default_segm;
}

template <typename captype, typename tcaptype, typename flowtype> inline int IBFSGraph<captype, tcaptype, flowtype>::add_node(int numNodes)
{
	int n = nNodes;
	nNodes += numNodes;
	return n;
}
#include <limits.h>
#include <cstdlib>
#define END_OF_ORPHANS   ( (node *) 1 )
#define PREVIOUSLY_ORPHAN ( (node *) 2 )
#define END_OF_LIST_NODE ( (node *) 1 )

#define PARENT_SRC_SINK ( (arc *) 1 )



//
// Orphan handling
//

#define ADD_ORPHAN_BACK(n)							\
if ((n)->nextOrphan == NULL)						\
{													\
	(n)->parent = (n)->firstArc;					\
}													\
if (orphanFirst != END_OF_ORPHANS)					\
{													\
	orphanLast = (orphanLast->nextOrphan = (n));	\
}													\
else												\
{													\
	orphanLast = (orphanFirst = (n));				\
}													\
(n)->nextOrphan = END_OF_ORPHANS





#define ADD_ORPHAN_FRONT(n)							\
if ((n)->nextOrphan == NULL)						\
{													\
	(n)->parent = (n)->firstArc;					\
}													\
if (orphanFirst == END_OF_ORPHANS)					\
{													\
	(n)->nextOrphan = END_OF_ORPHANS;				\
	orphanLast = (orphanFirst = (n));				\
}													\
else												\
{													\
	(n)->nextOrphan = orphanFirst;					\
	orphanFirst = (n);								\
}



//
// Active list handling
//

#define SET_ACTIVE(n)								\
if ((n)->nextActive == NULL)						\
{													\
	(n)->nextActive = END_OF_LIST_NODE;				\
	if (activeFirst1 == END_OF_LIST_NODE)			\
	{												\
		activeLast1 = (activeFirst1 = (n));			\
	}												\
	else											\
	{												\
		activeLast1 = (activeLast1->nextActive = (n));	\
	}												\
}




#define NODE_PTR_TO_INDEX(p)						\
( ((p) == NULL) ? (-1) : ((int)((p)-nodes)) )



#define NODE_INDEX_TO_PTR(ind)						\
( ((ind) == -1) ? (NULL) : (nodes+ind) )









template <typename captype, typename tcaptype, typename flowtype> IBFSGraph<captype, tcaptype, flowtype>::IBFSGraph(int numNodes, int numEdges)
{
	nNodes = 0;

	nodes = (node*) malloc((numNodes+1)*sizeof(node));
	arcs = (arc*) malloc((2*numEdges)*sizeof(arc));

	assert(nodes);
	assert(arcs);

	node *maxNode = nodes + numNodes - 1;
	for (nodeLast = nodes; nodeLast <= maxNode; nodeLast++)
	{
		nodeLast->firstArc = NULL;
		nodeLast->nextActive = NULL;
		nodeLast->label = 0;
		nodeLast->srcSinkCap = 0;
	}
	arcLast = arcs;
	flow = 0;
}

template <typename captype, typename tcaptype, typename flowtype> IBFSGraph<captype, tcaptype, flowtype>::~IBFSGraph()
{
	free(nodes);
	free(arcs);
}


template <typename captype, typename tcaptype, typename flowtype> void IBFSGraph<captype, tcaptype, flowtype>::augment(arc *bridge, AugmentationInfo* augInfo)
{
	node *x, *y;
	arc *a;
	captype bottleneck, flowCurr;

#ifdef STATS
	statsNumAugs++;
	statsNumPushes++;
	int augLen=1;
#endif

	//
	// Find Bottleneck In Source Tree
	//
	bottleneck = bridge->rCap;
	if (augInfo->remainingExcess == 0)
	{
		augInfo->remainingExcess = INT_MAX;
		for (x=bridge->sister->head; ; x=a->head)
		{
	#ifdef STATS
			augLen++;
			statsNumPushes++;
	#endif
			a = x->parent;
			if (a != PARENT_SRC_SINK)
			{
				if (augInfo->remainingExcess > a->sister->rCap)
				{
					augInfo->remainingExcess = a->sister->rCap;
				}
			}
			else
			{
				break;
			}
		}
		if (augInfo->remainingExcess > x->srcSinkCap)
		{
			augInfo->remainingExcess = x->srcSinkCap;
		}
	}
	if (bottleneck > augInfo->remainingExcess)
	{
		bottleneck = augInfo->remainingExcess;
	}

	//
	// Find Bottleneck In Sink Tree
	//
	if (augInfo->remainingDeficit == 0)
	{
		augInfo->remainingDeficit = INT_MAX;
		for (x=bridge->head; ; x=a->head)
		{
	#ifdef STATS
			augLen++;
			statsNumPushes++;
	#endif
			a = x->parent;
			if (a != PARENT_SRC_SINK)
			{
				if (augInfo->remainingDeficit > a->rCap)
				{
					augInfo->remainingDeficit = a->rCap;
				}
			}
			else
			{
				break;
			}
		}
		if (augInfo->remainingDeficit > (-x->srcSinkCap))
		{
			augInfo->remainingDeficit = (-x->srcSinkCap);
		}
	}
	if (bottleneck > augInfo->remainingDeficit)
	{
		bottleneck = augInfo->remainingDeficit;
	}

#ifdef STATS
	if (augLenMin > augLen)
	{
		augLenMin = augLen;
	}
	if (augLenMax < augLen)
	{
		augLenMax = augLen;
	}
#endif

	//
	// Augment Sink Tree
	//
	augInfo->remainingDeficit -= bottleneck;
	augInfo->flowDeficit += bottleneck;
	if (augInfo->remainingDeficit == 0)
	{
		flowCurr = augInfo->flowDeficit;
		augInfo->flowDeficit = 0;
		for (x=bridge->head; ; x=a->head)
		{
			a = x->parent;
			if (a != PARENT_SRC_SINK)
			{
				a->sister->rCap += flowCurr;
				a->sister_rCap = 1;
				a->rCap -= flowCurr;
				if (a->rCap == 0)
				{
					a->sister->sister_rCap = 0;
					y=x->parent->head->firstSon;
					if (y == x)
					{
						x->parent->head->firstSon = NODE_INDEX_TO_PTR(x->nextSibling);
					}
					else
					{
						for (; y->nextSibling != (x-nodes); y = (nodes+y->nextSibling));
						y->nextSibling = x->nextSibling;
					}
					ADD_ORPHAN_FRONT(x);
				}
			}
			else
			{
				break;
			}
		}
		x->srcSinkCap += flowCurr;
		if (x->srcSinkCap == 0)
		{
			ADD_ORPHAN_FRONT(x);
		}
		if (orphanFirst != END_OF_ORPHANS)
		{
			adoptionSink();
		}
	}

	//
	// Augment Source Tree
	//
	bridge->sister->rCap += bottleneck;
	bridge->sister_rCap = 1;
	bridge->rCap -= bottleneck;
	if (bridge->rCap == 0)
	{
		bridge->sister->sister_rCap = 0;
	}
	augInfo->remainingExcess -= bottleneck;
	augInfo->flowExcess += bottleneck;
	if (augInfo->remainingExcess == 0)
	{
		flowCurr = augInfo->flowExcess;
		augInfo->flowExcess = 0;
		for (x=bridge->sister->head; ; x=a->head)
		{
			a = x->parent;
			if (a != PARENT_SRC_SINK)
			{
				a->rCap += flowCurr;
				a->sister->sister_rCap = 1;
				a->sister->rCap -= flowCurr;
				if (a->sister->rCap == 0)
				{
					a->sister_rCap = 0;
					y=x->parent->head->firstSon;
					if (y == x)
					{
						x->parent->head->firstSon = NODE_INDEX_TO_PTR(x->nextSibling);
					}
					else
					{
						for (; y->nextSibling != (x-nodes); y = (nodes+y->nextSibling));
						y->nextSibling = x->nextSibling;
					}
					ADD_ORPHAN_FRONT(x);
				}
			}
			else
			{
				break;
			}
		}
		x->srcSinkCap -= flowCurr;
		if (x->srcSinkCap == 0)
		{
			ADD_ORPHAN_FRONT(x);
		}
		if (orphanFirst != END_OF_ORPHANS)
		{
			adoptionSrc();
		}
	}

	flow += bottleneck;
}




template <typename captype, typename tcaptype, typename flowtype> void IBFSGraph<captype, tcaptype, flowtype>::adoptionSrc()
{
	node *x, *y;
	arc *a, *aEnd;
	arc aTmp;
	int minLabel;

	while (orphanFirst != END_OF_ORPHANS)
	{
		x = orphanFirst;
		orphanFirst = x->nextOrphan;

#ifdef STATS
		statsNumOrphans++;
#endif

		// we need PREVIOUSLY_ORPHAN vs NULL
		// in order to establish whether the node
		// has already started a "new parent" scan
		// while in this level or not (used in ADD_ORPHAN)
		x->nextOrphan = PREVIOUSLY_ORPHAN;
		a = x->parent;
		x->parent = NULL;
		aEnd = (x+1)->firstArc;

		// check for rehook
		if (x->label != 1)
		{
			minLabel = x->label - 1;
			for (; a != aEnd; a++)
			{
#ifdef STATS
				orphanArcs1++;
#endif
				y = a->head;
				if (a->sister_rCap != 0 &&
					y->parent != NULL &&
					y->label == minLabel)
				{
					x->parent = a;
					x->nextSibling = NODE_PTR_TO_INDEX(y->firstSon);
					y->firstSon = x;
					break;
				}
			}
		}

		// give up on node - relabel it!
		if (x->parent == NULL)
		{
			minLabel = activeLevel + 1;
			for (a=x->firstArc; a != aEnd; a++)
			{
#ifdef STATS
				orphanArcs2++;
#endif
				y = a->head;
				if (y->parent != NULL &&
					y->label > 0 &&
					y->label < minLabel &&
					a->sister_rCap != 0)
				{
					minLabel = y->label;
					x->parent = a;
					if (minLabel == x->label) break;
				}
			}

			// create orphan sons
			for (y=x->firstSon; y; y=NODE_INDEX_TO_PTR(y->nextSibling))
			{
#ifdef STATS
				orphanArcs3++;
#endif
				if (minLabel == x->label &&
					y->parent != y->firstArc)
				{
					aTmp = *(y->parent);
					*(y->parent) = *(y->firstArc);
					*(y->firstArc) = aTmp;
					y->parent->sister->sister = y->parent;
					y->firstArc->sister->sister = y->firstArc;
				}
				ADD_ORPHAN_BACK(y);
			}
			x->firstSon = NULL;
			if (x->parent == NULL)
			{
				x->nextOrphan = NULL;
			}
			else
			{
				x->label = (minLabel+1);
				x->nextSibling = NODE_PTR_TO_INDEX(x->parent->head->firstSon);
				x->parent->head->firstSon = x;
				if (minLabel == activeLevel)
				{
					SET_ACTIVE(x);
				}
			}
		}
	}
}


template <typename captype, typename tcaptype, typename flowtype> void IBFSGraph<captype, tcaptype, flowtype>::adoptionSink()
{
	node *x, *y;
	arc *a, *aEnd;
	arc aTmp;
	int minLabel;

	while (orphanFirst != END_OF_ORPHANS)
	{
		x = orphanFirst;
		orphanFirst = x->nextOrphan;

		// we need PREVIOUSLY_ORPHAN vs NULL
		// in order to establish whether the node
		// has already started a "new parent" scan
		// while in this level or not (used in ADD_ORPHAN)
		x->nextOrphan = PREVIOUSLY_ORPHAN;
		a = x->parent;
		x->parent = NULL;
		aEnd = (x+1)->firstArc;

		// check for rehook
		if (x->label != -1)
		{
			minLabel = x->label + 1;
			for (; a != aEnd; a++)
			{
#ifdef STATS
				orphanArcs1++;
#endif
				y = a->head;
				if (a->rCap != 0 &&
					y->parent != NULL &&
					y->label == minLabel)
				{
					x->parent = a;
					x->nextSibling = NODE_PTR_TO_INDEX(y->firstSon);
					y->firstSon = x;
					break;
				}
			}
		}

		// give up on node - relabel it!
		if (x->parent == NULL)
		{
			minLabel = -(activeLevel+1);
			for (a=x->firstArc; a != aEnd; a++)
			{
#ifdef STATS
				orphanArcs2++;
#endif
				y = a->head;
				if (a->rCap != 0 &&
					y->parent != NULL &&
					y->label < 0 &&
					y->label > minLabel)
				{
					minLabel = y->label;
					x->parent = a;
					if (minLabel == x->label) break;
				}
			}

			// create orphan sons
			for (y=x->firstSon; y; y=NODE_INDEX_TO_PTR(y->nextSibling))
			{
#ifdef STATS
				orphanArcs3++;
#endif
				if (minLabel == x->label &&
					y->parent != y->firstArc)
				{
					aTmp = *(y->parent);
					*(y->parent) = *(y->firstArc);
					*(y->firstArc) = aTmp;
					y->parent->sister->sister = y->parent;
					y->firstArc->sister->sister = y->firstArc;
				}
				ADD_ORPHAN_BACK(y);
			}

			x->firstSon = NULL;
			if (x->parent == NULL)
			{
				x->nextOrphan = NULL;
			}
			else
			{
				x->label = (minLabel-1);
				x->nextSibling = NODE_PTR_TO_INDEX(x->parent->head->firstSon);
				x->parent->head->firstSon = x;
				if (minLabel == -activeLevel)
				{
					SET_ACTIVE(x);
				}
			}
		}
	}
}


template <typename captype, typename tcaptype, typename flowtype> void IBFSGraph<captype, tcaptype, flowtype>::prepareGraph()
{
	node *x, *y;
	arc *a, aTmp;

	//printf("c sizeof(ptr) = %ld \n", sizeof(node*));
	//printf("c sizeof(node) = %ld \n", sizeof(node));
	//printf("c sizeof(arc) = %ld \n", sizeof(arc));
	//printf("c #nodes = %ld \n", nodeLast-nodes);
	//printf("c #arcs = %ld \n", (arcLast-arcs) + (nodeLast-nodes));
	//printf("c #grid_arcs = %ld \n", arcLast-arcs);
	//printf("c trivial_flow = %ld \n", flow);

	// calculate start arc pointers for every node
	for (x=nodes; x<nodeLast; x++)
	{
		if (x > nodes)
		{
			x->label += (x-1)->label;
		}
	}
	for (x=nodeLast; x>=nodes; x--)
	{
		if (x > nodes)
		{
			x->label = (x-1)->label;
		}
		else
		{
			x->label = 0;
		}
		x->firstArc = arcs + x->label;
	}

	// swap arcs
	for (x=nodes; x<nodeLast; x++)
	{
		for (; x->firstArc != (arcs+((x+1)->label)); x->firstArc++)
		{
			for (y = x->firstArc->sister->head; y != x; y = x->firstArc->sister->head)
			{
				// get and advance last arc fwd in proper node
				a = y->firstArc;
				y->firstArc++;

				// prepare sister pointers
				if (a->sister == x->firstArc)
				{
					x->firstArc->sister = x->firstArc;
					a->sister = a;
				}
				else
				{
					a->sister->sister = x->firstArc;
					x->firstArc->sister->sister = a;
				}

				// swap
				aTmp = (*(x->firstArc));
				(*(x->firstArc)) = (*a);
				(*a) = aTmp;
			}
		}
	}

	// reset first arc pointers
	// and sister_rCap
	for (x=nodes; x <= nodeLast; x++)
	{
		if (x != nodeLast)
		{
			x->firstArc = arcs + x->label;
			x->label = 0;
		}
		if (x != nodes)
		{
			for (a=(x-1)->firstArc; a != x->firstArc; a++)
			{
				if (a->sister->rCap == 0)
				{
					a->sister_rCap = 0;
				}
				else
				{
					a->sister_rCap = 1;
				}
			}
		}
	}

	// check consistency
#ifdef DEBUG_GRAPH
	for (i=nodes; i<nodeLast; i++)
	{
		for (a=i->firstArc; a !=(i+1)->firstArc; a++)
		{
			if (a->sister->head != i ||
				a->sister->sister != a)
			{
				exit(1);
			}
		}
	}
#endif
}

template <typename captype, typename tcaptype, typename flowtype> flowtype IBFSGraph<captype, tcaptype, flowtype>::maxflowClean()
{
	node *x, *y, *xTmp, *prevTarget;
	arc *a, *aEnd, *aTmp;
	AugmentationInfo augInfo;

#ifdef STATS
	numAugs = 0;
	numOrphans = 0;
	grownSinkTree = 0;
	grownSourceTree = 0;

	numPushes = 0;
	orphanArcs1 = 0;
	orphanArcs2 = 0;
	orphanArcs3 = 0;
	growthArcs = 0;
	augLenMin = 999999;
	augLenMax = 0;
#endif

	//
	// init
	//
	orphanFirst = END_OF_ORPHANS;
	activeFirst1 = END_OF_LIST_NODE;
	activeLevel = 1;
	for (x=nodes; x<nodeLast; x++)
	{
		x->nextActive = NULL;
		x->firstSon = NULL;

		if (x->srcSinkCap == 0)
		{
			x->parent = NULL;
		}
		else
		{
			x->parent = PARENT_SRC_SINK;
			if (x->srcSinkCap > 0)
			{
				x->label = 1;
				SET_ACTIVE(x);
			}
			else
			{
				x->label = -1;
				SET_ACTIVE(x);
			}
		}
	}
	activeFirst0 = activeFirst1;
	activeFirst1 = END_OF_LIST_NODE;

	//
	// IBFS
	//
	prevTarget = NULL;
	augInfo.flowDeficit = 0;
	augInfo.flowExcess = 0;
	while (activeFirst0 != END_OF_LIST_NODE)
	{
		//
		// BFS level
		//
		while (activeFirst0 != END_OF_LIST_NODE)
		{
			x = activeFirst0;
			activeFirst0 = x->nextActive;
			x->nextActive = NULL;
			if (x->parent == NULL)
			{
				continue;
			}

			if (x->label > 0)
			{
				//
				// Source Tree
				//
				if (x->label != activeLevel)
				{
					SET_ACTIVE(x);
					continue;
				}

#ifdef STATS
				grownSourceTree++;
#endif
				aEnd = (x+1)->firstArc;
				for (a=x->firstArc; a != aEnd; a++)
				{
#ifdef STATS
					growthArcs++;
#endif
					if (a->rCap != 0)
					{
						y = a->head;
						if (y->parent == NULL)
						{
							y->label = x->label+1;
							y->parent = a->sister;
							y->nextSibling = NODE_PTR_TO_INDEX(x->firstSon);
							x->firstSon = y;
							SET_ACTIVE(y);
						}
						else if (y->label < 0)
						{
							x->nextActive = activeFirst0;
							activeFirst0 = x;
							if (prevTarget != y)
							{
								// clear deficit
								augInfo.remainingDeficit = 0;
								if (augInfo.flowDeficit != 0)
								{
									xTmp = prevTarget;
									for (; ; xTmp=aTmp->head)
									{
										aTmp = xTmp->parent;
										if (aTmp == PARENT_SRC_SINK) break;
										aTmp->sister->rCap += augInfo.flowDeficit;
										aTmp->sister_rCap = 1;
										aTmp->rCap -= augInfo.flowDeficit;
									}
									xTmp->srcSinkCap += augInfo.flowDeficit;
									augInfo.flowDeficit = 0;
								}
							}
							augment(a, &augInfo);
							prevTarget = y;
							if (x->parent == NULL || x->label != activeLevel)
							{
								break;
							}
							activeFirst0 = activeFirst0->nextActive;
							x->nextActive = NULL;
							a = ((x->firstArc)-1);
						}
					}
				}

				// clear excess
				augInfo.remainingExcess = 0;
				if (augInfo.flowExcess != 0)
				{
					xTmp = x;
					for (; ; xTmp=aTmp->head)
					{
						aTmp = xTmp->parent;
						if (aTmp == PARENT_SRC_SINK) break;
						aTmp->rCap += augInfo.flowExcess;
						aTmp->sister->sister_rCap = 1;
						aTmp->sister->rCap -= augInfo.flowExcess;
					}
					xTmp->srcSinkCap -= augInfo.flowExcess;
					augInfo.flowExcess = 0;
				}

				// clear deficit
				augInfo.remainingDeficit = 0;
				if (augInfo.flowDeficit != 0)
				{
					xTmp = prevTarget;
					for (; ; xTmp=aTmp->head)
					{
						aTmp = xTmp->parent;
						if (aTmp == PARENT_SRC_SINK) break;
						aTmp->sister->rCap += augInfo.flowDeficit;
						aTmp->sister_rCap = 1;
						aTmp->rCap -= augInfo.flowDeficit;
					}
					xTmp->srcSinkCap += augInfo.flowDeficit;
					augInfo.flowDeficit = 0;
				}
				prevTarget = NULL;
			}
			else
			{
				//
				// GROWTH SINK
				//
				if (-(x->label) != activeLevel)
				{
					SET_ACTIVE(x);
					continue;
				}

#ifdef STATS
				grownSinkTree++;
#endif
				aEnd = (x+1)->firstArc;
				for (a=x->firstArc; a != aEnd; a++)
				{
#ifdef STATS
					growthArcs++;
#endif
					if (a->sister_rCap != 0)
					{
						y = a->head;
						if (y->parent == NULL)
						{
							y->label = x->label-1;
							y->parent = a->sister;
							y->nextSibling = NODE_PTR_TO_INDEX(x->firstSon);
							x->firstSon = y;
							SET_ACTIVE(y);
						}
						else if (y->label > 0)
						{
							x->nextActive = activeFirst0;
							activeFirst0 = x;
//							if (prevTarget != j)
//							{
//								// clear excess
//								augInfo.remainingExcess = 0;
//								if (augInfo.flowExcess != 0)
//								{
//									i_tmp = prevTarget;
//									for (; ; i_tmp=a_tmp->head)
//									{
//										a_tmp = i_tmp->parent;
//										if (a_tmp == PARENT_SRC_SINK) break;
//										a_tmp->rCap += augInfo.flowExcess;
//										a_tmp->sister->sister_rCap = 1;
//										a_tmp->sister->rCap -= augInfo.flowExcess;
//									}
//									i_tmp->srcSinkCap -= augInfo.flowExcess;
//									augInfo.flowExcess = 0;
//								}
//							}
							augment(a->sister, &augInfo);
							prevTarget = y;
							if (x->parent == NULL || x->label != activeLevel)
							{
								break;
							}
							activeFirst0 = activeFirst0->nextActive;
							x->nextActive = NULL;
							a = ((x->firstArc)-1);
						}
					}
				}

				// clear excess
				augInfo.remainingExcess = 0;
				if (augInfo.flowExcess != 0)
				{
					xTmp = prevTarget;
					for (; ; xTmp=aTmp->head)
					{
						aTmp = xTmp->parent;
						if (aTmp == PARENT_SRC_SINK) break;
						aTmp->rCap += augInfo.flowExcess;
						aTmp->sister->sister_rCap = 1;
						aTmp->sister->rCap -= augInfo.flowExcess;
					}
					xTmp->srcSinkCap -= augInfo.flowExcess;
					augInfo.flowExcess = 0;
				}
				prevTarget = NULL;

				// clear deficit
				augInfo.remainingDeficit = 0;
				if (augInfo.flowDeficit != 0)
				{
					xTmp = x;
					for (; ; xTmp=aTmp->head)
					{
						aTmp = xTmp->parent;
						if (aTmp == PARENT_SRC_SINK) break;
						aTmp->sister->rCap += augInfo.flowDeficit;
						aTmp->sister_rCap = 1;
						aTmp->rCap -= augInfo.flowDeficit;
					}
					xTmp->srcSinkCap += augInfo.flowDeficit;
					augInfo.flowDeficit = 0;
				}
			}
		}

		//
		// switch to next level
		//
		activeFirst0 = activeFirst1;
		activeFirst1 = END_OF_LIST_NODE;
		activeLevel++;
	}

	return flow;
}

#endif


