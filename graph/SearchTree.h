/*
 * SearchTree.h
 *
 *  Created on: Mar 29, 2014
 *      Author: sam
 */
#ifndef SEARCH_TREE_H
#define SEARCH_TREE_H
template<typename Node>
class SearchTree{

public:

/*	class Node{
		virtual BST::Node* getParent()=0;
		virtual BST::Node* getLeft()=0;
		virtual BST::Node* getRight()=0;
	};*/

	SearchTree(){

	}

	virtual ~SearchTree(){

	}



	virtual void insertAfter(Node* insertAt,Node* toInsert) =0;
	virtual void insertBefore(Node* insertAt,Node* toInsert) =0;
	virtual Node* findRoot(Node* of)=0;
	virtual Node* findMin(Node* root)=0;

	virtual Node* findMax(Node* root)=0;

	virtual void remove(Node * toRemove)=0;
	virtual Node*splitAfter(Node * splitAt)=0;
	virtual Node*splitBefore(Node * splitAt)=0;
	//All nodes in right must be (strictly) greater than all nodes in left
	virtual Node*concat(Node * left, Node*right)=0;
};


template<typename Node>
class AugmentedSearchTree:public SearchTree<Node>{

public:

	AugmentedSearchTree(){

	}

	virtual ~AugmentedSearchTree(){

	}

	//Returns the number of nodes in the subtree rooted at n (including n)
	virtual int size(Node * n)=0;
	virtual void addToIncident(Node * x, int n_incident)=0;

	virtual void setIncident(Node * x, int n_incident)=0;
	virtual  int getIncident(Node * x)=0;
	virtual int getSubtreeIncident(Node * of)=0;

	virtual void insertAfter(Node* insertAt,Node* toInsert) =0;
	virtual void insertBefore(Node* insertAt,Node* toInsert) =0;
	virtual Node* findRoot(Node* of)=0;
	virtual Node* findMin(Node* root)=0;

	virtual Node* findMax(Node* root)=0;

	virtual int depth(Node * a)=0;

	//Return -1 if a< b, 1 if a>b, 0 if either they are the same node, or not in the same tree.
	virtual int compare(Node * a, Node * b)=0;
	virtual void remove(Node * toRemove)=0;
	virtual Node*splitAfter(Node * splitAt)=0;
	virtual Node*splitBefore(Node * splitAt)=0;
	//All nodes in right must be (strictly) greater than all nodes in left
	virtual Node*concat(Node * left, Node*right)=0;
};

#endif
