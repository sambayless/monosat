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
	virtual Node*split(Node * splitAt)=0;
	//All nodes in right must be (strictly) greater than all nodes in left
	virtual Node*concat(Node * left, Node*right)=0;
};

#endif
