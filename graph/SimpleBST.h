/*
 * BST.h
 *
 *  Created on: Mar 25, 2014
 *      Author: sam
 */

#include <cassert>
#include "BST.h"
template<typename Value>
class SimpleBST:public BST<Value>{

	int key;
	Value & v;

	SimpleBST*left;
	SimpleBST*right;
	SimpleBST*parent;
public:
	virtual SimpleBST(Value & v, int key):v(v),key(key),left(nullptr),right(nullptr),parent(nullptr){

	}
	virtual ~SimpleBST(){

	}
	virtual SimpleBST*getLeft(){
		return left;
	}

	virtual SimpleBST*getRight(){
		return right;
	}

	virtual SimpleBST*getParent(){
		return parent;
	}

	virtual Value &getValue(){
		return v;
	}

	virtual int getKey(){
		return key;
	}

	virtual void insert(SimpleBST* node) {
	    if (node->key <= key) {
	      if (!left) {
	        left = node;
	        left->parent = this;
	      } else {
	        left->insert(node);
	      }
	    } else if (node->key > key) {
	      if (!right) {
	        right = node;
	        right->parent = this;
	      } else {
	    	  right->insert(node);
	      }
	    }
	  }

	virtual SimpleBST* findRoot(){
		SimpleBST *m = this;
		while(m->getParent()){
			m=m->getParent();
		}
		return m;
	}

	virtual SimpleBST* findMin(){
		SimpleBST *m = this;
		while(m->getLeft())
			m=m->getLeft();

		return m;
	}

	virtual SimpleBST* findMax(){
		SimpleBST *m = this;
		while(m->getRight())
			m=m->getRight();

		return m;
	}

	virtual void remove(){
		if(!left && !right){
			if(parent){
				if(parent->left==this){
					parent->left=nullptr;
				}else {
					assert(parent->right==this);
					parent->right=nullptr;
				}
				parent->checkTree();
			}else{
				//do nothing
			}
		}else if (left && ! right){
			if(parent){
				if(parent->left==this){
					parent->left=left;
				}else {
					assert(parent->right==this);
					parent->right=left;
				}
			}else{
				left->parent=nullptr;
				left->checkTree();
			}
		}else if (right && ! left){
			if(parent){
				if(parent->left==this){
					parent->left=right;
				}else {
					assert(parent->right==this);
					parent->right=right;
				}
				parent->checkTree();
			}else{
				right->parent=nullptr;
				right->checkTree();
			}

		}else{
			SimpleBST * next = right->findMin();
			this->key=next->key;
			assert(!(next->right || next->left));
			next->remove();
			parent=nullptr;
			checkTree();
		}
		parent=nullptr;
	}


	virtual SimpleBST*concat(SimpleBST * toConcat){
		SimpleBST  * a = findMax();
		SimpleBST * b = toConcat->findMin();
		SimpleBST tmp;
		assert(a->key<b->key);
		tmp.key=(a->key+b->key)/2;
		tmp.left = a;
		tmp.right=b;
		tmp.remove();

		SimpleBST* root = a->findRoot();
		root->checkTree();
		return root;
	}
protected:
	virtual void checkTree(){
#ifndef NDEBUG
		if(getLeft()){
			assert(getLeft()->getKey()<=getKey());
			getLeft()->checkTree();
		}
		if(getRight()){
			assert(getRight()->getKey()>getKey());
			getRight()->checkTree();
		}
#endif
	}


};
