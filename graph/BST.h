/*
 * BST.h
 *
 *  Created on: Mar 26, 2014
 *      Author: sam
 */

#ifndef BST_H_
#define BST_H_


template<typename Value>
class BST{


public:
	BST(){

	}

	virtual ~BST(){

	}

	virtual BST*getLeft()=0;

	virtual BST*getRight()=0;

	virtual BST*getParent()=0;

	virtual Value &getValue()=0;

	virtual int getKey()=0;

	virtual void insert(BST* node) =0;

	virtual BST* findRoot(){
		BST *m = this;
		while(m->getParent()){
			m=m->getParent();
		}
		return m;
	}

	virtual BST* findMin(){
		BST *m = this;
		while(m->getLeft())
			m=m->getLeft();

		return m;
	}

	virtual BST* findMax(){
		BST *m = this;
		while(m->getRight())
			m=m->getRight();

		return m;
	}

	virtual void remove()=0;

	virtual BST*concat(BST * toConcat)=0;
};



#endif /* BST_H_ */
