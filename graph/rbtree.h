/*
 * rbtrees.h
 *
 *  BSD Licensed, from http://code.google.com/p/rbtrees
 */

#ifndef RBTREES_H_
#define RBTREES_H_
#include <iostream>
#define COLOR(color) (color == 0) ? "red" : "black"

namespace ASoliman {

static const short BLACK_ = 1;
static const short RED_ = 0;
static const short LEFT_ = 100;
static const short RIGHT_ = 200;

template<class Data>
class RBTreeNode {
	public:
		RBTreeNode(RBTreeNode<Data> *parent, int32_t key, Data & data) {
			this->parent_ = parent;
			this->left_ = NULL;
			this->right_ = NULL;
			this->color_ = RED_;
			this->key_ = key;
			this->data_ = data;
		}

		short getColor() {
			return this->color_;
		}

		Data & getData() {
			return this->data_;
		}

		int32_t getKey() {
			return this->key_;
		}

		RBTreeNode<Data> *getLeft() {
			return this->left_;
		}

		RBTreeNode<Data> *getRight() {
			return this->right_;
		}

		RBTreeNode<Data> *getParent() {
			return this->parent_;
		}

		void setLeft(RBTreeNode<Data> *node) {
			this->left_ = node;
		}

		void setRight(RBTreeNode<Data> *node) {
			this->right_ = node;
		}

		void setParent(RBTreeNode<Data> *node) {
			this->parent_ = node;
		}

		void setKey(int key) {
			this->key_ = key;
		}

		void setData(Data & data) {
			this->data_ = data;
		}
		void setColor(short color) {
			this->color_ = color;
		}
		~RBTreeNode() {
			delete this->left_;
			delete this->right_;
		}
	private:
		RBTreeNode<Data> *parent_;
		RBTreeNode<Data> *left_;
		RBTreeNode<Data> *right_;
		Data data_;
		int32_t key_;
		unsigned short color_;
};


template<class Data>
class RBTree {
	public:
		//Basic Functions


		//rbtreePtr createRBTree(int, void *);
		//void destroyRBTree();
	    RBTreeNode<Data> *queryTree(int key) {
			return queryTreeRecur_(this->root_, key);
		}
		RBTreeNode<Data> *insertNode(int32_t key, Data &data) {
				return treeInsert_(this->root_, key, data);
			}
	    int32_t deleteNode(int key) {
	    			RBTreeNode<Data> *z = queryTree(key);
	    			RBTreeNode<Data> *y;
	    			RBTreeNode<Data> *x;
	    			if (z == NULL)
	    				return 0;
	    			if (z->getLeft() == NULL || z->getRight() == NULL) {
	    				y = z;
	    			} else {
	    				y = getSuccessor_(z);
	    			}
	    			if (y->getLeft() != NULL) {
	    				x = y->getLeft();
	    			} else {
	    				x = y->getRight();
	    			}
	    			if (x != NULL) {
	    				x->setParent(y->getParent());
	    			}
	    			if (y->getParent() == NULL) {
	    				this->root_ = x;

	    			} else {
	    				if (y == y->getParent()->getLeft()) {
	    					y->getParent()->setLeft(y);
	    				} else {
	    					y->getParent()->setRight(y);
	    				}
	    			}
	    			if (y != z) {
	    				z->setKey(y->getKey());
	    				z->setData(y->getData());
	    			}
	    			return 0;
	    		}

		short isValidRedBlackTree() {
			return isValidRedBlackTreeRecur_(this->root_);
		}



		int32_t size() {
				return count_(this->root_, 0);
			}

			int maxDepth() {
				return maxDepthRecur_(this->root_);
			}

			RBTreeNode<Data> *getMinimum() {
				return getMinimum_(this->root_);
			}

			RBTreeNode<Data> *getMaximum() {
				return getMaximum_(this->root_);
			}

		//Visualization Functions
		void printPaths() {
			RBTreeNode<Data> *path[1000];
			printPathsRecur_(this->root_, path, 0);
		}
		void printTree() {
					printTreeRecur_(this->root_);
				}
		RBTree() {
			root_ = NULL;
		}

		RBTree(int32_t key, void *data) {
			root_ = new RBTreeNode<Data>(NULL, key, data);
		}
		~RBTree() {
			delete this->root_;
		}


	private:
		RBTreeNode<Data> *root_;

		RBTreeNode<Data> *createNode_(RBTreeNode<Data> *parent, short loc, int32_t key,
				void *data) {
			RBTreeNode<Data> *tmp = new RBTreeNode<Data>(parent, key, data);
			if (loc == LEFT_)
				parent->setLeft(tmp);
			else
				parent->setRight(tmp);
			return tmp;
		}

		RBTreeNode<Data> *queryTreeRecur_(RBTreeNode<Data> *tree, int32_t key) {
			if (tree == NULL) {
				return NULL;
			}
			if (key < tree->getKey()) {
				return queryTreeRecur_(tree->getLeft(), key);
			} else if (key > tree->getKey()) {
				return queryTreeRecur_(tree->getRight(), key);
			} else {
				return tree;
			}

		}


		int32_t maxDepthRecur_(RBTreeNode<Data> *tree) {
			if (tree == NULL) {
				return 0;
			}
			int32_t leftDepth = maxDepthRecur_(tree->getLeft());
			int32_t rightDepth = maxDepthRecur_(tree->getRight());
			if (leftDepth > rightDepth) {
				return leftDepth + 1;
			} else {
				return rightDepth + 1;
			}
		}






		short isValidRedBlackTreeRecur_(RBTreeNode<Data> *tree) {
			if (tree == NULL)
				return 1;
			if (!isValidRedBlackTreeRecur_(tree->getLeft())) {
				return 0;
			}
			if (!isValidRedBlackTreeRecur_(tree->getRight())) {
				return 0;
			}
			if (tree->getParent() == NULL && tree->getColor() == BLACK_) {
				return 1;
			}
			if (tree->getColor() == RED_ && tree->getParent() != NULL && tree->getParent()->getColor()
					== BLACK_) {
				return 1;
			}
			if (tree->getColor() == BLACK_ && tree->getParent() != NULL && tree->getParent()->getColor()
					== BLACK_) {
				return 1;
			}
			return 0;
		}



		RBTreeNode<Data> *treeInsert_(RBTreeNode<Data> *tree, int32_t key, void *data) {
			RBTreeNode<Data> *inserted;
			if (key < tree->getKey()) {
				if (tree->getLeft() == NULL) {
					inserted = createNode_(tree, LEFT_, key, data);
					insertFix_(inserted);
				} else {
					inserted = treeInsert_(tree->getLeft(), key, data);
				}
			} else {
				if (tree->getRight() == NULL) {
					inserted = createNode_(tree, RIGHT_, key, data);
					//printf("Inserted %d On Right\n", key);
					insertFix_(inserted);
				} else {
					inserted = treeInsert_(tree->getRight(), key, data);
				}
			}
			//FIX Location
			return inserted;
		}

		void insertFix_(RBTreeNode<Data> *node) {
			RBTreeNode<Data> *z = node;
			RBTreeNode<Data> *y = NULL;
			while (z != NULL && z->getParent() != NULL && z->getParent()->getColor()
					== RED_) {
				if (z->getParent()->getParent() != NULL) {
					if (z->getParent() == z->getParent()->getParent()->getLeft()) {
						//case 1,2,3
						y = uncle_(z);
						if (y != NULL && y->getColor() == RED_) {
							//case 1
							//RECOLOR
							z->getParent()->setColor(BLACK_);
							y->setColor(BLACK_);
							z->getParent()->getParent()->setColor(RED_);
							z = z->getParent()->getParent();
							continue;
						} else {
							if (z == z->getParent()->getRight()) {
								//case 2
								z = z->getParent();
								rotateLeft_(z);
							}
							z->getParent()->setColor(BLACK_);
							z->getParent()->getParent()->setColor(RED_);
							rotateRight_(z->getParent()->getParent());
							//case 3
						}
					} else {
						//case 4,5,6
						y = uncle_(z);

						if (y != NULL && y->getColor() == RED_) {
							//case 4
							//RECOLOR
							z->getParent()->setColor(BLACK_);
							y->setColor(BLACK_);
							z->getParent()->getParent()->setColor(RED_);
							z = z->getParent()->getParent();
						} else {
							if (z == z->getParent()->getLeft()) {
								//case 5
								z = z->getParent();
								rotateRight_(z);
							}
							//case 6
							z->getParent()->setColor(BLACK_);
							z->getParent()->getParent()->setColor(RED_);
							rotateLeft_(z->getParent()->getParent());
						}
					}
				} else {
					break;
				}
			}
			this->root_->setColor(BLACK_);
		}

		RBTreeNode<Data> *uncle_(RBTreeNode<Data> *myNode) {
			RBTreeNode<Data> *myGPNode = grandparent_(myNode);
			if (myGPNode == NULL) {
				return NULL;
			}
			if (myNode->getParent() == myGPNode->getLeft()) {
				return myGPNode->getRight();
			} else {
				return myGPNode->getLeft();
			}
		}

		short rotateRight_(RBTreeNode<Data> *myNode) {
			//printf("Rotate Right Of %d\n", myNode->key);
			//cool operation...
			if (myNode->getLeft() == NULL)
				return 1;
			RBTreeNode<Data> *leftNode = myNode->getLeft();
			RBTreeNode<Data> *correctParent = myNode->getParent();
			//let's fix the parent's links first...
			if (correctParent != NULL) {
				if (myNode == correctParent->getLeft()) {
					correctParent->setLeft(leftNode);
				} else {
					correctParent->setRight(leftNode);
				}
			}
			myNode->setLeft(leftNode->getRight());
			leftNode->setRight(myNode);
			leftNode->setParent(myNode->getParent());
			myNode->setParent(leftNode);

			if (this->root_ == myNode) {
				this->root_ = leftNode;
				leftNode->setParent(NULL);
			}
			return 0;
		}

		short rotateLeft_(RBTreeNode<Data> *myNode) {
			//printf("Rotate Left Of %d\n", myNode->key);
			//cool operation...
			if (myNode->getRight() == NULL)
				return 1;
			RBTreeNode<Data> *rightNode = myNode->getRight();
			RBTreeNode<Data> *correctParent = myNode->getParent();
			//let's fix the parent's links first...
			if (correctParent != NULL) {
				if (myNode == correctParent->getLeft()) {
					correctParent->setLeft(rightNode);
				} else {
					correctParent->setRight(rightNode);
				}
			}
			myNode->setRight(rightNode->getLeft());
			rightNode->setLeft(myNode);
			rightNode->setParent(myNode->getParent());
			myNode->setParent(rightNode);
			if (this->root_ == myNode) {
				this->root_ = rightNode;
				rightNode->setParent(NULL);
			}
			return 0;
		}

		short isLeaf_(RBTreeNode<Data> *node) {
			if (node->getLeft() == NULL && node->getRight() == NULL)
				return 1;
			return 0;
		}

		int32_t count_(RBTreeNode<Data> *node, int num) {
			if (node == NULL) {
				return num;
			}
			return count_(node->getLeft(), count_(node->getRight(), ++num));
		}


		void printPathsRecur_(RBTreeNode<Data> *node, RBTreeNode<Data> **path, int pathLen) {
			if (node == NULL)
				return;

			path[pathLen++] = node;
			if (isLeaf_(node)) {
				int32_t i;
				for (i = 0; i < pathLen; i++) {
					printf("%d ", path[i]->getKey());
				}
				printf("\n");
				return;
			}
			printPathsRecur_(node->getLeft(), path, pathLen);
			printPathsRecur_(node->getRight(), path, pathLen);
		}



		void printTreeRecur_(RBTreeNode<Data> *node) {
			if (node == NULL) {
				return;
			}
			printTreeRecur_(node->getLeft());
			printf("%d ", node->getKey());
			printTreeRecur_(node->getRight());
		}

		RBTreeNode<Data> *getSuccessor_(RBTreeNode<Data> *node) {
			if (node->getRight() != NULL) {
				return getMinimum_(node->getRight());
			}
			RBTreeNode<Data> *x = node;
			RBTreeNode<Data> *y = node->getParent();
			while (y != NULL && x == y->getRight()) {
				x = y;
				y = y->getParent();
			}
			return y;
		}

		RBTreeNode<Data> *getPredecessor_(RBTreeNode<Data> *node) {
			if (node->getLeft() != NULL) {
				return getMaximum_(node->getLeft());
			}
			RBTreeNode<Data> *x = node;
			RBTreeNode<Data> *y = node->getParent();
			while (y != NULL && x == y->getLeft()) {
				x = y;
				y = y->getParent();
			}
			return y;
		}
		RBTreeNode<Data> *grandparent_(RBTreeNode<Data> *myNode) {
			if ((myNode != NULL) && (myNode->getParent() != NULL)) {
				return myNode->getParent()->getParent();
			} else {
				return NULL;
			}
		}

		RBTreeNode<Data> *getMinimum_(RBTreeNode<Data> *tree) {
			if (tree == NULL) {
				return NULL;
			}
			if (tree->getLeft() != NULL) {
				return getMinimum_(tree->getLeft());
			}
			return tree;
		}



		RBTreeNode<Data> *getMaximum_(RBTreeNode<Data> *tree) {
			if (tree == NULL) {
				return NULL;
			}
			if (tree->getRight() != NULL) {
				return getMaximum_(tree->getRight());
			}
			return tree;
		}



};

}

#endif /* RBTREES_H_ */
