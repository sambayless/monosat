/*
 * NFATypes.h
 *
 *  Created on: Dec 21, 2014
 *      Author: sam
 */

#ifndef NFATYPES_H_
#define NFATYPES_H_


struct NFATransition{
	int edgeID;
	int input;
	int output;
};
struct FSMNullStatus {
	void accepts(int string, int state,int edgeID,int label) {

	}
	void generates(int string, bool generates){

	}
};
static FSMNullStatus fsmNullStatus;
#endif /* NFATYPES_H_ */
