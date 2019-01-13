/**************************************************************************************************
 The MIT License (MIT)

 Copyright (c) 2006, Emil Stefanov

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
 OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

// Disjoint Set Data Structure
// Author: Emil Stefanov
// Date: 03/28/06
// Implementaton is as described in http://en.wikipedia.org/wiki/Disjoint-set_data_structure
// Released under the MIT license by Emil Stefanov (http://www.emilstefanov.net/Projects/DisjointSets.aspx)
#include <cassert>

#include "DisjointSets.h"

DisjointSets::DisjointSets(){
    m_numElements = 0;
    m_numSets = 0;
}

DisjointSets::DisjointSets(int count){
    m_numElements = 0;
    m_numSets = 0;
    AddElements(count);
}

DisjointSets::DisjointSets(const DisjointSets& s){
    this->m_numElements = s.m_numElements;
    this->m_numSets = s.m_numSets;

    m_nodes = s.m_nodes;
}

DisjointSets::~DisjointSets(){

}

// Note: some internal data is modified for optimization even though this method is consant.
int DisjointSets::FindSet(int elementId){
    assert(elementId < m_numElements);
    // Find the root element that represents the set which `elementId` belongs to
    int curnodeID = elementId;

    //Node * curNode  = &m_nodes[elementId];
    while(m_nodes[curnodeID].parent != -1)
        curnodeID = m_nodes[curnodeID].parent;
    int rootID = curnodeID;

    // Walk to the root, updating the parents of `elementId`. Make those elements the direct
    // children of `root`. This optimizes the tree for future FindSet invokations.
    curnodeID = elementId;
    while(curnodeID != rootID){
        int nextID = m_nodes[curnodeID].parent;
        m_nodes[curnodeID].parent = rootID;
        curnodeID = nextID;
    }

    return m_nodes[rootID].index;
}

void DisjointSets::UnionSets(int setId1, int setId2){
    assert(setId1 < m_numElements);
    assert(setId2 < m_numElements);
    assert(setId1 == FindSet(setId1));
    assert(setId2 == FindSet(setId2));
    if(setId1 == setId2)
        return; // already unioned

    Node& set1 = m_nodes[setId1];
    Node& set2 = m_nodes[setId2];

    // Determine which node representing a set has a higher rank. The node with the higher rank is
    // likely to have a bigger subtree so in order to better balance the tree representing the
    // union, the node with the higher rank is made the parent of the one with the lower rank and
    // not the other way around.
    if(set1.rank > set2.rank)
        set2.parent = setId1;
    else if(set1.rank < set2.rank)
        set1.parent = setId2;
    else // set1->rank == set2->rank
    {
        set2.parent = setId1;
        ++set1.rank; // update rank
    }

    // Since two sets have fused into one, there is now one less set so update the set count.
    --m_numSets;
    assert(FindSet(setId1) == FindSet(setId2));
}

void DisjointSets::UnionElements(int elementID1, int elementID2){
    assert(elementID1 < m_numElements);
    assert(elementID2 < m_numElements);
    int setId1 = FindSet(elementID1);
    int setId2 = FindSet(elementID2);
    UnionSets(setId1, setId2);
}

void DisjointSets::AddElements(int numToAdd){
    assert(numToAdd >= 0);

    // insert and initialize the specified number of element nodes to the end of the `m_nodes` array
    //m_nodes.insert(m_nodes.end(), numToAdd, (Node*)NULL);

    for(int i = m_numElements; i < m_numElements + numToAdd; ++i){
        m_nodes.push_back({});
        assert(i == m_nodes.size() - 1);
        //m_nodes[i] = new Node();
        m_nodes[i].parent = -1;
        m_nodes[i].index = i;
        m_nodes[i].rank = 0;
    }

    // update element and set counts
    m_numElements += numToAdd;
    m_numSets += numToAdd;
}

int DisjointSets::GetElement(int fromSet){
    if(m_numSets == 1)
        return 0;    //shortcut common case.
    if(elements.size() != m_numSets){
        elements.clear();
        seen.clear();
        seen.resize(m_numElements);
        for(int i = 0; i < m_numElements; i++){
            int s = FindSet(i);
            if(!seen[s]){
                seen[s] = true;
                elements.push_back(s);
            }
        }
    }
    return elements[fromSet];
}

int DisjointSets::NumElements() const{
    return m_numElements;
}

int DisjointSets::NumSets() const{
    return m_numSets;
}
