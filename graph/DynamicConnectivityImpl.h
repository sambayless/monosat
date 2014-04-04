

#ifndef DYNAMIC_CONNECTIVITY_IMPL_H_
#define DYNAMIC_CONNECTIVITY_IMPL_H_
class DynamicConnectivityImpl{
public:

DynamicConnectivityImpl(){

}
virtual ~DynamicConnectivityImpl(){

}
virtual int numComponents()=0;
virtual bool connected(int u, int v)=0;
virtual void addNode()=0;
virtual void addEdge(int edgeID, int from, int to)=0;
virtual bool edgeEnabled(int edgeid)const=0;

virtual void setEdgeEnabled(int edgeid, bool enabled)=0;
};
#endif
