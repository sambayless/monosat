#include "LinkCutCost.h"
#include <climits>
#include <limits>
#include <gmpxx.h>
template<>
const int LinkCutCost<int>::INF=INT_MAX/2;
template<>
const long LinkCutCost<long>::INF=LONG_MAX/2;
template<>
const mpq_class LinkCutCost<mpq_class>::INF=LONG_MAX;//fix this...
template<>
const double LinkCutCost<double>::INF=std::numeric_limits<double>::infinity();
template<>
const float LinkCutCost<float>::INF=std::numeric_limits<float>::infinity();
//template<typename Weight>
//const Weight LinkCutCost<Weight>::INF=std::numeric_limits<Weight>::infinity();
