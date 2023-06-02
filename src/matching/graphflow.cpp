#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <set>
#include <map>
#include "../utils/types.h"
#include "../utils/globals.h"
#include "../utils/utils.h"
#include "../graph/graph.h"
#include "graphflow.h"
#include "LocalIndex.h"

struct pairCompare {
    bool operator()(const std::pair<float, int>& p1, const std::pair<float, int>& p2) {
        if (p1.first == p2.first) {
            return p1.second > p2.second;
        }
        return p1.first < p2.first;
    }
};

bool matchRecordCmp (MatchRecord* d1, MatchRecord* d2) {
    if(d1->getDensity()!=d2->getDensity()) {
        return d1->getDensity() > d2->getDensity();
    }
    else if(d1->getTmin()!=d2->getTmin()){
        return d1->getTmin() > d2->getTmin();
    }
    else {
        return (*d1->getVetex())>(*d2->getVetex());
    }
}
bool ForwardNeighborcmp(ForwardNeighbor*f1,ForwardNeighbor*f2){
    return (*f1)>(*f2);
}
bool areSame(float a, float b) {
    return std::fabs(a - b) < std::numeric_limits<float>::epsilon();
}
bool tupleVertexIdCmp(const std::tuple<int, int, float>& a, const std::tuple<int, int, float>& b) {
    if (std::get<0>(a) != std::get<0>(b)) {
        return std::get<0>(a) >std::get<0>(b);
    } else {
        return std::get<1>(a) > std::get<1>(b);
    }
}

bool tupleVertexIdCmp2(const std::tuple<int, int, float>& a, const std::tuple<int, int, float>& b) {
    if (std::get<0>(a) != std::get<0>(b)) {
        return std::get<0>(a) <std::get<0>(b);
    } else {
        return std::get<1>(a) < std::get<1>(b);
    }
}


bool tupleSumWeightCmp(const std::tuple<int, int, float>& a, const std::tuple<int, int, float>& b) {
    if (std::get<2>(a) != std::get<2>(b)) {
        return std::get<2>(a) >std::get<2>(b);
    } else {
        return std::get<0>(a) > std::get<0>(b);
    }
}
bool tupleResultCmp(const std::tuple<std::vector<int>,int,float>&a,const std::tuple<std::vector<int>,int,float>&b){
    if (std::get<2>(a) != std::get<2>(b)) {
        return std::get<2>(a) >std::get<2>(b);
    }
    else if(std::get<1>(a) != std::get<1>(b))  {
        return std::get<1>(a) > std::get<1>(b);
    }
    else{
        return std::get<0>(a) > std::get<0>(b);
    }
}
Graphflow::Graphflow(Graph& query_graph, Graph& data_graph,
                     uint max_num_results,
                     bool print_prep,
                     bool print_enum,
                     bool homo)
        : matching(query_graph, data_graph, max_num_results,
                   print_prep, print_enum, homo)
        , order_vs_(query_.NumEdges())
        , order_csrs_(query_.NumEdges())
        , order_offs_(query_.NumEdges())
        ,order_vertex_index(query_.NumEdges())
        ,topKSet(0)
        ,allMatchRecords(0)
        ,suffixMax(query_graph.NumVertices(),0)
        ,isolatedMax(query_graph.NumVertices(),-1)
        ,rightNeighbor(query_graph.NumEdges())
{
    qForwardNeighbors.resize(query_.NumEdges());
    for (uint i = 0; i < query_.NumEdges(); ++i)
    {
        order_vs_[i].resize(query_.NumVertices());//节点个数
        order_csrs_[i].resize(query_.NumEdges() + 1);//边的个数+1，
        order_offs_[i].resize(query_.NumVertices(), 0);//节点个数，初始化为0
        qForwardNeighbors[i].resize(query_.NumVertices());
        order_vertex_index[i].resize(query_.NumVertices());
        rightNeighbor[i].resize(query_.NumVertices());
    }
}


void Graphflow::Preprocessing()//预处理过程
{
    this->data_.InitLabelIndex();
    GenerateMatchingOrder();
    this->query_.InitLabelIndex();
    this->query_.InitMatchOrderType(this->order_vs_,this->qForwardNeighbors, this->order_vertex_index);
    //todo make starindex
    CreateStarIndex();
    //createSameLabelVertex();
    std::cout<<"Preprocess end"<<endl;

}
void Graphflow::updateStarIndex(uint match_index,uint caddidate_v,uint candidate_u) {
    int vertex_index=order_vertex_index[match_index][candidate_u];
    StarGraph* s=qForwardNeighbors[match_index][vertex_index];
    std::vector<ForwardNeighbor*>&queryVetex=s->GetqueryVertex();
    std::vector<Neighbor>vN= this->data_.vNeighbors[caddidate_v];
    std::sort(queryVetex.begin(),queryVetex.end(), ForwardNeighborcmp);
    std::sort(vN.begin(),vN.end(),greater<Neighbor>());
    int leftvN=0;
    int rightqV=0;
    StarGraph *tmpgraph=new StarGraph();
    while(leftvN<vN.size()&&rightqV<queryVetex.size()){
        if(vN[leftvN].GetelabelAndVertexLabel()<queryVetex[rightqV]->GetelabelAndVertexLabel())
        {
            break;
        }
        int flag=1;
        while(vN[leftvN].GetelabelAndVertexLabel()>queryVetex[rightqV]->GetelabelAndVertexLabel())
        {
            leftvN++;
            if(leftvN>=vN.size())
            {
                flag=0;
                break;
            }
        }
        if(!flag)
            break;
        if(vN[leftvN].GetelabelAndVertexLabel()==queryVetex[rightqV]->GetelabelAndVertexLabel()){
            uint tovId=queryVetex[rightqV]->GetVetexId();
            pair<uint,uint>tmppair=vN[leftvN].GetelabelAndVertexLabel();
            float edgeweight=vN[leftvN].GetEdgeWeight();
            ForwardNeighbor* f=new ForwardNeighbor(tovId,tmppair.second,tmppair.first,edgeweight);
            tmpgraph->AddForwardNeighbor(f);
            rightqV++;
            leftvN++;
        }
    }
    tmpgraph->computeMaxWeight();
    //find tmpgraph
    if(tmpgraph->GetForwardNeighborNum()==queryVetex.size()){
        if(s->GetStarMaxWeight()==queryVetex.size()*mw||s->GetStarMaxWeight()<tmpgraph->GetStarMaxWeight())
        {
            StarGraph * t=s;
            qForwardNeighbors[match_index][vertex_index]=tmpgraph;
            delete t;
        }
        else {
            delete tmpgraph;
        }
    }
}
float Graphflow::GetBackWeight(uint order_index,uint depth) {
    float sum=0;
    uint n=query_.NumVertices();
    std::vector<uint>& matchOrder= this->order_vs_[order_index];
    for(int i=depth;i<n;i++){
        uint id=matchOrder[i];
        float tmp=queryLocalIndexs[id].getMaxWeight();
        if(tmp!=FLT_MAX)
            sum+=tmp;
        else
        sum+=qForwardNeighbors[order_index][i]->GetStarMaxWeight();
    }
    return sum;
}
void Graphflow::CreateStarIndex() {
    //create matchOrderIndexToMaxWeight;
    //按照密度排序
   /* int qlen=qForwardNeighbors.size();
    for(int i=0;i<qlen;i++){
        int k=qForwardNeighbors[i].size();
        for(int j=1;j<k;j++){
            StarGraph* s=qForwardNeighbors[i][j];
            std::vector<ForwardNeighbor*>&queryVetex=s->GetqueryVertex();
            std::sort(queryVetex.begin(),queryVetex.end(), ForwardNeighborcmp);
        }
    }*/
    uint n= this->data_.vEdge.size();
    for(int i=0;i<n;i++){
        Edge*e=&data_.vEdge[i];
        uint v1=e->GetV1();
        uint v2=e->GetV2();
        uint v1label=e->GetV1Label();
        uint v2label=e->GetV2Label();
        float eWeight=e->GeteWeight();
        vector<int>m=EdgeisInMatchOrder(e);

        if(m.size()==0)
            continue;
        //for each edge<v1,v2>matches u1-->u2 or u2-->u1
            for(int i=0;i<m.size();i++){
               uint u1=order_vs_[m[i]][0];
               uint u2=order_vs_[m[i]][1];

               //对所有匹配序中
               for(int j=0;j<query_.NumEdges();j++){
                   uint candidate_u=order_vertex_index[j][u1]>order_vertex_index[j][u2]?u1:u2;
                   if(query_.GetVertexLabel(candidate_u)==v1label)
                   {

                       updateStarIndex(j,v1,candidate_u);
                   }
                   if(query_.GetVertexLabel(candidate_u)==v2label)
                   {
                       updateStarIndex(j,v2,candidate_u);
                   }
               }
        }
    }

}
vector<int> Graphflow::EdgeisInMatchOrder(uint v1, uint v2, uint v1label, uint v2label,uint velabel) {
    vector<int>result;
    for(int i=0;i<order_vs_.size();i++){
        uint u1=order_vs_[i][0];
        uint u2=order_vs_[i][1];
        uint u1label=query_.GetVertexLabel(u1);
        uint u2label=query_.GetVertexLabel(u2);
        uint qlabel=std::get<2>(query_.GetEdgeLabel(u1,u2));
        if((v1label==u1label&&v2label==u2label&&velabel==qlabel)||(v1label==u2label&&v2label==u1label&&velabel==qlabel))
        {
            result.emplace_back(i);
        }
    }
    return result;
}
vector<int> Graphflow::EdgeisInMatchOrder(Edge *edge){
    uint v1=edge->GetV1();
    uint v2=edge->GetV2();
    uint v1label=edge->GetV1Label();
    uint v2label=edge->GetV2Label();
    uint velabel=edge->GeteLabel();
    vector<int>result;
    for(int i=0;i<order_vs_.size();i++){
        uint u1=order_vs_[i][0];
        uint u2=order_vs_[i][1];
        uint u1label=query_.GetVertexLabel(u1);
        uint u2label=query_.GetVertexLabel(u2);
        uint qlabel=std::get<2>(query_.GetEdgeLabel(u1,u2));
        if((v1label==u1label&&v2label==u2label&&velabel==qlabel)||(v1label==u2label&&v2label==u1label&&velabel==qlabel))
        {
            result.emplace_back(i);
        }
    }
    return result;
    //返回其在哪一个matchorder中被匹配
}
void Graphflow::GenerateMatchingOrder()
{
    // generate the initial matching order, order_*s_[0]
    std::vector<bool> visited(query_.NumVertices(), false);
    uint max_degree = 0u;
    //首先找到的是度最大的节点
    for (size_t i = 0; i < query_.NumVertices(); i++)
    {
        if (query_.GetDegree(i) > max_degree)
        {
            max_degree = query_.GetDegree(i);
            order_vs_[0][0] = i;
            order_vertex_index[0][i]=0;
        }
    }
    visited[order_vs_[0][0]] = true;

    // loop over all remaining positions of the order
    for (uint i = 1; i < query_.NumVertices(); ++i)
    {
        uint max_adjacent = 0;
        uint max_adjacent_u = NOT_EXIST;
        //找到不在序列中，但是在序列中的邻居数量最多的顶点，添加到排序中
        for (size_t j = 0; j < query_.NumVertices(); j++)
        {
            uint cur_adjacent = 0u;
            if (visited[j]) continue;

            auto& q_nbrs = query_.GetNeighbors(j);
            for (auto& other : q_nbrs)
                if (visited[other])
                    cur_adjacent++;

            if (!cur_adjacent) continue;
            if (
                    max_adjacent_u == NOT_EXIST ||
                    (cur_adjacent == max_adjacent &&
                     query_.GetDegree(j) > query_.GetDegree(max_adjacent_u)) ||
                    cur_adjacent > max_adjacent
                    ) {
                max_adjacent = cur_adjacent;
                max_adjacent_u = j;
            }
        }
        order_vs_[0][i] = max_adjacent_u;
        order_vertex_index[0][max_adjacent_u]=i;
        visited[max_adjacent_u] = true;
        order_offs_[0][i] = order_offs_[0][i - 1];
        auto& q_nbrs = query_.GetNeighbors(max_adjacent_u);
        StarGraph*s=new StarGraph();
        for (auto &other: q_nbrs)
        {
            if (visited[other])
            {
                uint qlabel= std::get<2>(this->query_.GetEdgeLabel(max_adjacent_u,other));
                ForwardNeighbor* forwardNeighbor=new ForwardNeighbor(other, this->query_.GetVertexLabel(other),qlabel);
                s->AddForwardNeighbor(forwardNeighbor);
                order_csrs_[0][order_offs_[0][i]++] = other;
            }
        }
        s->InitalmaxWeight();
        qForwardNeighbors[0][i]=s;

    }

    // generate other incremental matching orders
    for (uint i = 1; i < query_.NumEdges(); ++i)
    {
        std::vector<bool> visited(query_.NumVertices(), false);

        // get the first edge
        std::vector<uint>::iterator it = std::lower_bound(
                order_offs_[0].begin(), order_offs_[0].end(), i + 1
        );
        uint tmp= *(order_vs_[0].begin() + std::distance(order_offs_[0].begin(), it));
        order_vs_[i][0]=tmp;
        order_vertex_index[i][tmp]=0;
        order_vs_[i][1] = order_csrs_[0][i];
        order_vertex_index[i][order_csrs_[0][i]]=1;
        order_csrs_[i][0] = order_vs_[i][0];
        StarGraph*s=new StarGraph();
        uint qlabel=std::get<2>(this->query_.GetEdgeLabel(order_vs_[i][0],order_vs_[i][1]));
        ForwardNeighbor *forwardNeighbor=new ForwardNeighbor(order_vs_[i][0], this->query_.GetVertexLabel(order_vs_[i][0]),qlabel);
        s->AddForwardNeighbor(forwardNeighbor);
        s->InitalmaxWeight();
        qForwardNeighbors[i][1]=(s);

        visited[order_vs_[i][0]] = true;
        visited[order_vs_[i][1]] = true;

        order_offs_[i][2] = order_offs_[i][1] = 1;
        for (uint j = 2; j < query_.NumVertices(); ++j)
        {
            uint max_adjacent = 0;
            uint max_adjacent_u = NOT_EXIST;
            for (size_t k = 0; k < query_.NumVertices(); k++)
            {
                uint cur_adjacent = 0u;
                if (visited[k]) continue;

                auto& q_nbrs = query_.GetNeighbors(k);
                for (auto& other : q_nbrs)
                    if (visited[other])
                        cur_adjacent++;

                if (!cur_adjacent) continue;
                if (
                        max_adjacent_u == NOT_EXIST ||
                        (cur_adjacent == max_adjacent &&
                         query_.GetDegree(k) > query_.GetDegree(max_adjacent_u)) ||
                        cur_adjacent > max_adjacent
                        ) {
                    max_adjacent = cur_adjacent;
                    max_adjacent_u = k;
                }
            }
            order_vs_[i][j] = max_adjacent_u;
            order_vertex_index[i][max_adjacent_u]=j;
            visited[max_adjacent_u] = true;

            order_offs_[i][j] = order_offs_[i][j - 1];
            StarGraph*s=new StarGraph();
            auto& q_nbrs = query_.GetNeighbors(max_adjacent_u);
            for (auto &other: q_nbrs)
            {
                if (visited[other])
                {
                    order_csrs_[i][order_offs_[i][j]++] = other;
                    qlabel=std::get<2>(this->query_.GetEdgeLabel(max_adjacent_u,other));
                    ForwardNeighbor* forwardNeighbor=new ForwardNeighbor(other, this->query_.GetVertexLabel(other),qlabel);
                    s->AddForwardNeighbor(forwardNeighbor);
                }
            }
            s->InitalmaxWeight();
            qForwardNeighbors[i][j]=(s);
        }
    }
    //创建所有节点的前向邻居数组
    this->query_.forwardNeighbors.resize(this->query_.NumEdges());
    for(int i=0;i<query_.forwardNeighbors.size();i++){
        query_.forwardNeighbors[i].resize(query_.NumVertices());
    }
    if (print_preprocessing_results_)
    {
        std::cout << "matching order: " << std::endl;
        std::cout << "-vertex(backward neighbors)-\n";
        for (uint i = 0; i < query_.NumEdges(); ++i)
        {
            std::cout << "#" << i << ": ";
            for (uint j = 0; j < query_.NumVertices(); ++j)
            {
                std::cout << order_vs_[i][j];
                if (j == 0)
                {
                    this->query_.forwardNeighbors[i][j]={};
                    std::cout << "-";
                    continue;
                }
                std::vector<ForwardNeighbor>currentQueryNeighbors;
                for (uint k = order_offs_[i][j - 1]; k < order_offs_[i][j]; k++)
                {
                    uint toVertexId=order_csrs_[i][k];
                    uint toVertexLabel=query_.GetVertexLabel(order_csrs_[i][k]);
                    uint edgelabel=std::get<2>(query_.GetEdgeLabel(order_vs_[i][j],toVertexId));
                    ForwardNeighbor f(toVertexId,toVertexLabel,edgelabel);
                    currentQueryNeighbors.push_back(f);
                    rightNeighbor[i][toVertexId].emplace_back(order_vs_[i][j]);
                    if (k == order_offs_[i][j - 1]) std::cout << "(";
                    std::cout << order_csrs_[i][k];
                    if (k != order_offs_[i][j] - 1) std::cout << ",";
                    else std::cout << ")";
                }
                this->query_.forwardNeighbors[i][j]=currentQueryNeighbors;
                if (j != query_.NumVertices() - 1)
                    std::cout << "-";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}
void Graphflow::InitialTopK(const std::string &path) {

    if (!io::file_exists(path.c_str()))
    {
        std::fstream fp(path,std::ios::out);
       for(auto t:topKSet){
           fp << t->printMatchRecord();
       }
       fp.close();
    }
#ifdef LOG_TRACK
    stringstream _ss;
    _ss<<"Initial Top k"<<std::endl;
    for(auto d:topKSet){
        if(d!=NULL){
            _ss<<"address "<<d;
            _ss<<" density:"<<d->getDensity()<<" tmin:"<<d->getTmin()
               <<" vetexs:";
            std::vector<uint> *vs=d->getVetex();
            for(int j=0;j<(*vs).size();j++){
                _ss<<(*vs)[j]<<" ";
            }
            _ss<<std::endl;
            Log::track1(_ss);
            _ss.clear();
            _ss.str("");
        }

    }
#endif

#ifdef RESULT_TRACK
    stringstream _ss1;
    _ss1<<"Initial Top k"<<std::endl;
   // std::sort(topKSet.begin(),topKSet.end(), matchRecordCmp);
    for(auto d:topKSet){
        if(d!=NULL){
            _ss1<<d->toString();
            Log::track2(_ss1);
            _ss1.clear();
            _ss1.str("");

        }
    }
#endif
}

void Graphflow::updateTopK(uint num) {
#ifdef PRINT_DEBUG
    /*for(auto it=edgeFlags.begin();it!=edgeFlags.end();it++){
        std::cout<<"edgeFlags["<<it->first.first<<","<<it->first.second<<"] "<<edgeFlags[std::make_pair(it->first.first,it->first.second)]<<std::endl;
    }*/
#endif
    //查看是否更新top k  s t是此次add edge的边的两个节点
 /*   if(num==0)
        return;
    int n=allMatchRecords.size();
    if(allMatchRecords.size()<k)
    {
        topKSet.resize(n);
        copy(allMatchRecords.begin(),allMatchRecords.end(),topKSet.begin());
    }
    else{
        topKSet.resize(k);
        copy(allMatchRecords.begin(),allMatchRecords.begin()+k,topKSet.begin());
    }*/
    stringstream _ss;
#ifdef LOG_TRACK
    for(auto d:topKSet){
        _ss<<"address "<<d;
        _ss<<" density:"<<d->getDensity()<<" tmin:"<<d->getTmin()
           <<" vetexs:";
        std::vector<uint>*vs=d->getVetex();
        for(int j=0;j<(*vs).size();j++){
            _ss<<(*vs)[j]<<" ";
        }
        _ss<<std::endl;
        Log::track1(_ss);
        _ss.clear();
        _ss.str("");
    }
#endif

#ifdef RESULT_TRACK
    _ss<<"after insert "<<std::endl;
    for(auto d:topKSet){
        _ss<<d->toString();
        Log::track2(_ss);
        _ss.clear();
        _ss.str("");
    }
#endif

}

void Graphflow::searchMatches(uint matchorderindex, searchType flag)  {
    std::chrono::high_resolution_clock::time_point start;
    /*if(this->visited_[61070]==true){
        std::cout<<"depth:"<<match.size()<<endl;
    }*/
    //1.找到前向邻居，找到候选解
    start=Get_Time();
    uint depth=this->match.size();
    std::vector<uint>& matchOrder= this->order_vs_[matchorderindex];
    uint queryVertex=matchOrder[depth];
    uint queryVertexLabel=this->query_.GetVertexLabel(queryVertex);
    const auto & fneighbors=this->query_.forwardNeighbors[matchorderindex][depth];
    vertexType currentSearchVertextype = this->query_.GetVertexType(matchorderindex,depth);
    std::vector<tuple<int,int,float>> singleVertexCandidate;
//    std::vector<uint>LDVertexs;
    uint tmin=UINT_MAX;
    float density=0;
    float maxweight=0;
    if(fneighbors.size()!=0){
       // uint u_min=UINT_MAX;
        uint u_min_index=UINT_MAX;
        uint u_min_neighbors_size=UINT_MAX;
        uint u_min_elabel=UINT_MAX;
        //find umin
       for (uint i=0;i<fneighbors.size();i++){
           ForwardNeighbor neighbor=fneighbors[i];
           uint neighbor_vertex=neighbor.GetVetexId();
           uint neighbor_index=order_vertex_index[matchorderindex][neighbor_vertex];
           uint cur_can_size=data_.GetNeighbors(neighbor_vertex).size();
           if(cur_can_size<u_min_neighbors_size){
               u_min_neighbors_size=cur_can_size;
               u_min_index=neighbor_index;
               u_min_elabel=neighbor.GetElabel();
           }
       }
       if(u_min_index!=UINT_MAX){
           uint u_min=std::get<0>(this->match[u_min_index]);
           const auto &u_min_nbrs=this->data_.GetNeighbors(u_min);
           const auto &u_min_nbrs_labels= this->data_.GetNeighborLabels(u_min);
           float tmp;
           for(int i=0;i<u_min_nbrs.size();i++){
               const uint v=u_min_nbrs[i];
               tmp=0;
               tmin=UINT_MAX;
               density=0;
               if(this->data_.GetVertexLabel(v)!=queryVertexLabel||u_min_nbrs_labels[i]!=u_min_elabel)
                   continue;
               if(!homomorphism_&& this->visited_[v]==true)
                   continue;
               tmp+=data_.GetEdgeWeight(u_min,v);
               tmin=std::min(tmin,data_.GetEdgeTime(u_min,v));
               bool joinable=true;
               for(int k=0;k<fneighbors.size();k++){
                   ForwardNeighbor neighbor=fneighbors[k];
                   uint neigh_vertex=neighbor.GetVetexId();
                   uint data_v_index=order_vertex_index[matchorderindex][neigh_vertex];
                   if(data_v_index==u_min_index)
                       continue;
                   uint data_v=std::get<0>(this->match[data_v_index]);
                   uint elabel=neighbor.GetElabel();
                   const auto & dataNeighbors=this->data_.GetNeighbors(data_v);
                   auto it =std::lower_bound(dataNeighbors.begin(),dataNeighbors.end(),v);
                   uint dis=std::distance(dataNeighbors.begin(),it);
                   if(it==dataNeighbors.end()||this->data_.GetNeighborLabels(data_v)[dis]!=elabel||*it!=v){
                       joinable= false;
                       break;
                   }
                   tmp+=data_.GetEdgeWeight(data_v,v);
                   tmin=std::min(tmin,data_.GetEdgeTime(data_v,v));
               }
               if(!joinable)
                   continue;
               density+=tmp;
               singleVertexCandidate.emplace_back(v,tmin, density);
               if(density>maxweight)
                   maxweight=density;
       }
       }
    }
    Print_Time2("findCandidate ",start);
        //若不是孤立节点 那么从候选解决中加入此节点，继续递归
        if(singleVertexCandidate.size()==0){
            return;
        }
    start=Get_Time();
    densityFilter(matchorderindex,depth,singleVertexCandidate);
    if(singleVertexCandidate.size()==0)
        return;
    Print_Time2("densityFilter ",start);
        if(currentSearchVertextype==freeVertex){
          //  start=Get_Time();
            //密度剪枝
            for(int i=0;i<singleVertexCandidate.size();i++){
                uint dataV=std::get<0>(singleVertexCandidate[i]);
                //递归
                matchVertex(dataV,i,depth,singleVertexCandidate);
                //updateweight;
                uint id=order_vs_[matchorderindex][match.size()-1];
                std::vector<uint>u1_rn=rightNeighbor[matchorderindex][id];
                updateDensityCandidate(matchorderindex,id,dataV,u1_rn);
                searchMatches(matchorderindex,flag);
                 clearDensityCandidate(u1_rn);
                popVertex(dataV,i,depth);

            }
        }
        else{
            isolatedMax[depth]=maxweight;
            this->matchVertex(singleVertexCandidate);
            // 若递归到终点，释放所有的孤立节点
            if(depth==this->query_.NumVertices()-1){
                //todo
                start=Get_Time();
                addMatchResult(matchorderindex,flag);
                //addMatchResultWithHeap(matchorderindex,flag);
                Print_Time2("addMatchResult ",start);
            }
            else{
                searchMatches(matchorderindex,flag);
            }
            this->popVertex();
        }
}


//flag==0为initial flag=1为update
void Graphflow::FindMatches(uint flag,uint order_index, uint depth, std::vector<uint> m, size_t &num_results, float density_s,
                            uint tmin) {
    if (reach_time_limit) return;
#ifdef PRINT_DEBUG
    if(density_s<0){
            std::cout<<"density <0 u:"<<order_vs_[order_index][depth]<<" depth: "<<depth<<"density: "<<density_s<<" match: ";
            for(int i=0;i<m.size();i++){
                std::cout<<m[i]<<" ";
            }
            std::cout<<std::endl;
        }
#endif
    if(flag==1) {
        float back_max_result = GetBackWeight(order_index, depth);
        uint n = query_.NumVertices();
        if (topKSet.size() == k) {
            float weight = topKSet.back()->getDensity();
            uint tmpdensity = density_s + back_max_result;
            if (tmpdensity / (sqrt(n) * (n - 1)) < weight)
                return;
        }
    }
    uint u = order_vs_[order_index][depth];
    uint u_min = NOT_EXIST;
    uint u_min_label = NOT_EXIST;
    uint u_min_size = UINT_MAX;
    uint tmptmin = UINT_MAX;
    uint copytmin = tmin;

    // find u_min
    const auto &q_nbrs = query_.GetNeighbors(u);
    const auto &q_nbr_labels = query_.GetNeighborLabels(u);
    for (uint i = 0u; i < q_nbrs.size(); i++) {
        //q_nbrs=0,3  u_other=0
        const uint u_other = q_nbrs[i];
        const uint u_other_label = q_nbr_labels[i];

        if (m[u_other] == UNMATCHED) continue;

        const uint cur_can_size = data_.GetNeighbors(m[u_other]).size();
        if (cur_can_size < u_min_size) {
            u_min_size = cur_can_size;
            u_min = u_other;
            u_min_label = u_other_label;
        }
    }



    //首先找到和当前需要匹配的查询点u有邻接关系的已经匹配的节点u_other，选择度最小的u_other为u_min
    //因此u_min已经匹配过了，其u和u_min是邻居关系。那么在数据图中的u_min的匹配的点邻居中，一定有标签和u相同并且边标签相同的
    //如果没有则不能匹配
    const auto &u_min_nbrs = data_.GetNeighbors(m[u_min]);
    const auto &u_min_nbr_labels = data_.GetNeighborLabels(m[u_min]);

    float tmp;
    bool candidate_empty = true;
    for (uint i = 0u; i < u_min_nbrs.size(); i++) {
        const uint v = u_min_nbrs[i];
        tmp = 0;
        tmin=copytmin;
        // 1. check labels
        num_intermediate_results_before_index_check_++;
        if (
                data_.GetVertexLabel(v) != query_.GetVertexLabel(u) ||
                u_min_nbr_labels[i] != u_min_label
                )
            continue;
        num_intermediate_results_after_index_check_++;
        //代表u_min和v之前的权值之和
        tmp += data_.GetEdgeWeight(m[u_min], v);
        tmptmin = data_.GetEdgeTime(m[u_min], v);
        tmin = std::min(tmin, tmptmin);


        //考虑在已经匹配的查询图u的邻居中，除了u_min以外的其他邻居，在数据图中都必须与u的匹配点v有边相连
        //因此若u_other在数据图中没有邻居是v。则也不匹配
        // 2. check if joinable
        bool joinable = true;
        for (uint j = 0u; j < q_nbrs.size(); j++) {
            const uint u_other = q_nbrs[j];
            const uint u_other_labels = q_nbr_labels[j];

            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto it = std::lower_bound(data_.GetNeighbors(m[u_other]).begin(), data_.GetNeighbors(m[u_other]).end(),
                                       v);
            uint dis = std::distance(data_.GetNeighbors(m[u_other]).begin(), it);
            if (
                    it == data_.GetNeighbors(m[u_other]).end() ||
                    *it != v ||
                    data_.GetNeighborLabels(m[u_other])[dis] != u_other_labels
                    ) {
                joinable = false;
                break;
            }
            tmp += data_.GetEdgeWeight(m[u_other], v);
            tmptmin = data_.GetEdgeTime(m[u_other], v);
            tmin = std::min(tmin, tmptmin);
        }
        if (!joinable) continue;
        num_intermediate_results_after_joinability_check_++;

        candidate_empty = false;

        // 3. check if visited
        if (!homomorphism_ && visited_[v]) continue;
        num_intermediate_results_after_visit_check_++;

        // 4. add a vertex mapping

            m[u] = v;
            visited_[v] = true;
            density_s += tmp;

        if (depth == query_.NumVertices() - 1) {


            //打印m中的所有匹配的节点
            /*   MatchRecord *matchrecord = new MatchRecord;

               for (auto j: m) {
                   matchrecord->AddVetex(j);
               }*/
//                density_s = density_s / (sqrt(m.size()) * (m.size() - 1));

                float lastds = density_s / (sqrt(m.size()) * (m.size() - 1));
                //sort(m.begin(),m.end());
              num_results++;
                MatchRecord *r = new MatchRecord(lastds, tmin, m);
                addMatchRecords(r);


            if (print_enumeration_results_) {
                for (auto j: m) {

                    std::cout << j << " ";
                }
            }

        } else {
            size_t num_results_before_recursion = num_results;
            FindMatches(flag,order_index, depth + 1, m, num_results, density_s, tmin);
            if (num_results == num_results_before_recursion) {
                num_intermediate_results_without_results_++;
            }

        }


            visited_[v] = false;
            m[u] = UNMATCHED;


#ifdef PRINT_DEBUG
            if (density_s < tmp) {
                std::cout << "u:" << order_vs_[order_index][depth] << " depth: " << depth << "density: " << density_s
                          << " tmp:" << tmp << " match: ";
                for (int i = 0; i < m.size(); i++) {
                    std::cout << m[i] << " ";
                }
                std::cout << std::endl;
            }
#endif
            density_s -= tmp;

//            tmin = copytmin;//回溯将tmin转为初始状态
        if (num_results >= max_num_results_) return;
        if (reach_time_limit) return;
    }
    if (candidate_empty) num_intermediate_results_with_empty_candidate_set_++;
}


bool Graphflow:: addMatchRecords(MatchRecord *r) {
    //sort(r->getVetex()->begin(),r->getVetex()->end());
    int n=topKSet.size();
    if(n<k){
        for(int j=n-1;j>=0;j--){
            if((*topKSet[j])>(*r)){
                topKSet.insert(topKSet.begin()+j+1,r);
                break;
            }
            if(j==0){
                topKSet.insert(topKSet.begin(),r);
            }
        }
        if(n==0){
            topKSet.insert(topKSet.begin(),r);
        }
        return true;
    }
    else{
        MatchRecord* d = topKSet.back();
        if((*r)>(*d))
        {
            delete d;
            topKSet.pop_back();
            int m=topKSet.size();
            for(int j=m-1;j>=0;j--){
                if((*topKSet[j])>(*r)){
                    topKSet.insert(topKSet.begin()+j+1,r);
                    break;
                }
                if(j==0){
                    topKSet.insert(topKSet.begin(),r);
                }
            }
#ifdef PRINT_DEBUG
            stringstream _ss;
            _ss<<"insert record: ";
            _ss<<r->toString();
            _ss<<"after insert top k "<<std::endl;
            for(auto d:topKSet){
                _ss<<d->toString();
                Log::track1(_ss);
                _ss.clear();
                _ss.str("");
            }
#endif
            return true;
        }
        else{
            delete r;
            return false;
        }
        return true;
    }
   /* int n=allMatchRecords.size();
    for(int j=n-1;j>=0;j--){
        if((*allMatchRecords[j])>(*r)){
            allMatchRecords.insert(allMatchRecords.begin()+j+1,r);
            break;
        }
        if(j==0){
            allMatchRecords.insert(allMatchRecords.begin(),r);
        }
    }
    if(n==0){
        allMatchRecords.insert(allMatchRecords.begin(),r);
    }*/
}

void Graphflow::InitialMatching(const std::string &path) {

    if (!io::file_exists(path.c_str()))
    {
        std::cout << "the file not exit " << path << std::endl;
        std::vector<uint> m(query_.NumVertices(), UNMATCHED);
        uint flag=0;
        float density_s=0;
        uint tmin=INT_MAX;
        uint order_index=0;
        uint depth=1;
        //#pragma omp parallel for num_threads(10) firstprivate(m) firstprivate(visited_) firstprivate(density_s) firstprivate(flag) firstprivate(tmin) firstprivate(order_index) firstprivate(depth)
        for (size_t i = 0; i < data_.NumVertices(); i++)
        {
            //std::cout<<"thread id"<<omp_get_thread_num<<endl;
            if (data_.GetVertexLabel(i) != NOT_EXIST)
            {
#ifdef PRINT_DEBUG
                stringstream _ss;
                _ss<<"vertex id:"<<i<<std::endl;
                Log::track1(_ss);
#endif
                if (query_.GetVertexLabel(order_vs_[0][0]) == data_.GetVertexLabel(i))
                {
                    m[order_vs_[0][0]] = i;
                    visited_[i] = true;

                    FindMatches(flag,order_index, depth, m, num_initial_results_,density_s,tmin);
                    visited_[i] = false;
                    m[order_vs_[0][0]] = UNMATCHED;
                }

            }

        }
    }
    else{
        std::ifstream ifs2(path);
        std::cout<<"load topk from file...."<<std::endl;
        char type;
        while (ifs2 >> type){
            if(type == 't'){
                float density;
                uint tmin,tmp;
                std::vector<uint>m;
                ifs2>>density>>tmin;
                for(int i=0;i<query_.NumVertices();i++){
                    ifs2>>tmp;
                    m.emplace_back(tmp);
                }
                MatchRecord* matchRecord=new MatchRecord(density,tmin,m);
                addMatchRecords(matchRecord);
            }
        }
    }
}

//动态的加边操作
void Graphflow::AddEdge(uint v1, uint v2, uint label, float weight, uint timestamp) {
    std::chrono::high_resolution_clock::time_point start;
    data_.AddEdge(v1, v2, label, weight, timestamp, 1);

    std::vector<uint> m(query_.NumVertices(), UNMATCHED);
    uint s = std::min(v1, v2);
    uint t = std::max(v1, v2);
    uint v1label=data_.GetVertexLabel(v1);
    uint v2label=data_.GetVertexLabel(v2);

    //update the index
    start=Get_Time();
    vector<int>match=EdgeisInMatchOrder(v1,v2,v1label,v2label,label);
    //for each edge<v1,v2>matches u1-->u2 or u2-->u1
    for(int i=0;i<match.size();i++){
        uint u1=order_vs_[match[i]][0];
        uint u2=order_vs_[match[i]][1];
        //对其他的边
        for(int j=0;j<query_.NumEdges();j++){
            uint candidate_u=order_vertex_index[j][u1]>order_vertex_index[j][u2]?u1:u2;
            if(query_.GetVertexLabel(candidate_u)==v1label)
            {
                updateStarIndex(j,v1,candidate_u);
            }
            if(query_.GetVertexLabel(candidate_u)==v2label)
            {
                updateStarIndex(j,v2,candidate_u);
            }
        }
    }

    Print_Time2("UpdateIndex ", start);
    size_t num_results = 0ul;

    start=Get_Time();
    this->data_.UpdateLabelIndex(v1,v2,label,1);
    std::vector<int>mCandidate;
    for(int i=0;i<match.size();i++){
        uint u1=order_vs_[match[i]][0];
        uint u2=order_vs_[match[i]][1];
        if(v1label==v2label){
            if((this->LabelFilter(v1,u1)&&this->LabelFilter(v2,u2))||(this->LabelFilter(v1,u2)&& this->LabelFilter(v2,u1)))
            {
                mCandidate.emplace_back(match[i]);
            }
        }
        else{
            if(v1label == query_.GetVertexLabel(u1)){
                if(this->LabelFilter(v1, u1) &&this->LabelFilter(v2, u2)) {
                    mCandidate.emplace_back(match[i]);
                }
            }
            else{
                if(this->LabelFilter(v2, u1) &&this->LabelFilter(v1, u2)) {
                    mCandidate.emplace_back(match[i]);
                }
            }
        }
    }
    Print_Time2("Labelfilter ",start);

    start=Get_Time();
    for(auto m:mCandidate){
        uint u1=order_vs_[m][0];
        uint u2=order_vs_[m][1];
        uint u1label= this->query_.GetVertexLabel(u1);
        uint u2label=this->query_.GetVertexLabel(u2);
        uint v1label=this->data_.GetVertexLabel(v1);
        uint v2label=this->data_.GetVertexLabel(v2);
        uint tmin=this->data_.GetEdgeTime(v1,v2);
        float weight=this->data_.GetEdgeWeight(v1,v2);
        const auto & matchOrder=order_vs_[m];
        if(v1label!=v2label){
            if(v1label!=u1label)
            {
                swap(v1,v2);
            }
            //todo
            this->matchCandidate.clear();
            this->match.clear();
            this->matchVertex(v1,INT_MAX,float(0));
            this->matchVertex(v2,tmin,weight);
            //compute suffix array
            for(int i=1;i<query_.NumVertices();i++){
                for(int j=i;j<query_.NumVertices();j++){
                    suffixMax[i]+=qForwardNeighbors[m][j]->GetStarMaxWeight();
                }
            }
            queryLocalIndexs.resize(query_.NumVertices());
            createQueryCandidate(m,v1,v2);
            searchMatches(m,positive);
            std::fill(suffixMax.begin(), suffixMax.end(), 0.0);
            std::fill(isolatedMax.begin(),isolatedMax.end(),-1);
            queryLocalIndexs.resize(0);
            //search()递归
            this->popVertex(v1);
            this->popVertex(v2);
        }
        else{
            for(int i = 0; i < 2; i++){
                    //todo
                    this->matchCandidate.clear();
                    this->matchVertex(v1,INT_MAX,float(0));
                    this->matchVertex(v2,tmin,weight);
                    //compute suffix array
                   /* for(int i=1;i<query_.NumVertices();i++){
                        for(int j=i;j<query_.NumVertices();j++){
                            suffixMax[i]+=qForwardNeighbors[m][j]->GetStarMaxWeight();
                        }
                    }*/
                    queryLocalIndexs.resize(query_.NumVertices());
                    createQueryCandidate(m,v1,v2);
                    searchMatches(m,positive);
                    std::fill(suffixMax.begin(), suffixMax.end(), 0.0);
                    std::fill(isolatedMax.begin(),isolatedMax.end(),-1);
                    queryLocalIndexs.resize(0);
                    this->popVertex(v2);
                   this->popVertex(v1);
              //  }
                std::swap(v1,v2);// round 2 need
            }
        }
    }
    Print_Time2("SearchMatches ", start);

  /*  start=Get_Time();
    if (max_num_results_ == 0) return;
    for (uint i = 0; i < query_.NumEdges(); i++) {
        uint u1 = order_vs_[i][0], u2 = order_vs_[i][1];
        auto temp_q_labels = query_.GetEdgeLabel(u1, u2);
        float density_s = weight;
        // check if any query edge match (v1 --> v2)
        if (
                std::get<0>(temp_q_labels) == data_.GetVertexLabel(v1) &&
                std::get<1>(temp_q_labels) == data_.GetVertexLabel(v2) &&
                std::get<2>(temp_q_labels) == label
                ) {
            m[u1] = v1;
            m[u2] = v2;
            visited_[v1] = true;
            visited_[v2] = true;
            uint tm = data_.GetEdgeTime(s, t);//最小时间戳首先等于插入的时间戳
            FindMatches(0,i, 2, m, num_results, density_s, tm);
            visited_[v1] = false;
            visited_[v2] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;
            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) goto END_ENUMERATION;
        }
        // check if any query edge match (v2 --> v1)
        if (
                std::get<0>(temp_q_labels) == data_.GetVertexLabel(v2) &&
                std::get<1>(temp_q_labels) == data_.GetVertexLabel(v1) &&
                std::get<2>(temp_q_labels) == label
                ) {
            m[u1] = v2;
            m[u2] = v1;
            visited_[v2] = true;
            visited_[v1] = true;

            uint tm = data_.GetEdgeTime(s, t);
            FindMatches(0,i, 2, m, num_results, density_s, tm);
            visited_[v2] = false;
            visited_[v1] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;

            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) goto END_ENUMERATION;
        }
    }
    END_ENUMERATION:
    num_positive_results_ += num_results;
    Print_Time2("SearchMatches ", start);*/
    start=Get_Time();
    updateTopK(num_results);
    Print_Time2("PrintTopk ", start);
}

//删除边 删除allmatch中tmin=td的记录
void Graphflow::deleteEdge(uint v1, uint v2) {
    std::chrono::high_resolution_clock::time_point start;
    start=Get_Time();
    //首先拿到v1 v2之间的tmin
    uint tmin=data_.GetEdgeTime(v1,v2);
    uint delete_num=0;
    uint cnt=0;
    stringstream _ss;
#ifdef LOG_TRACK
    /*    _ss<<"delete_num:"<<delete_num<<", delete edgeMaps["<<s<<","<<t<<"]"<<" tmin:"<<tmin<<"\n";
    Log::track(_ss);
    _ss.clear();
    _ss.str("");
    for(auto it:edgeMaps[std::make_pair(s,t)]){
        _ss<<"address:"<<it;
       _ss<<" density: "<<it->getDensity()<<" tmin:"<<it->getTmin()<<" vetexs:";
       for(auto i:it->getVetex()){
           _ss<<i<<" ";
       }
       _ss<<'\n';
       Log::track(_ss);
        _ss.clear();
        _ss.str("");
    }*/

#endif
    //删除对象 删除的时候边诱导的记录集合一定是空的
    /*for(auto it:edgeMaps[std::make_pair(s,t)]){
        delete it;
    }
    edgeMaps.erase(std::make_pair(s,t));*/
    //所有的edgeMaps里面tmin=（v1.v2)的tmin的都要被删除

    /*_ss<<"delete in topKSet when tmin =td"<<"\n";
    Log::track1(_ss);
    _ss.clear();
    _ss.str("");*/
    //统计topk 中包含的tmin=td的个数
    for(auto it=topKSet.begin();it!=topKSet.end();it++){
        if((*it)->getTmin()==tmin){
           /* _ss<<" topkSet address:"<<(*it);
            _ss<<"density: "<<(*it)->getDensity()<<" tmin:"<<(*it)->getTmin()<<" vetexs:";
            vector<uint>*ms=(*it)->getVetex();
            for(auto i:*ms){
                _ss<<i<<" ";
            }
            _ss<<'\n';*/
            cnt++;
            /*Log::track1(_ss);
            _ss.clear();
            _ss.str("");*/
        }
    }

    //删除所有记录中过期的匹配
    for(auto it=allMatchRecords.begin();it!=allMatchRecords.end();){
        if((*it)->getTmin()==tmin){
            /*_ss<<"allMatchRecords address:"<<(*it);
            _ss<<"density: "<<(*it)->getDensity()<<" tmin:"<<(*it)->getTmin()<<" vetexs:";
            vector<uint>*ms=(*it)->getVetex();
            for(auto i:*ms){
                _ss<<i<<" ";
            }
            _ss<<'\n';*/
            delete (*it);
            it=allMatchRecords.erase(it);
            delete_num++;
          /*  Log::track1(_ss);
            _ss.clear();
            _ss.str("");*/
        }
        else{
            it++;
        }
    }
    num_negative_results_+=delete_num;
    //std::cout<<"num_negative_results_"<<num_negative_results_<<endl;
    data_.RemoveEdge(v1, v2);

    Print_Time("deleteEdge  ", start);
    if(cnt ==0)
        return;
    //根据cnt的个数需要填补cnt个到Top k 否则就重新截取top k;
    start=Get_Time();
    deleteUpdateTopK();
    Print_Time("deleteUpdateTopK ", start);
}

//删除问题 需要对比补的是不是在topk中
void Graphflow::deleteUpdateTopK() {
    //直接将allmatch中的前k个赋值给topkSet；
    topKSet.clear();
    int n=allMatchRecords.size();
    if(n<k)
    {
        topKSet.resize(n);
        copy(allMatchRecords.begin(),allMatchRecords.end(),topKSet.begin());
    }else{
        topKSet.resize(k);
        copy(allMatchRecords.begin(),allMatchRecords.begin()+k,topKSet.begin());
    }

#ifdef LOG_TRACK
    stringstream _ss;
    _ss<<"Top k"<<std::endl;
    for(auto d:topKSet) {
        _ss << "address " << d;
        _ss << " density:" << d->getDensity() << " tmin:"
            << d->getTmin()
            << " vetexs:";
        std::vector<uint>* vs = d->getVetex();
        for (int j = 0; j < (*vs).size(); j++) {
            _ss << (*vs)[j] << " ";
        }
        _ss << std::endl;
        Log::track1(_ss);
        _ss.clear();
        _ss.str("");
    }
#endif

#ifdef RESULT_TRACK
  stringstream _ss;
 _ss<<"after delete"<<std::endl;
 //std::sort(topKSet.begin(),topKSet.end(), matchRecordCmp);
    for(auto d:topKSet){
        _ss<<d->toString();
        Log::track2(_ss);
        _ss.clear();
        _ss.str("");
    }
#endif

}

//动态地减边操作
void Graphflow::RemoveEdge(uint v1, uint v2) {
    std::vector<uint> m(query_.NumVertices(), UNMATCHED);

    std::tuple<uint, uint, uint> labels = data_.GetEdgeLabel(v1, v2);
    float weight = data_.GetEdgeWeight(v1, v2);

    size_t num_results = 0ul;
    if (max_num_results_ == 0) goto END_ENUMERATION;

    for (uint i = 0; i < query_.NumEdges(); i++) {
        uint u1 = order_vs_[i][1], u2 = order_csrs_[i][0];
        auto temp_q_labels = query_.GetEdgeLabel(u1, u2);
        float density_s;
        if (
                std::get<0>(temp_q_labels) == std::get<0>(labels) &&
                std::get<1>(temp_q_labels) == std::get<1>(labels) &&
                std::get<2>(temp_q_labels) == std::get<2>(labels)
                ) {
            m[u1] = v1;
            m[u2] = v2;
            visited_[v1] = true;
            visited_[v2] = true;
            uint s = std::min(v1, v2);
            uint t = std::max(v1, v2);
            uint tmin = data_.GetEdgeTime(s, t);
            FindMatches(1,i, 2, m, num_results, weight, tmin);

            visited_[v1] = false;
            visited_[v2] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;

            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) goto END_ENUMERATION;
        }
        if (
                std::get<1>(temp_q_labels) == std::get<0>(labels) &&
                std::get<0>(temp_q_labels) == std::get<1>(labels) &&
                std::get<2>(temp_q_labels) == std::get<2>(labels)
                ) {
            m[u1] = v2;
            m[u2] = v1;
            visited_[v2] = true;
            visited_[v1] = true;
            uint s = std::min(v1, v2);
            uint t = std::max(v1, v2);
            uint tmin = data_.GetEdgeTime(s, t);
            FindMatches(1,i, 2, m, num_results, weight, tmin);

            visited_[v2] = false;
            visited_[v1] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;

            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) goto END_ENUMERATION;
        }
    }
    END_ENUMERATION:

    num_negative_results_ += num_results;
    data_.RemoveEdge(v1, v2);
}

void Graphflow::AddVertex(uint id, uint label) {
    data_.AddVertex(id, label);

    visited_.resize(id + 1, false);
}

void Graphflow::RemoveVertex(uint id) {
    data_.RemoveVertex(id);
}

void Graphflow::GetMemoryCost(size_t &num_edges, size_t &num_vertices) {
    num_edges = 0ul;
    num_vertices = 0ul;
}
bool Graphflow::LabelFilter(uint data_v, uint query_v) {
    uint dataVlabelNum= this->data_.NumVLabels();
    uint dataElabelNum=this->data_.NumELabels();
    uint queryVlabelNum=this->query_.NumVLabels();
    uint queryElabelNum=this->query_.NumELabels();
    const auto & dataLabelIndex= this->data_.labelIndex[data_v];
    const auto & queryLabelIndex=this->query_.labelIndex[query_v];

    for(int i=0;i<queryVlabelNum;i++){
        if(i<dataVlabelNum){
            if(dataLabelIndex[i]<queryLabelIndex[i])
            {
                return false;
            }
        }
        else{
            if(queryLabelIndex[i]>0){
                return false;
            }
        }
    }
    for(int i=0;i<queryElabelNum;i++){
        if(i<dataElabelNum){
            if(dataLabelIndex[dataVlabelNum+i]<queryLabelIndex[queryVlabelNum+i]){
                return false;
            }
        }
        else{
            if(queryLabelIndex[queryVlabelNum+i]>0)
                return false;
        }
    }
    return true;
}


void Graphflow::matchVertex(uint data_v,int tmin,float density) {
    this->match.emplace_back(std::make_tuple(data_v,tmin,density));
    this->matchCandidate.push_back({std::make_tuple(data_v,tmin,density)});
    //this->matchVertexCandidate.push_back({});
    this->visited_[data_v]= true;
}


void Graphflow::matchVertex(uint data_v,int index,uint depth,std::vector<tuple<int,int,float>>&singleVertexCandidate) {
    int len=match.size();
    int preindex=0;
    for(int i=len-1;i>=0;i--){
        if(std::get<0>(match[i])!=-1){
            preindex=i;
            break;
        }
    }
    float weight=std::get<2>(match[preindex])+std::get<2>(singleVertexCandidate[index]);
    uint tmin=std::min(std::get<1>(match[preindex]),std::get<1>(singleVertexCandidate[index]));
    this->match.emplace_back(data_v,tmin,weight);
    this->matchCandidate.emplace_back(singleVertexCandidate);
    //this->matchVertexCandidate.push_back({});
    this->visited_[data_v]= true;
}

void Graphflow::matchVertex(std::vector<tuple<int,int,float>>&singleVertexCandidate) {
    this->match.emplace_back(-1,-1,-1);
    this->matchCandidate.emplace_back(singleVertexCandidate);
   // this->matchVertexCandidate.push_back(vertexCandidate);
}

void Graphflow::popVertex(uint data_v) {
    this->match.pop_back();
    //this->matchVertexCandidate.pop_back();
    this->matchCandidate.pop_back();
    this->visited_[data_v]=false;
}
void Graphflow::popVertex(uint data_v, int index, uint depth) {
    this->match.pop_back();
    //this->matchVertexCandidate.pop_back();
    this->matchCandidate.pop_back();
    this->visited_[data_v]=false;
}

void Graphflow::popVertex() {
    this->match.pop_back();
    //this->matchVertexCandidate.pop_back();
    this->matchCandidate.pop_back();
}
void Graphflow::densityFilter(uint matchorder_index,uint depth,std::vector<tuple<int,int,float>>&singleVertexCandidate) {
    uint n=query_.NumVertices();
    if(topKSet.size()<k){
        return;
    }
    float kw=topKSet.back()->getDensity();
    float sumWeight=0;
    //2.利用singleVertex剪枝
    bool flag= false;
    for(int i=depth-1;i>=0;i--) {
        if (!flag&&std::get<0>(match[i]) != -1) {
            sumWeight+=std::get<2>(this->match[i]);
            flag= true;
        }
        else if(std::get<0>(match[i]) == -1){
#ifndef PRINT_DEBUG
            if(isolatedMax[i]==-1){
             std::cout<<"isolatedMax is -1 "<<i<<endl;
            }
#endif
            sumWeight+=isolatedMax[i];
        }
    }
/*
    float isolateWeight=0;
    for(int i=0;i<isolatedVertexs.size();i++){
        isolateWeight+=std::get<2>(matchCandidate[isolatedVertexs[i]].front());
    }*/
    float backWeight=GetBackWeight(matchorder_index,depth+1);
  // sumWeight+=suffixMax[depth+1];
   sumWeight+=backWeight;
    int cnt=0;
    auto iter=singleVertexCandidate.begin();
    while(iter!=singleVertexCandidate.end()){
        float tmpweight=(sumWeight+std::get<2>(*iter))/(sqrt(n)*(n-1));
        if(tmpweight<kw){
            std::get<0>(*iter)=-1;
            cnt++;
        }
        iter++;
    }
    if(cnt==0)
        return;
    int newLen=singleVertexCandidate.size()-cnt;
    int svcLen=singleVertexCandidate.size();
    if(newLen==0)
    {
        singleVertexCandidate.resize(0);
        return;
    }
    else{
        int i=0;
        int j=0;
        while(j<svcLen){
            if(std::get<0>(singleVertexCandidate[j])!=-1){
                singleVertexCandidate[i]=singleVertexCandidate[j];
                i++;
            }
                j++;

        }
        singleVertexCandidate.resize(newLen);
    }
}

void Graphflow::combination_helper(std::vector<std::vector<int>>& result, std::vector<int>& current, const std::vector<std::vector<uint>>& nums, int k) {
    if (k == nums.size()) { // 如果已经处理完了所有的数组，将当前的组合加入结果中
        result.push_back(current);
        return;
    }
    for (int i = 0; i < nums[k].size(); i++) { // 从第 k+1 个数组中选出一个数
        current.push_back(nums[k][i]);
        combination_helper(result, current, nums, k + 1); // 递归处理后面的数组
        current.pop_back();
    }
}

std::vector<std::vector<int>>  Graphflow::combination(const std::vector<std::vector<uint>>& nums) {
    std::chrono::high_resolution_clock::time_point start;
   // start=Get_Time();
    std::vector<std::vector<int>> result;
    std::vector<int> current;
    combination_helper(result, current, nums, 0);
    return result;
   // Print_Time2("combination ",start);
}
bool Graphflow::LDVertexCandidateCheck(uint vertex, uint queryVertexLabel, const std::vector<uint> &needToIntersection,
                                       std::vector<uint> &intersectresult) {
    if(this->visited_[vertex]==true){
        return false;
    }
    std::vector<uint>result;
    const auto &vNeighbors=this->data_.GetNeighbors(vertex);
    if(!needToIntersection.empty()){
        std::set_intersection(vNeighbors.begin(),vNeighbors.end(),needToIntersection.begin(),needToIntersection.end(),std::insert_iterator<std::vector<uint>>(result,result.begin()));
        if(result.empty()){
            return false;
        }
    }
    else{
        for(int i = 0; i < vNeighbors.size(); i++){
            if(this->visited_[vNeighbors[i]] == false && this->data_.GetVertexLabel(vNeighbors[i]) == queryVertexLabel){
                result.push_back(vNeighbors[i]);
            }
        }
        if(result.empty())
            return false;
    }
    if(intersectresult.size()==0){
        intersectresult=result;
        return true;
    }
        std::vector<uint>intersectionCopy;
        intersectionCopy.reserve(intersectresult.size());
        std::swap(intersectionCopy,intersectresult);
        std::set_intersection(intersectionCopy.begin(),intersectionCopy.end(),result.begin(),result.end(),std::insert_iterator<std::vector<uint>>(intersectresult,intersectresult.begin()));
        if(intersectresult.empty())
        {
            return false;
        }
        return true;
}

void Graphflow::setSingleVertexByIntersectionResult(std::vector<tuple<int, int, float>> &singleVertex,
                                                    std::vector<uint> &intersectresult,std::vector<int>&r) {
    //intersectresult候选解，r是选择的LD顶点的候选解，singleVertex是对应的原来candidate候选边
    //利用intersectresult同步更新singleVertex
    sychronizeSingleVertexAndCandidate(singleVertex,intersectresult);
    //1.增加权值
    for(int i=0;i<singleVertex.size();i++){
        uint v1=std::get<0>(singleVertex[i]);
        float sumweight=0;
        int tmin=std::get<1>(singleVertex[i]);
        for(int j=0;j<r.size();j++){
            uint v2=r[j];
            sumweight+=this->data_.GetEdgeWeight(v1,v2);
            tmin=std::min(tmin,(int)this->data_.GetEdgeTime(v1,v2));
        }
        std::get<2>(singleVertex[i])=std::get<2>(singleVertex[i])+sumweight;
        std::get<1>(singleVertex[i])=tmin;
    }

}
void Graphflow::setLDVertexMatchResult(std::vector<int>&r,std::vector<uint>&LDVertexs) {
     for(int i=0;i<LDVertexs.size();i++){
         const auto candidate=this->matchCandidate[LDVertexs[i]];
         for(auto item:candidate){
             if(std::get<0>(item)==r[i]){
                 uint pre_index=LDVertexs[i]-1;
                 while(pre_index>0){
                     auto preItem=this->match[pre_index];
                     if(std::get<0>(preItem)!=-1){
                         int tmin=std::min(std::get<1>(preItem),std::get<1>(item));
                         float density=std::get<2>(item)+std::get<2>(preItem);
                         this->match[LDVertexs[i]]= make_tuple(r[i],tmin,density);
                         for(int k=LDVertexs[i]+1;k<match.size();k++){
                             //更新LD节点之后的密度和
                             if(std::get<0>(this->match[k])!=-1){
                                 std::get<2>(this->match[k])=std::get<2>(this->match[k])+std::get<2>(item);
                             }
                         }
                         break;
                     }else{
                         pre_index--;
//                         this->match[LDVertexs[i]]= make_tuple(r[i],std::get<1>(item),std::get<2>(item));
                     }
                 }
                 break;
             }
         }

     }
}
void Graphflow::setIsolateVertexMatchResult(std::vector<int> &r, std::vector<int> &isolateVertex,uint tmin,float density) {
    for(int i=0;i<isolateVertex.size();i++){
        uint index=isolateVertex[i];

        std::get<0>(this->match[index])=r[i];
        }
    int depth=this->match.size();
    std::get<1>(this->match[depth-1])=tmin;
    std::get<2>(this->match[depth-1])=density;

    }
void Graphflow::setIsolateVertexMatchResult(std::vector<pair<int, int>> &r, std::vector<int> &isolateVertex, uint tmin,
                                            float density) {
    for(int i=0;i<isolateVertex.size();i++){
        uint index=isolateVertex[i];

        std::get<0>(this->match[index])=r[i].first;
    }
    int depth=this->match.size();
    std::get<1>(this->match[depth-1])=tmin;
    std::get<2>(this->match[depth-1])=density;
}

void Graphflow::setBatchVisited(std::vector<int> &r,bool flag) {
    for(auto item:r){
      /*  if(item==61070){
            std::cout<<"setBatchVisited 61070: "<<flag<<std::endl;
        }*/
        this->visited_[item]= flag;
    }
}
void Graphflow::recoverLDVertexMatchResult(std::vector<int>&r,std::vector<uint> &LDVertexs) {
    for(int i=0;i<LDVertexs.size();i++)
    this->match[LDVertexs[i]]= make_tuple(-1,-1,-1);

    for(int i=0;i<LDVertexs.size();i++){
        const auto candidate=this->matchCandidate[LDVertexs[i]];
        for(auto item:candidate){
            if(std::get<0>(item)==r[i]){
                for(int k=LDVertexs[i]+1;k<match.size();k++){
                            //更新LD节点之后的密度和
                            if(std::get<0>(this->match[k])!=-1){
                                std::get<2>(this->match[k])=std::get<2>(this->match[k])-std::get<2>(item);
                            }
                        }
                break;
                }

            }
        }

    }

void Graphflow::recoverIsolateVertexMatchResult(std::vector<int> &IsolateVertexs) {
    for(int i=0;i<IsolateVertexs.size();i++)
        this->match[IsolateVertexs[i]]= make_tuple(-1,-1,-1);
}
void Graphflow::sychronizeSingleVertexAndCandidate(std::vector<tuple<int, int, float>> &singleVertex,
                                                   std::vector<uint> &intersectresult) {
    std::sort(singleVertex.begin(),singleVertex.end(), tupleVertexIdCmp2);
    if(singleVertex.size()==0){
        //todo
        for(auto item:intersectresult){
            singleVertex.push_back(make_tuple(item,INT_MAX,0));
        }
        return;
    }
    auto iter1=singleVertex.begin();
    auto iter2=intersectresult.begin();
    while(iter1!=singleVertex.end()){
        if(std::get<0>((*iter1))==(*iter2))
        {
            iter1++;
            iter2++;
        }else{
            iter1=singleVertex.erase(iter1);
        }
    }

}
void Graphflow::addMatchResult(uint matchorderindex, searchType type) {
    std::chrono::high_resolution_clock::time_point starttime;
    int n=query_.NumVertices();
    auto &isolateVertexs=query_.isolatedRecord[matchorderindex];
    std::vector<std::vector<tuple<int,int,float>>>combinezIsolateVertexs;
    combinezIsolateVertexs.reserve(isolateVertexs.size());
    for(int i=0;i<isolateVertexs.size();i++){
        std::vector<tuple<int,int,float>>singleVertex=matchCandidate[isolateVertexs[i]];
        std::sort(singleVertex.begin(),singleVertex.end(), tupleSumWeightCmp);
        combinezIsolateVertexs.emplace_back(singleVertex);
    }

    int len=combinezIsolateVertexs.size();
    //数组根据标签删减数据
   /* float sum=1;
    std::vector<int>vertexNum;
    for(int i=0;i<len;i++){
        vertexNum.push_back(combinezIsolateVertexs[i].size());
        sum*=combinezIsolateVertexs[i].size();
    }
    if(sum>k){
        for(int i=0;i<len;i++)
        {
            sum=1;
            bool isCut= true;
            std::vector<int>samlabelwithId;
            int id=order_vs_[matchorderindex][isolateVertexs[i]];
            auto &samlabelVertex=sameLabelVertex[id];
            for(int j=0;j<len;j++){
                if(j!=i){
                    int otherId=order_vs_[matchorderindex][isolateVertexs[j]];
                    if(std::find(samlabelVertex.begin(), samlabelVertex.end(),otherId)==samlabelVertex.end()){
                        sum*=vertexNum[j];
                    }
                    else{
                        samlabelwithId.push_back(vertexNum[j]);
                    }
                }
            }
            std::sort(samlabelwithId.begin(),samlabelwithId.end());
            int len=samlabelwithId.size();
            for(int k=0;k<len;k++){
                int tmp=samlabelwithId[k]-(k+1);
                if(tmp<0)
                {
                    isCut= false;
                    break;
                }
                sum*=(samlabelwithId[k]-(k+1));
            }
            if(!isCut){
                continue;
            }
            int num_bound= ceil(k*1.0/sum);
            if(num_bound<combinezIsolateVertexs[i].size())
            {
                combinezIsolateVertexs[i].resize(num_bound);
            }
        }
    }*/
#ifdef PRINT_DEBUG
    std::vector<int>testm(this->query_.NumVertices(),-1);
    for(int i=0;i<n;i++){
        if(std::get<0>(match[i])!=-1){
            testm[order_vs_[matchorderindex][i]]=std::get<0>(match[i]);
        }
    }
    std::cout<<"matchorderindex "<<matchorderindex<<endl;
    for(int i=0;i<n;i++){
        if(testm[i]!=-1){
           std::cout<<" m["<<i<<"]= "<<testm[i]<<" ";
        }
    }
    std::cout<<"Candidate"<<endl;
    for(int i:isolateVertexs){
        int id=order_vs_[matchorderindex][i];
        std::cout<<id<<" candidate: ";
        for(int j=0;j<matchCandidate[i].size();j++){
            int candi=std::get<0>(matchCandidate[i][j]);
            std::cout<<candi<<" ";
        }
        std::cout<<endl;
    }
   std::cout<<endl;
#endif
    auto wt=findWeightAndTminBeforeIsolated();
    int*pmax=new int[len]{};//最大pmax在组中的列号
    std::unordered_map<int,int>visitedId;
    float Tmax=wt.second;
    int Tmin=wt.first;
    bool TmaxisAllFirst=true;
   for(int i=0;i<len;i++){
        int index=0;
        int id=std::get<0>(combinezIsolateVertexs[i][index]);
        bool flag= true;
        while(visited_[id]){
            combinezIsolateVertexs[i].erase(combinezIsolateVertexs[i].begin()+index);
            if(index==combinezIsolateVertexs[i].size())
            {
                flag= false;
                break;
            }
            id=std::get<0>(combinezIsolateVertexs[i][index]);
        }
        if(!flag)
        {
            return;
        }
        if(!visitedId.count(id)){
            visitedId[id]=i;
        }
        else{
            TmaxisAllFirst= false;
            int pre_group_id=visitedId[id];
            float gap1=std::get<2>(combinezIsolateVertexs[pre_group_id][pmax[pre_group_id]])-std::get<2>(combinezIsolateVertexs[pre_group_id][pmax[pre_group_id]+1]);
            float gap2=std::get<2>(combinezIsolateVertexs[i][pmax[i]])-std::get<2>(combinezIsolateVertexs[i][pmax[i]+1]);
            if(gap1<gap2){
                pmax[pre_group_id]++;
                visitedId[id]=i;
                int preId=std::get<0>(combinezIsolateVertexs[pre_group_id][pmax[pre_group_id]]);
                visitedId[preId]=pre_group_id;
            }
            else{
                pmax[i]++;
            }
            int curId=std::get<0>(combinezIsolateVertexs[i][pmax[i]]);
            visitedId[curId]=i;
        }
    }
    float globalTmax=Tmax;
    std::vector<int>r(isolateVertexs.size());
    for(int i=0;i<len;i++){
        r[i]=std::get<0>(combinezIsolateVertexs[i][pmax[i]]);
        Tmax+=std::get<2>(combinezIsolateVertexs[i][pmax[i]]);
        Tmin=std::min(Tmin,std::get<1>(combinezIsolateVertexs[i][pmax[i]]));
    }
   if(TmaxisAllFirst){
       //1.
            globalTmax=Tmax;
           float density=Tmax/(sqrt(n)*(n-1));
           setIsolateVertexMatchResult(r,isolateVertexs,Tmin,Tmax);
           std::vector<uint>m(n);
           for(int i=0;i<match.size();i++){
               m[order_vs_[matchorderindex][i]]=std::get<0>(match[i]);
           }
           MatchRecord *record=new MatchRecord(density,Tmin,m);
           bool matchResult= addMatchRecords(record);
           if(matchResult){
               if(type==positive)
                   num_positive_results_++;
               else
                   num_negative_results_++;
           }
           recoverIsolateVertexMatchResult(isolateVertexs);
           if(topKSet.size()==k){
               float cur_density=topKSet.back()->getDensity();
               float density=Tmax/(sqrt(n)*(n-1));
               if(cur_density>density)
               {
                   delete[]pmax;
                   return;
               }
           }
   }
   else{
       for(int i=0;i<len;i++){
           globalTmax+=std::get<2>(combinezIsolateVertexs[i][0]);
       }
   }

   //2.选择下一个迭代的点
    int*hash=new int[len]{};
    float *Tbound=new float[len];
    std::fill(Tbound,Tbound+len,Tmax);

    int next_vertex_group=0;
    bool isover=false;
    int *noscan=new int[len]{};
    while(!isover){
        int tmpcnt=0;
        next_vertex_group= findTboundMaxIndex(Tbound,hash,noscan,combinezIsolateVertexs,len);

        if(isnoNextVertex(noscan,len))
        {
          isover= true;
            break;
        }

        hash[next_vertex_group]++;
        float pre=std::get<2>(combinezIsolateVertexs[next_vertex_group][pmax[next_vertex_group]]);
        float cur=std::get<2>(combinezIsolateVertexs[next_vertex_group][hash[next_vertex_group]]);
        Tbound[next_vertex_group]=globalTmax-pre+cur;
        float TboundNext=Tbound[next_vertex_group];
        if(TboundNext>Tmax)
            Tmax=Tbound[next_vertex_group];
        auto next_item=combinezIsolateVertexs[next_vertex_group][hash[next_vertex_group]];

        //3. 递归
        int id=std::get<0>(next_item);
        CatesianProductWithIndex(matchorderindex,type,next_vertex_group,0,len,hash,combinezIsolateVertexs,isolateVertexs,wt.first,wt.second);
        if(topKSet.size()==k){
            float cur_density=topKSet.back()->getDensity();
            float density=Tmax/(sqrt(n)*(n-1));
            if(cur_density>density)
            {
               isover= true;
                break;
            }
            float TboundDensity=TboundNext/ (sqrt(n)*(n-1));
            if(cur_density>TboundDensity){
                noscan[next_vertex_group]=1;
                for(int i=0;i<len;i++)
                    std::cout<<noscan[i]<<" ";
                  std::cout<<endl;
            }
        }

    }
       delete[]pmax;
       delete[]hash;
       delete[]Tbound;
       delete[]noscan;



    //starttime=Get_Time();

    //初始化
   /* auto wt=findWeightAndTminBeforeIsolated();
    std::vector<pair<int,int>>matchresult;
    for(int j=0;j<combinezIsolateVertexs.size();j++){
        matchresult.push_back(make_pair(std::get<0>(combinezIsolateVertexs[j][0]),0));//id号，以及在该列中的索引位置
    }
    //说明没有重复
    std::queue<std::vector<pair<int,int>>>que;
     que.push(matchresult);
        while(!que.empty()) {
            float weight = 0;
            int tmin = INT_MAX;
            auto ms = que.front();
            que.pop();
            bool isRepeat= false;
            //判断是否有已经访问过节点
            for(int i=0;i<ms.size();i++){
                int id=ms[i].first;
                int index=ms[i].second;
                if(visited_[id]==true){
                    if(index+1==combinezIsolateVertexs[i].size())
                        return;
                    int newId=std::get<0>(combinezIsolateVertexs[i][index+1]);
                    ms[i].first=newId;
                    ms[i].second=index+1;
                    isRepeat= true;
                }
            }
            if(isRepeat){
                que.push(ms);
                continue;
            }
            std::set<int> isrepeat;
            for (auto item: ms) {
                isrepeat.insert(item.first);
            }
            if (ms.size() == isrepeat.size()) {
                for (int j = 0; j < ms.size(); j++) {
                    int index = ms[j].second;
                    weight += std::get<2>(combinezIsolateVertexs[j][index]);
                    tmin = std::min(tmin, std::get<1>(combinezIsolateVertexs[j][index]));
                }
                tmin = std::min(tmin, wt.first);
                float t = (sqrt(n) * (n - 1));
                float density = (weight + wt.second) / t;
                setIsolateVertexMatchResult(ms, isolateVertexs, tmin, density);
                //加入结果集
                std::vector<uint> m(n);
                for (int i = 0; i < match.size(); i++) {
                    m[order_vs_[matchorderindex][i]] = std::get<0>(match[i]);
                }
                MatchRecord *r = new MatchRecord(density, tmin, m);
                if(!addMatchRecords(r))
                {
                    recoverIsolateVertexMatchResult(isolateVertexs);
                    return;
                }
                else{
                    num_positive_results_++;
                }
                recoverIsolateVertexMatchResult(isolateVertexs);
                //找到下一行中的差值最小的替换点
                std::map<float, vector<std::pair<int, int>>> minWeight;//列号，索引位置
                float diffmin = INT_MAX;
                for (int i = 0; i < ms.size(); i++) {
                    auto index = ms[i].second;
                    if (index == combinezIsolateVertexs[i].size() - 1) {
                        continue;
                    }
                    float d = std::get<2>(combinezIsolateVertexs[i][index]) -
                              std::get<2>(combinezIsolateVertexs[i][index+1]);
                    diffmin = std::min(diffmin, d);
                    minWeight[d].emplace_back(i, index + 1);
                }
                for (auto item: minWeight[diffmin]) {
                    std::vector<pair<int, int>> mscopy = ms;
                    int id = std::get<0>(combinezIsolateVertexs[item.first][item.second]);
                    ms[item.first] = make_pair(id, item.second);
                    que.push(ms);
                    ms = mscopy;
                }
            }
            else {
                //有重复的，首先找到各个重复的key，并将diff差值较小的替换，较大的留下，然后进行下一次循环
                std::map<float, vector<std::pair<int, int>>> minWeight;//列号，索引位置
                std::unordered_map<int,vector<pair<int,int>>>sumId;//id号，（列号，索引号）
                std::unordered_map<int,float>idMaxDiff;//id号，最大的差值
                for(int i=0;i<ms.size();i++){
                    int id=ms[i].first;
                    int index=ms[i].second;
                    sumId[id].emplace_back(i,index);
                    if(idMaxDiff.find(id)!=idMaxDiff.end()){
                        float d = std::get<2>(combinezIsolateVertexs[i][index]) -
                                  std::get<2>(combinezIsolateVertexs[i][index+1]);
                        idMaxDiff[id]=std::max(d,idMaxDiff[id]);
                    }
                    else{
                        if (index == combinezIsolateVertexs[i].size() - 1) {
                          idMaxDiff[id]=-1;
                        }
                        else{
                            float d = std::get<2>(combinezIsolateVertexs[i][index]) -
                                      std::get<2>(combinezIsolateVertexs[i][index+1]);
                            idMaxDiff[id]=d;
                        }

                    }
                }
                for(auto iter=sumId.begin();iter!=sumId.end();iter++){
                    if((*iter).second.size()>1){
                        int vertexId=(*iter).first;
                        int maxdiff=idMaxDiff[vertexId];
                        if(maxdiff==-1){
                            return;
                        }
                        auto ids=(*iter).second;
                        for(auto item:ids){
                            float d=std::get<2>(combinezIsolateVertexs[item.first][item.second]) -
                                    std::get<2>(combinezIsolateVertexs[item.first][item.second+1]);
                            if(d<maxdiff){
                             ms[item.first]= make_pair(vertexId,item.second+1);
                            }
                        }

                    }
                }
                que.push(ms);
            }

        }*/


    //先剪枝再笛卡尔积
    /*starttime=Get_Time();

    std::vector<float>colmax;
    if(topKSet.size()==k){
        float topkdensity=topKSet.back()->getDensity();
        for(int i=0;i<combinezIsolateVertexs.size();i++){
            int max=0;
            for(int j=0;j<combinezIsolateVertexs.size();j++){
                if(j==i)
                    continue;
                max+=std::get<2>(combinezIsolateVertexs[j][0]);
            }
            colmax.emplace_back(max);
        }
        for(int i=0;i<combinezIsolateVertexs.size();i++){
            for(int j=1;j<combinezIsolateVertexs[i].size();j++){
                float weight=wt.second+std::get<2>(combinezIsolateVertexs[i][j])+colmax[i];
                float density=weight/(sqrt(n)*(n-1));
                if(density<topkdensity)
                    combinezIsolateVertexs[i].erase(combinezIsolateVertexs[i].begin()+j,combinezIsolateVertexs[i].end());
            }
        }
    }

    std::vector<tuple<std::vector<int>,int,float>>result= combinationMatchResult(combinezIsolateVertexs);
    sort(result.begin(),result.end(), tupleResultCmp);
    Print_Time2("combinationMatchResult ",starttime);
    bool flag=true;
    bool isVisted= false;
    int i;
    while(flag){
        for(i=0;i<result.size();i++){
            for(int j=0;j<isolateVertexs.size();j++){
                uint v=std::get<0>(result[i])[j];
                if(visited_[v]==true){
                    isVisted= true;
                    break;
                }
            }
            if(isVisted){
                continue;
            }
            //result中的元素可能重复
            std::vector<int>tmpresult=std::get<0>(result[i]);
            std::set<int>isrepeat (tmpresult.begin(),tmpresult.end());
            if(isrepeat.size()!=tmpresult.size())
            {
                continue;
            }

            //更新match record
            //find weight/tmin before

            uint tmin=std::min(std::get<1>(result[i]),wt.first);
            float weight=std::get<2>(result[i]);
            float t=(sqrt(n)*(n-1));
            float density=(weight+wt.second)/ (sqrt(n)*(n-1));
            setIsolateVertexMatchResult(std::get<0>(result[i]),isolateVertexs,tmin,density);
            //设置为已访问
            setBatchVisited(std::get<0>(result[i]), true);
            //加入结果集
            std::vector<uint>m(n);
            for(int i=0;i<match.size();i++){
                m[order_vs_[matchorderindex][i]]=std::get<0>(match[i]);
            }
            MatchRecord *r=new MatchRecord(density,tmin,m);
            bool matchResult= addMatchRecords(r);
            if(!matchResult){
                flag= false;
                recoverIsolateVertexMatchResult(isolateVertexs);
                setBatchVisited(std::get<0>(result[i]), false);
                break;
            }
            else{
                if(type==positive){
                    num_positive_results_++;
                }
                else if(type==negative){
                    num_negative_results_++;
                }
            }
            //还原isolatevertex和Visited
            recoverIsolateVertexMatchResult(isolateVertexs);
            setBatchVisited(std::get<0>(result[i]), false);
        }
        if(i==result.size()){
            flag= false;
        }
    }*/

}
std::vector<tuple<std::vector<int>,int,float>> Graphflow::combinationMatchResult(std::vector<std::vector<tuple<int, int, float>>> combinezIsolateVertexs) {
    std::chrono::high_resolution_clock::time_point start;
   // start=Get_Time();
    std::vector<tuple<std::vector<int>, int, float>> result;
    std::vector<int>current;
    combinationMatchResultHelp(result,current,combinezIsolateVertexs,0,INT_MAX,0);
  //  Print_Time2("combinationMatchResult ",start);
    return result;
}
void Graphflow::combinationMatchResultHelp(std::vector<tuple<std::vector<int>, int, float>> &result,
                                           std::vector<int> &current,
                                           std::vector<std::vector<tuple<int, int, float>>> &combinezIsolateVertexs,
                                           int k,int tmin,float density)
                                           {
    if(k==combinezIsolateVertexs.size()){
        result.push_back(make_tuple(current,tmin,density));
        return;
    }
    for(int i=0;i<combinezIsolateVertexs[k].size();i++){

        int copytmin=tmin;
        float copydensity=density;
        tmin=std::min(std::get<1>(combinezIsolateVertexs[k][i]),tmin);
        density+=std::get<2>(combinezIsolateVertexs[k][i]);
        current.push_back(std::get<0>(combinezIsolateVertexs[k][i]));
        combinationMatchResultHelp(result,current,combinezIsolateVertexs,k+1,tmin,density);
        current.pop_back();
        tmin=copytmin;
        density=copydensity;
    }
}
std::pair<int,float>Graphflow::findWeightAndTminBeforeIsolated() {
    int depth=match.size()-1;
    for(int i=depth;i>=0;i--){
        if(std::get<0>(match[i])==-1){
            continue;
        }
        else{
            float weight=std::get<2>(match[i]);
            int tmin=std::get<1>(match[i]);
            return make_pair(tmin,weight);
        }
    }
}
void Graphflow::CatesianProductWithIndex(int matchorderindex,searchType type,int curIndex,int depth,int len,int*hash,std::vector<std::vector<tuple<int,int,float>>>&combinezIsolateVertexs,std::vector<int>&isolateVertexs,int &tmin,float &weight) {
    if(depth==len){
        int n=query_.NumVertices();
        std::vector<uint>m(n);
        for(int i=0;i<match.size();i++){
            m[order_vs_[matchorderindex][i]]=std::get<0>(match[i]);
        }
        float density=weight/ (sqrt(n)*(n-1));
        MatchRecord *record=new MatchRecord(density,tmin,m);
        bool matchResult= addMatchRecords(record);
        if(matchResult){
            if(type==positive)
                num_positive_results_++;
            else
                num_negative_results_++;
        }
        return;
    }
    if(depth!=curIndex){
        int right=hash[depth];
        for(int i=0;i<=right;i++){
            int id=std::get<0>(combinezIsolateVertexs[depth][i]);
            if(visited_[id])
                continue;
            visited_[id]=true;
            int copytmin=tmin;
            float copyweight=weight;
            std::get<0>(this->match[isolateVertexs[depth]])=id;
            tmin=std::min(tmin,std::get<1>(combinezIsolateVertexs[depth][i]));
            weight+=std::get<2>(combinezIsolateVertexs[depth][i]);
            CatesianProductWithIndex(matchorderindex,type,curIndex,depth+1,len,hash,combinezIsolateVertexs,isolateVertexs,tmin,weight);
            visited_[id]=false;
            std::get<0>(this->match[isolateVertexs[depth]])=-1;
            tmin=copytmin;
            weight=copyweight;
        }
    }else{
        std::cout<<"hash["<<curIndex<<"]"<<hash[curIndex]<<" depth["<<depth<<"]"<<std::get<0>(combinezIsolateVertexs[depth][hash[curIndex]])<<endl;
        int id=std::get<0>(combinezIsolateVertexs[depth][hash[curIndex]]);
        if(visited_[id]){
            return;
        }
        visited_[id]=true;
        int copytmin=tmin;
        float copyweight=weight;
        std::get<0>(this->match[isolateVertexs[depth]])=id;
        tmin=std::min(tmin,std::get<1>(combinezIsolateVertexs[depth][hash[curIndex]]));
        weight+=std::get<2>(combinezIsolateVertexs[depth][hash[curIndex]]);
        CatesianProductWithIndex(matchorderindex, type,curIndex,depth+1,len,hash,combinezIsolateVertexs,isolateVertexs,tmin,weight);
        visited_[id]=false;
        std::get<0>(this->match[isolateVertexs[depth]])=-1;
        tmin=copytmin;
        weight=copyweight;
    }

}
int Graphflow::findTboundMaxIndex(float *Tbound,int*hash,int* noscan,std::vector<std::vector<tuple<int,int,float>>>&combinezIsolateVertexs,int len) {
    float max=0;
    int maxindex=0;
    for(int i=0;i<len;i++){
        int index=hash[i]+1;
        if(noscan[i]==1)
        {
           /* if(i==0)
            {
                if(i+1!=len){
                    max=Tbound[i+1];
                    maxindex=i+1;
                }

            }
            else{
                if(i+1!=len){
                    if(Tbound[i+1]>max)
                    {
                        max=Tbound[i+1];
                        maxindex=i+1;
                    }
                }
            }*/
           /*if(i+1!=len){
               max=Tbound[i+1];
               maxindex=i+1;
           }
*/
            continue;
        }
        if(index>=combinezIsolateVertexs[i].size())
        {
            noscan[i]=1;
            /*if(i==0)
            {
                if(i+1!=len){
                    max=Tbound[i+1];
                    maxindex=i+1;
                }

            }
            else{
                if(i+1!=len){
                    if(Tbound[i+1]>max)
                    {
                        max=Tbound[i+1];
                        maxindex=i+1;
                    }
                }
            }*/
            /*if(i+1!=len){
                max=Tbound[i+1];
                maxindex=i+1;
            }
*/
            continue;
        }
        int id=std::get<0>(combinezIsolateVertexs[i][index]);
        bool flag= true;
        while(visited_[id])
        {
            combinezIsolateVertexs[i].erase(combinezIsolateVertexs[i].begin()+index);
            if(index>=combinezIsolateVertexs[i].size()){
                noscan[i]=1;
                flag= false;
                break;
            }
            else{
                id=std::get<0>(combinezIsolateVertexs[i][index]);
            }
        }
        if(!flag)
            continue;
       /* while(visitedId[id]){
            hash[i]++;
            index=hash[i]+1;
            if(index>=combinezIsolateVertexs[i].size()){
                noscan[i]=1;
                flag= false;
                break;
            }
            else{
                id=std::get<0>(combinezIsolateVertexs[i][index]);
            }
        }
        if(!flag)
            continue;*/
        if(Tbound[i]>max)
        {
            max=Tbound[i];
            maxindex=i;
        }
    }
    return maxindex;
}
bool Graphflow::isnoNextVertex(int *noscan,int len) {
    for(int i=0;i<len;i++){
        if(noscan[i]==0)
            return false;
    }
    return true;
}

void Graphflow::addMatchResultWithHeap(uint matchorderindex, searchType type) {
    std::chrono::high_resolution_clock::time_point starttime;
    int n=query_.NumVertices();
    auto &isolateVertexs=query_.isolatedRecord[matchorderindex];
    std::vector<std::vector<tuple<int,int,float>>>combinezIsolateVertexs;
    combinezIsolateVertexs.reserve(isolateVertexs.size());
    for(int i=0;i<isolateVertexs.size();i++){
        std::vector<tuple<int,int,float>>singleVertex=matchCandidate[isolateVertexs[i]];
       // std::sort(singleVertex.begin(),singleVertex.end(), tupleSumWeightCmp);
        combinezIsolateVertexs.push_back(singleVertex);
    }

    int len=combinezIsolateVertexs.size();
#ifdef PRINT_DEBUG
    std::vector<int>testm(this->query_.NumVertices(),-1);
    for(int i=0;i<n;i++){
        if(std::get<0>(match[i])!=-1){
            testm[order_vs_[matchorderindex][i]]=std::get<0>(match[i]);
        }
    }
    std::cout<<"matchorderindex "<<matchorderindex<<endl;
    for(int i=0;i<n;i++){
        if(testm[i]!=-1){
            std::cout<<" m["<<i<<"]= "<<testm[i]<<" ";
        }
    }
    std::cout<<"Candidate"<<endl;
    for(int i:isolateVertexs){
        int id=order_vs_[matchorderindex][i];
        std::cout<<id<<" candidate: ";
        for(int j=0;j<matchCandidate[i].size();j++){
            int candi=std::get<0>(matchCandidate[i][j]);
            std::cout<<candi<<" ";
        }
        std::cout<<endl;
    }
    std::cout<<endl;
#endif
    auto wt=findWeightAndTminBeforeIsolated();
    int*hash=new int[len]{};//索引指针
    int*noscan=new int[len]{};//记录是否停止
    std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int>>, pairCompare> maxHeap;
    for(int i=0;i<combinezIsolateVertexs.size();i++){
        float weight=std::get<2>(combinezIsolateVertexs[i][0]);
        maxHeap.push(std::make_pair(weight,i));
    }
    //计算全局Tmax//大于 则直接return
    float Tmax=wt.second;
    std::unordered_map<int,int>visitedId;
    int*pmax=new int[len]{};
    for(int i=0;i<len;i++){
        int index=0;
        int id=std::get<0>(combinezIsolateVertexs[i][index]);
        bool flag= true;
        while(visited_[id]){
            combinezIsolateVertexs[i].erase(combinezIsolateVertexs[i].begin()+index);
            if(index==combinezIsolateVertexs[i].size())
            {
                flag= false;
                break;
            }
            id=std::get<0>(combinezIsolateVertexs[i][index]);
        }
        if(!flag)
        {
            return;
        }
        if(!visitedId.count(id)){
            visitedId[id]=i;
        }
        else{
            int pre_group_id=visitedId[id];
            float gap1=std::get<2>(combinezIsolateVertexs[pre_group_id][pmax[pre_group_id]])-std::get<2>(combinezIsolateVertexs[pre_group_id][pmax[pre_group_id]+1]);
            float gap2=std::get<2>(combinezIsolateVertexs[i][pmax[i]])-std::get<2>(combinezIsolateVertexs[i][pmax[i]+1]);
            if(gap1<gap2){
                pmax[pre_group_id]++;
                visitedId[id]=i;
                int preId=std::get<0>(combinezIsolateVertexs[pre_group_id][pmax[pre_group_id]]);
                visitedId[preId]=pre_group_id;
            }
            else{
                pmax[i]++;
            }
            int curId=std::get<0>(combinezIsolateVertexs[i][pmax[i]]);
            visitedId[curId]=i;
        }
    }
    for(int i=0;i<len;i++){
        Tmax+=std::get<2>(combinezIsolateVertexs[i][pmax[i]]);
    }
    if(topKSet.size()==k){
          float density=Tmax/(sqrt(n)*(n-1));
        if(topKSet.back()->getDensity()>density)
            return;
    }
    while(!maxHeap.empty()&&!isnoNextVertex(noscan,len)){
        auto item=maxHeap.top();
        maxHeap.pop();
        int group=item.second;
        int pre= hash[group];
        hash[group]++;
        if(pre==combinezIsolateVertexs[group].size()-1){
            noscan[group]=1;
        }
        else{
            float w=std::get<2>(combinezIsolateVertexs[group][hash[group]]);
            maxHeap.push(std::make_pair(w,group));
        }
        int id=std::get<0>(combinezIsolateVertexs[group][pre]);
        if(visited_[id]){
            combinezIsolateVertexs[group].erase(combinezIsolateVertexs[group].begin()+pre);
        }
        else{
            //笛卡尔积 从最小的seen number 开始
            int minGroup=0;
            int minNum=INT_MAX;
            bool flag=true;
            for(int i=0;i<len;i++){
                if(i==group)
                    continue;
                if(hash[i]==0)
                {
                    flag= false;
                    break;
                }
                else if(hash[i]<minNum){
                    minGroup=i;
                    minNum=hash[group];
                }
            }

            if(!flag)
                continue;
            int depthLen=len-1;
            std::vector<int>copyisolateIndex;
            copyisolateIndex.emplace_back(minGroup);
            for(int i=0;i<len;i++){
                if(i==minGroup||i==group)
                    continue;
                copyisolateIndex.emplace_back(i);
            }
            visited_[id]=true;
            std::get<0>(this->match[isolateVertexs[group]])=id;
            float weight=wt.second+std::get<2>(combinezIsolateVertexs[group][pre]);
            int tmin=std::min(wt.first,std::get<1>(combinezIsolateVertexs[group][pre]));
            float Tmax=0;
            //递归
            CatesianProductWithHeap(matchorderindex,type,0,depthLen,hash,combinezIsolateVertexs,isolateVertexs,copyisolateIndex,tmin,weight,Tmax);
            visited_[id]=false;
            std::get<0>(this->match[isolateVertexs[group]])=-1;
            float curDensity=topKSet.back()->getDensity();
            if(Tmax==0)
                continue;
            if(curDensity>Tmax){
                noscan[group]=1;
            }

        }
    }
    delete[]hash;
    delete[]noscan;
    delete[]pmax;
}
void Graphflow::CatesianProductWithHeap(int matchorderindex, searchType type, int depth, int len, int *hash,
                                        std::vector<std::vector<tuple<int, int, float>>> &combinezIsolateVertexs,
                                        std::vector<int> &isolateVertexs, std::vector<int> &isolatedIndex, int &tmin, float &weight,float &Tmax) {
    if(depth==len){
        int n=query_.NumVertices();
        std::vector<uint>m(n);
        for(int i=0;i<match.size();i++){
            m[order_vs_[matchorderindex][i]]=std::get<0>(match[i]);
        }
        float density=weight/ (sqrt(n)*(n-1));
        if(density>Tmax)
        {
            Tmax=density;
        }
        MatchRecord *record=new MatchRecord(density,tmin,m);
        bool matchResult= addMatchRecords(record);
        if(matchResult){
            if(type==positive)
                num_positive_results_++;
            else
                num_negative_results_++;
        }
        return;
    }
        int group=isolatedIndex[depth];
        int right=hash[group];
        for(int i=0;i<right;i++){
            int id=std::get<0>(combinezIsolateVertexs[group][i]);
            if(visited_[id])
                continue;
            visited_[id]=true;
            int copytmin=tmin;
            float copyweight=weight;
            std::get<0>(this->match[isolateVertexs[group]])=id;
            tmin=std::min(tmin,std::get<1>(combinezIsolateVertexs[group][i]));
            weight+=std::get<2>(combinezIsolateVertexs[group][i]);
            CatesianProductWithHeap(matchorderindex,type,depth+1,len,hash,combinezIsolateVertexs,isolateVertexs,isolatedIndex,tmin,weight,Tmax);
            visited_[id]=false;
            std::get<0>(this->match[isolateVertexs[group]])=-1;
            tmin=copytmin;
            weight=copyweight;
        }
}
void Graphflow::createQueryCandidate(int matchorderindex, uint v1, uint v2) {
    std::vector<uint>& matchOrder= this->order_vs_[matchorderindex];
    uint u1=matchOrder[0];
    uint u2=matchOrder[1];
    std::vector<uint>u1_rn=rightNeighbor[matchorderindex][u1];
    std::vector<uint>u2_rn=rightNeighbor[matchorderindex][u2];
    int len=1;
    //u1_rn
    updateDensityCandidate(matchorderindex,u1,v1,u1_rn);
    updateDensityCandidate(matchorderindex,u2,v2,u2_rn);

}
void Graphflow::updateDensityCandidate(int matchorderindex, uint u1, uint v1, std::vector<uint> &u1_rn) {
    const std::vector<Neighbor> &v1N=data_.vNeighbors[v1];
    for(uint u1_neighbor:u1_rn){
        LocalIndex &ql=queryLocalIndexs[u1_neighbor];
        uint u1_neighbor_index = order_vertex_index[matchorderindex][u1_neighbor];
        auto &candiates= ql.getCandidate();
        if(ql.getIsFirst()){
            ql.setIsFirst(false);
            for(auto neighbor:v1N){
                int neighbor_id=neighbor.getVertexId();
                std::vector<Neighbor>vN=data_.vNeighbors[neighbor_id];
                    if(neighbor.GetEdgelabel() == std::get<2>(query_.GetEdgeLabel(u1,u1_neighbor))&&neighbor.getVertexLabel()==query_.GetVertexLabel(u1_neighbor)){
                        StarGraph* s = qForwardNeighbors[matchorderindex][u1_neighbor_index];
                        std::vector<ForwardNeighbor*>& queryVetex=s->GetqueryVertex();
                       // std::sort(queryVetex.begin(),queryVetex.end(),ForwardNeighborcmp);
                        std::sort(vN.begin(),vN.end(),greater<Neighbor>());
                        int leftvN=0;
                        int rightqV=0;
                        int sumweight=0;
                        int cnt=0;
                        while(leftvN<vN.size()&&rightqV<queryVetex.size()){
                            if(vN[leftvN].GetelabelAndVertexLabel()<queryVetex[rightqV]->GetelabelAndVertexLabel())
                            {
                                break;
                            }
                            int flag=1;
                            while(vN[leftvN].GetelabelAndVertexLabel()>queryVetex[rightqV]->GetelabelAndVertexLabel())
                            {
                                leftvN++;
                                if(leftvN>=vN.size())
                                {
                                    flag=0;
                                    break;
                                }
                            }
                            if(!flag)
                                break;
                            if(vN[leftvN].GetelabelAndVertexLabel()==queryVetex[rightqV]->GetelabelAndVertexLabel()){
                                pair<uint,uint>tmppair=vN[leftvN].GetelabelAndVertexLabel();
                                float edgeweight=vN[leftvN].GetEdgeWeight();
                                sumweight+=edgeweight;
                                rightqV++;
                                leftvN++;
                                cnt++;
                            }
                        }
                        if(cnt==queryVetex.size()){
                            ql.insertCandidate(neighbor_id,sumweight);
                            if(ql.getMaxWeight()==FLT_MAX||sumweight>ql.getMaxWeight())
                            {
                                ql.setMaxWeght(sumweight);
                                ql.setmaxId(neighbor_id);
                            }
                        }
                    }
            }
        }
        else{
            auto candidate=ql.getCandidate();
            //检查maxIndex是否在vN中，在直接返回
            int maxId=ql.getMaxId();
            //不在则判断candidate是否在v1N中，取交集
            /*for(int i=0;i<v1N.size();i++){
                if(v1N[i].getVertexId()==maxId)
                    return;
            }*/
            int maxWeight=0;
            for(int i=0;i<candiates.size();i++){
                int id=candiates[i].first;
                if(id==maxId)
                    return;
                float weight=candiates[i].second;
                bool flag= false;
                if(id==-1)
                {
                    continue;
                }
                for(int i=0;i<v1N.size();i++){
                    if(v1N[i].getVertexId()==id)
                    {
                        flag= true;
                        break;
                    }
                }
                if(flag){
                    if(weight>maxWeight){
                        maxWeight=weight;
                        maxId=id;
                    }
                }
                else{
                    candiates[i].first=-1;
                }
            }
        }


    }
}
void Graphflow::clearDensityCandidate(std::vector<uint> &u1_rn) {
    for(uint u1_neighbor:u1_rn){
        queryLocalIndexs[u1_neighbor].clearCandidate();

    }
}