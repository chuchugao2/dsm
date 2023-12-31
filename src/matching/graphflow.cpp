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


bool ForwardNeighborcmp(ForwardNeighbor*f1,ForwardNeighbor*f2){
    return (*f1)>(*f2);
}
bool areSame(float a, float b) {
    return std::fabs(a - b) < std::numeric_limits<float>::epsilon();
}



bool tupleVertexIdCmp2(const std::tuple<int, int, float>& a, const std::tuple<int, int, float>& b) {
    if (std::get<0>(a) != std::get<0>(b)) {
        return std::get<0>(a) <std::get<0>(b);
    } else {
        return std::get<1>(a) < std::get<1>(b);
    }
}



Graphflow::Graphflow(Graph& query_graph, Graph& data_graph,Subgraph& global_subgraph,
                     uint max_num_results,
                     bool print_prep,
                     bool print_enum,
                     bool homo)
        : matching(query_graph, data_graph, global_subgraph,max_num_results,
                   print_prep, print_enum, homo)
        , order_vs_(query_.NumEdges())
        , order_csrs_(query_.NumEdges())
        , order_offs_(query_.NumEdges())
        ,order_vertex_index(query_.NumEdges())
        ,topKSet(0)
        ,suffixMax(query_graph.NumVertices(),0)
        ,isolatedMax(query_graph.NumVertices(),-1)
        ,rightNeighbor(query_graph.NumEdges())
        ,matchCandidate(query_graph.NumVertices())
        ,match(query_graph.NumVertices())
        ,labelToQueryVertex(query_graph.NumVLabels())
        ,globalVkMatchUk(data_graph.NumVertices())
        ,globalStarIndex(query_.NumEdges())
        ,queryVertexIndexInlabel(query_.NumVertices())
        ,LocalStarIndex(query_.NumVertices())
        ,matchLeftNeighborSum(query_.NumEdges())
        ,leftNeighborIdSum(query_.NumEdges())
        ,matchVetexLeftNeighbor(query_.NumVertices())
        , matchVetexSumweight(query_.NumVertices())
{
//    globalStarIndex.resize(query_.NumEdges());
    for (uint i = 0; i < query_.NumEdges(); ++i)
    {
        order_vs_[i].resize(query_.NumVertices());//节点个数
        order_csrs_[i].resize(query_.NumEdges() + 1);//边的个数+1，
        order_offs_[i].resize(query_.NumVertices(), 0);//节点个数，初始化为0
        globalStarIndex[i].resize(query_.NumVertices());
        order_vertex_index[i].resize(query_.NumVertices());
        rightNeighbor[i].resize(query_.NumVertices());
        matchLeftNeighborSum[i].resize(query_.NumVertices());
        leftNeighborIdSum[i].resize(query_.NumVertices());
    }
}
Graphflow::~Graphflow() noexcept {
    for(MatchRecord* item:topKSet){
        delete item;
        item= nullptr;
    }
    for(int i=0;i<query_.NumEdges();i++){
        for(StarGraph*s:globalStarIndex[i]){
            delete s;
        }
    }
}

void Graphflow::Preprocessing()//预处理过程
{
    this->data_.InitLabelIndex();
    GenerateMatchingOrder();
    this->query_.InitLabelIndex();
    //create globalsubgraph;
    createGlobalSubgraph();

#ifdef LOG_TRACK
    stringstream _ss;
    for(int i=0;i<query_.NumVertices();i++){
        _ss<<i<<"candidate:"<<endl;
        for(auto m:globalsubgraph_.matchCandidate[i])
        {
            _ss<<m<<" ";
        }
        _ss<<endl;
    }
    Log::track1(_ss);
#endif
    this->query_.InitMatchOrderType(this->order_vs_,this->rightNeighbor);
    createLabelToQueryVertex();
    //createStarIndex with globalSubgraph?
    CreateStarIndex();
    std::cout<<"Preprocess end"<<endl;
}
void Graphflow::updateStarIndex(uint match_index,uint caddidate_v,const std::vector<uint>&canditeQueryVertexs) {
    std::vector<int>&result=globalVkMatchUk[caddidate_v][match_index];
    result.resize(canditeQueryVertexs.size());
    std::vector<Neighbor>&vN= this->data_.vNeighbors[caddidate_v];
    std::sort(vN.begin(),vN.end(),greater<Neighbor>());
    //std::sort(queryVetex.begin(),queryVetex.end(), ForwardNeighborcmp);
    float sumWeight=0;
    for(int i=0;i<canditeQueryVertexs.size();i++){
        uint candidate_u=canditeQueryVertexs[i];
        sumWeight=0;
        int vertex_index=order_vertex_index[match_index][candidate_u];
        if(vertex_index==0)
            continue;
        StarGraph* s=globalStarIndex[match_index][vertex_index];
        std::vector<ForwardNeighbor*>&queryVetex=s->GetqueryVertex();
        /* if(candidate_u==4&&caddidate_v==309819){
             for(int i=0;i<queryVetex.size();i++){
                 auto item=queryVetex[i]->GetelabelAndVertexLabel();
                 std::cout<<"elabel:"<<item.first<<" vlabel:"<<item.second<<endl;
             }
             std::cout<<caddidate_v<<"neighbors"<<endl;
             for(auto v:vN){
                 auto item=v.GetelabelAndVertexLabel();
                 std::cout<<"elabel:"<<item.first<<" vlabel:"<<item.second<<" toVertexId:"<<v.getVertexId()<<endl;
             }
             exit(-1);
         }*/

        int leftvN=0;
        int rightqV=0;
        //StarGraph *tmpgraph=new StarGraph();
        int flag=1;
        int qvSize=queryVetex.size();
        int vNsize=vN.size();
        while(leftvN<vNsize&&rightqV<qvSize){
            if(vN[leftvN].GetelabelAndVertexLabel()<queryVetex[rightqV]->GetelabelAndVertexLabel())
            {
                flag=0;
                break;
            }
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
//                uint tovId=queryVetex[rightqV]->GetVetexId();
//                pair<uint,uint>tmppair=vN[leftvN].GetelabelAndVertexLabel();
                float edgeweight=vN[leftvN].GetEdgeWeight();
                sumWeight+=edgeweight;
                //ForwardNeighbor* f=new ForwardNeighbor(tovId,tmppair.second,tmppair.first,edgeweight);
                //tmpgraph->AddForwardNeighbor(f);
                rightqV++;
                leftvN++;
            }
        }
        if(!flag)
        {
            result[i]=s->getStarMaxWeight();
        }
        if(rightqV==qvSize){
            result[i]=sumWeight;
            if(s->getStarMaxWeight()==queryVetex.size()*mw||s->getStarMaxWeight()<sumWeight)
            {
                s->setStarMaxWeight(sumWeight);
            }
        }
        //tmpgraph->computeMaxWeight();
        //find tmpgraph
        /*  if(tmpgraph->GetForwardNeighborNum()==queryVetex.size()){
              if(s->GetStarMaxWeight()==queryVetex.size()*mw||s->GetStarMaxWeight()<tmpgraph->GetStarMaxWeight())
              {
                  StarGraph * t=s;
                  globalStarIndex[match_index][vertex_index]=tmpgraph;
                  delete t;
              }
              else {
                  delete tmpgraph;
              }
          }*/
    }

}
void Graphflow::updateStarIndex(uint match_index, uint caddidate_v, uint candidate_u,int candidate_v_index) {
    std::vector<int>&result=globalVkMatchUk[caddidate_v][match_index];
    int vertex_index=order_vertex_index[match_index][candidate_u];
    if(vertex_index==0)
    {
        return;
    }
    StarGraph* s=globalStarIndex[match_index][vertex_index];
    const std::vector<ForwardNeighbor*>&queryVetex=s->GetqueryVertex();
    std::vector<Neighbor>&vN= this->globalsubgraph_.vNeighbors[caddidate_v];
   // std::sort(vN.begin(),vN.end(),greater<Neighbor>());
    int leftvN=0;
    int rightqV=0;
    int vNsize=vN.size();
    int qVsize=queryVetex.size();
    int flag=1;
    float sumweight=0;
   // std::vector<uint>MatchId;
    while(leftvN<vNsize&&rightqV<qVsize){
        if(vN[leftvN].GetelabelAndVertexLabel()<queryVetex[rightqV]->GetelabelAndVertexLabel())
        {
            flag=0;
            break;
        }

        while(vN[leftvN].GetelabelAndVertexLabel()>queryVetex[rightqV]->GetelabelAndVertexLabel()||vN[leftvN].getMatchQueryVertexId()!=queryVetex[rightqV]->GetVetexId()||vN[leftvN].getfromVertexId()!=candidate_u)
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
           // MatchId.emplace_back(vN[leftvN].getVertexId());
            float edgeweight=vN[leftvN].GetEdgeWeight();
            rightqV++;
            leftvN++;
            sumweight+=edgeweight;
        }
    }
    if(!flag)
    {
        return;
    }
    else if(rightqV==qVsize){
        //globalVkMatchUk更新
        result[candidate_v_index]=sumweight;
        // globalStarIndex更新
        if(s->getStarMaxWeight()==queryVetex.size()*mw||s->getStarMaxWeight()<sumweight)
        {
            s->setStarMaxWeight(sumweight);
            s->setMatchDataVertexId(caddidate_v);
        }
    }
}
float Graphflow::GetBackWeight(uint order_index,uint depth) {
    float sum=0;
    uint n=query_.NumVertices();
    std::vector<uint>& matchOrder= this->order_vs_[order_index];
    for(int i=depth;i<n;i++){
        sum+=LocalStarIndex[i];
    }
    return sum;
}
void Graphflow::CreateStarIndex() {
    int m=query_.NumEdges();
    int n=query_.NumVertices();
    for(int i=0;i<n;i++){
        const std::vector<uint>&candidate=globalsubgraph_.matchCandidate[i];
        const std::vector<uint>&candidate_u=labelToQueryVertex[query_.GetVertexLabel(i)];
        for(uint v:candidate){
            globalVkMatchUk[v].resize(m);
            for(int j=0;j<m;j++){
                globalVkMatchUk[v][j].resize(candidate_u.size());
                if(i==order_vs_[j][0])
                    continue;
                int candidate_index=queryVertexIndexInlabel[i];
                updateStarIndex(j,v,i,candidate_index);
            }

        }
        }
    }
  /*  for(int i=0;i<n;i++){
        int label=data_.GetVertexLabel(i);
        if(label==-1)
            continue;
        const std::vector<uint>&candidate_u=labelToQueryVertex[label];
        if(candidate_u.size()==0)
            continue;
        globalVkMatchUk[i].resize(m);
        for(int j=0;j<m;j++){
            updateStarIndex(j,i,candidate_u);
        }
}*/
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
                //  std::cout<<"globalStarIndex[0]["<<i<<"] other:"<<other<<" other label"<< this->query_.GetVertexLabel(other)<<endl;
                ForwardNeighbor* forwardNeighbor=new ForwardNeighbor(other, this->query_.GetVertexLabel(other),qlabel);
                s->AddForwardNeighbor(forwardNeighbor);
                order_csrs_[0][order_offs_[0][i]++] = other;
            }
        }
        s->InitalmaxWeight();
        globalStarIndex[0][i]=s;

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
        globalStarIndex[i][1]=(s);

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
                    // std::cout<<"globalStarIndex["<<i<<"]"<<"["<<j<<"] "<<"other:"<<other<<" other label"<< this->query_.GetVertexLabel(other)<<endl;
                    order_csrs_[i][order_offs_[i][j]++] = other;
                    qlabel=std::get<2>(this->query_.GetEdgeLabel(max_adjacent_u,other));
                    ForwardNeighbor* forwardNeighbor=new ForwardNeighbor(other, this->query_.GetVertexLabel(other),qlabel);
                    s->AddForwardNeighbor(forwardNeighbor);
                }
            }
            s->InitalmaxWeight();
            globalStarIndex[i][j]=(s);
        }
    }
    //对globalStarIndex中的所有匹配序的所有节点的前向邻居按照<el,vl>排序
    for(int i=0;i<query_.NumEdges();i++){
        for(int j=1;j<query_.NumVertices();j++) {
            StarGraph *s = globalStarIndex[i][j];
            std::vector<ForwardNeighbor*> &globalIndex = s->GetqueryVertex();
            std::sort(globalIndex.begin(), globalIndex.end(), ForwardNeighborcmp);
        }
    }


    size_t sumId=0;
    //创建所有节点的右邻居数组
    if (print_preprocessing_results_)
    {
        std::cout << "matching order: " << std::endl;
        std::cout << "-vertex(backward neighbors)-\n";
        for (uint i = 0; i < query_.NumEdges(); ++i)
        {
            std::cout << "#" << i << ": ";
            for (uint j = 0; j < query_.NumVertices(); ++j)
            {
                sumId=0;
                std::cout << order_vs_[i][j];
                if (j == 0)
                {
                    // this->query_.forwardNeighbors[i][j]={};
                    std::cout << "-";
                    continue;
                }
                std::vector<ForwardNeighbor>currentQueryNeighbors;
                matchLeftNeighborSum[i][j]=order_offs_[i][j]-order_offs_[i][j-1];
                for (uint k = order_offs_[i][j - 1]; k < order_offs_[i][j]; k++)
                {
                    uint toVertexId=order_csrs_[i][k];
                    uint toVertexIndex=order_vertex_index[i][toVertexId];
                    uint toVertexLabel=query_.GetVertexLabel(order_csrs_[i][k]);
                    uint edgelabel=std::get<2>(query_.GetEdgeLabel(order_vs_[i][j],toVertexId));
                    ForwardNeighbor f(toVertexIndex,toVertexId,toVertexLabel,edgelabel);
                    currentQueryNeighbors.push_back(f);
                    rightNeighbor[i][toVertexId].emplace_back(order_vs_[i][j]);
                    if (k == order_offs_[i][j - 1]) std::cout << "(";
                    std::cout << order_csrs_[i][k];
                    if (k != order_offs_[i][j] - 1) std::cout << ",";
                    else std::cout << ")";
                }
                leftNeighborIdSum[i][j]=sumId;
                //  this->query_.forwardNeighbors[i][j]=currentQueryNeighbors;
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

void Graphflow::updateTopK() {
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
    if(!isUpdateIntopkset){
        return;
    }
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

void Graphflow::searchMatches(int depth,uint matchorderindex, searchType flag)  {
    //1.找到前向邻居，找到候选解
    std::vector<uint>& matchOrder= this->order_vs_[matchorderindex];
    uint queryVertex=matchOrder[depth];
    uint queryVertexLabel=this->query_.GetVertexLabel(queryVertex);
    // const auto & fneighbors=this->query_.forwardNeighbors[matchorderindex][depth];
    //vertexType currentSearchVertextype = this->query_.GetVertexType(matchorderindex,depth);
    std::vector<SingleCandidate>& singleVertexCandidate=this->matchCandidate[depth];
    std::vector<SingleCandidate>copySingleVertexCandidate=this->matchCandidate[depth];
    getIntersetSingleCandidate(singleVertexCandidate,matchorderindex,depth);
    if(singleVertexCandidate.size()==0){
        this->matchCandidate[depth]=copySingleVertexCandidate;
        return;
    }

    total_densityFilter_time.StartTimer();
    densityFilter(matchorderindex,depth,singleVertexCandidate);
    total_densityFilter_time.StopTimer();
    if(singleVertexCandidate.size()==0)
    {
        this->matchCandidate[depth]=copySingleVertexCandidate;
        return;
    }
    if(isInsert)
        IsearchSpace+=singleVertexCandidate.size();
    else
        DsearchSpace+=singleVertexCandidate.size();
    //Print_Time2("densityFilter ",start);
    //顺序扩展
    if(depth==query_.NumVertices()-1){
        //add matchresult;
        std::sort(singleVertexCandidate.begin(),singleVertexCandidate.end());
        for(const SingleCandidate &single:singleVertexCandidate){
            uint dataV=single.getVertexId();
            if(visited_[dataV])
                continue;
            float sumWeight= this->match[depth-1].getSumWeight();
            this->match[depth].setVertexId(dataV);
            sumWeight+=single.getSumWeight();
            this->match[depth].setSumWeight(sumWeight);
            int n=query_.NumVertices();
            std::vector<uint>m(n);
            for(int i=0;i<match.size();i++){
                m[order_vs_[matchorderindex][i]]=match[i].getVertexId();
            }
            float density=sumWeight/ (sqrt(n)*(n-1));

            MatchRecord *record=new MatchRecord(density,m);
            int matchResult= addMatchRecords(record);
            allMatchFind++;
            if(matchResult==1){
                if(flag==positive)
                {
                    num_positive_results_++;
                    numAddTopk++;
                }
            }
            else if(matchResult==3){
                this->matchCandidate[depth]=copySingleVertexCandidate;
                this->match[depth].clearSingleCandidate();
                return;
            }
        }
        //clear candidate;
        this->matchCandidate[depth]=copySingleVertexCandidate;
        this->match[depth].clearSingleCandidate();
        return;
    }
    else{
        for(int i=0;i<singleVertexCandidate.size();i++){
            uint dataV=singleVertexCandidate[i].getVertexId();
            if(visited_[dataV])
                continue;
            float weight=singleVertexCandidate[i].getSumWeight();
            //递归
            matchVertex(0,depth,dataV,weight);

            const std::vector<uint>&uk_neighbor=rightNeighbor[matchorderindex][queryVertex];
            std::vector<std::vector<SingleCandidate>>copyCandidate(uk_neighbor.size());
            std::vector<int>copyLocalStarIndex(query_.NumVertices());
            for(int i=0;i<uk_neighbor.size();i++){
                int uk_neighbor_index=order_vertex_index[matchorderindex][uk_neighbor[i]];
                copyCandidate[i]=matchCandidate[uk_neighbor_index];
                copyLocalStarIndex[uk_neighbor_index]=LocalStarIndex[uk_neighbor_index];
            }
            // std::cout<<"depth :"<<depth<<" data:"<<dataV<<endl;
            total_updaterightNeighborCandidate_time.StartTimer();
            bool isNull=updaterightNeighborCandidate(matchorderindex,queryVertex,0, false,dataV,uk_neighbor);
            total_updaterightNeighborCandidate_time.StopTimer();
            if(isNull)
            {
                for(int i=0;i<uk_neighbor.size();i++){
                    int uk_neighbor_index=order_vertex_index[matchorderindex][uk_neighbor[i]];
                    matchCandidate[uk_neighbor_index]= copyCandidate[i];
                    LocalStarIndex[uk_neighbor_index]=copyLocalStarIndex[uk_neighbor_index];
                }
                this->visited_[dataV]= false;
                continue;
            }
            //copy SingleCandidate
            //updateweight;
            searchMatches(depth+1,matchorderindex,flag);
            //返回candidate状态
            for(int i=0;i<uk_neighbor.size();i++){
                int uk_neighbor_index=order_vertex_index[matchorderindex][uk_neighbor[i]];
                matchCandidate[uk_neighbor_index]= copyCandidate[i];
                LocalStarIndex[uk_neighbor_index]=copyLocalStarIndex[uk_neighbor_index];
            }
            this->visited_[dataV]= false;
        }
        this->matchCandidate[depth]=copySingleVertexCandidate;
        this->match[depth].clearSingleCandidate();
    }


    //带孤立节点
   /* if(currentSearchVertextype==freeVertex){
        //  start=Get_Time();
        //密度剪枝
        for(int i=0;i<singleVertexCandidate.size();i++){
            uint dataV=singleVertexCandidate[i].getVertexId();
            if(visited_[dataV])
                continue;
            float weight=singleVertexCandidate[i].getSumWeight();
            //递归
            matchVertex(0,depth,dataV,weight);

            const std::vector<uint>&uk_neighbor=rightNeighbor[matchorderindex][queryVertex];
            std::vector<std::vector<SingleCandidate>>copyCandidate(uk_neighbor.size());
            std::vector<int>copyLocalStarIndex(query_.NumVertices());
            for(int i=0;i<uk_neighbor.size();i++){
                int uk_neighbor_index=order_vertex_index[matchorderindex][uk_neighbor[i]];
                copyCandidate[i]=matchCandidate[uk_neighbor_index];
                copyLocalStarIndex[uk_neighbor_index]=LocalStarIndex[uk_neighbor_index];
            }
            // std::cout<<"depth :"<<depth<<" data:"<<dataV<<endl;
            std::chrono::high_resolution_clock::time_point start;
            start=Get_Time();
            bool isNull=updaterightNeighborCandidate(matchorderindex,queryVertex,0, false,dataV,uk_neighbor);
            total_update_localIndex_time+= Duration2(start);
            if(isNull)
            {
                for(int i=0;i<uk_neighbor.size();i++){
                    int uk_neighbor_index=order_vertex_index[matchorderindex][uk_neighbor[i]];
                    matchCandidate[uk_neighbor_index]= copyCandidate[i];
                    LocalStarIndex[uk_neighbor_index]=copyLocalStarIndex[uk_neighbor_index];
                }
                this->visited_[dataV]= false;
                continue;
            }
            //copy SingleCandidate
            //updateweight;
            searchMatches(depth+1,matchorderindex,flag);
            //返回candidate状态
            for(int i=0;i<uk_neighbor.size();i++){
                int uk_neighbor_index=order_vertex_index[matchorderindex][uk_neighbor[i]];
                matchCandidate[uk_neighbor_index]= copyCandidate[i];
                LocalStarIndex[uk_neighbor_index]=copyLocalStarIndex[uk_neighbor_index];
            }
            this->visited_[dataV]= false;
        }
        this->matchCandidate[depth]=copySingleVertexCandidate;
        this->popVertex(depth);
    }
    else{
        isolatedMax[depth]=maxWeight;
        this->matchVertex(depth);
        // 若递归到终点，释放所有的孤立节点
        if(depth==this->query_.NumVertices()-1){
            //todo
            start=Get_Time();
            //addMatchResult(matchorderindex,flag);
            addMatchResultWithHeap(matchorderindex,flag);
            //Print_Time2("addMatchResult ",start);
            total_addMatchResult_time+= Duration2(start);
        }
        else{
            searchMatches(depth+1,matchorderindex,flag);
        }
        this->matchCandidate[depth]=copySingleVertexCandidate;
        this->popVertex(depth);

    }*/
}



//flag==0为initial flag=1为update
void Graphflow::FindMatches(uint flag,uint order_index, uint depth, std::vector<uint> m, size_t &num_results, float density_s) {
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
            float lastds = density_s / (sqrt(m.size()) * (m.size() - 1));
            //sort(m.begin(),m.end());
            num_results++;
            MatchRecord *r = new MatchRecord(lastds, m);
            addMatchRecords(r);


            if (print_enumeration_results_) {
                for (auto j: m) {

                    std::cout << j << " ";
                }
            }

        } else {
            size_t num_results_before_recursion = num_results;
            FindMatches(flag,order_index, depth + 1, m, num_results, density_s);
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


int Graphflow:: addMatchRecords(MatchRecord *r) {
    //sort(r->getVetex()->begin(),r->getVetex()->end());
    int n=topKSet.size();
    if(n<k){
        for(int j=n-1;j>=0;j--){
            if((*topKSet[j])==(*r)){
                delete r;
                return 2;
            }
        }
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
        isUpdateIntopkset= true;
        return 1;
    }
    else{
        for(int j=n-1;j>=0;j--){
            if((*topKSet[j])==(*r)){
                delete r;
                return 2;
            }
        }
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
            isUpdateIntopkset=true;
            return 1;
        }
        else{
            delete r;
            return 3;
        }
        //return true;
    }
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

                    FindMatches(flag,order_index, depth, m, num_initial_results_,density_s);
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
                uint tmp;
                std::vector<uint>m;
                ifs2>>density;
                for(int i=0;i<query_.NumVertices();i++){
                    ifs2>>tmp;
                    m.emplace_back(tmp);
                }
                MatchRecord* matchRecord=new MatchRecord(density,m);
                addMatchRecords(matchRecord);
            }
        }
    }
}

//动态的加边操作
void Graphflow::AddEdge(uint v1, uint v2, uint label, float weight) {
#ifdef LOG_TRACK
     stringstream _ss;

    for(auto n:globalsubgraph_.vNeighbors[63117]){
        if(n.getVertexId()==1709&&n.getMatchQueryVertexId()==0){
            _ss<<"find 1709"<<endl;
        }
    }
    Log::track1(_ss);
#endif
    total_update_globalIndex_time.StartTimer();
    int numAddTopk=0;
    allMatchFind=0;
    bool flag= true;
    data_.AddEdge(v1, v2, label, weight,  1);
    this->data_.UpdateLabelIndex(v1,v2,label,1);
    uint v1label=data_.GetVertexLabel(v1);
    uint v2label=data_.GetVertexLabel(v2);
     vector<int>match=EdgeisInMatchOrder(v1,v2,v1label,v2label,label);
    if(match.size()==0){
        total_update_globalIndex_time.StartTimer();
        return;
    }
    //update globalsubgraph and starIndex
    bool isInGlobalSubgraph=updateGlobalSubgraph(v1,v2,label,weight,match);;
    if(!isInGlobalSubgraph){
        total_update_globalIndex_time.StopTimer();
        return;
    }
    total_update_globalIndex_time.StopTimer();
#ifdef LOG_TRACK
    stringstream _ss1;
    for(auto n:globalsubgraph_.matchCandidate[1]){
        if(n==62501)
        {
            _ss<<"find 62501"<<endl;
        }
    }
#endif
    total_search_time.StartTimer();
    isUpdateIntopkset= false;
    for(auto m:match){
        uint u1=order_vs_[m][0];
        uint u2=order_vs_[m][1];
        uint u1label= this->query_.GetVertexLabel(u1);
        uint u2label=this->query_.GetVertexLabel(u2);
        uint v1label=this->data_.GetVertexLabel(v1);
        uint v2label=this->data_.GetVertexLabel(v2);
        float weight=this->data_.GetEdgeWeight(v1,v2);
        InitialLocalIndex(m);
        if(v1label!=v2label){
            if(v1label!=u1label)
            {
                swap(v1,v2);
            }
            //todo
           flag=SearchMatchesWithEdge(m,v1,v2,weight,u1,u2,positive);
            if(flag)
                continue;
        }
        else{
            for(int i = 0; i < 2; i++) {
              flag= SearchMatchesWithEdge(m,v1,v2,weight,u1,u2,positive);
                if(flag)
                {
                    std::swap(v1, v2);
                    continue;
                }
                std::swap(v1, v2);
            }
        }
    }
    total_search_time.StopTimer();
    //Print_Time2("SearchMatches ", start);
    END_ENUMERATION:
    total_print_time.StartTimer();
    sumAllMatchFind+=allMatchFind;
    std::cout<<"num add top k:"<<numAddTopk<<endl;
    std::cout<<"all match find:"<<allMatchFind<<endl;
    updateTopK();
    total_print_time.StopTimer();
    //Print_Time2("PrintTopk ", start);
}

//删除边 删除allmatch中tmin=td的记录
void Graphflow::deleteEdge(uint v1, uint v2) {
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
    /* //删除对象 删除的时候边诱导的记录集合一定是空的
     //统计topk 中包含的tmin=td的个数
     for(auto it=topKSet.begin();it!=topKSet.end();it++){
         if((*it)->getTmin()==tmin){

         }
     }

     //删除所有记录中过期的匹配
     for(auto it=allMatchRecords.begin();it!=allMatchRecords.end();){
         if((*it)->getTmin()==tmin){

             delete (*it);
             it=allMatchRecords.erase(it);
             delete_num++;
         }
         else{
             it++;
         }
     }*/
    num_negative_results_+=delete_num;
    //std::cout<<"num_negative_results_"<<num_negative_results_<<endl;
    data_.RemoveEdge(1,v1, v2);

    //Print_Time("deleteEdge  ", start);
    if(cnt ==0)
        return;
    //根据cnt的个数需要填补cnt个到Top k 否则就重新截取top k;
    deleteUpdateTopK();
    //Print_Time("deleteUpdateTopK ", start);
}

//删除问题 需要对比补的是不是在topk中
void Graphflow::deleteUpdateTopK() {

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
void Graphflow::RemoveEdge(uint v1, uint v2,uint label) {
    //1.更新data、更新nlf标签
    total_delete_update_time.StartTimer();
    allMatchFind=0;
    uint v1label=data_.GetVertexLabel(v1);
    uint v2label=data_.GetVertexLabel(v2);
    float weight=data_.GetEdgeWeight(v1,v2);
    data_.RemoveEdge(0,v1, v2);
    data_.UpdateLabelIndex(v1,v2,label,0);
    bool flag;
    //2.更新subgraph中的点和边
    vector<int>match=EdgeisInMatchOrder(v1,v2,v1label,v2label,label);
    //3.删除边并且更新索引
    deleteGlobalSubgraph(v1,v2,label,weight,match);
    total_delete_update_time.StopTimer();
#ifdef LOG_TRACK
 /*   stringstream _ss;
    if(globalStarIndex[0][7]->getStarMaxWeight()==1){
        _ss<<"max=1"<<endl;
    }
    Log::track1(_ss);*/

#endif
    //4.判断topk中是否包含删除边
    total_delete_time.StartTimer();
     flag=deleteMatchRecordWithEdge(v1,v1label,v2,v2label,label,match);
    if(!flag)
    {
        total_delete_time.StartTimer();
        return;
    }
    //5 从subgraph中的Edge进行重搜
    //从候选最少的节点，以及它的邻居候选次少的节点进行重搜
    int n=query_.NumVertices();
    uint minVertexSize=UINT_MAX;
    uint minVertex=0;
    for(int i=0;i<n;i++){
        int s=globalsubgraph_.matchCandidate[i].size();
        if(s<minVertexSize){
            minVertexSize=s;
            minVertex=i;
        }
    }
    const std::vector<uint>&minNeighbor=query_.GetNeighbors(minVertex);
    uint minVertexNeighorSize=UINT_MAX;
    uint minVertexNeighor=0;
    for(uint u:minNeighbor){
        int us=globalsubgraph_.matchCandidate[u].size();
        if(us<minVertexNeighorSize){
            minVertexNeighorSize=us;
            minVertexNeighor=u;
        }
    }
    int m=query_.NumEdges();
    int matchIndex=0;
    for(int i=0;i<m;i++){
        if((order_vs_[i][0]==minVertex&&order_vs_[i][1]==minVertexNeighor)||((order_vs_[i][1]==minVertex&&order_vs_[i][0]==minVertexNeighor))){
            matchIndex=i;
            break;
        }
    }
    const std::vector<uint>&minVertexCandidate=globalsubgraph_.matchCandidate[minVertex];
    for(uint mv:minVertexCandidate){
        const std::vector<Neighbor>&neighbors=globalsubgraph_.GetVNeighbors(mv);
        for(const Neighbor&n:neighbors){
            if(n.getMatchQueryVertexId()==minVertexNeighor&&n.getfromVertexId()==minVertex){
                //search
                uint u1=order_vs_[matchIndex][0];
                uint u2=order_vs_[matchIndex][1];
                uint u1label= this->query_.GetVertexLabel(u1);
                uint u2label=this->query_.GetVertexLabel(u2);
                uint v1=mv;
                uint v2=n.getVertexId();
                uint v1label=this->data_.GetVertexLabel(mv);
                uint v2label=n.getVertexLabel();
                float weight=n.GetEdgeWeight();
                InitialLocalIndex(matchIndex);

                if(v1label!=v2label){
                    if(v1label!=u1label)
                    {
                        swap(v1,v2);
                    }
                    //todo
                    SearchMatchesWithEdge(matchIndex,v1,v2,weight,u1,u2,negative);

                }
                else{
                    for(int i = 0; i < 2; i++) {
                        SearchMatchesWithEdge(matchIndex,v1,v2,weight,u1,u2,negative);
                        std::swap(v1, v2);
                    }
                }
            }
        }
    }
    total_delete_time.StopTimer();

    sumDeleteallMatchFind+=allMatchFind;
    std::cout<<"delete research matches:"<<allMatchFind<<endl;
    deleteUpdateTopK();
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




void Graphflow::matchVertex(bool isFirstEdge,uint depth,uint data_v,float w) {
    if(isFirstEdge){
        this->match[depth].setVertexId(data_v);
        this->match[depth].setSumWeight(w);
        this->matchCandidate[depth].emplace_back(SingleCandidate(data_v,w));
        //this->matchVertexCandidate.push_back({});

        this->visited_[data_v]= true;
    }
    else{
        int preindex=0;
        for(int i=depth-1;i>=0;i--){
            if(match[i].getVertexId()!=-1){
                preindex=i;
                break;
            }
        }
        float weight=match[preindex].getSumWeight()+w;
        this->match[depth].setVertexId(data_v);
        this->match[depth].setSumWeight(weight);
        //this->matchCandidate.emplace_back(singleVertexCandidate);
        //this->matchVertexCandidate.push_back({});

        this->visited_[data_v]= true;
    }

}

void Graphflow::matchVertex(int depth) {
    this->match[depth].setIsolateSingleCandate();
    //this->matchCandidate.emplace_back(singleVertexCandidate);
    // this->matchVertexCandidate.push_back(vertexCandidate);
}

void Graphflow::popVertex(uint depth,uint data_v) {
    this->match[depth].clearSingleCandidate();
    //this->matchVertexCandidate.pop_back();
    this->matchCandidate[depth].clear();

    this->visited_[data_v]=false;
}
void Graphflow::popVertex(uint data_v, uint matchorderindex, uint depth, const std::vector<uint>&uk_neighbor) {
    this->match[depth].clearSingleCandidate();
    const int n=uk_neighbor.size();
    for(int u_id:uk_neighbor){
        int query_order_index=order_vertex_index[matchorderindex][u_id];
        matchCandidate[query_order_index].clear();
    }
    //this->matchVertexCandidate.pop_back();
    this->matchCandidate[depth].clear();

    this->visited_[data_v]=false;
}

void Graphflow::popVertex(uint depth) {
    this->match[depth].clearSingleCandidate();
    // this->matchCandidate[i].clear();
}
void Graphflow::densityFilter(uint matchorder_index,uint depth,std::vector<SingleCandidate>&singleVertexCandidate) {
    uint n=query_.NumVertices();
    if(topKSet.size()<k){
        auto iter=singleVertexCandidate.begin();
        while(iter!=singleVertexCandidate.end()) {
            float md = (*iter).getSumWeight();
            iter++;
        }
        return;
    }
    float kw=topKSet.back()->getDensity();
    float sumWeight=0;
    //2.利用singleVertex剪枝
    bool flag= false;
    //顺序扩展
    sumWeight+=this->match[depth-1].getSumWeight();
    //延迟扩展
/*    for(int i=depth-1;i>=0;i--) {
        if (!flag&&(match[i]).getVertexId() != -1) {
            sumWeight+=this->match[i].getSumWeight();
            flag= true;
        }
        else if(match[i].getVertexId() == -1){
#ifdef PRINT_DEBUG
            if(isolatedMax[i]==-1){
             std::cout<<"isolatedMax is -1 "<<i<<endl;
            }
#endif
            sumWeight+=isolatedMax[i];
        }
    }*/

    float backWeight=GetBackWeight(matchorder_index,depth+1);
    sumWeight+=backWeight;
    /* sumWeight+=suffixMax[depth+1];*/
    int cnt=0;

    // std::cout<<"sumWeight: "<<sumWeight<<endl;
    auto iter=singleVertexCandidate.begin();
    while(iter!=singleVertexCandidate.end()){
        float md=(*iter).getSumWeight();
        float tmpweight=(sumWeight+md)/(sqrt(n)*(n-1));
        if(tmpweight<kw){
            (*iter).setVertexId(-1);
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
            if(singleVertexCandidate[j].getVertexId()!=-1){
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

void Graphflow::setIsolateVertexMatchResult(std::vector<int> &r, std::vector<int> &isolateVertex,float density) {
    for(int i=0;i<isolateVertex.size();i++){
        uint index=isolateVertex[i];
        this->match[index].setVertexId(r[i]);
    }
    int depth=this->match.size();
    this->match[depth-1].setSumWeight(density);
}


void Graphflow::setBatchVisited(std::vector<int> &r,bool flag) {
    for(auto item:r){
        /*  if(item==61070){
              std::cout<<"setBatchVisited 61070: "<<flag<<std::endl;
          }*/
        this->visited_[item]= flag;
    }
}


void Graphflow::recoverIsolateVertexMatchResult(std::vector<int> &IsolateVertexs) {
    for(int i=0;i<IsolateVertexs.size();i++)
        this->match[IsolateVertexs[i]].setIsolateSingleCandate();
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
    std::vector<std::vector<SingleCandidate>>combinezIsolateVertexs;
    combinezIsolateVertexs.reserve(isolateVertexs.size());
    for(int i=0;i<isolateVertexs.size();i++){
        std::vector<SingleCandidate>&singleVertex=matchCandidate[isolateVertexs[i]];
        std::sort(singleVertex.begin(),singleVertex.end());
        combinezIsolateVertexs.emplace_back(singleVertex);
    }

    int len=combinezIsolateVertexs.size();
#ifdef PRINT_DEBUG
    std::vector<int>testm(this->query_.NumVertices(),-1);
    for(int i=0;i<n;i++){
        if(match[i].getVertexId()!=-1){
            testm[order_vs_[matchorderindex][i]]=match[i].getVertexId();
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
            int candi=matchCandidate[i][j].getVertexId();
            std::cout<<candi<<" ";
        }
        std::cout<<endl;
    }
   std::cout<<endl;
#endif
    auto wt=findWeightBeforeIsolated();
    int*pmax=new int[len]{};//最大pmax在组中的列号
    std::unordered_map<int,int>visitedId;
    float Tmax=wt;
    bool TmaxisAllFirst=true;
    for(int i=0;i<len;i++){
        int index=0;
        int id=combinezIsolateVertexs[i][index].getVertexId();
        bool flag= true;
        while(visited_[id]){
            combinezIsolateVertexs[i].erase(combinezIsolateVertexs[i].begin()+index);
            if(index==combinezIsolateVertexs[i].size())
            {
                flag= false;
                break;
            }
            id=combinezIsolateVertexs[i][index].getVertexId();
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
            float gap1=combinezIsolateVertexs[pre_group_id][pmax[pre_group_id]].getSumWeight()-combinezIsolateVertexs[pre_group_id][pmax[pre_group_id]+1].getSumWeight();
            float gap2=combinezIsolateVertexs[i][pmax[i]].getSumWeight()-combinezIsolateVertexs[i][pmax[i]+1].getSumWeight();
            if(gap1<gap2){
                pmax[pre_group_id]++;
                visitedId[id]=i;
                int preId=combinezIsolateVertexs[pre_group_id][pmax[pre_group_id]].getVertexId();
                visitedId[preId]=pre_group_id;
            }
            else{
                pmax[i]++;
            }
            int curId=combinezIsolateVertexs[i][pmax[i]].getVertexId();
            visitedId[curId]=i;
        }
    }
    float globalTmax=Tmax;
    std::vector<int>r(isolateVertexs.size());
    for(int i=0;i<len;i++){
        r[i]=combinezIsolateVertexs[i][pmax[i]].getVertexId();
        Tmax+=combinezIsolateVertexs[i][pmax[i]].getSumWeight();
    }
    if(TmaxisAllFirst){
        //1.
        globalTmax=Tmax;
        float density=Tmax/(sqrt(n)*(n-1));
        setIsolateVertexMatchResult(r,isolateVertexs,Tmax);
        std::vector<uint>m(n);
        for(int i=0;i<match.size();i++){
            m[order_vs_[matchorderindex][i]]=match[i].getVertexId();
        }
        MatchRecord *record=new MatchRecord(density,m);
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
            globalTmax+=combinezIsolateVertexs[i][0].getSumWeight();
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
        float pre=combinezIsolateVertexs[next_vertex_group][pmax[next_vertex_group]].getSumWeight();
        float cur=combinezIsolateVertexs[next_vertex_group][hash[next_vertex_group]].getSumWeight();
        Tbound[next_vertex_group]=globalTmax-pre+cur;
        float TboundNext=Tbound[next_vertex_group];
        if(TboundNext>Tmax)
            Tmax=Tbound[next_vertex_group];
        auto next_item=combinezIsolateVertexs[next_vertex_group][hash[next_vertex_group]];

        //3. 递归
        int id=next_item.getVertexId();
        CatesianProductWithIndex(matchorderindex,type,next_vertex_group,0,len,hash,combinezIsolateVertexs,isolateVertexs,wt);
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
                /* for(int i=0;i<len;i++)
                     std::cout<<noscan[i]<<" ";
                   std::cout<<endl;*/
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
        float copydensity=density;
        tmin=std::min(std::get<1>(combinezIsolateVertexs[k][i]),tmin);
        density+=std::get<2>(combinezIsolateVertexs[k][i]);
        current.push_back(std::get<0>(combinezIsolateVertexs[k][i]));
        combinationMatchResultHelp(result,current,combinezIsolateVertexs,k+1,tmin,density);
        current.pop_back();
        density=copydensity;
    }
}
float Graphflow::findWeightBeforeIsolated() {
    int depth=match.size()-1;
    for(int i=depth;i>=0;i--){
        if(match[i].getVertexId()==-1){
            continue;
        }
        else{
            float weight=match[i].getSumWeight();
            return weight;
        }
    }
}
void Graphflow::CatesianProductWithIndex(int matchorderindex,searchType type,int curIndex,int depth,int len,int*hash,std::vector<std::vector<SingleCandidate>>&combinezIsolateVertexs,std::vector<int>&isolateVertexs,float &weight) {
    if(depth==len){
        int n=query_.NumVertices();
        std::vector<uint>m(n);
        for(int i=0;i<match.size();i++){
            m[order_vs_[matchorderindex][i]]=match[i].getVertexId();
        }
        float density=weight/ (sqrt(n)*(n-1));
        MatchRecord *record=new MatchRecord(density,m);
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
            int id=combinezIsolateVertexs[depth][i].getVertexId();
            if(visited_[id])
                continue;
            visited_[id]=true;
            float copyweight=weight;
            this->match[isolateVertexs[depth]].setVertexId(id);
            weight+=combinezIsolateVertexs[depth][i].getSumWeight();
            CatesianProductWithIndex(matchorderindex,type,curIndex,depth+1,len,hash,combinezIsolateVertexs,isolateVertexs,weight);
            visited_[id]=false;
            this->match[isolateVertexs[depth]].setVertexId(-1);
            weight=copyweight;
        }
    }else{
        //  std::cout<<"hash["<<curIndex<<"]"<<hash[curIndex]<<" depth["<<depth<<"]"<<combinezIsolateVertexs[depth][hash[curIndex]].getVertexId()<<endl;
        const int &id=combinezIsolateVertexs[depth][hash[curIndex]].getVertexId();
        if(visited_[id]){
            return;
        }
        visited_[id]=true;
        float copyweight=weight;
        this->match[isolateVertexs[depth]].setVertexId(id);
        weight+=combinezIsolateVertexs[depth][hash[curIndex]].getSumWeight();
        CatesianProductWithIndex(matchorderindex, type,curIndex,depth+1,len,hash,combinezIsolateVertexs,isolateVertexs,weight);
        visited_[id]=false;
        this->match[isolateVertexs[depth]].setVertexId(-1);
        weight=copyweight;
    }

}
int Graphflow::findTboundMaxIndex(float *Tbound,int*hash,int* noscan,std::vector<std::vector<SingleCandidate>>&combinezIsolateVertexs,int len) {
    float max=0;
    int maxindex=0;
    for(int i=0;i<len;i++){
        int index=hash[i]+1;
        if(noscan[i]==1)
        {
            continue;
        }
        if(index>=combinezIsolateVertexs[i].size())
        {
            noscan[i]=1;
            continue;
        }
        int id=combinezIsolateVertexs[i][index].getVertexId();
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
                id=combinezIsolateVertexs[i][index].getVertexId();
            }
        }
        if(!flag)
            continue;
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
    std::vector<std::vector<SingleCandidate>>combinezIsolateVertexs;
    combinezIsolateVertexs.reserve(isolateVertexs.size());
    for(int i=0;i<isolateVertexs.size();i++){
        std::vector<SingleCandidate>&singleVertex=matchCandidate[isolateVertexs[i]];
        std::sort(singleVertex.begin(),singleVertex.end());
        combinezIsolateVertexs.emplace_back(singleVertex);
    }

    int len=combinezIsolateVertexs.size();
#ifdef PRINT_DEBUG
    std::vector<int>testm(this->query_.NumVertices(),-1);
    for(int i=0;i<n;i++){
        if(match[i].getVertexId()!=-1){
            testm[order_vs_[matchorderindex][i]]=match[i].getVertexId();
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
            int candi=matchCandidate[i][j].getVertexId();
            std::cout<<candi<<" ";
        }
        std::cout<<endl;
    }
    std::cout<<endl;
#endif
    auto wt=findWeightBeforeIsolated();
    int*hash=new int[len]{};//索引指针
    int*noscan=new int[len]{};//记录是否停止
    std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int>>, pairCompare> maxHeap;
    for(int i=0;i<combinezIsolateVertexs.size();i++){
        float weight=combinezIsolateVertexs[i][0].getSumWeight();
        maxHeap.push(std::make_pair(weight,i));
    }
    //计算全局Tmax//大于 则直接return
/*    float Tmax=wt;
    std::unordered_map<int,int>visitedId;
    int*pmax=new int[len]{};
    for(int i=0;i<len;i++){
        int index=0;
        int id=combinezIsolateVertexs[i][index].getVertexId();
        bool flag= true;
        while(visited_[id]){
            if(index==combinezIsolateVertexs[i].size())
            {
                flag= false;
                break;
            }
            combinezIsolateVertexs[i][index].setVertexId(-1);
           // combinezIsolateVertexs[i].erase(combinezIsolateVertexs[i].begin()+index);
            index++;
            id=combinezIsolateVertexs[i][index].getVertexId();
        }
        pmax[i]=index;
        if(!flag)
        {
            return;
        }
        if(!visitedId.count(id)){
            visitedId[id]=i;
        }
        else{
            int pre_group_id=visitedId[id];
            float gap1=combinezIsolateVertexs[pre_group_id][pmax[pre_group_id]].getSumWeight()-combinezIsolateVertexs[pre_group_id][pmax[pre_group_id]+1].getSumWeight();
            float gap2=combinezIsolateVertexs[i][pmax[i]].getSumWeight()-combinezIsolateVertexs[i][pmax[i]+1].getSumWeight();
            if(gap1<gap2){
                if(pmax[pre_group_id]<combinezIsolateVertexs[pre_group_id].size()-1)
                {
                    pmax[pre_group_id]++;
                    visitedId[id]=i;
                    int preId=combinezIsolateVertexs[pre_group_id][pmax[pre_group_id]].getVertexId();
                    visitedId[preId]=pre_group_id;
                }
                else{
                    pmax[i]++;
                }
            }
            else{
                if(pmax[i]<combinezIsolateVertexs[i].size()-1)
                {
                    pmax[i]++;
                }
                else{
                    pmax[pre_group_id]++;
                    visitedId[id]=i;
                    int preId=combinezIsolateVertexs[pre_group_id][pmax[pre_group_id]].getVertexId();
                    visitedId[preId]=pre_group_id;
                }
            }
            int curId=combinezIsolateVertexs[i][pmax[i]].getVertexId();
            visitedId[curId]=i;
        }
    }
    for(int i=0;i<len;i++){
        Tmax+=combinezIsolateVertexs[i][pmax[i]].getSumWeight();
    }
    if(topKSet.size()==k){
          float density=Tmax/(sqrt(n)*(n-1));
        if(topKSet.back()->getDensity()>density)
            return;
    }*/
    float globalTmax=wt;
    for(int i=0;i<len;i++)
        globalTmax+=combinezIsolateVertexs[i][0].getSumWeight();
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
            float w=combinezIsolateVertexs[group][hash[group]].getSumWeight();
            maxHeap.push(std::make_pair(w,group));
        }
        int id=combinezIsolateVertexs[group][pre].getVertexId();
        if(visited_[id]){
            combinezIsolateVertexs[group].erase(combinezIsolateVertexs[group].begin()+pre);
            hash[group]--;
            continue;
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
            this->match[isolateVertexs[group]].setVertexId(id);
            float weight=wt+combinezIsolateVertexs[group][pre].getSumWeight();
            float Tmax=0;
            //递归
            CatesianProductWithHeap(matchorderindex,type,0,depthLen,hash,combinezIsolateVertexs,isolateVertexs,copyisolateIndex,weight);
            visited_[id]=false;
            this->match[isolateVertexs[group]].setVertexId(-1);
            if(topKSet.size()==k)
            {
                float curDensity=topKSet.back()->getDensity();
                Tmax=globalTmax-combinezIsolateVertexs[group][0].getSumWeight()+combinezIsolateVertexs[group][pre].getSumWeight();
                float Tmaxdensity=Tmax/ (sqrt(n)*(n-1));
                if(curDensity>Tmaxdensity){
                    noscan[group]=1;
                }
            }
        }
    }
    delete[]hash;
    delete[]noscan;
}
void Graphflow::CatesianProductWithHeap(int matchorderindex, searchType type, int depth, int len, int *hash,
                                        std::vector<std::vector<SingleCandidate>> &combinezIsolateVertexs,
                                        std::vector<int> &isolateVertexs, std::vector<int> &isolatedIndex,  float &weight) {
    if(depth==len){
        int n=query_.NumVertices();
        std::vector<uint>m(n);
        for(int i=0;i<match.size();i++){
            m[order_vs_[matchorderindex][i]]=match[i].getVertexId();
        }
        float density=weight/ (sqrt(n)*(n-1));

        MatchRecord *record=new MatchRecord(density,m);
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
        int id=combinezIsolateVertexs[group][i].getVertexId();
        if(visited_[id])
            continue;
        visited_[id]=true;
        float copyweight=weight;
        this->match[isolateVertexs[group]].setVertexId(id);
        weight+=combinezIsolateVertexs[group][i].getSumWeight();
        CatesianProductWithHeap(matchorderindex,type,depth+1,len,hash,combinezIsolateVertexs,isolateVertexs,isolatedIndex,weight);
        visited_[id]=false;
        this->match[isolateVertexs[group]].setVertexId(-1);
        weight=copyweight;
    }
}


void Graphflow::createLabelToQueryVertex() {
    for(int i=0;i<query_.NumVertices();i++){
        uint label=query_.GetVertexLabel(i);
        labelToQueryVertex[label].emplace_back(i);
        queryVertexIndexInlabel[i]=labelToQueryVertex[label].size()-1;
    }
}

bool Graphflow::updaterightNeighborCandidate(int matchorderindex,uint uk,uint uk_neigh,bool isFirstEdge, uint vk,const std::vector<uint>&uk_neighbor) {
    const  std::vector<Neighbor>&vN= this->globalsubgraph_.vNeighbors[vk];
#ifdef LOG_TRACK
#endif
    const int n=uk_neighbor.size();
    //1.对于所有的右邻居，找其候选解
    for(int i=0;i<n;i++) {
        uint query_id = uk_neighbor[i];
        if(isFirstEdge){
            if(query_id==uk_neigh)
            {
                isFirstEdge= false;
                continue;
            }
        }
        uint query_vertex_label = query_.GetVertexLabel(query_id);
        int query_order_index=order_vertex_index[matchorderindex][query_id];
        uint query_elabel=std::get<2>(query_.GetEdgeLabel(uk,query_id));
        StarGraph*s=globalStarIndex[matchorderindex][query_order_index];
        bool isFirstVertex= true;
        bool isCandidateFirstNull= true;
        float maxweight=0;
        float curWeight=0;
        if(matchCandidate[query_order_index].size()!=0)
            isCandidateFirstNull= false;
        //对于vk的每个邻居neighbor
        for (const auto &neighbor: vN) {
            uint neighbor_id = neighbor.getVertexId();
            const uint neighbor_match_id=neighbor.getMatchQueryVertexId();
            const uint neighbor_from_id=neighbor.getfromVertexId();
            if(visited_[neighbor_id]||neighbor_match_id!=query_id||neighbor_from_id!=uk)
                continue;
            uint  v_elabel=neighbor.GetEdgelabel();
            if (neighbor.getVertexLabel() == query_vertex_label&&query_elabel==v_elabel) {
                maxweight = globalVkMatchUk[neighbor_id][matchorderindex][queryVertexIndexInlabel[query_id]];
                //update LocalStarIndex
                //add candidate
                curWeight=neighbor.GetEdgeWeight();
                if (isCandidateFirstNull) {
                    if(isFirstVertex){
                        isFirstVertex= false;
                        LocalStarIndex[query_order_index]=maxweight;
                    }
                    else{
                        if(LocalStarIndex[query_order_index]<maxweight){
                            LocalStarIndex[query_order_index]=maxweight;
                        }
                    }
                    matchCandidate[query_order_index].emplace_back( neighbor_id,curWeight);
                }
                else{
                    std::vector<SingleCandidate>&neigh_candidate=matchCandidate[query_order_index];
                    for(SingleCandidate&s:neigh_candidate){
                        if(s.getVertexId()==neighbor_id){
                            if(isFirstVertex){
                                isFirstVertex= false;
                                LocalStarIndex[query_order_index]=maxweight;
                            }
                            else{
                                if(LocalStarIndex[query_order_index]<maxweight){
                                    LocalStarIndex[query_order_index]=maxweight;
                                }
                            }
                            s.addFlag();
                            s.addSumWeight(curWeight);
                            break;
                        }
                    }
                }
            }
            if(isInsert)
                IdeterminCandite++;
            else
                DdeterminCandite++;
        }
        if(matchCandidate[query_order_index].size()==0) {
            //isFirst恢复true;
            for (int i = 0; i < n; i++) {
                uint query_id = uk_neighbor[i];
                int query_order_index = order_vertex_index[matchorderindex][query_id];
                matchCandidate[query_order_index].clear();
            }
            return true;
        }
    }
    return false;
}
void Graphflow::InitialLocalIndex(int matchorderindex) {
    const std::vector<StarGraph*> & gs=globalStarIndex[matchorderindex];
    int n=gs.size();
    for(int i=1;i<n;i++){
        LocalStarIndex[i]=gs[i]->getStarMaxWeight();
    }
}
void Graphflow::getIntersetSingleCandidate(std::vector<SingleCandidate> &singleVertexCandidate,int matchorderindex,int depth) {
    int i=0;
    int j=0;
    int csize=singleVertexCandidate.size();
    int len=0;
    int flag=matchLeftNeighborSum[matchorderindex][depth];
    if(flag==1)
        return;
    while(j<csize){
        if(singleVertexCandidate[j].getFlag()==flag){
            singleVertexCandidate[i]=singleVertexCandidate[j];
            i++;
            len++;
        }
        j++;

    }
    singleVertexCandidate.resize(len);
}
void Graphflow::PrintAverageTime(int len) {
     int ilen=10000-len;
    std::cout <<"average query graph degree:"<< std::fixed << std::setprecision(2)<<query_.NumEdges()*2.0/query_.NumVertices()<<endl;
    std::cout<<"average data graph degree:"<<std::fixed << std::setprecision(2)<<data_.NumEdges()*2.0/data_.NumVertices()<<endl;
    std::cout << "average serach time: " << std::fixed << std::setprecision(2)
              << total_search_time.GetTimer() * 1.0 / ilen << " microseconds" << endl;
    std::cout << "average update global index time: " << std::fixed << std::setprecision(2)
              << total_update_globalIndex_time.GetTimer() * 1.0 / ilen << " microseconds" << endl;
    std::cout << "average update time " << std::fixed << std::setprecision(2)
              << (total_update_globalIndex_time.GetTimer() * 1.0 / ilen + total_search_time.GetTimer() * 1.0 / ilen)
              << " microseconds" << endl;
    std::cout << "average insert density filter time: " << std::fixed << std::setprecision(2)
              << Itotal_densityfilter_time * 1.0 / ilen << " microseconds" << endl;
    std::cout << "average insert updaterightNeighborCandidate time " << std::fixed << std::setprecision(2)
              <<  (Itotal_updaterightNeighborCandidate_time* 1.0 / ilen)<< " microseconds" << endl;
    std::cout << "average print time: " << std::fixed << std::setprecision(2) << total_print_time.GetTimer() * 1.0 / ilen
              << " microseconds" << endl;
    std::cout << "average delete search time:" << std::fixed << std::setprecision(2)
              << total_delete_time.GetTimer() * 1.0 / len << " microseconds" << endl;
    std::cout << "average delete update global subgraph time:" << std::fixed << std::setprecision(2)
              << (total_delete_update_time.GetTimer() * 1.0) / len << " microseconds" << endl;
    std::cout << "average delete update time:" << std::fixed << std::setprecision(2)
              << (total_delete_time.GetTimer() * 1.0 / len + total_delete_update_time.GetTimer() * 1.0 / len)
              << " microseconds" << endl;
    std::cout << "average delete density filter time: " << std::fixed << std::setprecision(2)
              << total_densityFilter_time.GetTimer() * 1.0 / len << " microseconds" << endl;
    std::cout << "average delete updaterightNeighborCandidate time: " << std::fixed << std::setprecision(2)
              << total_updaterightNeighborCandidate_time.GetTimer() * 1.0 / len << " microseconds" << endl;
    std::cout << "num add edge update:" << numupdatestar / query_.NumEdges() << endl;
    std::cout << "num add matchresult:" <<  sumAllMatchFind << endl;
    std::cout << "num delete matchresult:" <<  sumDeleteallMatchFind << endl;
    std::cout << "num insert searchspace:" <<  IsearchSpace << endl;
    std::cout << "num delete searchspace:" <<  DsearchSpace << endl;
    std::cout<<"num insert check neighbor:"<<IdeterminCandite<<endl;
    std::cout<<"num delete check neighbor:"<<DdeterminCandite<<endl;
#ifdef COMPUTE_TRACK
stringstream _ss;
    _ss<< (total_update_globalIndex_time.GetTimer() * 1.0 / ilen + total_search_time.GetTimer() * 1.0 / ilen)<<","
       <<total_update_globalIndex_time.GetTimer() * 1.0 / ilen<<","
       <<total_search_time.GetTimer() * 1.0 / ilen<<","
       <<(total_delete_time.GetTimer() * 1.0 / len + total_delete_update_time.GetTimer() * 1.0 / len)<<","
       <<(total_delete_update_time.GetTimer() * 1.0) / len<<","
       <<total_delete_time.GetTimer() * 1.0 / len<<","
       <<sumAllMatchFind<<","
       <<sumDeleteallMatchFind<<","
       <<query_.NumEdges()*2.0/query_.NumVertices()<<","
       <<data_.NumEdges()*2.0/data_.NumVertices()<<","
       <<Itotal_densityfilter_time * 1.0 / ilen<<","
       <<(Itotal_updaterightNeighborCandidate_time* 1.0 / ilen)<<","
       <<total_densityFilter_time.GetTimer() * 1.0 / len<<","
       <<total_updaterightNeighborCandidate_time.GetTimer() * 1.0 / len<<","
       << IsearchSpace<<","
       <<DsearchSpace<<","
       <<IdeterminCandite<<","
       <<DdeterminCandite<<","
       <<space_cost<<endl;
    Log::track3(_ss);
#endif
}
void Graphflow::createGlobalSubgraph() {
    //对q于每个节点用nlf检查
    int n = data_.NumVertices();
    int m = query_.NumVertices();
    for (int i = 0; i < n; i++) {
        uint vlabel = data_.GetVertexLabel(i);
        for (int j = 0; j < m; j++) {
            uint qlabel = query_.GetVertexLabel(j);
            if (vlabel == qlabel) {
                if (LabelFilter(i, j)) {
                    globalsubgraph_.addQueryVertexCandidate(j, i);
                }
            }
        }
    }
    //addEdge
    for(int i=0;i< this->query_.NumVertices();i++){
        std::vector<uint>& irightneighbors=rightNeighbor[0][i];
        const std::vector<uint>&i_candidate=globalsubgraph_.matchCandidate[i];
        for(uint ic:i_candidate){
            const std::vector<uint>&i_nbrs=data_.GetNeighbors(ic);
            for(uint inbr:i_nbrs){
                for(uint j:irightneighbors){
                    const std::vector<uint>&j_candidate=globalsubgraph_.matchCandidate[j];
                    if(std::binary_search(j_candidate.begin(),j_candidate.end(), inbr)){
                        uint v1label=data_.GetVertexLabel(ic);
                        uint v2label=data_.GetVertexLabel(inbr);
                        uint edgelabel=std::get<2>(data_.GetEdgeLabel(ic,inbr));
                        float weight=data_.GetEdgeWeight(ic,inbr);
                        globalsubgraph_.AddEdge(i,j,ic,data_.GetVertexLabel(ic),inbr,data_.GetVertexLabel(inbr),edgelabel,weight);
                    }
                }
            }
        }
    }


}
bool Graphflow::updateGlobalGraphHelp(int m, uint u1, uint u2, uint u1label, uint u2label, uint v1,uint v2,uint v1label, uint v2label,
                                      uint elabel, const std::vector<std::vector<uint>>&mcandidate,bool &flag) {
    bool isMatch= false;
    uint n=query_.NumEdges();
    bool isNew1= false;
    bool isNew2= false;
    bool isContain1= true;
    bool isContain2= true;
    if(v1label!=u1label)
    {
        swap(v1,v2);
    }
    if(!std::binary_search(mcandidate[u1].begin(),mcandidate[u1].end(), v1)){
        if(this->LabelFilter(v1,u1)){
            isNew1= true;
        }
        else{
            isContain1= false;
        }
    }
    if(!std::binary_search(mcandidate[u2].begin(),mcandidate[u2].end(), v2)){
        if(this->LabelFilter(v2,u2)){
            isNew2= true;
        }
        else{
            isContain2= false;
        }
    }
    if(isContain1&&isContain2){
        flag= true;
        isMatch= true;
        if(isNew1&&isNew2){
            updateglobalVertexStarIndex(u1,v1,u1label,elabel,n,mcandidate);
            updateglobalVertexStarIndex(u2,v2,u2label,elabel,n,mcandidate);
        }
        else if(isNew1&&!isNew2){
            updateglobalVertexStarIndex(u1,v1,u1label,elabel,n,mcandidate);
        }
        else if(!isNew1&&isNew2){
            updateglobalVertexStarIndex(u2,v2,u2label,elabel,n,mcandidate);
        }
        else{
            float w=data_.GetEdgeWeight(v1,v2);
            globalsubgraph_.AddEdge(u1,u2,v1,u1label,v2,u2label,elabel,w);
            int candidate_index= queryVertexIndexInlabel[u1];
            for(int j=0;j<query_.NumEdges();j++) {
                updateStarIndex(j,v1,u1,candidate_index);
            }
            int candidate_index2= queryVertexIndexInlabel[u2];
            for(int j=0;j<query_.NumEdges();j++) {
                updateStarIndex(j,v2,u2,candidate_index2);
            }
            numupdatestar+=2;
        }
    }
    else if(isContain1&&isNew1){
        updateglobalVertexStarIndex(u1,v1,u1label,elabel,n,mcandidate);
    }
    else if(isContain2&&isNew2){
        updateglobalVertexStarIndex(u2,v2,u2label,elabel,n,mcandidate);
    }
    return isMatch;
 }
bool Graphflow::updateGlobalSubgraph(uint v1, uint v2, uint label, float weight,std::vector<int>&match) {
    bool flag= false;
        //查看是否满足nlf条件，如果满足，判断是否已经在候选集中，如果是，直接加边，如果不是更新候选集
        //如果不满足nlf条件，跳过
        uint n=query_.NumEdges();
        const auto &mcandidate=globalsubgraph_.matchCandidate;
        for(auto it=match.begin();it!=match.end();){
            bool isMatch= false;
            uint u1=order_vs_[*it][0];
            uint u2=order_vs_[*it][1];
            uint u1label=query_.GetVertexLabel(u1);
            uint u2label=query_.GetVertexLabel(u2);
            uint v1label=data_.GetVertexLabel(v1);
            uint v2label=data_.GetVertexLabel(v2);
            uint elabel=std::get<2>(query_.GetEdgeLabel(u1,u2));
            if(v1label!=v2label) {
                if(updateGlobalGraphHelp(*it, u1, u2, u1label, u2label, v1, v2, v1label, v2label, elabel,
                                      globalsubgraph_.matchCandidate, flag)){
                    isMatch= true;
                }
            }
            else{
                for(int i=0;i<2;i++){
                    if(updateGlobalGraphHelp(*it, u1, u2, u1label, u2label, v1, v2, v1label, v2label, elabel,
                                          globalsubgraph_.matchCandidate, flag)){
                        isMatch= true;
                    }
                    else{
                        if(isMatch){
                            continue;
                        }
                    }
                    std::swap(v1,v2);
                }
            }
           if(isMatch== false){
               it=match.erase(it);
           }
           else{
               ++it;
           }
        }
#ifdef LOG_TRACK
        stringstream _ss;
        for(auto i:globalsubgraph_.matchCandidate[1]){
         if(i==63117){
            _ss<<"first 63117"<<endl;
         }
        }
    for(auto n:globalsubgraph_.vNeighbors[63117]){
        if(n.getVertexId()==1709&&n.getMatchQueryVertexId()==0){
            _ss<<"find 1709"<<endl;
        }
    }
        Log::track1(_ss);
#endif
    return flag;
}
void Graphflow::updateglobalVertexStarIndex(uint u1,uint v1,uint u1label,uint elabel,uint n, const std::vector<std::vector<uint>>&mcandidate) {
    globalsubgraph_.addQueryVertexCandidate(u1,v1);
    int v1size=labelToQueryVertex[query_.GetVertexLabel(u1)].size();
    globalVkMatchUk[v1].resize(n);
    for(int w=0;w<n;w++){
        globalVkMatchUk[v1][w].resize(v1size);
    }
    const std::vector<uint>&neighbors1=query_.GetNeighbors(u1);
    const std::vector<uint>&v1_nbrs=data_.GetNeighbors(v1);
    bool flagAdd;
    for(uint v1_nbr:v1_nbrs){
        for(uint n1:neighbors1){
            uint n1label=query_.GetVertexLabel(n1);
            if(std::binary_search(mcandidate[n1].begin(), mcandidate[n1].end(),v1_nbr)){
                float weight=data_.GetEdgeWeight(v1,v1_nbr);
                flagAdd=globalsubgraph_.AddEdge(u1,n1,v1,u1label,v1_nbr,n1label,elabel,weight);
                if(flagAdd) {
                    int candidate_index = queryVertexIndexInlabel[n1];
                    numupdatestar++;
                    for (int j = 0; j < query_.NumEdges(); j++) {
                        updateStarIndex( j, v1_nbr, n1, candidate_index);
                    }

                }
            }
        }
    }
    if(flagAdd) {
        int candidate_index = queryVertexIndexInlabel[u1];
        numupdatestar++;
        for (int j = 0; j < query_.NumEdges(); j++) {
            updateStarIndex( j, v1, u1, candidate_index);
        }
    }
 }
bool Graphflow::deleteMatchRecordWithEdge(uint v1, uint v1label,uint v2, uint v2label,uint label,std::vector<int>&match) {
    bool flag = false;
    for (auto it = topKSet.begin(); it != topKSet.end();) {
        MatchRecord *record = *it;
        const std::vector<uint> &m = record->getVetex();
        bool iterflag= true;
        for (int mindex: match) {
            uint u1 = order_vs_[mindex][0];
            uint u2 = order_vs_[mindex][1];
            if ((m[u1] == v1 && m[u2] == v2) || (m[u2] == v1 && m[u1] == v2)) {
                delete record;
                record= nullptr;
                it = topKSet.erase(it);
                iterflag= false;
                flag = true;
                num_negative_results_++;
                break;
            }
        }
        if(iterflag){
            ++it;
        }
    }
    return flag;
}
bool Graphflow::SearchMatchesWithEdge(uint m,uint v1,uint v2,float weight,uint u1,uint u2,searchType type){
    this->matchVertex(true, 0, v1, float(0));
    this->matchVertex(true, 1, v2, weight);
    if(isInsert)
        IsearchSpace+=2;
    else
        DsearchSpace+=2;
    bool isNull;
    const std::vector<uint>&uk_neighbor1=rightNeighbor[m][u1];
    total_updaterightNeighborCandidate_time.StartTimer();
    isNull=updaterightNeighborCandidate(m, u1, u2,true,v1, uk_neighbor1);
    total_updaterightNeighborCandidate_time.StopTimer();
    if(isNull)
    {
        this->popVertex(1, v2);
        this->popVertex(0, v1);
        return true;
    }
    const std::vector<uint>&uk_neighbor2=rightNeighbor[m][u2];
    total_updaterightNeighborCandidate_time.StartTimer();
    isNull=updaterightNeighborCandidate(m, u2, u1, true,v2, uk_neighbor2);
    total_updaterightNeighborCandidate_time.StopTimer();
    if(isNull)
    {
        for(int u_id:uk_neighbor1){
            int query_order_index=order_vertex_index[m][u_id];
            matchCandidate[query_order_index].clear();
        }
        this->popVertex(1, v2);
        this->popVertex(0, v1);
        return true;
    }
    searchMatches(2, m, type);
    this->popVertex( v2,m,1,uk_neighbor1);
    this->popVertex(v1,m,0,uk_neighbor2);
    return false;
 }

void Graphflow::deleteGlobalSubgraph(uint v1, uint v2,uint elabel,float weight, std::vector<int> &match) {
    //查看是否满足nlf条件，如果满足，判断是否已经在候选集中，如果是，直接加边，如果不是更新候选集
    //如果不满足nlf条件，跳过
    uint n=query_.NumEdges();
    bool isMatch= false;
    auto &mcandidate=globalsubgraph_.matchCandidate;
    for(auto it=match.begin();it!=match.end();it++){
        uint u1=order_vs_[*it][0];
        uint u2=order_vs_[*it][1];
        uint u1label=query_.GetVertexLabel(u1);
        uint u2label=query_.GetVertexLabel(u2);
        uint v1label=data_.GetVertexLabel(v1);
        uint v2label=data_.GetVertexLabel(v2);
        uint elabel=std::get<2>(query_.GetEdgeLabel(u1,u2));
        if(v1label!=v2label) {
            if(v1label!=u1label){
                std::swap(v1,v2);
            }
//            std::chrono::high_resolution_clock::time_point start;
            globalsubgraph_.RemoveEdge(v1,u1label,v2,u2label,u1,u2,elabel,weight);
            isMatch=deleteGlobalSubgraphHelp(*it,u1,u2,u1label,u2label,v1,v2,v1label,v2label,elabel,weight,mcandidate);
        }
        else{
            for(int i=0;i<2;i++){
                globalsubgraph_.RemoveEdge(v1,v1label,v2,v2label,u1,u2,elabel,weight);
                isMatch=deleteGlobalSubgraphHelp(*it,u1,u2,u1label,u2label,v1,v2,v1label,v2label,elabel,weight,mcandidate);
                std::swap(v1,v2);
            }
        }
       /* if(isMatch== false){
            it=match.erase(it);
        }
        else{
            ++it;
        }*/
    }
    //total_delete_update_time+= Duration2(start);
 }
 void Graphflow::deleteUpdateglobalVertexStarIndex(uint u1,uint v1,uint v2,uint n) {
     int candidate_index= queryVertexIndexInlabel[u1];
     bool isContain;
     const std::vector<uint>&candidate=globalsubgraph_.matchCandidate[u1];
     //更新u1位置的索引
     for(int j=0;j<n;j++){
        // isContain= false;
         int vertex_index=order_vertex_index[j][u1];
         if(vertex_index==0)
             continue;
         StarGraph *s=globalStarIndex[j][vertex_index];
         if(s->getMatchDataVertexId()==v1){
             s->setStarMaxWeight(s->GetForwardNeighborNum()*mw);
             updateStarIndex(j,v1,u1,candidate_index);
             for(uint v:candidate){
                 if(globalVkMatchUk[v][j][candidate_index]>s->getStarMaxWeight()||s->getStarMaxWeight()==s->GetForwardNeighborNum()*mw)
                 {
                     s->setStarMaxWeight(globalVkMatchUk[v][j][candidate_index]);
                     s->setMatchDataVertexId(v);
                 }
             }
         }
     }
 }

 bool Graphflow::deleteGlobalSubgraphHelp(int m,uint u1,uint u2,uint u1label,uint u2label, uint v1,uint v2,uint v1label,uint v2label,
                                          uint elabel,float weight, std::vector<std::vector<uint>>&mcandidate) {
     bool isMatch= false;
     uint n=query_.NumEdges();
     bool isContain1= true;
     bool isContain2= true;
     if(v1label!=u1label)
     {
         swap(v1,v2);
     }
     if(!std::binary_search(mcandidate[u1].begin(),mcandidate[u1].end(), v1)){
             isContain1= false;
     }
     if(!std::binary_search(mcandidate[u2].begin(),mcandidate[u2].end(), v2)){
             isContain2= false;
     }
     //对其他的边
     //判断两个点是否还符合条件
     if(isContain1&&isContain2){
         isMatch= true;
         bool flag1= LabelFilter(v1,u1);
         bool flag2= LabelFilter(v2,u2);
         uint v1label=data_.GetVertexLabel(v1);
         uint v2label=data_.GetVertexLabel(v2);
         //globalsubgraph_.RemoveEdge(v1,v1label,v2,v2label,u1,u2,elabel,weight);
         if(flag1&&flag2){
             //仍在候选节点中，删除边，更新v1,v2 局部索引
             deleteUpdateglobalVertexStarIndex(u1,v1,v2,n);
             deleteUpdateglobalVertexStarIndex(u2,v2,v1,n);
         }
         else if(!flag1&&flag2){
             //删除v1并更新v1_nbr的局部索引
             deleteGlobalGraphCandidateEdges(m,u1,v1,mcandidate);
         }
         else if(flag1&&!flag2){
             //删除v2并更新v2_nbr的局部索引
             deleteGlobalGraphCandidateEdges(m,u2,v2,mcandidate);
         }
         else{
             deleteGlobalGraphCandidateEdges(m,u1,v1,mcandidate);
             deleteGlobalGraphCandidateEdges(m,u2,v2,mcandidate);
         }
     }
     return isMatch;
 }



 void Graphflow::deleteGlobalGraphCandidateEdges(uint m,uint u1,uint v1,std::vector<std::vector<uint>>&mcandidate) {
     const std::vector<uint> &neighbors1 = query_.GetNeighbors(u1);
     globalsubgraph_.deleteQueryVertexCandidate(u1, v1);
     const std::vector<uint> &v1_nbrs = data_.GetNeighbors(v1);
     bool flagDel;
     for (uint v1_nbr: v1_nbrs) {
         for (uint n1: neighbors1) {
             uint n1label = query_.GetVertexLabel(n1);
             if (std::binary_search(mcandidate[n1].begin(), mcandidate[n1].end(), v1_nbr)) {
                 uint v1label = data_.GetVertexLabel(v1);
                 uint v1_nbr_label = data_.GetVertexLabel(v1_nbr);
                 uint elabel = std::get<2>(data_.GetEdgeLabel(v1, v1_nbr));
                 uint weight = data_.GetEdgeWeight(v1, v1_nbr);
                 flagDel = globalsubgraph_.RemoveEdge(v1, v1label, v1_nbr, v1_nbr_label, u1, n1, elabel, weight);
                 if (flagDel) {
                     //判断是否对其最大权值有影响
                     int candidate_index = queryVertexIndexInlabel[n1];
                     // bool isContain= false;
                     const std::vector<uint> &candidate = globalsubgraph_.matchCandidate[n1];
                     for (int j = 0; j < query_.NumEdges(); j++) {
                         //isContain= false;
                         int vertex_index = order_vertex_index[j][n1];
                         if (vertex_index == 0)
                             continue;
                         StarGraph *s = globalStarIndex[j][vertex_index];
                         if (s->getMatchDataVertexId() == n1) {
                             int vertex_index = order_vs_[j][n1];
                             if (vertex_index == 0)
                                 continue;
                             StarGraph *s = globalStarIndex[j][vertex_index];
                             s->setStarMaxWeight(s->GetForwardNeighborNum() * mw);
                             updateStarIndex(j, v1_nbr, n1, candidate_index);
                             for (uint v: candidate) {
                                 if (globalVkMatchUk[v][j][candidate_index] > s->getStarMaxWeight() ||
                                     s->getStarMaxWeight() == s->GetForwardNeighborNum() * mw) {
                                     s->setStarMaxWeight(globalVkMatchUk[v][j][candidate_index]);
                                 }
                             }

                         }

                     }
                 }
             }
         }
     }
 }
