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
    bool operator()(const std::pair<float, int> &p1, const std::pair<float, int> &p2) {
        if (p1.first == p2.first) {
            return p1.second > p2.second;
        }
        return p1.first < p2.first;
    }
};


bool ForwardNeighborcmp(ForwardNeighbor *f1, ForwardNeighbor *f2) {
    return (*f1) > (*f2);
}

bool areSame(float a, float b) {
    return std::fabs(a - b) < std::numeric_limits<float>::epsilon();
}

bool CompareNeighbors2(const Neighbor& a, const Neighbor& b) {
    if (a.GetEdgelabel() != b.GetEdgelabel()) {
        return a.GetEdgelabel() >b.GetEdgelabel();
    } else {
        return a.getVertexLabel() > b.getVertexLabel();
    }
}
bool compareMatchRecord( MatchRecord* record1,  MatchRecord* record2) {
    return (*record1)>(*record2);
}
bool tupleVertexIdCmp2(const std::tuple<int, int, float> &a, const std::tuple<int, int, float> &b) {
    if (std::get<0>(a) != std::get<0>(b)) {
        return std::get<0>(a) < std::get<0>(b);
    } else {
        return std::get<1>(a) < std::get<1>(b);
    }
}


Graphflow::Graphflow(Graph &query_graph, Graph &data_graph,
                     uint max_num_results,
                     bool print_prep,
                     bool print_enum,
                     bool homo)
        : matching(query_graph, data_graph, max_num_results,
                   print_prep, print_enum, homo), order_vs_(query_.NumEdges()), order_csrs_(query_.NumEdges()),
          order_offs_(query_.NumEdges()), order_vertex_index(query_.NumEdges()), topKSet(0),
          suffixMax(query_graph.NumVertices(), 0), isolatedMax(query_graph.NumVertices(), -1),
          rightNeighbor(query_graph.NumEdges()), matchCandidate(query_graph.NumVertices()),
          match(query_graph.NumVertices()), labelToQueryVertex(data_graph.NumVLabels()),
          globalVkMatchUk(data_graph.NumVertices()), globalStarIndex(query_.NumEdges()),
          queryVertexIndexInlabel(query_.NumVertices()), LocalStarIndex(query_.NumVertices()),
          matchLeftNeighborSum(query_.NumEdges()) {
//    globalStarIndex.resize(query_.NumEdges());
    for (uint i = 0; i < query_.NumEdges(); ++i) {
        order_vs_[i].resize(query_.NumVertices());//节点个数
        order_csrs_[i].resize(query_.NumEdges() + 1);//边的个数+1，
        order_offs_[i].resize(query_.NumVertices(), 0);//节点个数，初始化为0
        globalStarIndex[i].resize(query_.NumVertices());
        order_vertex_index[i].resize(query_.NumVertices());
        rightNeighbor[i].resize(query_.NumVertices());
        matchLeftNeighborSum[i].resize(query_.NumVertices());
    }
    for(uint i=0;i<query_.NumVertices();i++){
        matchCandidate[i].reserve(10000);
    }
}

Graphflow::~Graphflow() noexcept {
    for (MatchRecord *item: topKSet) {
        delete item;
        item = nullptr;
    }
    for (int i = 0; i < query_.NumEdges(); i++) {
        for (StarGraph *s: globalStarIndex[i]) {
            delete s;
        }
    }
}

void Graphflow::Preprocessing()//预处理过程
{
    this->data_.InitLabelIndex();
    GenerateMatchingOrder();
    this->query_.InitLabelIndex();
#ifdef ISOLATE
    this->query_.InitMatchOrderType(this->order_vs_, this->rightNeighbor);
#endif
    createLabelToQueryVertex();
    CreateStarIndex();
    std::cout << "Preprocess end" << endl;
}

void Graphflow::updateStarIndex(uint match_index, uint caddidate_v, const std::vector<uint> &canditeQueryVertexs) {
    std::vector<int> &result = globalVkMatchUk[caddidate_v][match_index];
    result.resize(canditeQueryVertexs.size());
    std::vector<Neighbor> &vN = this->data_.vNeighbors[caddidate_v];
    float sumWeight = 0;
    for (int i = 0; i < canditeQueryVertexs.size(); i++) {
        uint candidate_u = canditeQueryVertexs[i];
        sumWeight = 0;
        int vertex_index = order_vertex_index[match_index][candidate_u];
        if (vertex_index == 0)
            continue;
        StarGraph *s = globalStarIndex[match_index][vertex_index];
        std::vector<ForwardNeighbor *> &queryVetex = s->GetqueryVertex();
        int leftvN = 0;
        int rightqV = 0;
        int flag = 1;
        int qvSize = queryVetex.size();
        int vNsize = vN.size();
        while (leftvN < vNsize && rightqV < qvSize) {
            if (vN[leftvN].GetelabelAndVertexLabel() < queryVetex[rightqV]->GetelabelAndVertexLabel()) {
                flag = 0;
                break;
            }
            while (vN[leftvN].GetelabelAndVertexLabel() > queryVetex[rightqV]->GetelabelAndVertexLabel()) {
                leftvN++;
                if (leftvN >= vN.size()) {
                    flag = 0;
                    break;
                }
            }
            if (!flag)
                break;
            if (vN[leftvN].GetelabelAndVertexLabel() == queryVetex[rightqV]->GetelabelAndVertexLabel()) {
                uint tovId = queryVetex[rightqV]->GetVetexId();
                pair<uint, uint> tmppair = vN[leftvN].GetelabelAndVertexLabel();
                float edgeweight = vN[leftvN].GetEdgeWeight();
                sumWeight += edgeweight;
                rightqV++;
                leftvN++;
            }
        }
        if (!flag) {
            result[i] = s->getStarMaxWeight();
        }
        if (rightqV == qvSize) {
            result[i] = sumWeight;
            if (s->getStarMaxWeight() == queryVetex.size() * mw || s->getStarMaxWeight() < sumWeight) {
                s->setStarMaxWeight(sumWeight);
                s->setMatchDataVertexId(caddidate_v);
            }
        }
    }
}

void Graphflow::updateStarIndex(uint match_index, uint caddidate_v, uint candidate_u, int candidate_v_index) {
    std::vector<int> &result = globalVkMatchUk[caddidate_v][match_index];
    int vertex_index = order_vertex_index[match_index][candidate_u];
    StarGraph *s = globalStarIndex[match_index][vertex_index];
    const std::vector<ForwardNeighbor *> &queryVetex = s->GetqueryVertex();
    std::vector<Neighbor> &vN = this->data_.vNeighbors[caddidate_v];
    int leftvN = 0;
    int rightqV = 0;
    int vNsize = vN.size();
    int qVsize = queryVetex.size();
    int flag = 1;
    float sumweight = 0;
    while (leftvN < vNsize && rightqV < qVsize) {
        if (vN[leftvN].GetelabelAndVertexLabel() < queryVetex[rightqV]->GetelabelAndVertexLabel()) {
            flag = 0;
            break;
        }

        while (vN[leftvN].GetelabelAndVertexLabel() > queryVetex[rightqV]->GetelabelAndVertexLabel()) {
            leftvN++;
            if (leftvN >= vN.size()) {
                flag = 0;
                break;
            }
        }
        if (!flag)
            break;
        if (vN[leftvN].GetelabelAndVertexLabel() == queryVetex[rightqV]->GetelabelAndVertexLabel()) {
            float edgeweight = vN[leftvN].GetEdgeWeight();
            rightqV++;
            leftvN++;
            sumweight += edgeweight;
        }
    }
    if (!flag) {
        return;
    } else if (rightqV == qVsize) {
        //globalVkMatchUk更新
        result[candidate_v_index] = sumweight;
        // globalStarIndex更新
        if (s->getStarMaxWeight() == queryVetex.size() * mw || s->getStarMaxWeight() < sumweight) {
            s->setStarMaxWeight(sumweight);
            s->setMatchDataVertexId(caddidate_v);
        }
    }
}

void Graphflow::updateStarIndex(uint match_index, uint caddidate_v, uint candidate_u) {
    int vertex_index = order_vertex_index[match_index][candidate_u];
    StarGraph *s = globalStarIndex[match_index][vertex_index];
    const std::vector<ForwardNeighbor *> &queryVetex = s->GetqueryVertex();
    std::vector<Neighbor> &vN = this->data_.vNeighbors[caddidate_v];
    int leftvN = 0;
    int rightqV = 0;
    float sumweight = 0;
    int vNsize=vN.size();
    int qVSize=queryVetex.size();
    int flag=1;
    while (leftvN < vNsize && rightqV < qVSize) {
        if (vN[leftvN].GetelabelAndVertexLabel() < queryVetex[rightqV]->GetelabelAndVertexLabel()) {
           flag=0;
            break;
        }
        while (vN[leftvN].GetelabelAndVertexLabel() > queryVetex[rightqV]->GetelabelAndVertexLabel()) {
            leftvN++;
            if (leftvN >= vN.size()) {
                flag = 0;
                break;
            }
        }
        if (!flag)
            break;
        if (vN[leftvN].GetelabelAndVertexLabel() == queryVetex[rightqV]->GetelabelAndVertexLabel()) {
            float edgeweight = vN[leftvN].GetEdgeWeight();
            rightqV++;
            leftvN++;
            sumweight += edgeweight;
        }
    }
    if(!flag){
        return;
    }
    else if (rightqV == qVSize) {
        if (s->getStarMaxWeight() == queryVetex.size() * mw || s->getStarMaxWeight() < sumweight) {
            s->setStarMaxWeight(sumweight);
            s->setMatchDataVertexId(caddidate_v);
        }
    }
}

float Graphflow::GetBackWeight(uint order_index, uint depth) {
    float sum = 0;
    uint n = query_.NumVertices();
    std::vector<uint> &matchOrder = this->order_vs_[order_index];
    for (int i = depth; i < n; i++) {
        sum += LocalStarIndex[i];
    }
    return sum;
}

void Graphflow::CreateStarIndex() {
    int n = data_.NumVertices();
    int m = query_.NumEdges();
    for (int i = 0; i < n; i++) {
        int label = data_.GetVertexLabel(i);
        if (label == -1)
            continue;
        const std::vector<uint> &candidate_us = labelToQueryVertex[label];
        if (candidate_us.size() == 0)
            continue;
#ifdef LOCAL
        globalVkMatchUk[i].resize(m);
          for(int j=0;j<m;j++){
            updateStarIndex(j,i,candidate_us);
        }
#endif
#ifdef GLOBAL
        for (int j = 0; j < m; j++) {
            for (auto candidate_u: candidate_us) {
                int vertex_index = order_vertex_index[j][candidate_u];
                if(vertex_index==0)
                    continue;
                updateStarIndex(j, i, candidate_u);
            }

        }
#endif

    }
}

vector<int> Graphflow::EdgeisInMatchOrder(uint v1, uint v2, uint v1label, uint v2label, uint velabel) {
    vector<int> result;
    for (int i = 0; i < order_vs_.size(); i++) {
        uint u1 = order_vs_[i][0];
        uint u2 = order_vs_[i][1];
        uint u1label = query_.GetVertexLabel(u1);
        uint u2label = query_.GetVertexLabel(u2);
        uint qlabel = std::get<2>(query_.GetEdgeLabel(u1, u2));
        if ((v1label == u1label && v2label == u2label && velabel == qlabel) ||
            (v1label == u2label && v2label == u1label && velabel == qlabel)) {
            result.emplace_back(i);
        }
    }
    return result;
}

vector<int> Graphflow::EdgeisInMatchOrder(Edge *edge) {
    uint v1 = edge->GetV1();
    uint v2 = edge->GetV2();
    uint v1label = edge->GetV1Label();
    uint v2label = edge->GetV2Label();
    uint velabel = edge->GeteLabel();
    vector<int> result;
    for (int i = 0; i < order_vs_.size(); i++) {
        uint u1 = order_vs_[i][0];
        uint u2 = order_vs_[i][1];
        uint u1label = query_.GetVertexLabel(u1);
        uint u2label = query_.GetVertexLabel(u2);
        uint qlabel = std::get<2>(query_.GetEdgeLabel(u1, u2));
        if ((v1label == u1label && v2label == u2label && velabel == qlabel) ||
            (v1label == u2label && v2label == u1label && velabel == qlabel)) {
            result.emplace_back(i);
        }
    }
    return result;
    //返回其在哪一个matchorder中被匹配
}

void Graphflow::GenerateMatchingOrder() {
    // generate the initial matching order, order_*s_[0]
    std::vector<bool> visited(query_.NumVertices(), false);
    uint max_degree = 0u;
    //首先找到的是度最大的节点
    for (size_t i = 0; i < query_.NumVertices(); i++) {
        if (query_.GetDegree(i) > max_degree) {
            max_degree = query_.GetDegree(i);
            order_vs_[0][0] = i;
            order_vertex_index[0][i] = 0;
        }
    }
    visited[order_vs_[0][0]] = true;

    // loop over all remaining positions of the order
    for (uint i = 1; i < query_.NumVertices(); ++i) {
        uint max_adjacent = 0;
        uint max_adjacent_u = NOT_EXIST;
        //找到不在序列中，但是在序列中的邻居数量最多的顶点，添加到排序中
        for (size_t j = 0; j < query_.NumVertices(); j++) {
            uint cur_adjacent = 0u;
            if (visited[j]) continue;

            auto &q_nbrs = query_.GetNeighbors(j);
            for (auto &other: q_nbrs)
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
        order_vertex_index[0][max_adjacent_u] = i;
        visited[max_adjacent_u] = true;
        order_offs_[0][i] = order_offs_[0][i - 1];
        auto &q_nbrs = query_.GetNeighbors(max_adjacent_u);
        StarGraph *s = new StarGraph();
        for (auto &other: q_nbrs) {
            if (visited[other]) {
                uint qlabel = std::get<2>(this->query_.GetEdgeLabel(max_adjacent_u, other));
                //  std::cout<<"globalStarIndex[0]["<<i<<"] other:"<<other<<" other label"<< this->query_.GetVertexLabel(other)<<endl;
                ForwardNeighbor *forwardNeighbor = new ForwardNeighbor(other, this->query_.GetVertexLabel(other),
                                                                       qlabel);
                s->AddForwardNeighbor(forwardNeighbor);
                order_csrs_[0][order_offs_[0][i]++] = other;
            }
        }
        s->InitalmaxWeight();
        globalStarIndex[0][i] = s;

    }

    // generate other incremental matching orders
    for (uint i = 1; i < query_.NumEdges(); ++i) {
        std::vector<bool> visited(query_.NumVertices(), false);

        // get the first edge
        std::vector<uint>::iterator it = std::lower_bound(
                order_offs_[0].begin(), order_offs_[0].end(), i + 1
        );
        uint tmp = *(order_vs_[0].begin() + std::distance(order_offs_[0].begin(), it));
        order_vs_[i][0] = tmp;
        order_vertex_index[i][tmp] = 0;
        order_vs_[i][1] = order_csrs_[0][i];
        order_vertex_index[i][order_csrs_[0][i]] = 1;
        order_csrs_[i][0] = order_vs_[i][0];
        StarGraph *s = new StarGraph();
        uint qlabel = std::get<2>(this->query_.GetEdgeLabel(order_vs_[i][0], order_vs_[i][1]));

        ForwardNeighbor *forwardNeighbor = new ForwardNeighbor(order_vs_[i][0],
                                                               this->query_.GetVertexLabel(order_vs_[i][0]), qlabel);
        s->AddForwardNeighbor(forwardNeighbor);
        s->InitalmaxWeight();
        globalStarIndex[i][1] = (s);

        visited[order_vs_[i][0]] = true;
        visited[order_vs_[i][1]] = true;

        order_offs_[i][2] = order_offs_[i][1] = 1;
        for (uint j = 2; j < query_.NumVertices(); ++j) {
            uint max_adjacent = 0;
            uint max_adjacent_u = NOT_EXIST;
            for (size_t k = 0; k < query_.NumVertices(); k++) {
                uint cur_adjacent = 0u;
                if (visited[k]) continue;

                auto &q_nbrs = query_.GetNeighbors(k);
                for (auto &other: q_nbrs)
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
            order_vertex_index[i][max_adjacent_u] = j;
            visited[max_adjacent_u] = true;

            order_offs_[i][j] = order_offs_[i][j - 1];
            StarGraph *s = new StarGraph();
            auto &q_nbrs = query_.GetNeighbors(max_adjacent_u);
            for (auto &other: q_nbrs) {
                if (visited[other]) {
                    // std::cout<<"globalStarIndex["<<i<<"]"<<"["<<j<<"] "<<"other:"<<other<<" other label"<< this->query_.GetVertexLabel(other)<<endl;
                    order_csrs_[i][order_offs_[i][j]++] = other;
                    qlabel = std::get<2>(this->query_.GetEdgeLabel(max_adjacent_u, other));
                    ForwardNeighbor *forwardNeighbor = new ForwardNeighbor(other, this->query_.GetVertexLabel(other),
                                                                           qlabel);
                    s->AddForwardNeighbor(forwardNeighbor);
                }
            }
            s->InitalmaxWeight();
            globalStarIndex[i][j] = (s);
        }
    }
    //对globalStarIndex中的所有匹配序的所有节点的前向邻居按照<el,vl>排序
    for (int i = 0; i < query_.NumEdges(); i++) {
        for (int j = 1; j < query_.NumVertices(); j++) {
            StarGraph *s = globalStarIndex[i][j];
            std::vector<ForwardNeighbor *> &globalIndex = s->GetqueryVertex();
            std::sort(globalIndex.begin(), globalIndex.end(), ForwardNeighborcmp);
        }
    }



    //创建所有节点的右邻居数组
    if (print_preprocessing_results_) {
        std::cout << "matching order: " << std::endl;
        std::cout << "-vertex(backward neighbors)-\n";
        for (uint i = 0; i < query_.NumEdges(); ++i) {
            std::cout << "#" << i << ": ";
            for (uint j = 0; j < query_.NumVertices(); ++j) {
                std::cout << order_vs_[i][j];
                if (j == 0) {
                    // this->query_.forwardNeighbors[i][j]={};
                    std::cout << "-";
                    continue;
                }
                std::vector<ForwardNeighbor> currentQueryNeighbors;
                matchLeftNeighborSum[i][j] = order_offs_[i][j] - order_offs_[i][j - 1];
                for (uint k = order_offs_[i][j - 1]; k < order_offs_[i][j]; k++) {
                    uint toVertexId = order_csrs_[i][k];
                    uint toVertexIndex = order_vertex_index[i][toVertexId];
                    uint toVertexLabel = query_.GetVertexLabel(order_csrs_[i][k]);
                    uint edgelabel = std::get<2>(query_.GetEdgeLabel(order_vs_[i][j], toVertexId));
                    ForwardNeighbor f(toVertexIndex, toVertexId, toVertexLabel, edgelabel);
                    currentQueryNeighbors.push_back(f);
                    uint rightNeighborLabel=query_.GetVertexLabel(order_vs_[i][j]);
                    Neighbor rNeighbor(order_vs_[i][j],rightNeighborLabel,edgelabel);
                    auto lower = std::lower_bound(rightNeighbor[i][toVertexId].begin(), rightNeighbor[i][toVertexId].end(), rNeighbor, CompareNeighbors);
                    rightNeighbor[i][toVertexId].insert(lower,rNeighbor);
                   // rightNeighbor[i][toVertexId].emplace_back(order_vs_[i][j]);
                    if (k == order_offs_[i][j - 1]) std::cout << "(";
                    std::cout << order_csrs_[i][k];
                    if (k != order_offs_[i][j] - 1) std::cout << ",";
                    else std::cout << ")";
                }
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

    if (!io::file_exists(path.c_str())) {
        std::fstream fp(path, std::ios::out);
        for (auto t: topKSet) {
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
    _ss1 << "Initial Top k" << std::endl;
    // std::sort(topKSet.begin(),topKSet.end(), matchRecordCmp);
    for (auto d: topKSet) {
        if (d != NULL) {
            _ss1 << d->toString();
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
    if (!isUpdateIntopkset) {
        return;
    }
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
    _ss << "after insert " << std::endl;
    for (auto d: topKSet) {
        _ss << d->toString();
        Log::track2(_ss);
        _ss.clear();
        _ss.str("");
    }
#endif
}

void Graphflow::searchMatches(int depth, uint matchorderindex, searchType flag) {
    // start=Get_Time();
    //1.找到前向邻居，找到候选解
    std::vector<uint> &matchOrder = this->order_vs_[matchorderindex];
    uint queryVertex = matchOrder[depth];
    uint queryVertexLabel = this->query_.GetVertexLabel(queryVertex);
    // const auto & fneighbors=this->query_.forwardNeighbors[matchorderindex][depth];
    std::vector<SingleCandidate> &singleVertexCandidate = this->matchCandidate[depth];
    std::vector<SingleCandidate> copySingleVertexCandidate = this->matchCandidate[depth];
    /* singleVertexCandidate.clear();
     float maxweight=0;
     uint u_min=NOT_EXIST;
     uint u_min_label=NOT_EXIST;
     uint u_min_size=UINT_MAX;
     uint u_min_index=NOT_EXIST;

     //find u_min
     const auto &q_nbrs=query_.GetNeighbors(queryVertex);
     const auto &q_nbr_labels=query_.GetNeighborLabels(queryVertex);
     for(uint i=0u;i<q_nbrs.size();i++){
         const uint u_other=q_nbrs[i];
         const uint u_other_index=order_vertex_index[matchorderindex][u_other];
         const uint u_other_label=q_nbr_labels[i];
         if(u_other_index>=depth)
             continue;
         const uint cur_can_size=data_.GetNeighbors(match[u_other_index].getVertexId()).size();
         if(cur_can_size<u_min_size){
             u_min_size=cur_can_size;
             u_min=u_other;
             u_min_label=u_other_label;
             u_min_index=order_vertex_index[matchorderindex][u_min];
         }
     }
     uint u_min_id=match[u_min_index].getVertexId();
     const auto &u_min_nbrs = data_.GetNeighbors(u_min_id);
     const auto &u_min_nbr_labels = data_.GetNeighborLabels(u_min_id);

     float tmp;
     for (uint i = 0u; i < u_min_nbrs.size(); i++) {
         const uint v = u_min_nbrs[i];
         tmp = 0;
         // 1. check labels
         if (
                 data_.GetVertexLabel(v) != query_.GetVertexLabel(queryVertex) ||
                 u_min_nbr_labels[i] != u_min_label
                 )
             continue;
         // 2. check if visited
         if (!homomorphism_ && visited_[v]) continue;
         //代表u_min和v之前的权值之和
         tmp += data_.GetEdgeWeight(u_min_id, v);

         //考虑在已经匹配的查询图u的邻居中，除了u_min以外的其他邻居，在数据图中都必须与u的匹配点v有边相连
         //因此若u_other在数据图中没有邻居是v。则也不匹配
         // 3. check if joinable
         bool joinable = true;
         for (uint j = 0u; j < q_nbrs.size(); j++) {
             const uint u_other = q_nbrs[j];
             const uint u_other_labels = q_nbr_labels[j];
             const uint u_other_index = order_vertex_index[matchorderindex][u_other];
             const uint u_other_id = match[u_other_index].getVertexId();
             if (u_other_index >= depth || u_other == u_min) continue;
             auto it = std::lower_bound(data_.GetNeighbors(u_other_id).begin(), data_.GetNeighbors(u_other_id).end(),
                                        v);
             uint dis = std::distance(data_.GetNeighbors(u_other_id).begin(), it);
             if (
                     it == data_.GetNeighbors(u_other_id).end() ||
                     *it != v ||
                     data_.GetNeighborLabels(u_other_id)[dis] != u_other_labels
                     ) {
                 joinable = false;
                 break;
             }
             tmp += data_.GetEdgeWeight(u_other_id, v);
         }
         if (!joinable) continue;
         singleVertexCandidate.emplace_back(v,tmp);
         if(tmp>maxweight){
             maxweight=tmp;
         }
     }*/
    //取singleCandidate的交集部分
    /*getIntersetSingleCandidate(singleVertexCandidate, matchorderindex, depth);
    if (singleVertexCandidate.size() == 0) {
        this->matchCandidate[depth] = copySingleVertexCandidate;
        return;
    }*/
    //Print_Time2("findCandidate ",start);
    //若不是孤立节点 那么从候选解决中加入此节点，继续递归
    /* if(singleVertexCandidate.size()==0){
         return;
     }*/
    total_densityFilter_time.StartTimer();
    densityFilter(matchorderindex, depth, singleVertexCandidate);
    total_densityFilter_time.StopTimer();
    if (singleVertexCandidate.size() == 0) {
        this->matchCandidate[depth] = copySingleVertexCandidate;
        return;
    }
    if(isInsert)
        IsearchSpace+=singleVertexCandidate.size();
    else
        DsearchSpace+=singleVertexCandidate.size();

    if (depth == query_.NumVertices() - 1) {
        //add matchresult;
        std::sort(singleVertexCandidate.begin(), singleVertexCandidate.end());
        for (const SingleCandidate &single: singleVertexCandidate) {
            uint dataV = single.getVertexId();
            if (visited_[dataV])
                continue;
            float sumWeight = this->match[depth - 1].getSumWeight();
            this->match[depth].setVertexId(dataV);
            sumWeight += single.getSumWeight();
            this->match[depth].setSumWeight(sumWeight);
            int n = query_.NumVertices();
            std::vector<uint> m(n);
            for (int i = 0; i < match.size(); i++) {
                m[order_vs_[matchorderindex][i]] = match[i].getVertexId();
            }
            float density = sumWeight / (sqrt(n) * (n - 1));
            MatchRecord *record = new MatchRecord(density, m);
#ifdef DEBUG
            stringstream _ss;
            _ss<<record->toString();
            Log::track1(_ss);
#endif
            allMatchFind++;
            int matchResult = addMatchRecords(record);
            if (matchResult == 1) {
                if (flag == positive) {
                    num_positive_results_++;
                    numAddTopk++;
                }
            } else if (matchResult == 3) {
                this->matchCandidate[depth] = copySingleVertexCandidate;
                this->match[depth].clearSingleCandidate();
                return;
            }
        }
        //clear candidate;
        this->matchCandidate[depth] = copySingleVertexCandidate;
        this->match[depth].clearSingleCandidate();
        return;
    } else {
        for (int i = 0; i < singleVertexCandidate.size(); i++) {
            uint dataV = singleVertexCandidate[i].getVertexId();
            if (visited_[dataV])
                continue;
            float weight = singleVertexCandidate[i].getSumWeight();
            matchVertex(0, depth, dataV, weight);
            const std::vector<Neighbor> &uk_neighbor = rightNeighbor[matchorderindex][queryVertex];
            std::vector<std::vector<SingleCandidate>> copyCandidate(uk_neighbor.size());
            std::vector<int> copyLocalStarIndex(query_.NumVertices());
            for (int i = 0; i < uk_neighbor.size(); i++) {
                int uk_neighbor_index = order_vertex_index[matchorderindex][uk_neighbor[i].getVertexId()];
                copyCandidate[i] = matchCandidate[uk_neighbor_index];
                copyLocalStarIndex[uk_neighbor_index] = LocalStarIndex[uk_neighbor_index];
            }
            // std::cout<<"depth :"<<depth<<" data:"<<dataV<<endl;
            total_updaterightNeighborCandidate_time.StartTimer();
            bool isNull = updaterightNeighborCandidate(matchorderindex, queryVertex, 0, false, dataV, uk_neighbor);
            total_updaterightNeighborCandidate_time.StopTimer();
            if (isNull) {
                for (int i = 0; i < uk_neighbor.size(); i++) {
                    int uk_neighbor_index = order_vertex_index[matchorderindex][uk_neighbor[i].getVertexId()];
                    matchCandidate[uk_neighbor_index] = copyCandidate[i];
                    LocalStarIndex[uk_neighbor_index] = copyLocalStarIndex[uk_neighbor_index];
                    //todo传值
                }
                this->visited_[dataV] = false;
                continue;
            }
            //copy SingleCandidate
            //updateweight;
            searchMatches(depth + 1, matchorderindex, flag);
            for (int i = 0; i < uk_neighbor.size(); i++) {
                int uk_neighbor_index = order_vertex_index[matchorderindex][uk_neighbor[i].getVertexId()];
                matchCandidate[uk_neighbor_index] = copyCandidate[i];
                LocalStarIndex[uk_neighbor_index] = copyLocalStarIndex[uk_neighbor_index];
            }
            this->visited_[dataV] = false;
        }
        this->matchCandidate[depth] = copySingleVertexCandidate;
        this->match[depth].clearSingleCandidate();
    }
}


//flag==0为initial flag=1为update
void Graphflow::FindMatches(uint flag, uint order_index, uint depth, std::vector<uint> m, size_t &num_results,
                            float density_s, searchType type) {
    if (reach_time_limit) return;
    if(isInsert)
        IsearchSpace++;
    else
        DsearchSpace++;
    if (flag == 1) {
        float back_max_result = suffixMax[depth];
        uint n = query_.NumVertices();
        if (topKSet.size() == k) {
            float weight = topKSet.back()->getDensity();
            float tmpdensity = (density_s + back_max_result)/(sqrt(n) * (n - 1));
            if (tmpdensity < weight)
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
        total_updaterightNeighborCandidate_time.StartTimer();
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
            if(isInsert)
                IdeterminCandite++;
            else
                DdeterminCandite++;
            tmp += data_.GetEdgeWeight(m[u_other], v);
        }
        total_updaterightNeighborCandidate_time.StopTimer();
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
            MatchRecord *record = new MatchRecord(lastds, m);
            allMatchFind++;
            int matchResult = addMatchRecords(record);
            if (matchResult == 1) {
                if (type == positive) {
                    num_positive_results_++;
                    numAddTopk++;
                }
            }
            if (print_enumeration_results_) {
                for (auto j: m) {

                    std::cout << j << " ";
                }
            }

        } else {
            size_t num_results_before_recursion = num_results;
            FindMatches(flag, order_index, depth + 1, m, num_results, density_s, type);
            if (num_results == num_results_before_recursion) {
                num_intermediate_results_without_results_++;
            }

        }


        visited_[v] = false;
        m[u] = UNMATCHED;
        density_s -= tmp;

//            tmin = copytmin;//回溯将tmin转为初始状态
        if (num_results >= max_num_results_) return;
        if (reach_time_limit) return;
    }
    if (candidate_empty) num_intermediate_results_with_empty_candidate_set_++;
}


int Graphflow::addMatchRecords(MatchRecord *r) {

    int n = topKSet.size();
    if (n < k) {
/*
       if(n==0){
           topKSet.insert(topKSet.begin(),r);
           isUpdateIntopkset= true;
           return 1;
       }
        auto insertPosition = std::lower_bound(topKSet.begin(), topKSet.end(), r, compareMatchRecord);
       if(insertPosition!=topKSet.end()&&(**insertPosition)==(*r)){
            delete r;
            return 2;
        }
        else{
            topKSet.insert(insertPosition,r);
            isUpdateIntopkset= true;
            return 1;
        }*/

         for (int j = n - 1; j >= 0; j--) {
                 if ((*topKSet[j]) == (*r)) {
                     delete r;
                     return 2;
                 }
             }
        for (int j = n - 1; j >= 0; j--) {
            if ((*topKSet[j]) > (*r)) {
                topKSet.insert(topKSet.begin() + j + 1, r);
                break;
            }
            if (j == 0) {
                topKSet.insert(topKSet.begin(), r);
            }
        }
        if (n == 0) {
            topKSet.insert(topKSet.begin(), r);
        }
        isUpdateIntopkset = true;
        return 1;
    } else {
       /* auto insertPosition = std::lower_bound(topKSet.begin(), topKSet.end(), r, compareMatchRecord);
        if(insertPosition==topKSet.end()){
            delete r;
            return 3;
        }
        else if((**insertPosition)==(*r)){
            delete r;
            return 2;
        }
        else{
            topKSet.insert(insertPosition,r);
            MatchRecord* d=topKSet.back();
            delete d;
            topKSet.pop_back();
            isUpdateIntopkset= true;
            return 1;

        }*/
        for (int j = n - 1; j >= 0; j--) {
            if ((*topKSet[j]) == (*r)) {
                delete r;
                return 2;
            }
        }
        MatchRecord *d = topKSet.back();
        if ((*r) > (*d)) {
            delete d;
            topKSet.pop_back();
            int m = topKSet.size();
            for (int j = m - 1; j >= 0; j--) {
                if ((*topKSet[j]) > (*r)) {
                    topKSet.insert(topKSet.begin() + j + 1, r);
                    break;
                }
                if (j == 0) {
                    topKSet.insert(topKSet.begin(), r);
                }
            }
            isUpdateIntopkset = true;
            return 1;
        } else {
            delete r;
            return 3;
        }
        //return true;
    }
}

void Graphflow::InitialMatching(const std::string &path) {

    if (!io::file_exists(path.c_str())) {
        std::cout << "the file not exit " << path << std::endl;
        std::vector<uint> m(query_.NumVertices(), UNMATCHED);
        uint flag = 0;
        float density_s = 0;
        uint tmin = INT_MAX;
        uint order_index = 0;
        uint depth = 1;
        //#pragma omp parallel for num_threads(10) firstprivate(m) firstprivate(visited_) firstprivate(density_s) firstprivate(flag) firstprivate(tmin) firstprivate(order_index) firstprivate(depth)
        for (size_t i = 0; i < data_.NumVertices(); i++) {
            //std::cout<<"thread id"<<omp_get_thread_num<<endl;
            if (data_.GetVertexLabel(i) != NOT_EXIST) {
#ifdef PRINT_DEBUG
                stringstream _ss;
                _ss<<"vertex id:"<<i<<std::endl;
                Log::track1(_ss);
#endif
                if (query_.GetVertexLabel(order_vs_[0][0]) == data_.GetVertexLabel(i)) {
                    m[order_vs_[0][0]] = i;
                    visited_[i] = true;

                    FindMatches(flag, order_index, depth, m, num_initial_results_, density_s, positive);
                    visited_[i] = false;
                    m[order_vs_[0][0]] = UNMATCHED;
                }

            }

        }
    } else {
        std::ifstream ifs2(path);
        std::cout << "load topk from file...." << std::endl;
        char type;
        while (ifs2 >> type) {
            if (type == 't') {
                float density;
                uint tmp;
                std::vector<uint> m;
                ifs2 >> density;
                for (int i = 0; i < query_.NumVertices(); i++) {
                    ifs2 >> tmp;
                    m.emplace_back(tmp);
                }
                MatchRecord *matchRecord = new MatchRecord(density, m);
                addMatchRecords(matchRecord);
            }
        }
    }
}

void Graphflow::AddEdgeWithGlobalIndex(uint v1, uint v2, uint label, float weight, uint timestamp) {
    total_update_globalIndex_time.StartTimer();
    numAddTopk = 0;
    allMatchFind = 0;
    data_.AddEdge(v1, v2, label, weight, timestamp, 1);
    data_.UpdateLabelIndex(v1, v2, label, 1);
    uint v1label = data_.GetVertexLabel(v1);
    uint v2label = data_.GetVertexLabel(v2);
    vector<int> match = EdgeisInMatchOrder(v1, v2, v1label, v2label, label);
    if(match.size()==0){
        return;
    }
    //for each edge<v1,v2>matches u1-->u2 or u2-->u1
    for (int i = 0; i < match.size(); i++) {
        uint u1 = order_vs_[match[i]][0];
        uint u2 = order_vs_[match[i]][1];
        //对其他的边
        for (int j = 0; j < query_.NumEdges(); j++) {
            uint candidate_u = order_vertex_index[j][u1] > order_vertex_index[j][u2] ? u1 : u2;
            if (query_.GetVertexLabel(candidate_u) == v1label) {
                numupdatestar++;
                updateStarIndex(j, v1, candidate_u);
            }
            if (query_.GetVertexLabel(candidate_u) == v2label) {
                numupdatestar++;
                updateStarIndex(j, v2, candidate_u);
            }
        }
    }
    total_update_globalIndex_time.StopTimer();

    total_search_time.StartTimer();
    std::vector<uint> m(query_.NumVertices(), UNMATCHED);
    isUpdateIntopkset = false;
    if (max_num_results_ == 0) return;
    size_t num_results = 0ul;
    for (auto ma:match) {
        uint u1 = order_vs_[ma][0], u2 = order_vs_[ma][1];
        auto temp_q_labels = query_.GetEdgeLabel(u1, u2);
        float density_s = weight;
        //compute suffix array
        for (int j = query_.NumVertices() - 1; j > 0; j--) {
            if (j == query_.NumVertices() - 1) {
                suffixMax[j] = globalStarIndex[ma][j]->getStarMaxWeight();
            } else {
                suffixMax[j] = suffixMax[j + 1] + globalStarIndex[ma][j]->getStarMaxWeight();
            }
        }
        // check if any query edge match (v1 --> v2)
        if (
                std::get<0>(temp_q_labels) == data_.GetVertexLabel(v1) &&
                std::get<1>(temp_q_labels) == data_.GetVertexLabel(v2) &&
                std::get<2>(temp_q_labels) == label
                ) {
            m[u1] = v1;
            m[u2] = v2;
            if(isInsert)
            IsearchSpace+=2;
            visited_[v1] = true;
            visited_[v2] = true;
            FindMatches(1, ma, 2, m, num_results, density_s, positive);
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
            if(isInsert)
            IsearchSpace+=2;
            visited_[v2] = true;
            visited_[v1] = true;
            FindMatches(1, ma, 2, m, num_results, density_s, positive);
            visited_[v2] = false;
            visited_[v1] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;

            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) goto END_ENUMERATION;
        }
    }
    END_ENUMERATION:
//    num_positive_results_ += num_results;
    total_search_time.StopTimer();
    total_print_time.StartTimer();
    sumAllMatchFind+=allMatchFind;
    std::cout << "num add top k:" << numAddTopk << endl;
    std::cout << "all match find:" << num_results << endl;
    updateTopK();
    total_print_time.StopTimer();
}

//动态的加边操作
void Graphflow::AddEdge(uint v1, uint v2, uint label, float weight, uint timestamp) {
#ifdef DEBUG
    stringstream _ss;
    for(int i=0;i<query_.NumVertices();i++){
        if(matchCandidate[i].size()>0){
            _ss<<"3 matchCandidate["<<i<<"] not null"<<endl;
            Log::track1(_ss);
            _ss.clear();
            _ss.str("");
        }
    }
#endif
    total_update_globalIndex_time.StartTimer();
    numAddTopk = 0;
    allMatchFind = 0;
    data_.AddEdge(v1, v2, label, weight, timestamp, 1);
    data_.UpdateLabelIndex(v1, v2, label, 1);
    uint v1label = data_.GetVertexLabel(v1);
    uint v2label = data_.GetVertexLabel(v2);
    vector<int> match = EdgeisInMatchOrder(v1, v2, v1label, v2label, label);
    if (match.size() == 0) {
        return;
    }
    for (int i = 0; i < match.size(); i++) {
        uint u1 = order_vs_[match[i]][0];
        uint u2 = order_vs_[match[i]][1];
        //对其他的边
        for (int j = 0; j < query_.NumEdges(); j++) {
            uint candidate_u = order_vertex_index[j][u1] > order_vertex_index[j][u2] ? u1 : u2;
            if (query_.GetVertexLabel(candidate_u) == v1label) {
                numupdatestar++;
                int candidate_index = queryVertexIndexInlabel[candidate_u];
                updateStarIndex(j, v1, candidate_u, candidate_index);
            }
            if (query_.GetVertexLabel(candidate_u) == v2label) {
                numupdatestar++;
                int candidate_index = queryVertexIndexInlabel[candidate_u];
                updateStarIndex(j, v2, candidate_u, candidate_index);
            }
        }
    }

    total_update_globalIndex_time.StopTimer();
    //Print_Time2("UpdateIndex ", start);
    size_t num_results = 0ul;

    total_search_time.StartTimer();
    isUpdateIntopkset = false;
    for (auto m: match) {
#ifdef DEBUG
        stringstream _ss;
        for(int i=0;i<query_.NumVertices();i++){
            if(matchCandidate[i].size()>0){
                _ss<<"2 matchCandidate["<<i<<"] not null orderindex:"<<m<<endl;
                Log::track1(_ss);
                _ss.clear();
                _ss.str("");
            }
        }
#endif
        uint u1 = order_vs_[m][0];
        uint u2 = order_vs_[m][1];
        uint u1label = this->query_.GetVertexLabel(u1);
        uint u2label = this->query_.GetVertexLabel(u2);
        uint v1label = this->data_.GetVertexLabel(v1);
        uint v2label = this->data_.GetVertexLabel(v2);
        float weight = this->data_.GetEdgeWeight(v1, v2);
        const auto &matchOrder = order_vs_[m];
        InitialLocalIndex(m);
        if (v1label != v2label) {
            if (v1label != u1label) {
                swap(v1, v2);
            }
            //todo
            if (this->LabelFilter(v1, u1) && this->LabelFilter(v2, u2)) {
                bool flag=SearchMatchesWithEdge(m,v1,v2,weight,u1,u2,positive);
                if(flag)
                    continue;
            }
        } else {
            for (int i = 0; i < 2; i++) {
                if (this->LabelFilter(v1, u1) && this->LabelFilter(v2, u2)) {
                    bool flag=SearchMatchesWithEdge(m,v1,v2,weight,u1,u2,positive);
                    if(flag)
                    {
                        std::swap(v1, v2);
                        continue;
                    }
                }
                std::swap(v1, v2);
            }
        }
#ifdef DEBUG
        for(int i=0;i<query_.NumVertices();i++){
            if(matchCandidate[i].size()>0){
                _ss<<"2 matchCandidate["<<i<<"] not null orderindex:"<<m<<endl;
                Log::track1(_ss);
                _ss.clear();
                _ss.str("");
            }
        }
#endif
    }
    total_search_time.StopTimer();
    //Print_Time2("SearchMatches ", start);

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
    END_ENUMERATION:
    total_print_time.StartTimer();
    sumAllMatchFind+=allMatchFind;
    std::cout << "num add top k:" << numAddTopk << endl;
    std::cout << "all match find:" << allMatchFind << endl;
    updateTopK();
    total_print_time.StopTimer();
    //Print_Time2("PrintTopk ", start);
}

//删除边 删除allmatch中tmin=td的记录
void Graphflow::deleteEdge(uint v1, uint v2) {
    std::chrono::high_resolution_clock::time_point start;
    start = Get_Time();
    //首先拿到v1 v2之间的tmin
    uint tmin = data_.GetEdgeTime(v1, v2);
    uint delete_num = 0;
    uint cnt = 0;
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
    num_negative_results_ += delete_num;
    //std::cout<<"num_negative_results_"<<num_negative_results_<<endl;
    data_.RemoveEdge(v1, v2);

    //Print_Time("deleteEdge  ", start);
    if (cnt == 0)
        return;
    //根据cnt的个数需要填补cnt个到Top k 否则就重新截取top k;
    start = Get_Time();
    deleteUpdateTopK();
    //Print_Time("deleteUpdateTopK ", start);
}

//删除问题 需要对比补的是不是在topk中
void Graphflow::deleteUpdateTopK() {
#ifdef RESULT_TRACK
    stringstream _ss;
    _ss << "after delete" << std::endl;
    //std::sort(topKSet.begin(),topKSet.end(), matchRecordCmp);
    for (auto d: topKSet) {
        _ss << d->toString();
        Log::track2(_ss);
        _ss.clear();
        _ss.str("");
    }
#endif
}

//动态地减边操作
void Graphflow::RemoveEdge(uint v1, uint v2, uint label) {
    //删除边
    total_delete_update_time.StartTimer();
    allMatchFind = 0;
    uint v1label = data_.GetVertexLabel(v1);
    uint v2label = data_.GetVertexLabel(v2);
    data_.RemoveEdge(v1, v2);
    data_.UpdateLabelIndex(v1, v2, label, 0);
    uint n = query_.NumEdges();
    vector<int> match = EdgeisInMatchOrder(v1, v2, v1label, v2label, label);
    deleteUpdateStarIndex(v1, v2, match);
    total_delete_update_time.StopTimer();
    total_delete_time.StartTimer();
    bool flag = deleteMatchRecordWithEdge(v1, v1label, v2, v2label, label, match);
    if (!flag) {
        return;
    }
#ifdef GLOBAL
    std::vector<uint> mm(query_.NumVertices(), UNMATCHED);
    size_t num_results = 0ul;
#endif
    for (int i = 0; i < data_.NumVertices(); i++) {
        const std::vector<uint> &neighbors = data_.GetNeighbors(i);
        auto lower = std::lower_bound(neighbors.begin(), neighbors.end(), i);
        while (lower != neighbors.end()) {
            uint v1 = i;
            uint v2 = (*lower);
            uint v1label = data_.GetVertexLabel(v1);
            uint v2label = data_.GetVertexLabel(v2);
            uint label = std::get<2>(data_.GetEdgeLabel(v1, v2));
            vector<int> match = EdgeisInMatchOrder(v1, v2, v1label, v2label, label);
            for (auto m: match) {
                uint u1 = order_vs_[m][0];
                uint u2 = order_vs_[m][1];
                uint u1label = this->query_.GetVertexLabel(u1);
                uint u2label = this->query_.GetVertexLabel(u2);
                uint v1label = this->data_.GetVertexLabel(v1);
                uint v2label = this->data_.GetVertexLabel(v2);
                float weight = this->data_.GetEdgeWeight(v1, v2);
#ifdef LOCAL
                InitialLocalIndex(m);
                if (v1label != v2label) {
                    if (v1label != u1label) {
                        swap(v1, v2);
                    }
                    //todo
                    if (this->LabelFilter(v1, u1) && this->LabelFilter(v2, u2)) {
                        bool flag=SearchMatchesWithEdge(m,v1,v2,weight,u1,u2,negative);
                        if(flag)
                            continue;
                    }
                } else {
                    for (int i = 0; i < 2; i++) {
                        if (this->LabelFilter(v1, u1) && this->LabelFilter(v2, u2)) {
                            //todo
                            //this->matchCandidate.clear();
                            bool flag=SearchMatchesWithEdge(m,v1,v2,weight,u1,u2,negative);
                            if(flag)
                            {
                                std::swap(v1, v2);
                                continue;
                            }
                        }
                        std::swap(v1, v2);
                    }
                }
#endif
#ifdef GLOBAL
                for (int j = query_.NumVertices() - 1; j > 0; j--) {
                    if (j == query_.NumVertices() - 1) {
                        suffixMax[j] = globalStarIndex[m][j]->getStarMaxWeight();
                    } else {
                        suffixMax[j] = suffixMax[j + 1] + globalStarIndex[m][j]->getStarMaxWeight();
                    }
                }
                if (v1label != v2label) {
                    if (v1label != u1label) {
                        swap(v1, v2);
                    }
                    //todo
                    if (this->LabelFilter(v1, u1) && this->LabelFilter(v2, u2)) {
                        //compute suffix array
                        mm[u1] = v1;
                        mm[u2] = v2;
                        visited_[v1] = true;
                        visited_[v2] = true;
                        if(!isInsert)
                            DsearchSpace+=2;
                        //total_test.StartTimer();
                        FindMatches(1, m, 2, mm, num_results, weight, negative);
                        //total_test.StopTimer();
                        visited_[v1] = false;
                        visited_[v2] = false;
                        mm[u1] = UNMATCHED;
                        mm[u2] = UNMATCHED;
                    }
                } else {
                    for (int i = 0; i < 2; i++) {
                        if (this->LabelFilter(v1, u1) && this->LabelFilter(v2, u2)) {
                            mm[u1] = v1;
                            mm[u2] = v2;
                            visited_[v1] = true;
                            visited_[v2] = true;
                            if(!isInsert)
                                DsearchSpace+=2;
                            //total_test.StartTimer();
                            FindMatches(1, m, 2, mm, num_results, weight, negative);
                            //total_test.StopTimer();
                            visited_[v1] = false;
                            visited_[v2] = false;
                            mm[u1] = UNMATCHED;
                            mm[u2] = UNMATCHED;
                        }

                        std::swap(v1, v2);
                    }
                }
#endif
            }
            lower++;
        }
    }
    total_delete_time.StopTimer();
    sumDeleteallMatchFind += allMatchFind;
    std::cout << "delete research matches:" << allMatchFind << endl;
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
    uint dataVlabelNum = this->data_.NumVLabels();
    uint dataElabelNum = this->data_.NumELabels();
    uint queryVlabelNum = this->query_.NumVLabels();
    uint queryElabelNum = this->query_.NumELabels();
    const auto &dataLabelIndex = this->data_.labelIndex[data_v];
    const auto &queryLabelIndex = this->query_.labelIndex[query_v];

    for (int i = 0; i < queryVlabelNum; i++) {
        if (i < dataVlabelNum) {
            if (dataLabelIndex[i] < queryLabelIndex[i]) {
                return false;
            }
        } else {
            if (queryLabelIndex[i] > 0) {
                return false;
            }
        }
    }
    for (int i = 0; i < queryElabelNum; i++) {
        if (i < dataElabelNum) {
            if (dataLabelIndex[dataVlabelNum + i] < queryLabelIndex[queryVlabelNum + i]) {
                return false;
            }
        } else {
            if (queryLabelIndex[queryVlabelNum + i] > 0)
                return false;
        }
    }
    return true;
}


void Graphflow::matchVertex(bool isFirstEdge, uint depth, uint data_v, float w) {
    if (isFirstEdge) {
        this->match[depth].setVertexId(data_v);
        this->match[depth].setSumWeight(w);
        this->matchCandidate[depth].emplace_back(SingleCandidate(data_v, w));
        //this->matchVertexCandidate.push_back({});

        this->visited_[data_v] = true;
    } else {
        float weight = match[depth-1].getSumWeight() + w;
        this->match[depth].setVertexId(data_v);
        this->match[depth].setSumWeight(weight);
        //this->matchCandidate.emplace_back(singleVertexCandidate);
        //this->matchVertexCandidate.push_back({});

        this->visited_[data_v] = true;
    }

}

void Graphflow::matchVertex(int depth) {
    this->match[depth].setIsolateSingleCandate();
    //this->matchCandidate.emplace_back(singleVertexCandidate);
    // this->matchVertexCandidate.push_back(vertexCandidate);
}

void Graphflow::popVertex(uint depth, uint data_v) {
    this->match[depth].clearSingleCandidate();
    //this->matchVertexCandidate.pop_back();
    this->matchCandidate[depth].clear();

    this->visited_[data_v] = false;
}

void Graphflow::popVertex(uint data_v, uint matchorderindex, uint depth, const std::vector<Neighbor> &uk_neighbor) {
    this->match[depth].clearSingleCandidate();
    const int n = uk_neighbor.size();
    for (auto u: uk_neighbor) {
        int query_order_index = order_vertex_index[matchorderindex][u.getVertexId()];
        matchCandidate[query_order_index].resize(0);
    }
    //this->matchVertexCandidate.pop_back();
    this->matchCandidate[depth].resize(0);

    this->visited_[data_v] = false;
}

void Graphflow::popVertex(uint depth) {
    this->match[depth].clearSingleCandidate();
    // this->matchCandidate[i].clear();
}

void Graphflow::densityFilter(uint matchorder_index, uint depth, std::vector<SingleCandidate> &singleVertexCandidate) {
    uint n = query_.NumVertices();
    if (topKSet.size() < k) {
        auto iter = singleVertexCandidate.begin();
        while (iter != singleVertexCandidate.end()) {
            float md = (*iter).getSumWeight();
            iter++;
        }
        return;
    }
    float kw = topKSet.back()->getDensity();
    float sumWeight = 0;
    sumWeight += this->match[depth - 1].getSumWeight();

    float backWeight = GetBackWeight(matchorder_index, depth + 1);
    sumWeight += backWeight;
    int cnt = 0;

    auto iter = singleVertexCandidate.begin();
    while (iter != singleVertexCandidate.end()) {
        float md = (*iter).getSumWeight();
        float tmpweight = (sumWeight + md) / (sqrt(n) * (n - 1));
        if (tmpweight < kw) {
            (*iter).setVertexId(-1);
            cnt++;
        }
        iter++;
    }
    if (cnt == 0)
        return;
    int newLen = singleVertexCandidate.size() - cnt;
    int svcLen = singleVertexCandidate.size();
    if (newLen == 0) {
        singleVertexCandidate.resize(0);
        return;
    } else {
        int i = 0;
        int j = 0;
        while (j < svcLen) {
            if (singleVertexCandidate[j].getVertexId() != -1) {
                singleVertexCandidate[i] = singleVertexCandidate[j];
                i++;
            }
            j++;

        }
        singleVertexCandidate.resize(newLen);
    }
}

void Graphflow::combination_helper(std::vector<std::vector<int>> &result, std::vector<int> &current,
                                   const std::vector<std::vector<uint>> &nums, int k) {
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

std::vector<std::vector<int>> Graphflow::combination(const std::vector<std::vector<uint>> &nums) {
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
    if (this->visited_[vertex] == true) {
        return false;
    }
    std::vector<uint> result;
    const auto &vNeighbors = this->data_.GetNeighbors(vertex);
    if (!needToIntersection.empty()) {
        std::set_intersection(vNeighbors.begin(), vNeighbors.end(), needToIntersection.begin(),
                              needToIntersection.end(),
                              std::insert_iterator<std::vector<uint>>(result, result.begin()));
        if (result.empty()) {
            return false;
        }
    } else {
        for (int i = 0; i < vNeighbors.size(); i++) {
            if (this->visited_[vNeighbors[i]] == false &&
                this->data_.GetVertexLabel(vNeighbors[i]) == queryVertexLabel) {
                result.push_back(vNeighbors[i]);
            }
        }
        if (result.empty())
            return false;
    }
    if (intersectresult.size() == 0) {
        intersectresult = result;
        return true;
    }
    std::vector<uint> intersectionCopy;
    intersectionCopy.reserve(intersectresult.size());
    std::swap(intersectionCopy, intersectresult);
    std::set_intersection(intersectionCopy.begin(), intersectionCopy.end(), result.begin(), result.end(),
                          std::insert_iterator<std::vector<uint>>(intersectresult, intersectresult.begin()));
    if (intersectresult.empty()) {
        return false;
    }
    return true;
}

void Graphflow::setSingleVertexByIntersectionResult(std::vector<tuple<int, int, float>> &singleVertex,
                                                    std::vector<uint> &intersectresult, std::vector<int> &r) {
    //intersectresult候选解，r是选择的LD顶点的候选解，singleVertex是对应的原来candidate候选边
    //利用intersectresult同步更新singleVertex
    sychronizeSingleVertexAndCandidate(singleVertex, intersectresult);
    //1.增加权值
    for (int i = 0; i < singleVertex.size(); i++) {
        uint v1 = std::get<0>(singleVertex[i]);
        float sumweight = 0;
        int tmin = std::get<1>(singleVertex[i]);
        for (int j = 0; j < r.size(); j++) {
            uint v2 = r[j];
            sumweight += this->data_.GetEdgeWeight(v1, v2);
            tmin = std::min(tmin, (int) this->data_.GetEdgeTime(v1, v2));
        }
        std::get<2>(singleVertex[i]) = std::get<2>(singleVertex[i]) + sumweight;
        std::get<1>(singleVertex[i]) = tmin;
    }

}

void Graphflow::setIsolateVertexMatchResult(std::vector<int> &r, std::vector<int> &isolateVertex, float density) {
    for (int i = 0; i < isolateVertex.size(); i++) {
        uint index = isolateVertex[i];
        this->match[index].setVertexId(r[i]);
    }
    int depth = this->match.size();
    this->match[depth - 1].setSumWeight(density);
}


void Graphflow::setBatchVisited(std::vector<int> &r, bool flag) {
    for (auto item: r) {
        /*  if(item==61070){
              std::cout<<"setBatchVisited 61070: "<<flag<<std::endl;
          }*/
        this->visited_[item] = flag;
    }
}


void Graphflow::recoverIsolateVertexMatchResult(std::vector<int> &IsolateVertexs) {
    for (int i = 0; i < IsolateVertexs.size(); i++)
        this->match[IsolateVertexs[i]].setIsolateSingleCandate();
}

void Graphflow::sychronizeSingleVertexAndCandidate(std::vector<tuple<int, int, float>> &singleVertex,
                                                   std::vector<uint> &intersectresult) {
    std::sort(singleVertex.begin(), singleVertex.end(), tupleVertexIdCmp2);
    if (singleVertex.size() == 0) {
        //todo
        for (auto item: intersectresult) {
            singleVertex.push_back(make_tuple(item, INT_MAX, 0));
        }
        return;
    }
    auto iter1 = singleVertex.begin();
    auto iter2 = intersectresult.begin();
    while (iter1 != singleVertex.end()) {
        if (std::get<0>((*iter1)) == (*iter2)) {
            iter1++;
            iter2++;
        } else {
            iter1 = singleVertex.erase(iter1);
        }
    }

}

void Graphflow::addMatchResult(uint matchorderindex, searchType type) {
    std::chrono::high_resolution_clock::time_point starttime;
    int n = query_.NumVertices();
    auto &isolateVertexs = query_.isolatedRecord[matchorderindex];
    std::vector<std::vector<SingleCandidate>> combinezIsolateVertexs;
    combinezIsolateVertexs.reserve(isolateVertexs.size());
    for (int i = 0; i < isolateVertexs.size(); i++) {
        std::vector<SingleCandidate> &singleVertex = matchCandidate[isolateVertexs[i]];
        std::sort(singleVertex.begin(), singleVertex.end());
        combinezIsolateVertexs.emplace_back(singleVertex);
    }

    int len = combinezIsolateVertexs.size();
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
    auto wt = findWeightBeforeIsolated();
    int *pmax = new int[len]{};//最大pmax在组中的列号
    std::unordered_map<int, int> visitedId;
    float Tmax = wt;
    bool TmaxisAllFirst = true;
    for (int i = 0; i < len; i++) {
        int index = 0;
        int id = combinezIsolateVertexs[i][index].getVertexId();
        bool flag = true;
        while (visited_[id]) {
            combinezIsolateVertexs[i].erase(combinezIsolateVertexs[i].begin() + index);
            if (index == combinezIsolateVertexs[i].size()) {
                flag = false;
                break;
            }
            id = combinezIsolateVertexs[i][index].getVertexId();
        }
        if (!flag) {
            return;
        }
        if (!visitedId.count(id)) {
            visitedId[id] = i;
        } else {
            TmaxisAllFirst = false;
            int pre_group_id = visitedId[id];
            float gap1 = combinezIsolateVertexs[pre_group_id][pmax[pre_group_id]].getSumWeight() -
                         combinezIsolateVertexs[pre_group_id][pmax[pre_group_id] + 1].getSumWeight();
            float gap2 = combinezIsolateVertexs[i][pmax[i]].getSumWeight() -
                         combinezIsolateVertexs[i][pmax[i] + 1].getSumWeight();
            if (gap1 < gap2) {
                pmax[pre_group_id]++;
                visitedId[id] = i;
                int preId = combinezIsolateVertexs[pre_group_id][pmax[pre_group_id]].getVertexId();
                visitedId[preId] = pre_group_id;
            } else {
                pmax[i]++;
            }
            int curId = combinezIsolateVertexs[i][pmax[i]].getVertexId();
            visitedId[curId] = i;
        }
    }
    float globalTmax = Tmax;
    std::vector<int> r(isolateVertexs.size());
    for (int i = 0; i < len; i++) {
        r[i] = combinezIsolateVertexs[i][pmax[i]].getVertexId();
        Tmax += combinezIsolateVertexs[i][pmax[i]].getSumWeight();
    }
    if (TmaxisAllFirst) {
        //1.
        globalTmax = Tmax;
        float density = Tmax / (sqrt(n) * (n - 1));
        setIsolateVertexMatchResult(r, isolateVertexs, Tmax);
        std::vector<uint> m(n);
        for (int i = 0; i < match.size(); i++) {
            m[order_vs_[matchorderindex][i]] = match[i].getVertexId();
        }
        MatchRecord *record = new MatchRecord(density, m);
        bool matchResult = addMatchRecords(record);
        if (matchResult) {
            if (type == positive)
                num_positive_results_++;
            else
                num_negative_results_++;
        }
        recoverIsolateVertexMatchResult(isolateVertexs);
        if (topKSet.size() == k) {
            float cur_density = topKSet.back()->getDensity();
            float density = Tmax / (sqrt(n) * (n - 1));
            if (cur_density > density) {
                delete[]pmax;
                return;
            }
        }
    } else {
        for (int i = 0; i < len; i++) {
            globalTmax += combinezIsolateVertexs[i][0].getSumWeight();
        }
    }

    //2.选择下一个迭代的点
    int *hash = new int[len]{};
    float *Tbound = new float[len];
    std::fill(Tbound, Tbound + len, Tmax);

    int next_vertex_group = 0;
    bool isover = false;
    int *noscan = new int[len]{};
    while (!isover) {
        int tmpcnt = 0;
        next_vertex_group = findTboundMaxIndex(Tbound, hash, noscan, combinezIsolateVertexs, len);

        if (isnoNextVertex(noscan, len)) {
            isover = true;
            break;
        }

        hash[next_vertex_group]++;
        float pre = combinezIsolateVertexs[next_vertex_group][pmax[next_vertex_group]].getSumWeight();
        float cur = combinezIsolateVertexs[next_vertex_group][hash[next_vertex_group]].getSumWeight();
        Tbound[next_vertex_group] = globalTmax - pre + cur;
        float TboundNext = Tbound[next_vertex_group];
        if (TboundNext > Tmax)
            Tmax = Tbound[next_vertex_group];
        auto next_item = combinezIsolateVertexs[next_vertex_group][hash[next_vertex_group]];

        //3. 递归
        int id = next_item.getVertexId();
        CatesianProductWithIndex(matchorderindex, type, next_vertex_group, 0, len, hash, combinezIsolateVertexs,
                                 isolateVertexs, wt);
        if (topKSet.size() == k) {
            float cur_density = topKSet.back()->getDensity();
            float density = Tmax / (sqrt(n) * (n - 1));
            if (cur_density > density) {
                isover = true;
                break;
            }
            float TboundDensity = TboundNext / (sqrt(n) * (n - 1));
            if (cur_density > TboundDensity) {
                noscan[next_vertex_group] = 1;
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

std::vector<tuple<std::vector<int>, int, float>>
Graphflow::combinationMatchResult(std::vector<std::vector<tuple<int, int, float>>> combinezIsolateVertexs) {
    std::chrono::high_resolution_clock::time_point start;
    // start=Get_Time();
    std::vector<tuple<std::vector<int>, int, float>> result;
    std::vector<int> current;
    combinationMatchResultHelp(result, current, combinezIsolateVertexs, 0, INT_MAX, 0);
    //  Print_Time2("combinationMatchResult ",start);
    return result;
}

void Graphflow::combinationMatchResultHelp(std::vector<tuple<std::vector<int>, int, float>> &result,
                                           std::vector<int> &current,
                                           std::vector<std::vector<tuple<int, int, float>>> &combinezIsolateVertexs,
                                           int k, int tmin, float density) {
    if (k == combinezIsolateVertexs.size()) {
        result.push_back(make_tuple(current, tmin, density));
        return;
    }
    for (int i = 0; i < combinezIsolateVertexs[k].size(); i++) {
        float copydensity = density;
        tmin = std::min(std::get<1>(combinezIsolateVertexs[k][i]), tmin);
        density += std::get<2>(combinezIsolateVertexs[k][i]);
        current.push_back(std::get<0>(combinezIsolateVertexs[k][i]));
        combinationMatchResultHelp(result, current, combinezIsolateVertexs, k + 1, tmin, density);
        current.pop_back();
        density = copydensity;
    }
}

float Graphflow::findWeightBeforeIsolated() {
    int depth = match.size() - 1;
    for (int i = depth; i >= 0; i--) {
        if (match[i].getVertexId() == -1) {
            continue;
        } else {
            float weight = match[i].getSumWeight();
            return weight;
        }
    }
}

void
Graphflow::CatesianProductWithIndex(int matchorderindex, searchType type, int curIndex, int depth, int len, int *hash,
                                    std::vector<std::vector<SingleCandidate>> &combinezIsolateVertexs,
                                    std::vector<int> &isolateVertexs, float &weight) {
    if (depth == len) {
        int n = query_.NumVertices();
        std::vector<uint> m(n);
        for (int i = 0; i < match.size(); i++) {
            m[order_vs_[matchorderindex][i]] = match[i].getVertexId();
        }
        float density = weight / (sqrt(n) * (n - 1));
        MatchRecord *record = new MatchRecord(density, m);
        bool matchResult = addMatchRecords(record);
        if (matchResult) {
            if (type == positive)
                num_positive_results_++;
            else
                num_negative_results_++;
        }
        return;
    }
    if (depth != curIndex) {
        int right = hash[depth];
        for (int i = 0; i <= right; i++) {
            int id = combinezIsolateVertexs[depth][i].getVertexId();
            if (visited_[id])
                continue;
            visited_[id] = true;
            float copyweight = weight;
            this->match[isolateVertexs[depth]].setVertexId(id);
            weight += combinezIsolateVertexs[depth][i].getSumWeight();
            CatesianProductWithIndex(matchorderindex, type, curIndex, depth + 1, len, hash, combinezIsolateVertexs,
                                     isolateVertexs, weight);
            visited_[id] = false;
            this->match[isolateVertexs[depth]].setVertexId(-1);
            weight = copyweight;
        }
    } else {
        //  std::cout<<"hash["<<curIndex<<"]"<<hash[curIndex]<<" depth["<<depth<<"]"<<combinezIsolateVertexs[depth][hash[curIndex]].getVertexId()<<endl;
        const int &id = combinezIsolateVertexs[depth][hash[curIndex]].getVertexId();
        if (visited_[id]) {
            return;
        }
        visited_[id] = true;
        float copyweight = weight;
        this->match[isolateVertexs[depth]].setVertexId(id);
        weight += combinezIsolateVertexs[depth][hash[curIndex]].getSumWeight();
        CatesianProductWithIndex(matchorderindex, type, curIndex, depth + 1, len, hash, combinezIsolateVertexs,
                                 isolateVertexs, weight);
        visited_[id] = false;
        this->match[isolateVertexs[depth]].setVertexId(-1);
        weight = copyweight;
    }

}

int Graphflow::findTboundMaxIndex(float *Tbound, int *hash, int *noscan,
                                  std::vector<std::vector<SingleCandidate>> &combinezIsolateVertexs, int len) {
    float max = 0;
    int maxindex = 0;
    for (int i = 0; i < len; i++) {
        int index = hash[i] + 1;
        if (noscan[i] == 1) {
            continue;
        }
        if (index >= combinezIsolateVertexs[i].size()) {
            noscan[i] = 1;
            continue;
        }
        int id = combinezIsolateVertexs[i][index].getVertexId();
        bool flag = true;
        while (visited_[id]) {
            combinezIsolateVertexs[i].erase(combinezIsolateVertexs[i].begin() + index);
            if (index >= combinezIsolateVertexs[i].size()) {
                noscan[i] = 1;
                flag = false;
                break;
            } else {
                id = combinezIsolateVertexs[i][index].getVertexId();
            }
        }
        if (!flag)
            continue;
        if (Tbound[i] > max) {
            max = Tbound[i];
            maxindex = i;
        }
    }
    return maxindex;
}

bool Graphflow::isnoNextVertex(int *noscan, int len) {
    for (int i = 0; i < len; i++) {
        if (noscan[i] == 0)
            return false;
    }
    return true;
}

void Graphflow::addMatchResultWithHeap(uint matchorderindex, searchType type) {
    std::chrono::high_resolution_clock::time_point starttime;
    int n = query_.NumVertices();
    auto &isolateVertexs = query_.isolatedRecord[matchorderindex];
    std::vector<std::vector<SingleCandidate>> combinezIsolateVertexs;
    combinezIsolateVertexs.reserve(isolateVertexs.size());
    for (int i = 0; i < isolateVertexs.size(); i++) {
        std::vector<SingleCandidate> &singleVertex = matchCandidate[isolateVertexs[i]];
        std::sort(singleVertex.begin(), singleVertex.end());
        combinezIsolateVertexs.emplace_back(singleVertex);
    }

    int len = combinezIsolateVertexs.size();
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
    auto wt = findWeightBeforeIsolated();
    int *hash = new int[len]{};//索引指针
    int *noscan = new int[len]{};//记录是否停止
    std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int>>, pairCompare> maxHeap;
    for (int i = 0; i < combinezIsolateVertexs.size(); i++) {
        float weight = combinezIsolateVertexs[i][0].getSumWeight();
        maxHeap.push(std::make_pair(weight, i));
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
    float globalTmax = wt;
    for (int i = 0; i < len; i++)
        globalTmax += combinezIsolateVertexs[i][0].getSumWeight();
    while (!maxHeap.empty() && !isnoNextVertex(noscan, len)) {
        auto item = maxHeap.top();
        maxHeap.pop();
        int group = item.second;
        int pre = hash[group];
        hash[group]++;
        if (pre == combinezIsolateVertexs[group].size() - 1) {
            noscan[group] = 1;
        } else {
            float w = combinezIsolateVertexs[group][hash[group]].getSumWeight();
            maxHeap.push(std::make_pair(w, group));
        }
        int id = combinezIsolateVertexs[group][pre].getVertexId();
        if (visited_[id]) {
            combinezIsolateVertexs[group].erase(combinezIsolateVertexs[group].begin() + pre);
            hash[group]--;
            continue;
        } else {
            //笛卡尔积 从最小的seen number 开始
            int minGroup = 0;
            int minNum = INT_MAX;
            bool flag = true;
            for (int i = 0; i < len; i++) {
                if (i == group)
                    continue;
                if (hash[i] == 0) {
                    flag = false;
                    break;
                } else if (hash[i] < minNum) {
                    minGroup = i;
                    minNum = hash[group];
                }
            }

            if (!flag)
                continue;
            int depthLen = len - 1;
            std::vector<int> copyisolateIndex;
            copyisolateIndex.emplace_back(minGroup);
            for (int i = 0; i < len; i++) {
                if (i == minGroup || i == group)
                    continue;
                copyisolateIndex.emplace_back(i);
            }
            visited_[id] = true;
            this->match[isolateVertexs[group]].setVertexId(id);
            float weight = wt + combinezIsolateVertexs[group][pre].getSumWeight();
            float Tmax = 0;
            //递归
            CatesianProductWithHeap(matchorderindex, type, 0, depthLen, hash, combinezIsolateVertexs, isolateVertexs,
                                    copyisolateIndex, weight);
            visited_[id] = false;
            this->match[isolateVertexs[group]].setVertexId(-1);
            if (topKSet.size() == k) {
                float curDensity = topKSet.back()->getDensity();
                Tmax = globalTmax - combinezIsolateVertexs[group][0].getSumWeight() +
                       combinezIsolateVertexs[group][pre].getSumWeight();
                float Tmaxdensity = Tmax / (sqrt(n) * (n - 1));
                if (curDensity > Tmaxdensity) {
                    noscan[group] = 1;
                }
            }
        }
    }
    delete[]hash;
    delete[]noscan;
}

void Graphflow::CatesianProductWithHeap(int matchorderindex, searchType type, int depth, int len, int *hash,
                                        std::vector<std::vector<SingleCandidate>> &combinezIsolateVertexs,
                                        std::vector<int> &isolateVertexs, std::vector<int> &isolatedIndex,
                                        float &weight) {
    if (depth == len) {
        int n = query_.NumVertices();
        std::vector<uint> m(n);
        for (int i = 0; i < match.size(); i++) {
            m[order_vs_[matchorderindex][i]] = match[i].getVertexId();
        }
        float density = weight / (sqrt(n) * (n - 1));

        MatchRecord *record = new MatchRecord(density, m);
        bool matchResult = addMatchRecords(record);
        if (matchResult) {
            if (type == positive)
                num_positive_results_++;
            else
                num_negative_results_++;
        }
        return;
    }
    int group = isolatedIndex[depth];
    int right = hash[group];
    for (int i = 0; i < right; i++) {
        int id = combinezIsolateVertexs[group][i].getVertexId();
        if (visited_[id])
            continue;
        visited_[id] = true;
        float copyweight = weight;
        this->match[isolateVertexs[group]].setVertexId(id);
        weight += combinezIsolateVertexs[group][i].getSumWeight();
        CatesianProductWithHeap(matchorderindex, type, depth + 1, len, hash, combinezIsolateVertexs, isolateVertexs,
                                isolatedIndex, weight);
        visited_[id] = false;
        this->match[isolateVertexs[group]].setVertexId(-1);
        weight = copyweight;
    }
}


void Graphflow::createLabelToQueryVertex() {
    for (int i = 0; i < query_.NumVertices(); i++) {
        uint label = query_.GetVertexLabel(i);
        labelToQueryVertex[label].emplace_back(i);
        queryVertexIndexInlabel[i] = labelToQueryVertex[label].size() - 1;
    }
}

bool Graphflow::updaterightNeighborCandidate(int matchorderindex, uint uk, uint uk_neigh, bool isFirstEdge, uint vk,
                                             const std::vector<Neighbor> &uk_neighbor) {
    const std::vector<Neighbor> &vN = this->data_.vNeighbors2[vk];
    const int n = uk_neighbor.size();
    std::vector<Neighbor>::const_iterator iter1,iter2,start,lower;
    iter1=iter2=start=vN.begin();
    //1.对于所有的右邻居，找其候选解
    for (int i = 0; i < n; i++) {
        uint query_id = uk_neighbor[i].getVertexId();
        if (isFirstEdge) {
            if (query_id == uk_neigh) {
                isFirstEdge = false;
                continue;
            }
        }
        uint query_vertex_label = query_.GetVertexLabel(query_id);
        int query_order_index = order_vertex_index[matchorderindex][query_id];
        uint query_elabel = std::get<2>(query_.GetEdgeLabel(uk, query_id));
        StarGraph *s = globalStarIndex[matchorderindex][query_order_index];
        bool isFirstVertex = true;
        bool isCandidateFirstNull = true;
        float maxweight = 0;
        float curWeight = 0;
        Neighbor neighbor(query_vertex_label, query_elabel);
        if(i>0){
            lower=uk_neighbor[i].GetelabelAndVertexLabel()!=uk_neighbor[i-1].GetelabelAndVertexLabel()?iter2:iter1;
        }
        else{
            lower = std::lower_bound(start, vN.end(), neighbor, CompareNeighbors2);
        }

        iter1=lower;
        std::pair<uint, uint> evl = make_pair(query_elabel, query_vertex_label);
        int flag=1;
       // auto lower = std::lower_bound(vN.begin(), vN.end(), neighbor, CompareNeighbors2);
        //auto upper = std::upper_bound(vN.begin(), vN.end(), neighbor, CompareNeighbors2);
        if(matchCandidate[query_order_index].size()==0) {
            //auto lower = std::lower_bound(vN.begin(), vN.end(), neighbor, CompareNeighbors2);
            if (lower==vN.end()) {
                return true;
            } else {
                while(lower<vN.end()){
                    if((*lower).GetelabelAndVertexLabel()<evl){
                        break;
                    }
                    while((*lower).GetelabelAndVertexLabel()>evl){
                        lower++;
                        if(lower==vN.end()){
                            flag=0;
                            break;
                        }
                    }
                    if(!flag){
                        break;
                    }
                    if((*lower).GetelabelAndVertexLabel()==evl){
                        uint neighbor_id = (*lower).getVertexId();
                        if (visited_[neighbor_id]) {
                            lower++;
                            continue;
                        }
                        maxweight = globalVkMatchUk[neighbor_id][matchorderindex][queryVertexIndexInlabel[query_id]];
                        //update LocalStarIndex
                        //add candidate
                        curWeight = (*lower).GetEdgeWeight();
                        if (isFirstVertex) {
                            isFirstVertex = false;
                            LocalStarIndex[query_order_index] = maxweight;
                        } else {
                            if (LocalStarIndex[query_order_index] < maxweight) {
                                LocalStarIndex[query_order_index] = maxweight;
                            }
                        }
                        matchCandidate[query_order_index].emplace_back(neighbor_id, curWeight);
                        lower++;
                    }
                }
                iter2=lower;
            }
        }
        else{
            int qSize=matchCandidate[query_order_index].size();
            int rightV=0;
            int flag=1;
            while (lower<vN.end()&&rightV<qSize) {
                uint neighbor_id = (*lower).getVertexId();
                //uint curmatchCandidateId=matchCandidate[query_order_index][rightV].getVertexId();
                if (visited_[neighbor_id]) {
                    lower++;
                    continue;
                }
                if((*lower).GetelabelAndVertexLabel()<evl){
                    break;
                }
                while((*lower).GetelabelAndVertexLabel()>evl){
                    lower++;
                    if(lower==vN.end()){
                        flag=0;
                        break;
                    }
                }
                if(!flag)
                    break;
                if((*lower).GetelabelAndVertexLabel()==evl)
                {  //双指针取交集
                    while(lower<vN.end()&&rightV<qSize){
                        if(lower->GetelabelAndVertexLabel()!=evl){
                            flag=0;
                            break;
                        }
                        else if((*lower).getVertexId()<matchCandidate[query_order_index][rightV].getVertexId()){
                            rightV++;
                        }
                        else if((*lower).getVertexId()>matchCandidate[query_order_index][rightV].getVertexId()){
                            lower++;
                            if(lower==vN.end()){
                                flag=0;
                                break;
                            }
                        }
                        else{
                            matchCandidate[query_order_index][rightV].setFlag(1);
                            curWeight=(*lower).GetEdgeWeight();
                            neighbor_id=(*lower).getVertexId();
                            maxweight = globalVkMatchUk[neighbor_id][matchorderindex][queryVertexIndexInlabel[query_id]];
                            if(isFirstVertex){
                                isFirstVertex= false;
                                LocalStarIndex[query_order_index]=maxweight;
                            }
                            else{
                                if(LocalStarIndex[query_order_index]<maxweight){
                                    LocalStarIndex[query_order_index]=maxweight;
                                }
                            }
                            matchCandidate[query_order_index][rightV].addSumWeight(curWeight);
                            lower++;
                            rightV++;
                        }
                    }
                    if(!flag)
                        break;
                }
            }
            iter2=lower;
            int k = 0;
            int j = 0;
            int len=0;
            int csize = matchCandidate[query_order_index].size();
            while (j < csize) {
                if (matchCandidate[query_order_index][j].getFlag() == 1) {
                    matchCandidate[query_order_index][j].setFlag(0);
                    matchCandidate[query_order_index][k] = matchCandidate[query_order_index][j];
                    k++;
                    len++;
                }
                j++;
            }
            matchCandidate[query_order_index].resize(len);
            if(isInsert)
                IdeterminCandite++;
            else
                DdeterminCandite++;
        }
        if(matchCandidate[query_order_index].size()==0){
            for (int j = 0; j < i; j++) {
                uint query_id = uk_neighbor[j].getVertexId();
                int query_order_index = order_vertex_index[matchorderindex][query_id];
                matchCandidate[query_order_index].resize(0);
            }
            return true;
        }
    }
    return false;
}

void Graphflow::InitialLocalIndex(int matchorderindex) {
    const std::vector<StarGraph *> &gs = globalStarIndex[matchorderindex];
    int n = gs.size();
    for (int i = 1; i < n; i++) {
        LocalStarIndex[i] = gs[i]->getStarMaxWeight();
    }
}

void Graphflow::getIntersetSingleCandidate(std::vector<SingleCandidate> &singleVertexCandidate, int matchorderindex,
                                           int depth) {
    int i = 0;
    int j = 0;
    int csize = singleVertexCandidate.size();
    int len = 0;
    int flag = matchLeftNeighborSum[matchorderindex][depth];
    if (flag == 1)
        return;
    while (j < csize) {
        if (singleVertexCandidate[j].getFlag() == flag) {
            singleVertexCandidate[i] = singleVertexCandidate[j];
            i++;
            len++;
        }
        j++;

    }
    singleVertexCandidate.resize(len);
}

void Graphflow::PrintAverageTime(int len) {
    std::cout <<"average query graph degree:"<< std::fixed << std::setprecision(2)<<query_.NumEdges()*2.0/query_.NumVertices()<<endl;
    std::cout<<"average data graph degree:"<<std::fixed << std::setprecision(2)<<data_.NumEdges()*2.0/data_.NumVertices()<<endl;
    std::cout << "average serach time: " << std::fixed << std::setprecision(2)
              << total_search_time.GetTimer() * 1.0 / len << " microseconds" << endl;
    std::cout << "average update global index time: " << std::fixed << std::setprecision(2)
              << total_update_globalIndex_time.GetTimer() * 1.0 / len << " microseconds" << endl;
    std::cout << "average update time " << std::fixed << std::setprecision(2)
              << (total_update_globalIndex_time.GetTimer() * 1.0 / len + total_search_time.GetTimer() * 1.0 / len)
              << " microseconds" << endl;
    std::cout << "average insert density filter time: " << std::fixed << std::setprecision(2)
              << Itotal_densityfilter_time * 1.0 / len << " microseconds" << endl;
    std::cout << "average insert updaterightNeighborCandidate time " << std::fixed << std::setprecision(2)
              <<  (Itotal_updaterightNeighborCandidate_time* 1.0 / len)<< " microseconds" << endl;
    std::cout << "average print time: " << std::fixed << std::setprecision(2) << total_print_time.GetTimer() * 1.0 / len
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
    std::cout<<"test"<<total_test.GetTimer()*1.0/len<<endl;
#ifdef COMPUTE_TIME
    stringstream _ss;
    _ss<< (total_update_globalIndex_time.GetTimer() * 1.0 / len + total_search_time.GetTimer() * 1.0 / len)<<","
       <<total_update_globalIndex_time.GetTimer() * 1.0 / len<<","
       <<total_search_time.GetTimer() * 1.0 / len<<","
       <<(total_delete_time.GetTimer() * 1.0 / len + total_delete_update_time.GetTimer() * 1.0 / len)<<","
       <<(total_delete_update_time.GetTimer() * 1.0) / len<<","
       <<total_delete_time.GetTimer() * 1.0 / len<<","
       <<sumAllMatchFind<<","
       <<sumDeleteallMatchFind<<","
       <<query_.NumEdges()*2.0/query_.NumVertices()<<","
       <<data_.NumEdges()*2.0/data_.NumVertices()<<","
       <<Itotal_densityfilter_time * 1.0 / len<<","
       <<(Itotal_updaterightNeighborCandidate_time* 1.0 / len)<<","
       <<total_densityFilter_time.GetTimer() * 1.0 / len<<","
       <<total_updaterightNeighborCandidate_time.GetTimer() * 1.0 / len<<","
       << IsearchSpace<<","
       <<DsearchSpace<<","
       <<IdeterminCandite<<","
       <<DdeterminCandite<<endl;

    Log::track3(_ss);
#endif
}

void Graphflow::deleteUpdateglobalVertexStarIndex(uint u1, uint v1, uint n) {
#ifdef LOCAL
    int candidate_index= queryVertexIndexInlabel[u1];
#endif
    uint u1label = query_.GetVertexLabel(u1);
    for (int j = 0; j < n; j++) {
        int vertex_index = order_vertex_index[j][u1];
        if (vertex_index == 0)
            continue;
        StarGraph *s = globalStarIndex[j][vertex_index];
        if (s->getMatchDataVertexId() == v1) {
            StarGraph *s = globalStarIndex[j][vertex_index];
            s->setStarMaxWeight(s->GetForwardNeighborNum() * mw);
#ifdef LOCAL
            updateStarIndex(j,v1,u1,candidate_index);
#endif
#ifdef GLOBAL
            updateStarIndex(j, v1, u1);
#endif
            for (int i = 0; i < data_.NumVertices(); i++) {
                if (i == v1)
                    continue;
                if (data_.GetVertexLabel(i) == u1label) {
#ifdef LOCAL
                    if(LabelFilter(i,u1)){
                        if(globalVkMatchUk[i][j][candidate_index]>s->getStarMaxWeight()||s->getStarMaxWeight()==s->GetForwardNeighborNum()*mw){
                            s->setStarMaxWeight(globalVkMatchUk[i][j][candidate_index]);
                            s->setMatchDataVertexId(i);
                        }
                    }
#endif
#ifdef GLOBAL
                    if (LabelFilter(i, u1)) {
                        updateStarIndex(j, i, u1);
                    }
#endif
                }
            }
        }
    }
}

void Graphflow::deleteUpdateStarIndex(uint v1, uint v2, std::vector<int> &match) {
    uint n = query_.NumEdges();
    for (auto it = match.begin(); it != match.end(); it++) {
        uint u1 = order_vs_[*it][0];
        uint u2 = order_vs_[*it][1];
        uint u1label = query_.GetVertexLabel(u1);
        uint u2label = query_.GetVertexLabel(u2);
        uint v1label = data_.GetVertexLabel(v1);
        uint v2label = data_.GetVertexLabel(v2);
        uint elabel = std::get<2>(query_.GetEdgeLabel(u1, u2));
        if (v1label != v2label) {
            deleteUpdateglobalVertexStarIndex(u1, v1, n);
            deleteUpdateglobalVertexStarIndex(u2, v2, n);
        } else {
            for (int i = 0; i < 2; i++) {
                deleteUpdateglobalVertexStarIndex(u1, v1, n);
                deleteUpdateglobalVertexStarIndex(u2, v2, n);
                std::swap(v1, v2);
            }
        }
    }
}

bool Graphflow::deleteMatchRecordWithEdge(uint v1, uint v1label, uint v2, uint v2label, uint label,
                                          std::vector<int> &match) {
    bool flag = false;
    for (auto it = topKSet.begin(); it != topKSet.end();) {
        MatchRecord *record = *it;
        const std::vector<uint> &m = record->getVetex();
        bool iterflag = true;
        for (int mindex: match) {
            uint u1 = order_vs_[mindex][0];
            uint u2 = order_vs_[mindex][1];
            if ((m[u1] == v1 && m[u2] == v2) || (m[u2] == v1 && m[u1] == v2)) {
                delete record;
                record = nullptr;
                it = topKSet.erase(it);
                iterflag = false;
                flag = true;
                num_negative_results_++;
                break;
            }
        }
        if (iterflag) {
            ++it;
        }
    }
    return flag;
}
bool Graphflow::SearchMatchesWithEdge(uint m,uint v1,uint v2,uint weight,uint u1,uint u2,searchType type){
    this->matchVertex(true, 0, v1, float(0));
    this->matchVertex(true, 1, v2, weight);
    if(isInsert)
        IsearchSpace+=2;
    else
        DsearchSpace+=2;
    bool isNull;
    const std::vector<Neighbor> &uk_neighbor1 = rightNeighbor[m][u1];
    total_updaterightNeighborCandidate_time.StartTimer();
    isNull = updaterightNeighborCandidate(m, u1, u2, true, v1, uk_neighbor1);
    total_updaterightNeighborCandidate_time.StopTimer();
    if (isNull) {
        this->popVertex(1, v2);
        this->popVertex(0, v1);
        return true;
    }
    const std::vector<Neighbor> &uk_neighbor2 = rightNeighbor[m][u2];
    total_updaterightNeighborCandidate_time.StartTimer();
    isNull = updaterightNeighborCandidate(m, u2, u1, true, v2, uk_neighbor2);
    total_updaterightNeighborCandidate_time.StopTimer();
    if (isNull) {
        for (auto u: uk_neighbor1) {
            int query_order_index = order_vertex_index[m][u.getVertexId()];
            matchCandidate[query_order_index].resize(0);
        }
        this->popVertex(1, v2);
        this->popVertex(0, v1);
        return true;
    }
    searchMatches(2, m, type);
    //search()递归
    this->popVertex(v2, m, 1, uk_neighbor1);
    this->popVertex(v1, m, 0, uk_neighbor2);
    return false;
#ifdef DEBUG
    stringstream _ss;
    for(int i=0;i<query_.NumVertices();i++){
        if(matchCandidate[i].size()>0){
            _ss<<"matchCandidate["<<i<<"] not null orderindex:"<<m<<endl;
            Log::track1(_ss);
            _ss.clear();
            _ss.str("");
        }
    }
#endif
}