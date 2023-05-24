//
// Created by 高楚楚 on 2023/5/16.
//

#include <set>
#include <cfloat>
#include "Instopk.h"
bool edgesCmp ( const Edge&e1,const Edge&e2) {
   return e1.GeteWeight()>e2.GeteWeight();
}


Instopk::Instopk(Graph &query_graph, Graph &data_graph, uint max_num_results, bool print_prep, bool print_enum,
                 bool homo,uint d): matching(query_graph, data_graph, max_num_results,
               print_prep, print_enum, homo)
                , order_vs_(query_.NumEdges())
                , order_csrs_(query_.NumEdges())
                , order_offs_(query_.NumEdges())
                ,order_vertex_index(query_.NumEdges())
                ,topKSet(0)
                ,allMatchRecords(0)
                ,dist(d)
               {
                   for (uint i = 0; i < query_.NumEdges(); ++i)
                   {
                       order_vs_[i].resize(query_.NumVertices());//节点个数
                       order_csrs_[i].resize(query_.NumEdges() + 1);//边的个数+1，
                       order_offs_[i].resize(query_.NumVertices(), 0);//节点个数，初始化为0
                       order_vertex_index[i].resize(query_.NumVertices());
                   }
               }
void Instopk::AddVertex(uint id, uint label) {
    data_.AddVertex(id, label);
    visited_.resize(id + 1, false);
}
void Instopk::AddEdge(uint v1, uint v2, uint label, float weight, uint timestamp) {
    data_.AddEdge(v1, v2, label, weight, timestamp, 1);
}
void Instopk::GetMemoryCost(size_t &num_edges, size_t &num_vertices) {
    num_edges = 0ul;
    num_vertices = 0ul;
}
void Instopk::InitialTopK(const std::string &path) {

}
void Instopk::deleteUpdateTopK() {

}
void Instopk::deleteEdge(uint v1, uint v2) {

}
void Instopk::updateTopK(uint num) {

}
void Instopk::InitialMatching(const std::string &path) {

}
void Instopk::RemoveEdge(uint v1, uint v2) {

}
void Instopk::RemoveVertex(uint id) {

}
void Instopk::Preprocessing() {
    //create sortedlist
    CreateSortEdgeList();
    //todo generate matchorder
    GenerateMatchingOrder();
    //todo create Topology /NNW index
    CreateTopologyAndMNWIndex();
    //todo create queryTopology;
    CreateQueryTopology();

}
void Instopk::CreateSortEdgeList() {
    //根据查询图找到边的组合，然后数据图中加入组合，排序
    for (auto edge:this->data_.vEdge){
        uint w1=edge.GetV1Label();
        uint w2=edge.GetV2Label();
        if(w1>w2){
            std::swap(w1,w2);
        }
        std::string key=""+std::to_string(w1)+"#"+std::to_string(w2);
        sortEdgeList[key].push_back(edge);
    }
    for(auto & item:sortEdgeList){
        std::sort(item.second.begin(),item.second.end(), edgesCmp);
    }
    //创建每个点到Edgelist的指针
    for(auto it=sortEdgeList.begin();it!=sortEdgeList.end();it++){
        std::string key=it->first;
        std::vector<Edge>list=it->second;
        std::map<int,std::vector<int>>map;
        for(int i=0;i<list.size();i++){
            int v1=list[i].GetV1();
            int v2=list[i].GetV2();
            map[v1].push_back(i);
            map[v2].push_back(i);
        }
        node2EdgeListPointers[key]=map;
    }

}
void Instopk::CreateTopologyAndMNWIndex() {
    //initial MMW TopologyIndex
    int labelNum = data_.NumVLabels();
    int verticeNum = data_.NumVertices();

    //d=1
    for (int i = 0; i < labelNum; i++) {
        std::string str = std::to_string(i);

        dToOrderType[1].push_back(str);
        for (int j = 0; j < verticeNum; j++) {
            std::map<std::string, int> s;
            s[str] = 0;
            TopologyIndex[1][j][str]=0;
            MNW[1][j][str]=0;
        }
    }
    //d=others
    for (int d = 2; d <= dist; d++) {
        std::vector<std::string> strs = dToOrderType[d - 1];
        for(int k=0;k<verticeNum;k++){
            for (auto s: strs) {
                for (int j = 0; j < labelNum; j++) {
                    std::string copys = s+"#";
                    copys += std::to_string(j);
                    dToOrderType[d].push_back(copys);
                    std::map<std::string, int> smap;
                    smap[copys] = 0;
                    TopologyIndex[d][k][copys]=0;
                    MNW[d][k][copys]=0;
                }
            }
        }

    }
    //bfs createindex
    for (int i = 0; i < verticeNum; i++) {
        std::queue<std::tuple<int, std::string, int>> que;//id,type,weight
        std::queue<int>depth;//id depth
        depth.push((0));
        std::queue<std::set<int>> vertexs;
        std::set<int> verset;
        verset.insert(i);
        vertexs.push(verset);
        que.push(std::make_tuple(i, "", 0));
        while (!que.empty()) {
            auto item = que.front();
            que.pop();
            auto veritem = vertexs.front();
            vertexs.pop();
            int d=depth.front();
            depth.pop();
            int id = std::get<0>(item);
            std::string str = std::get<1>(item);
            int weight = std::get<2>(item);
            if (d == dist) {
                continue;
            }
            auto neighbors = data_.GetNeighbors(id);
            for (auto n: neighbors) {
                if (!veritem.count(n)) {
                    float newWeight = weight + data_.GetEdgeWeight(id, n);
                    std::string newStr;
                    if(str.size()!=0)
                         newStr= str +"#"+std::to_string(data_.GetVertexLabel(n));
                    else
                         newStr = str +std::to_string(data_.GetVertexLabel(n));
                    que.push(std::make_tuple(n, newStr, newWeight));
                    std::set<int> newset = veritem;
                    newset.insert(n);
                    vertexs.push(newset);
                    int newd =d+1;
                    depth.push(newd);
                    TopologyIndex[d+1][i][newStr]++;
                    MNW[d+1][i][newStr] = std::max(MNW[d+1][i][newStr], newWeight);
                }
            }

        }
    }
}
void Instopk::GenerateMatchingOrder() {
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
        for (auto &other: q_nbrs)
        {
            if (visited[other])
            {
                order_csrs_[0][order_offs_[0][i]++] = other;
            }
        }
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
            auto& q_nbrs = query_.GetNeighbors(max_adjacent_u);
            for (auto &other: q_nbrs)
            {
                if (visited[other])
                {
                    order_csrs_[i][order_offs_[i][j]++] = other;

                }
            }
        }
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
                    std::cout << "-";
                    continue;
                }

                for (uint k = order_offs_[i][j - 1]; k < order_offs_[i][j]; k++)
                {

                    if (k == order_offs_[i][j - 1]) std::cout << "(";
                    std::cout << order_csrs_[i][k];
                    if (k != order_offs_[i][j] - 1) std::cout << ",";
                    else std::cout << ")";
                }

                if (j != query_.NumVertices() - 1)
                    std::cout << "-";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}
void Instopk::CreateQueryTopology() {
    //initial TopologyIndex
    int labelNum = query_.NumVLabels();
    int verticeNum = query_.NumVertices();

    //d=1
    for (int i = 0; i < labelNum; i++) {
        std::string str = std::to_string(i);
        querydToOrderType[1].push_back(str);
        for (int j = 0; j < verticeNum; j++) {
            std::map<std::string, int> s;
            s[str] = 0;
            queryTopologyIndex[1][j][str]=0;
        }
    }
    //d=others
    for (int d = 2; d <= dist; d++) {
        std::vector<std::string> strs = querydToOrderType[d - 1];
        for(int k=0;k<verticeNum;k++){
            for (auto s: strs) {
                for (int j = 0; j < labelNum; j++) {
                    std::string copys = s+"#";
                    copys += std::to_string(j);
                    querydToOrderType[d].push_back(copys);
                    std::map<std::string, int> smap;
                    smap[copys] = 0;
                    queryTopologyIndex[d][k][copys]=0;
                }
            }
        }
    }
    //bfs createindex
    for (int i = 0; i < verticeNum; i++) {
        std::queue<std::tuple<int, std::string>> que;//id,type,weight
        std::queue<int>depth;//id depth
        depth.push((0));
        std::queue<std::set<int>> vertexs;
        std::set<int> verset;
        verset.insert(i);
        vertexs.push(verset);
        que.push(std::make_tuple(i, ""));
        while (!que.empty()) {
            auto item = que.front();
            que.pop();
            auto veritem = vertexs.front();
            vertexs.pop();
            int d=depth.front();
            depth.pop();
            int id = std::get<0>(item);
            std::string str = std::get<1>(item);
            if (d == dist) {
                continue;
            }
            auto neighbors = query_.GetNeighbors(id);
            for (auto n: neighbors) {
                if (!veritem.count(n)) {

                    std::string newStr;
                    if(str.size()!=0)
                        newStr= str +"#"+std::to_string(query_.GetVertexLabel(n));
                    else
                        newStr = str +std::to_string(query_.GetVertexLabel(n));
                    que.push(std::make_tuple(n, newStr));
                    std::set<int> newset = veritem;
                    newset.insert(n);
                    vertexs.push(newset);
                    int newd =d+1;
                    depth.push(newd);
                    queryTopologyIndex[d+1][i][newStr]++;
                }
            }

        }
    }
    //create actualQueryEdges=vEdge
    //create queryEdgeToIndex 每条查询边对应的index索引
    //queryEdgesToEdgeType 查询边id对应的查询边标签类型
    // queryEdgeTypes  所有查询边标签种类
    // queryEdgeType2Edges 某种查询边标签对应哪些查询边id
    for(int i=0;i<query_.vEdge.size();i++){
        Edge edge=query_.vEdge[i];
        int v1=edge.GetV1();
        int v2=edge.GetV2();
        int w1=edge.GetV1Label();
        int w2=edge.GetV2Label();
        if(w1>w2)
            std::swap(w1,w2);
        std::string str= std::to_string(v1)+"#"+ std::to_string(v2);
        std::string typestr=std::to_string(w1)+"#"+ std::to_string(w2);
        queryEdgeToIndex[str]=i;
        queryEdgesToEdgeType[str]=typestr;
        queryEdgeTypes.insert(typestr);
        queryEdgeType2Edges[typestr].push_back(edge);
    }
}
void Instopk::SearchMatchesWithSortedList() {

    //初始化指针
    InitialPointers();
    //初始化查询候选节点
    InitialqueryCandidate();

    //搜索 首边利用path剪枝
    while(true){
        bool end=0;
        for(auto p:pointers){
            std::string s=p.first;
            int index=p.second;
            if(index>sortEdgeList[s].size()-1){
               end=1;
                break;
            }
        }
        if(end){
            break;
        }
        //找到最大边
        std::string max="";
        float maxWeight=0;
        int maxIndex=0;
        for(auto p:pointers){
            std::string s=p.first;
            int index=p.second;
            float val=sortEdgeList[s][index].GeteWeight();
            if(val>maxWeight)
            {
                maxWeight=val;
                max=s;
                maxIndex=index;
            }
        }
        //        std::vector<Edge>edgesOfMaxTypes=queryEdgeType2Edges[max];


        //一次匹配一条边
        std::vector<Edge>edgesofMaxType=queryEdgeType2Edges[max];
        for(Edge edge:edgesofMaxType) {
            std::set<Edge> considerEdges;
            std::set<int> considerEdgesIndex;
            considerEdges.insert(edge);
            std::string str = std::to_string(edge.GetV1()) + "#" + std::to_string(edge.GetV2());
            int index = queryEdgeToIndex[str];
            considerEdgesIndex.insert(index);
            std::set<std::vector<std::string>> curCandidates;
            std::set<std::vector<std::string>> pcCurr;
            Edge maxEdge = sortEdgeList[max][maxIndex];
            if (edge.GetV1Label() == maxEdge.GetV1Label() && edge.GetV2Label() == maxEdge.GetV2Label()) {
                std::vector<std::string> list(query_.NumEdges(), "");
                std::string str = maxEdge.GetV1() + "#" + maxEdge.GetV2();
                list[index] = str;
                pcCurr.insert(list);
            } else if (edge.GetV1Label() == maxEdge.GetV2Label() && edge.GetV2Label() == maxEdge.GetV1Label()) {
                std::vector<std::string> list(query_.NumEdges(), "");
                std::string str = maxEdge.GetV2() + "#" + maxEdge.GetV1();
                list[index] = str;
                pcCurr.insert(list);
            }
            //第一条边计算上界
            float upScoreOfNonConsideredEdges1 = 0;
            for (auto edge: queryEdgeToIndex) {
                std::string str = edge.first;
                int index = edge.second;
                if (!considerEdgesIndex.count(index)) {
                    std::string tmp = queryEdgesToEdgeType[str];
                    upScoreOfNonConsideredEdges1 += sortEdgeList[tmp][pointers[tmp]].GeteWeight();
                }
            }
            for (std::vector<std::string> pc: pcCurr) {
                //计算已经实例化边的分数
                float actualScore = 0;
                for (int i: considerEdgesIndex) {
                    std::string edge = pc[i];
                    std::pair<int, int> vertex = splitString(edge);
                    actualScore += data_.GetEdgeWeight(vertex.first, vertex.second);
                }
                //利用路径计算未实例化的边
                float upScoreOfNonConsideredEdges2 = FLT_MAX;
                upScoreOfNonConsideredEdges2 = GetUperBoundWithPath(considerEdgesIndex, pc);
                float upScoreOfNonConsideredEdges = std::min(upScoreOfNonConsideredEdges1,
                                                             upScoreOfNonConsideredEdges2);
                float upBound = actualScore + upScoreOfNonConsideredEdges;
                if (topKSet.size() == k) {
                    int n = query_.NumVertices();
                    float weight = topKSet.back()->getDensity();
                    if (upBound / sqrt(n) * (n - 1) > weight) {
                        curCandidates.insert(pc);
                    }
                } else {
                    curCandidates.insert(pc);
                }
            }
            if (curCandidates.size() == 0)
                continue;

            while (considerEdges.size() != queryEdgeToIndex.size()) {
                std::set<std::vector<std::string>> newCandidates;
                std::set<int> verticesCoverd;
                //寻找下一扩展边
                for (Edge edge: considerEdges) {
                    verticesCoverd.insert(edge.GetV1());
                    verticesCoverd.insert(edge.GetV2());
                }
                std::vector<Edge> nextEdgeCandidates;
                for (auto e: query_.vEdge) {
                    int v1 = e.GetV1();
                    int v2 = e.GetV2();
                    if (verticesCoverd.count(v1) && !considerEdges.count(e) ||
                        verticesCoverd.count(v2) && !considerEdges.count(e))
                        nextEdgeCandidates.push_back(e);
                }
                if (nextEdgeCandidates.size() == 0) {
                    std::cout << "error in nextEdgeCandidates" << std::endl;
                }
                Edge nextEdge = nextEdgeCandidates[0];
                int n1 = nextEdge.GetV1();
                int n2 = nextEdge.GetV2();
                int e1 = -1;
                int e2 = -1;
                int pos1 = -1;
                int pos2 = -1;
                //检验下一查询边的两个顶点是一个或者两个都在已经实例化的边中
                for (Edge edge: considerEdges) {
                    int v1 = edge.GetV1();
                    int v2 = edge.GetV2();
                    std::string s = std::to_string(v1) + "#" + std::to_string(v2);
                    if (v1 == n1 || v2 == n1) {
                        e1 = queryEdgeToIndex[s];
                        if (v1 == n1) {
                            pos1 = 1;
                        } else
                            pos1 = 2;
                    }
                    if (v1 == n2 || v2 == n2) {
                        e2 = queryEdgeToIndex[s];
                        if (v1 == n2)
                            pos1 = 1;
                        else
                            pos1 = 2;
                    }
                }
                int t1 = nextEdge.GetV1Label();
                int t2 = nextEdge.GetV2Label();
                if (t1 > t2)
                    std::swap(t1, t2);
                std::string nextEdgeType = std::to_string(t1) + "#" + std::to_string(t2);
                std::string nextEdgestr = std::to_string(n1) + "#" + std::to_string(n2);
                int edgeIndex = queryEdgeToIndex[nextEdgestr];
                considerEdges.insert(nextEdge);
                considerEdgesIndex.insert(edgeIndex);

                //找到对应匹配nextEdge的数据边
                for (auto c: curCandidates) {
                    int node1 = -1;
                    int node2 = -1;
                    std::vector<int> edgeIDs1;
                    std::vector<int> edgeIDs2;
                    if (e1 != -1) {
                        if (pos1 == 1)
                            node1 = splitString(c[e1]).first;
                        else
                            node1 = splitString(c[e1]).second;
                    }
                    if (e2 != -1) {
                        if (pos2 == 1)
                            node2 = splitString(c[e2]).first;
                        else
                            node2 = splitString(c[e2]).second;
                    }
                    if (node1 != -1)
                        edgeIDs1 = node2EdgeListPointers[nextEdgeType][node1];
                    if (node2 != -1)
                        edgeIDs2 = node2EdgeListPointers[nextEdgeType][node2];
                    if (node1 == node2)
                        continue;
                    std::set<std::vector<std::string>> pontentialCandidates;
                    if (node1 != -1 && node2 != -1) {
                        std::vector<int> intersection;
                        if (edgeIDs1.size() != 0 && edgeIDs2.size() != 0) {
                            for (int k: edgeIDs1)
                                if (std::find(edgeIDs2.begin(), edgeIDs2.end(), k) != edgeIDs2.end())
                                    intersection.push_back(k);
                        }
                        for (int k: intersection) {
                            std::vector<std::string> newCandi = c;
                            int flag = 0;
                            if (node1 > node2)
                                std::swap(node1, node2);
                            std::string ee = std::to_string(node1) + "#" + std::to_string(node2);
                            std::string ree= std::to_string(node2) + "#" + std::to_string(node1);
                            for (std::string s: c) {
                                if (s == ee) {
                                    flag = 1;
                                    break;
                                }
                            }
                            if (flag == 1)
                                continue;
                            newCandi[edgeIndex] = ee;
                            pontentialCandidates.insert(newCandi);
                        }
                    } else if (node1 == -1 && node2 != -1) {
                        if (edgeIDs2.size() != 0) {
                            for (int k: edgeIDs2) {
                                std::vector<std::string> newCandi = c;
                                int flag = 0;
                                std::string ee;
                                std::string ree;
                                if(sortEdgeList[nextEdgeType][k].GetV2()!=node2)
                                {
                                    ee=sortEdgeList[nextEdgeType][k].GetV2()+"#"+sortEdgeList[nextEdgeType][k].GetV1();
                                    ree = sortEdgeList[nextEdgeType][k].GetV1() + "#" +sortEdgeList[nextEdgeType][k].GetV2();
                                }
                                else{
                                    ee=sortEdgeList[nextEdgeType][k].GetV1()+"#"+sortEdgeList[nextEdgeType][k].GetV2();
                                    ree = sortEdgeList[nextEdgeType][k].GetV2() + "#" +sortEdgeList[nextEdgeType][k].GetV1();
                                }

                                for (std::string s: c) {
                                    if (s == ee||s==ree) {
                                        flag = 1;
                                        break;
                                    }
                                }
                                if (flag)
                                    continue;
                                newCandi[edgeIndex] = ee;
                                pontentialCandidates.insert(newCandi);
                            }
                        }
                    } else if (node1 != -1 && node2 == -1) {
                        if (edgeIDs1.size() != 0) {
                            for (int k: edgeIDs1) {
                                std::vector<std::string> newCandi = c;
                                int flag = 0;
                                std::string ee;
                                std::string ree;
                                if(sortEdgeList[nextEdgeType][k].GetV1()!=node1)
                                {
                                    ee=sortEdgeList[nextEdgeType][k].GetV2()+"#"+sortEdgeList[nextEdgeType][k].GetV1();
                                    ree = sortEdgeList[nextEdgeType][k].GetV1() + "#" +sortEdgeList[nextEdgeType][k].GetV2();
                                }
                                else{
                                    ee=sortEdgeList[nextEdgeType][k].GetV1()+"#"+sortEdgeList[nextEdgeType][k].GetV2();
                                    ree = sortEdgeList[nextEdgeType][k].GetV2() + "#" +sortEdgeList[nextEdgeType][k].GetV1();
                                }
                                for (std::string s: c) {
                                    if (s == ee||s==ree) {
                                        flag = 1;
                                        break;
                                    }
                                }
                                if (flag)
                                    continue;
                                newCandi[edgeIndex] = ee;
                                pontentialCandidates.insert(newCandi);
                            }
                        }
                    }

                    //计算加入nextEdge之后对应的上界
                    upScoreOfNonConsideredEdges1 = 0;
                    for (auto edge: queryEdgeToIndex) {
                        std::string str = edge.first;
                        int index = edge.second;
                        if (!considerEdgesIndex.count(index)) {
                            std::string tmp = queryEdgesToEdgeType[str];
                            upScoreOfNonConsideredEdges1 += sortEdgeList[tmp][pointers[tmp]].GeteWeight();
                        }
                    }
                    for (std::vector<std::string> pc: pontentialCandidates) {
                        float actualScore = 0;
                        for (int i: considerEdgesIndex) {
                            std::string edge = pc[i];
                            std::pair<int, int> vertex = splitString(edge);
                            actualScore += data_.GetEdgeWeight(vertex.first, vertex.second);
                        }
                        float upSCoreOfNonConsideredEdges2 = FLT_MAX;
                        float upSCoreOfNonConsideredEdges = std::min(upScoreOfNonConsideredEdges1,
                                                                     upSCoreOfNonConsideredEdges2);
                        float upBound = actualScore + upSCoreOfNonConsideredEdges;
                        if (topKSet.size() == k) {
                            int n = query_.NumVertices();
                            float weight = topKSet.back()->getDensity();
                            if (upBound / sqrt(n) * (n - 1) > weight) {
                                curCandidates.insert(pc);
                            }
                        } else {
                            newCandidates.insert(pc);

                        }
                    }
                }
                curCandidates = newCandidates;
            }
            //使用最后得到的结果更新topk
            for (std::vector<std::string> c: curCandidates) {
                float actualSCore = 0;
                int n = query_.NumEdges();
                std::vector<uint> m(n);
                for (int i = 0; i < n; i++) {
                    int q1 = query_.vEdge[i].GetV1();
                    int q2 = query_.vEdge[i].GetV2();
                    std::pair<int, int> vertex = splitString(c[i]);
                    int v1 = vertex.first;
                    int v2 = vertex.second;
                    if(query_.vEdge[i].GetV1Label()>query_.vEdge[i].GetV2Label()){
                        m[q1] = v2;
                        m[q2] = v1;
                    }
                    else{
                        m[q1] = v1;
                        m[q2] = v2;
                    }
                }
                for (int i: considerEdgesIndex) {
                    std::string edge = c[i];
                    std::pair<int, int> vertex = splitString(edge);
                    actualSCore += data_.GetEdgeWeight(vertex.first, vertex.second);
                }
                float density = actualSCore / sqrt(n) * (n - 1);
                //是否会与topk结果重复？？
                MatchRecord *r = new MatchRecord(density, 0, m);
                if(!isContainMatchRecord(r))
                {
                    addMatchRecords(r);
                }
            }
        }
          //移动指针，下一次遍历查询边
          std::vector<Edge>list=sortEdgeList[max];
          std::vector<Edge> arr1=queryEdgeType2Edges[max];
          int old=pointers[max];
          for(int c=pointers[max]+1;c<list.size();c++){
              int n1=list[c].GetV1();
              int n2=list[c].GetV2();
              int flag=0;
              for(Edge ee:arr1){
                  int v1=ee.GetV1();
                  int v2=ee.GetV2();
                  bool v1flag1=std::find(matchVertexCandidate[v1].begin(),matchVertexCandidate[v1].end(),n1)!=matchVertexCandidate[v1].end()?1:0;
                  bool v1flag2=std::find(matchVertexCandidate[v2].begin(),matchVertexCandidate[v2].end(),n2)!=matchVertexCandidate[v2].end()?1:0;
                  bool v2flag1=std::find(matchVertexCandidate[v1].begin(),matchVertexCandidate[v1].end(),n2)!=matchVertexCandidate[v1].end()?1:0;
                  bool v2flag2=std::find(matchVertexCandidate[v2].begin(),matchVertexCandidate[v2].end(),n1)!=matchVertexCandidate[v2].end()?1:0;
                  if((v1flag1&&v1flag2)||(v2flag1&&v2flag2)){
                      flag=1;
                      break;
                  }

              }
              if(flag==1)
              {
                  pointers[max]=c;
              }
          }
          if(old==pointers[max]){
              pointers[max]=list.size();
          }
          if(topKSet.size()==k){
              float maxUpperBoud=0;
              for(int i=0;i<query_.vEdge.size();i++)
              {
                  std::string edge= std::to_string(query_.vEdge[i].GetV1())+std::to_string(query_.vEdge[i].GetV2());
                  std::string type=queryEdgesToEdgeType[edge];
                  if(sortEdgeList[type].size()>pointers[type])
                      maxUpperBoud+=sortEdgeList[type][pointers[type]].GeteWeight();
              }
              float weight=topKSet.back()->getDensity();
              maxUpperBoud=maxUpperBoud/ sqrt(query_.NumVertices())*(query_.NumVertices());
              if(maxUpperBoud<weight)
                  return;

        }
    }

}
void Instopk::InitialPointers() {
    for(std::string edgeType:queryEdgeTypes){
        int pos=edgeType.find("#");
        int t1= atoi(edgeType.substr(0,pos).c_str());
        int len=edgeType.length()-1-pos;
        int t2= atoi(edgeType.substr(pos+1,len).c_str());
        std::string orderedEdgeType=edgeType;
        auto list=sortEdgeList[orderedEdgeType];
        auto arr1=queryEdgeType2Edges[orderedEdgeType];
        for(int c=0;c<list.size();c++){
            int n1=list[c].GetV1();
            int n2=list[c].GetV2();
            int flag=0;
            for(Edge &e:arr1){
                int v1=e.GetV1();
                int v2=e.GetV2();
                bool v1flag1=std::find(matchVertexCandidate[v1].begin(),matchVertexCandidate[v1].end(),n1)!=matchVertexCandidate[v1].end()?1:0;
                bool v1flag2=std::find(matchVertexCandidate[v2].begin(),matchVertexCandidate[v2].end(),n2)!=matchVertexCandidate[v2].end()?1:0;
                bool v2flag1=std::find(matchVertexCandidate[v1].begin(),matchVertexCandidate[v1].end(),n2)!=matchVertexCandidate[v1].end()?1:0;
                bool v2flag2=std::find(matchVertexCandidate[v2].begin(),matchVertexCandidate[v2].end(),n1)!=matchVertexCandidate[v2].end()?1:0;
                if((v1flag1&&v1flag2)||(v2flag1&&v2flag2)){
                    flag=1;
                    break;
                }
            }
            if(flag==1){
                pointers[orderedEdgeType]=c;
                break;
            }
        }
        if(!pointers.count(orderedEdgeType)){
            pointers[orderedEdgeType]=list.size();
        }
    }
}
std::pair<int,int> Instopk::splitString(std::string s) {
    int pos=s.find("#");
    int v1= atoi(s.substr(0,pos).c_str());
    int len=s.size()-1-pos;
    int v2= atoi(s.substr(pos+1,len).c_str());
    return std::make_pair(v1,v2);
}
std::set<Path> Instopk::GetPaths(int i,  std::vector<std::pair<int,int>>coverEdges) {
        std::queue<std::pair<int,Path>> que;//id,type,weight
        std::queue<int>depth;//id depth
        depth.push((0));
        std::queue<std::set<int>> vertexs;
        std::set<int> verset;
        verset.insert(i);
        vertexs.push(verset);
        Path path;
        path.nodes.push_back(i);
        que.push(std::make_pair(i,path));
        std::set<Path>pathSet;
        while (!que.empty()) {
            auto item = que.front();
            que.pop();
            auto veritem = vertexs.front();
            vertexs.pop();
            int d=depth.front();
            depth.pop();
            int id = std::get<0>(item);
             Path p = std::get<1>(item);
            if (d == dist) {
                continue;
            }
            auto neighbors = data_.GetNeighbors(id);
            for (auto n: neighbors) {
                int v1=std::min(id,(int)n);
                int v2=std::max(id,(int)n);
                std::pair<int,int>e=std::make_pair(v1,v2);
                if(std::find(coverEdges.begin(),coverEdges.end(),e)!=coverEdges.end())
                    continue;
                if (!veritem.count(n)) {
                    p.nodes.push_back(n);
                    pathSet.insert(p);
                    que.push(std::make_pair(n, p));
                    std::set<int> newset = veritem;
                    newset.insert(n);
                    vertexs.push(newset);
                    int newd =d+1;
                    depth.push(newd);
                }
            }
        }
        return pathSet;
}
float Instopk::GetUperBoundWithPath(std::set<int> consideredEdgeIndex, std::vector<std::string> pc) {
    float score=0;
    std::set<int>instantiatedVertices;
    std::vector<std::pair<int,int>>coverEdges;
    std::map<int,int>queryIdtoDataId;
    for(int c:consideredEdgeIndex){
        Edge e=query_.vEdge[c];
        int n1=e.GetV1();
        int n2=e.GetV2();
        coverEdges.push_back(std::make_pair(n1,n2));
        instantiatedVertices.insert(n1);
        instantiatedVertices.insert(n2);
        std::pair<int,int>vertices= splitString(pc[c]);
        queryIdtoDataId[n1]=vertices.first;
        queryIdtoDataId[n2]=vertices.second;
    }
    std::set<Path>gloabalPathSet;
    for(int i:instantiatedVertices){
        std::set<Path>result=GetPaths(i,coverEdges);
        gloabalPathSet.insert(result.begin(),result.end());
    }
    //计算所有从实例化节点开始的path的最大权重值
    while(coverEdges.size()!=query_.NumEdges()&&gloabalPathSet.size()!=0){
        int maxLength=0;
        std::vector<int>maxPath;
        for(Path p:gloabalPathSet){
            if(p.nodes.size()>maxLength){
                maxLength=p.nodes.size()-1;
                maxPath=p.nodes;
            }
        }
        //将maxPath经过的边都加入coverEdges;
        for(int i=1;i<maxPath.size();i++){
            int v1=std::min(maxPath[i],maxPath[i-1]);
            int v2=std::max(maxPath[i],maxPath[i-1]);
            coverEdges.push_back(std::make_pair(v1,v2));
        }
       //计算maxPath的估计权重
       std::string typeStr="";
       for(int i=1;i<maxPath.size();i++)
       {
           if(i==maxPath.size()-1)
               typeStr= std::to_string(query_.GetVertexLabel(maxPath[i]));
           else
           typeStr= std::to_string(query_.GetVertexLabel(maxPath[i]))+"#";
       }
       score+=MNW[maxLength][queryIdtoDataId[maxPath[0]]][typeStr];
       //将globalPath中包含maxPath中的边的路径都删除
       std::set<Path>gloablPathSetNew;
       for(Path p:gloabalPathSet){
           int flag=0;
           for(int ii=1;ii<p.nodes.size();ii++)
           {
              for(int jj=1;jj<maxPath.size();jj++){
                  if(maxPath[jj]==p.nodes[ii]&&maxPath[jj-1]==p.nodes[ii-1]){
                      flag=1;
                      break;
                  }
              }
           }
           if(flag==0)
               gloablPathSetNew.insert(p);
       }
       gloabalPathSet=gloablPathSetNew;
    }
    //计算剩余第三部分，不在d跳路径内也不在实例化边的节点
    for(auto edge:query_.vEdge){
        int n1=edge.GetV1();
        int n2=edge.GetV2();
        std::pair<int,int>vertices= std::make_pair(n1,n2);
        if(std::find(coverEdges.begin(), coverEdges.end(),vertices)==coverEdges.end()){
            int w1=edge.GetV1Label();
            int w2=edge.GetV2Label();
            if(w1>w2)
                std::swap(w1,w2);
            std::string type= std::to_string(w1)+"#"+ std::to_string(w2);
            int index=pointers[type];
            score+=sortEdgeList[type][index].GeteWeight();
        }
    }
    return score;
}

void Instopk::InitialqueryCandidate() {
    //初始化每个查询点的候选顶点
    matchVertexCandidate.reserve(query_.NumVertices());
    for (int i = 0; i < query_.NumVertices(); i++) {
        int qlabel = query_.GetVertexLabel(i);
        for (int j = 0; j < data_.NumVertices(); j++) {
            if (data_.GetVertexLabel(j) == qlabel) {
                //filter
                bool flag = true;
                for (int d = 1; d <= dist; d++) {
                    auto it = TopologyIndex[d][i].begin();
                    auto it2 = queryTopologyIndex[d][i].begin();
                    flag = true;
                    while (it2 != queryTopologyIndex[d][i].end()) {
                        if (it2->second > it->second) {
                            flag = false;
                            break;
                        }
                        it++;
                        it2++;
                    }
                    if (!flag) {
                        break;
                    }
                }
                if (flag) {
                    matchVertexCandidate[i].emplace_back(j);
                }
            }
        }
    }
}
bool Instopk::addMatchRecords(MatchRecord *r) {
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

            return true;
        }
        else{
            delete r;
            return false;
        }
        return true;
    }
}
bool Instopk::isContainMatchRecord(MatchRecord*m) {
    for(auto it=topKSet.begin();it!=topKSet.begin();it++){
        if((*(*it))==(*m))
        return true;
    }
    return false;
}


