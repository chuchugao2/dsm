//
// Created by 高楚楚 on 2023/5/16.
//

#include <set>
#include <cfloat>
#include <fstream>
#include "Instopk.h"
#include "../utils/Log.h"

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
                ,dToOrderType(d+1)
                ,TopologyIndex(d+1)
                ,queryTopologyIndex(d+1)
                ,MNW(d+1)
               {
                    for(int i=0;i<d+1;i++){
                        TopologyIndex[i].resize(data_.NumVertices());
                        MNW[i].resize(data_.NumVertices());
                        queryTopologyIndex[i].resize(query_.NumVertices());
                    }
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
    total_search_time.StartTimer();
    isUpdateIntopkset= false;
    data_.AddEdge(v1, v2, label, weight, timestamp, 1);
    //1.update sortEdgeList and pointers
    updateSortEdgelist(v1,v2, true);
    //2.update MNWindex Topologyindex
    updateMNWIndexAndDataTopologyIndex(v1,v2,label,weight, true);
    //3.searchTopkResult
    //3.1 updateQueryCandidate
    updateQueryCandidate(v1,v2,label, true);
    InitialPointers();
    SearchMatchesWithSortedList();
    total_search_time.StopTimer();
    updateTopK(0);
}

void Instopk::GetMemoryCost(size_t &num_edges, size_t &num_vertices) {
    num_edges = 0ul;
    num_vertices = 0ul;
}
void Instopk::InitialTopK(const std::string &path) {
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
    std::stringstream _ss1;
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
void Instopk::deleteUpdateTopK() {
    stringstream _ss;
#ifdef RESULT_TRACK
    _ss<<"after delete "<<std::endl;
    for(auto d:topKSet){
        _ss<<d->toString();
        Log::track2(_ss);
        _ss.clear();
        _ss.str("");
    }
#endif
}
void Instopk::deleteEdge(uint v1, uint v2) {
//todo deletion

}
void Instopk::updateTopK(uint num) {
    stringstream _ss;
    if (!isUpdateIntopkset) {
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
void Instopk::InitialMatching(const std::string &path) {
    if (!io::file_exists(path.c_str()))
    {
        std::cout << "the file not exit " << path << std::endl;
        std::vector<uint> m(query_.NumVertices(), UNMATCHED);
        //初始化查询候选节点
        InitialqueryCandidate();
        //初始化指针
        InitialPointers();
        SearchMatchesWithSortedList();
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
void Instopk::RemoveEdge(uint v1, uint v2) {
    total_delete_time.StartTimer();
    uint label=std::get<2>(data_.GetEdgeLabel(v1,v2));
    float weight=data_.GetEdgeWeight(v1,v2);
    data_.RemoveEdge(v1, v2);
    updateSortEdgelist(v1,v2, false);
    updateMNWIndexAndDataTopologyIndex(v1,v2,label, weight,false);
    updateQueryCandidate(v1,v2,label, false);
    InitialPointers();
    SearchMatchesWithSortedList();
    deleteUpdateTopK();
    total_delete_time.StopTimer();
}
void Instopk::RemoveVertex(uint id) {

}
void Instopk::PrintAverageTime(int len) {
    int ilen=10000-len;
    std::cout <<"average query graph degree:"<< std::fixed << std::setprecision(2)<<query_.NumEdges()*2.0/query_.NumVertices()<<endl;
    std::cout<<"average data graph degree:"<<std::fixed << std::setprecision(2)<<data_.NumEdges()*2.0/data_.NumVertices()<<endl;
    std::cout << "average serach time: " << std::fixed << std::setprecision(2)
              << total_search_time.GetTimer() * 1.0 / ilen << " microseconds" << endl;
    std::cout << "average delete search time:" << std::fixed << std::setprecision(2)
              << total_delete_time.GetTimer() * 1.0 / len << " microseconds" << endl;
#ifdef COMPUTE_TIME
    stringstream _ss;
    _ss<<total_search_time.GetTimer() * 1.0 / ilen<<","
       <<InitialSpace<<","
       <<total_delete_time.GetTimer() * 1.0 / len<<endl;
    Log::track3(_ss);
#endif
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
        const std::vector<Edge>&list=it->second;
        std::vector<std::vector<int>>map(data_.NumVertices());
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

    //搜索 首边利用path剪枝
    while(true){
        bool end=0;
        for(const auto &p:pointers){
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
    /*   Edge test= sortEdgeList[max][maxIndex];
        if((test.GetV1()==237894&&test.GetV2()==342222)||(test.GetV2()==237894&&test.GetV1()==342222)){
            std::cout<<"237894 342222"<<endl;
            std::cout<<maxIndex<<endl;
            std::cout<<max<<endl;
        }
        if((test.GetV1()==237894&&test.GetV2()==402161)||(test.GetV2()==237894&&test.GetV1()==402161)){
            std::cout<<"237894 402161"<<endl;
            std::cout<<maxIndex<<endl;
            std::cout<<max<<endl;
        }
        if((test.GetV1()==237894&&test.GetV2()==244418)||(test.GetV2()==237894&&test.GetV1()==244418)){
            std::cout<<"237894 244418"<<endl;
            std::cout<<maxIndex<<endl;
            std::cout<<max<<endl;
        }
        if((test.GetV1()==402161&&test.GetV2()==244418)||(test.GetV2()==402161&&test.GetV1()==244418)){
            std::cout<<"402161 244418"<<endl;
            std::cout<<maxIndex<<endl;
            std::cout<<max<<endl;
        }
*/
        //一次匹配一条边
        std::vector<Edge>edgesofMaxType=queryEdgeType2Edges[max];
        for(Edge edge:edgesofMaxType) {
            std::vector<Edge> considerEdges;
            std::vector<int> considerEdgesIndex;
            considerEdges.push_back(edge);
            std::string str = std::to_string(edge.GetV1()) + "#" + std::to_string(edge.GetV2());
            int index = queryEdgeToIndex[str];
            considerEdgesIndex.emplace_back(index);
            std::vector<std::vector<std::string>> curCandidates;
            std::vector<std::vector<std::string>> pcCurr;
            Edge maxEdge = sortEdgeList[max][maxIndex];
            if (edge.GetV1Label() == maxEdge.GetV1Label() && edge.GetV2Label() == maxEdge.GetV2Label()) {
                std::vector<std::string> list(query_.NumEdges(), "");
                std::string str = std::to_string(maxEdge.GetV1()) + "#" + std::to_string(maxEdge.GetV2());
                list[index] = str;
                pcCurr.push_back(list);
            }
            if (edge.GetV1Label() == maxEdge.GetV2Label() && edge.GetV2Label() == maxEdge.GetV1Label()) {
                std::vector<std::string> list(query_.NumEdges(), "");
                std::string str = std::to_string(maxEdge.GetV2()) + "#" + std::to_string(maxEdge.GetV1());
                list[index] = str;
                pcCurr.push_back(list);
            }
            //第一条边计算上界
            float upScoreOfNonConsideredEdges1 = 0;
            for (const auto &edge: queryEdgeToIndex) {
                std::string str = edge.first;
                int index = edge.second;
                if (!count(considerEdgesIndex,index)) {
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
                if(upScoreOfNonConsideredEdges2==-1){
                    continue;
                }
                float upScoreOfNonConsideredEdges = std::min(upScoreOfNonConsideredEdges1,
                                                             upScoreOfNonConsideredEdges2);
                float upBound = actualScore + upScoreOfNonConsideredEdges;
                if (topKSet.size() == k) {
                    int n = query_.NumVertices();
                    float weight = topKSet.back()->getDensity();
                    upBound=upBound / (sqrt(n) * (n - 1));
                    if ( !(upBound<weight)) {
                        curCandidates.push_back(pc);
                    }
                } else {
                    curCandidates.push_back(pc);
                }
            }
            if (curCandidates.size() == 0)
                continue;

            //匹配其他边
            while (considerEdges.size() != queryEdgeToIndex.size()) {
                std::vector<std::vector<std::string>> newCandidates;
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
                 /*   std::cout<<verticesCoverd.count(v1)<<endl;
                    std::cout<<count(considerEdges,e)<<endl;
                    std::cout<<verticesCoverd.count(v2)<<endl;
                    std::cout<<count(considerEdges,e)<<endl;*/
                    if ((verticesCoverd.count(v1) && !count(considerEdges,e)) ||
                            (verticesCoverd.count(v2) && !count(considerEdges,e)))
                    {
                        nextEdgeCandidates.push_back(e);
                        break;
                    }
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
                            pos2 = 1;
                        else
                            pos2 = 2;
                    }
                }
                int t1 = nextEdge.GetV1Label();
                int t2 = nextEdge.GetV2Label();
                if (t1 > t2)
                    std::swap(t1, t2);
                std::string nextEdgeType = std::to_string(t1) + "#" + std::to_string(t2);
                std::string nextEdgestr = std::to_string(n1) + "#" + std::to_string(n2);
                int edgeIndex = queryEdgeToIndex[nextEdgestr];
                considerEdges.push_back(nextEdge);
                considerEdgesIndex.emplace_back(edgeIndex);

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
                    std::vector<std::vector<std::string>> pontentialCandidates;
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
                            pontentialCandidates.push_back(newCandi);
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
                                    ee=std::to_string(sortEdgeList[nextEdgeType][k].GetV2())+"#"+std::to_string(sortEdgeList[nextEdgeType][k].GetV1());
                                    ree = std::to_string(sortEdgeList[nextEdgeType][k].GetV1()) + "#" +std::to_string(sortEdgeList[nextEdgeType][k].GetV2());
                                }
                                else{
                                    ee=std::to_string(sortEdgeList[nextEdgeType][k].GetV1())+"#"+std::to_string(sortEdgeList[nextEdgeType][k].GetV2());
                                    ree = std::to_string(sortEdgeList[nextEdgeType][k].GetV2()) + "#" +std::to_string(sortEdgeList[nextEdgeType][k].GetV1());
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
                                pontentialCandidates.push_back(newCandi);
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
                                    ee=std::to_string(sortEdgeList[nextEdgeType][k].GetV2())+"#"+std::to_string(sortEdgeList[nextEdgeType][k].GetV1());
                                    ree = std::to_string(sortEdgeList[nextEdgeType][k].GetV1()) + "#" +std::to_string(sortEdgeList[nextEdgeType][k].GetV2());
                                }
                                else{
                                    ee=std::to_string(sortEdgeList[nextEdgeType][k].GetV1())+"#"+std::to_string(sortEdgeList[nextEdgeType][k].GetV2());
                                    ree = std::to_string(sortEdgeList[nextEdgeType][k].GetV2()) + "#" +std::to_string(sortEdgeList[nextEdgeType][k].GetV1());
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
                                pontentialCandidates.push_back(newCandi);
                            }
                        }
                    }
                    //计算加入nextEdge之后对应的上界
                    upScoreOfNonConsideredEdges1 = 0;
                    for (auto edge: queryEdgeToIndex) {
                        std::string str = edge.first;
                        int index = edge.second;
                        if (!count(considerEdgesIndex,index)) {
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
                            upBound=upBound / (sqrt(n) * (n - 1));
                            if (!(upBound <weight)) {
                                newCandidates.push_back(pc);
                            }
                        } else {
                            newCandidates.push_back(pc);
                        }
                    }
                }
                curCandidates = newCandidates;
            }
            //使用最后得到的结果更新topk
            int n=query_.NumVertices();
            std::vector<uint> m(n);
            for(int i=0;i<curCandidates.size();i++){
                float actualSCore=0;
                const auto & matchresult=curCandidates[i];
                for(int i=0;i<matchresult.size();i++){
                    int q1=query_.vEdge[i].GetV1();
                    int q2=query_.vEdge[i].GetV2();
                    std::pair<int, int> vertex = splitString(matchresult[i]);
                    int v1 = vertex.first;
                    int v2 = vertex.second;
                    m[q1]=v1;
                    m[q2]=v2;
                    actualSCore+=data_.GetEdgeWeight(vertex.first, vertex.second);
                }
                /*if(m[0]==237894&&m[1]==342222&&m[2]==402161&&m[3]==244418){
                    std::cout<<"find result"<<endl;
                }*/
                float density = actualSCore / (sqrt(n) * (n - 1));
                //是否会与topk结果重复？？
                MatchRecord *r = new MatchRecord(density, m);
                if(!isContainMatchRecord(r))
                {
                    int flag=addMatchRecords(r);
                    if(flag)
                        isUpdateIntopkset= true;
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
                  break;
              }
          }
          if(old==pointers[max]){
              pointers[max]=list.size();
          }
          if(topKSet.size()==k){
              float maxUpperBoud=0;
              for(int i=0;i<query_.vEdge.size();i++)
              {
                  std::string edge= std::to_string(query_.vEdge[i].GetV1())+"#"+std::to_string(query_.vEdge[i].GetV2());
                  std::string type=queryEdgesToEdgeType[edge];
                  if(sortEdgeList[type].size()>pointers[type])
                      maxUpperBoud+=sortEdgeList[type][pointers[type]].GeteWeight();
              }
              float weight=topKSet.back()->getDensity();
              int n=query_.NumVertices();
              maxUpperBoud=maxUpperBoud/ (sqrt(n)*(n-1));
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
        auto &list=sortEdgeList[orderedEdgeType];
        auto &arr1=queryEdgeType2Edges[orderedEdgeType];
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
std::set<Path> Instopk::GetPaths(int i,  const std::vector<std::pair<int,int>>&coverEdges) {
        std::queue<std::pair<int,Path>> que;//id,type,weight
        std::queue<int>depth;//id depth
        depth.push((0));
        std::queue<std::set<int>> vertexs;
        std::set<int> verset;
        const pair<int,int>p=coverEdges[0];
        verset.insert(p.first);
        verset.insert(p.second);
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
            const auto &neighbors = query_.GetNeighbors(id);
            for (auto n: neighbors) {
                int v1=std::min(id,(int)n);
                int v2=std::max(id,(int)n);
                std::pair<int,int>e=std::make_pair(v1,v2);
                if(std::find(coverEdges.begin(),coverEdges.end(),e)!=coverEdges.end())
                    continue;
                if (!veritem.count(n)) {
                    Path newp=p;
                    newp.nodes.push_back(n);
                    pathSet.insert(newp);
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
float Instopk::GetUperBoundWithPath(const std::vector<int> &consideredEdgeIndex, const std::vector<std::string>&pc) {
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
            if(p.nodes.size()-1>maxLength){
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
               typeStr+= std::to_string(query_.GetVertexLabel(maxPath[i]));
           else
           typeStr+= std::to_string(query_.GetVertexLabel(maxPath[i]))+"#";
       }
       if(MNW[maxLength][queryIdtoDataId[maxPath[0]]][typeStr]==0)
           return -1;
       score+=MNW[maxLength][queryIdtoDataId[maxPath[0]]][typeStr];
       //将globalPath中包含maxPath中的边的路径都删除
       std::set<Path>gloablPathSetNew;
       for(Path p:gloabalPathSet){
           int flag=0;
           if(maxPath==p.nodes)
           {
               continue;
           }
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
    matchVertexCandidate.resize(query_.NumVertices());
    for (int i = 0; i < query_.NumVertices(); i++) {
        int qlabel = query_.GetVertexLabel(i);
        for (int j = 0; j < data_.NumVertices(); j++) {
            if (data_.GetVertexLabel(j) == qlabel) {
                //filter
                bool flag = true;
                for (int d = 1; d <= dist; d++) {
                    auto it = TopologyIndex[d][j].begin();
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
    for(auto it=topKSet.begin();it!=topKSet.end();it++){
        if((*(*it))==(*m))
        return true;
    }
    return false;
}
bool Instopk::count(std::vector<int> &t, int e) {
    if(std::find(t.begin(), t.end(),e)==t.end())
        return false;
    else
        return true;
}
bool Instopk::count(std::vector<Edge> &t, Edge e) {
    if(std::find(t.begin(), t.end(),e)==t.end())
        return false;
    else
        return true;
}
void Instopk::updateSortEdgelist(uint v1, uint v2,bool isAdd) {
    int v1label=data_.GetVertexLabel(v1);
    int v2label=data_.GetVertexLabel(v2);
    if(v1label>v2label){
        swap(v1label,v2label);
    }
    std::string str=std::to_string(v1label)+"#"+std::to_string(v2label);
    if(isAdd){
        Edge edge=data_.vEdge[data_.vEdge.size()-1];
        insertEdgeToSortEdgelist(str,edge);
    }
    else{
        int elabel=std::get<2>(data_.GetEdgeLabel(v1,v2));
        if(v1>v2)
        {
            swap(v1,v2);
        }
        Edge edge(v1,v2,v1label,v2label,elabel);
        deleteEdgeInSortEdgelist(str, edge);
    }
}
void Instopk::updateMNWIndexAndDataTopologyIndex(uint v1,uint v2,uint label,float weight,bool isAdd) {
    if(isAdd){
        //update TopologyIndex and MNWindex
        uint v1label=data_.GetVertexLabel(v1);
        uint v2label=data_.GetVertexLabel(v2);
        int n=data_.NumVertices();
        for(int i=0;i<n;i++){
            if(i==v1){
                std::string str=std::to_string(v2label);
                TopologyIndex[1][i][str]++;
                MNW[1][i][str]=std::max( MNW[1][i][str],weight);
                const std::vector<uint>&neighbors=data_.GetNeighbors(v2);
                for(uint n:neighbors){
                    float w=data_.GetEdgeWeight(v2,n);
                    w=w+weight;
                    uint nlabel=data_.GetVertexLabel(n);
                    str=str+"#"+std::to_string(nlabel);
                    TopologyIndex[2][i][str]++;
                    MNW[2][i][str]=std::max(MNW[2][i][str],w);
                }
            }
            else if(i==v2){
                std::string str=std::to_string(v1label);
                TopologyIndex[1][i][str]++;
                MNW[1][i][str]=std::max( MNW[1][i][str],weight);
                const std::vector<uint>&neighbors=data_.GetNeighbors(v1);
                for(uint n:neighbors){
                    float w=data_.GetEdgeWeight(v1,n);
                    w=w+weight;
                    uint nlabel=data_.GetVertexLabel(n);
                    str=str+"#"+std::to_string(nlabel);
                    TopologyIndex[2][i][str]++;
                    MNW[2][i][str]=std::max(MNW[2][i][str],w);
                }
            }
            else{
                if(isNeighbor(i,v1)){
                    std::string str=std::to_string(v1label)+"#"+std::to_string(v2label);
                    TopologyIndex[2][i][str]++;
                    float w=data_.GetEdgeWeight(i,v1);
                    w+=weight;
                    MNW[2][i][str]=std::max(MNW[2][i][str],w);
                }
                else if(isNeighbor(i,v2)){
                    std::string str=std::to_string(v2label)+"#"+std::to_string(v1label);
                    TopologyIndex[2][i][str]++;
                    float w=data_.GetEdgeWeight(i,v2);
                    w+=weight;
                    MNW[2][i][str]=std::max(MNW[2][i][str],w);
                }
            }
        }
    }
    else{
        //update TopologyIndex and MNWindex
        uint v1label=data_.GetVertexLabel(v1);
        uint v2label=data_.GetVertexLabel(v2);
        int n=data_.NumVertices();
        std::set<uint>updatelabels;
        updatelabels.insert(v1);
        updatelabels.insert(v2);
        for(int i=0;i<n;i++){
            if(i==v1){
                std::string str=std::to_string(v2label);
                TopologyIndex[1][i][str]--;
                //update MNW[1][i][str]
                MNW[1][i][str]=0;
                //update MNW[2][i][str]
                const std::vector<uint>&neighbors2=data_.GetNeighbors(v2);
                for(uint n:neighbors2){
                    uint nlabel=data_.GetVertexLabel(n);
                    float w=data_.GetEdgeWeight(v2,n);
                    w=w+weight;
                    string tmpstr=str+"#"+std::to_string(nlabel);
                    TopologyIndex[2][i][tmpstr]--;
                    if(MNW[2][i][tmpstr]==w)
                        MNW[2][i][tmpstr]=0;
                }
            }
            else if(i==v2){
                std::string str=std::to_string(v2label);
                TopologyIndex[1][i][str]--;
                //update MNW[1][i][str]
                MNW[1][i][str]=0;
                //update MNW[2][i][str]
                const std::vector<uint>&neighbors2=data_.GetNeighbors(v2);
                for(uint n:neighbors2){
                    uint nlabel=data_.GetVertexLabel(n);
                    float w=data_.GetEdgeWeight(v2,n);
                    w=w+weight;
                    string tmpstr=str+"#"+std::to_string(nlabel);
                    TopologyIndex[2][i][tmpstr]--;
                    if(MNW[2][i][tmpstr]==w)
                        MNW[2][i][tmpstr]=0;
                }
            }
            else{
                if(isNeighbor(i,v1)){
                    std::string str=std::to_string(v1label)+"#"+std::to_string(v2label);
                    TopologyIndex[2][i][str]--;
                    //重算str=v1label$v2label
                    const std::vector<uint>&neighbors1=data_.GetNeighbors(i);
                    float w=data_.GetEdgeWeight(i,v1);
                    w+=weight;
                    if( MNW[2][i][str]==w){
                        updatelabels.insert(data_.GetVertexLabel(i));
                        MNW[2][i][str]=0;
                    }
                }
                else if(isNeighbor(i,v2)){
                    std::string str=std::to_string(v2label)+"#"+std::to_string(v1label);
                    TopologyIndex[2][i][str]--;
                    //重算str=v1label$v2label
                    const std::vector<uint>&neighbors1=data_.GetNeighbors(i);
                    float w=data_.GetEdgeWeight(i,v2);
                    w+=weight;
                    if( MNW[2][i][str]==w){
                        updatelabels.insert(data_.GetVertexLabel(i));
                        MNW[2][i][str]=0;
                    }
                }
            }
        }
        //update MNW[2][i][str]  start label in updatelabel
        for (int i = 0; i < data_.NumVertices(); i++) {
            uint ilabel=data_.GetVertexLabel(i);
            if(updatelabels.count(i))
                continue;
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
}
void Instopk::insertEdgeToSortEdgelist(std::string &str, Edge edge) {
    std::vector<Edge>& list=sortEdgeList[str];
    int n=list.size();
    for(int j=n-1;j>=0;j--){
        if(list[j].GeteWeight()>edge.GeteWeight()){
            list.insert(list.begin()+j+1,edge);
        }
        if(j==0){
         list.insert(list.begin(),edge);
        }
    }
    if(n==0){
        list.insert(list.begin(),edge);
    }
    //update node2EdgelistPointers in str
    const std::vector<Edge>&strlist=sortEdgeList[str];
    std::vector<std::vector<int>>map(data_.NumVertices());
    for(int i=0;i<strlist.size();i++){
        int v1=strlist[i].GetV1();
        int v2=strlist[i].GetV2();
        map[v1].push_back(i);
        map[v2].push_back(i);
    }
    node2EdgeListPointers[str]=map;
}
void Instopk::updateQueryCandidate(uint v1,uint v2,uint label,bool isAdd) {
    if(isAdd){
        uint v1label=data_.GetVertexLabel(v1);
        uint v2label=data_.GetVertexLabel(v2);
        std::string str=std::to_string(v1label)+"#"+std::to_string(v2label);
        int n=data_.NumVertices();
        int m=query_.NumVertices();
        for(int i=0;i<n;i++){
            if(i==v1){
                for(int j=0;j<m;j++){
                    if(query_.GetVertexLabel(j)==v1label){
                        if(count(matchVertexCandidate[j],v1)){
                            continue;
                        }
                        else{
                            bool flag= true;
                            for(int d=1;d<=dist;d++){
                                auto it=TopologyIndex[d][i].begin();
                                auto it2=queryTopologyIndex[d][j].begin();
                                flag= true;
                                while(it2!=queryTopologyIndex[d][j].end()){
                                    if(it2->second>it->second){
                                        flag= false;
                                        break;
                                    }
                                    it++;
                                    it2++;
                                }
                                if(!flag){
                                    break;
                                }
                            }
                            if (flag) {
                                matchVertexCandidate[j].emplace_back(i);
                            }
                        }
                    }
                }
            }
            else if(i==v2){
                for(int j=0;j<m;j++){
                    if(query_.GetVertexLabel(j)==v2label){
                        if(count(matchVertexCandidate[j],v2)){
                            continue;
                        }
                        else{
                            bool flag= true;
                            for(int d=1;d<=dist;d++){
                                auto it=TopologyIndex[d][i].begin();
                                auto it2=queryTopologyIndex[d][j].begin();
                                flag= true;
                                while(it2!=queryTopologyIndex[d][j].end()){
                                    if(it2->second>it->second){
                                        flag= false;
                                        break;
                                    }
                                    it++;
                                    it2++;
                                }
                                if(!flag){
                                    break;
                                }
                            }
                            if (flag) {
                                matchVertexCandidate[j].emplace_back(i);
                            }
                        }
                    }
                }
            }
            else{
                if(isNeighbor(i,v1)|| isNeighbor(i,v2)){
                    uint ilabel=data_.GetVertexLabel(i);
                    for(int j=0;j<m;j++){
                        if(query_.GetVertexLabel(j)==ilabel){
                            if(count(matchVertexCandidate[j],i)){
                                continue;
                            }
                            else{
                                auto it=TopologyIndex[2][i].begin();
                                auto it2=queryTopologyIndex[2][j].begin();
                                bool flag= true;
                                while(it2!=queryTopologyIndex[2][j].end()){
                                    if(it2->second>it->second){
                                        flag= false;
                                        break;
                                    }
                                    it++;
                                    it2++;
                                }
                                if(flag){
                                    matchVertexCandidate[j].emplace_back(i);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else{
        uint v1label=data_.GetVertexLabel(v1);
        uint v2label=data_.GetVertexLabel(v2);
        std::string str=std::to_string(v1label)+"#"+std::to_string(v2label);
        int n=data_.NumVertices();
        int m=query_.NumVertices();
        for(int i=0;i<n;i++){
            if(i==v1){
                for(int j=0;j<m;j++){
                    if(query_.GetVertexLabel(j)==v1label){
                        if(!count(matchVertexCandidate[j],v1)){
                            continue;
                        }
                        else{
                            bool flag= true;
                            for(int d=1;d<=dist;d++){
                                auto it=TopologyIndex[d][i].begin();
                                auto it2=queryTopologyIndex[d][j].begin();
                                flag= true;
                                while(it2!=queryTopologyIndex[d][j].end()){
                                    if(it2->second>it->second){
                                        flag= false;
                                        break;
                                    }
                                    it++;
                                    it2++;
                                }
                                if(!flag){
                                    break;
                                }
                            }
                            if (!flag) {
                                for(int m=0;m<matchVertexCandidate[j].size();m++){
                                    if(matchVertexCandidate[j][m]==i){
                                        matchVertexCandidate[j].erase(matchVertexCandidate[j].begin()+m);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else if(i==v2){
                for(int j=0;j<m;j++){
                    if(query_.GetVertexLabel(j)==v2label){
                        if(!count(matchVertexCandidate[j],v2)){
                            continue;
                        }
                        else{
                            bool flag= true;
                            for(int d=1;d<=dist;d++){
                                auto it=TopologyIndex[d][i].begin();
                                auto it2=queryTopologyIndex[d][j].begin();
                                flag= true;
                                while(it2!=queryTopologyIndex[d][j].end()){
                                    if(it2->second>it->second){
                                        flag= false;
                                        break;
                                    }
                                    it++;
                                    it2++;
                                }
                                if(!flag){
                                    break;
                                }
                            }
                            if (!flag) {
                                matchVertexCandidate[j].emplace_back(i);
                            }
                        }
                    }
                }
            }
            else{
                if(isNeighbor(i,v1)|| isNeighbor(i,v2)){
                    uint ilabel=data_.GetVertexLabel(i);
                    for(int j=0;j<m;j++){
                        if(query_.GetVertexLabel(j)==ilabel){
                            if(count(matchVertexCandidate[j],i)){
                                continue;
                            }
                            else{
                                auto it=TopologyIndex[2][i].begin();
                                auto it2=queryTopologyIndex[2][j].begin();
                                bool flag= true;
                                while(it2!=queryTopologyIndex[2][j].end()){
                                    if(it2->second>it->second){
                                        flag= false;
                                        break;
                                    }
                                    it++;
                                    it2++;
                                }
                                if(flag){
                                    for(int m=0;m<matchVertexCandidate[j].size();m++){
                                        if(matchVertexCandidate[j][m]==i){
                                            matchVertexCandidate[j].erase(matchVertexCandidate[j].begin()+m);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

}
bool Instopk::isNeighbor(uint curVertex, uint v1) {
    const std::vector<uint>&neighbors=data_.GetNeighbors(curVertex);
    for(uint n:neighbors){
        if(n==v1)
        {
            return true;
        }
    }
    return false;
}
void Instopk::deleteEdgeInSortEdgelist(std::string &str, Edge edge) {
    std::vector<Edge>& list=sortEdgeList[str];
    int n=list.size();
    for(int j=n-1;j>=0;j--){
        if(list[j]==edge){
            list.erase(list.begin()+j);
            break;
        }
    }
    //update node2EdgelistPointers in str
    const std::vector<Edge>&strlist=sortEdgeList[str];
    std::vector<std::vector<int>>map(data_.NumVertices());
    for(int i=0;i<strlist.size();i++){
        int v1=strlist[i].GetV1();
        int v2=strlist[i].GetV2();
        map[v1].push_back(i);
        map[v2].push_back(i);
    }
    node2EdgeListPointers[str]=map;
}

