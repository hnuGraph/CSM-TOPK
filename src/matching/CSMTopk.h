#ifndef MATCHING_CSMTOPK
#define MATCHING_CSMTOPK

#include <vector>
#include "../utils/types.h"
#include "../graph/graph.h"
#include "../matching/matching.h"
#include "unordered_map"
#include "sstream"
#include "../utils/Log.h"
#include "../utils/globals.h"
#include "../graph/MatchRecord.h"
#include "../graph/StarGraph.h"
#include "algorithm"
#include "cfloat"
#include "SingleCandidate.h"
#include "../utils/Timer.h"

class CSMTopk : public matching
{
public:
    // a list of matching orders starting from each query edge
    // the first matching order also applies to the initial matching
    std::vector<std::vector<uint> > order_vs_; //ƥ��ĵ�˳��
    std::vector<std::vector<uint> > order_csrs_;//ƥ��ڵ�ǰ���ھ�
    std::vector<std::vector<uint> > order_offs_;//ƥ��ڵ������
    std::vector<std::vector<uint>> order_vertex_index;//ÿ���ڵ��ڸ�ƥ�����е�λ�� order_vertex_index[0][u1]��ʾ��һ��ƥ�����У�u1������λ��
//        std::unordered_map<std::pair<uint,uint>,std::vector<MatchRecord*>,pair_hash> edgeMaps;//��Ӧ��ÿ���ߵĸ��������ṹ
    // uint s,t;//findMatchʱ���¼�ֵ
    //��¼top k��¼����
    std::vector<MatchRecord*> topKSet;
    //��¼����ƥ��Ľ��
    std::vector<std::vector<StarGraph*>>globalStarIndex;//��¼ÿ��ƥ����ÿ���ڵ����ھ��Լ����Ȩֵ
    std::vector<SingleCandidate>match;//ÿ���ڵ��ƥ����  vertex density
    std::vector<std::vector<SingleCandidate>>matchCandidate;//id density
    std::vector<float>suffixMax;
    std::vector<float>isolatedMax;
    std::vector<std::vector<std::vector<uint>>>rightNeighbor;//ƥ�������ţ�id��
    //std::vector<LocalIndex>queryLocalIndexs;
    std::vector<std::vector<std::vector<int>>>globalVkMatchUk;//<vk,ak,uk>
    std::vector<std::vector<uint>>labelToQueryVertex;//ÿ����ǩ��Ӧ�Ĳ�ѯ���ǩ
    std::vector<uint>queryVertexIndexInlabel;//ÿ����ѯ����label�����е�������
    std::vector<float>LocalStarIndex;//�ֲ�����
    std::vector<std::vector<int>>matchLeftNeighborSum;//���нڵ����ھӵĸ���
    std::vector<std::vector<size_t>>matchVetexLeftNeighbor;//����ƥ�����������ھ������
    std::vector<vector<StarGraph*>>matchVetexSumweight;//ÿ����ϸ��µõ������Ȩֵ
    std::vector<std::vector<size_t>>leftNeighborIdSum;//ÿ���ڵ����ھ�id��
    bool isUpdateIntopkset;
    int numAddTopk=0;
    int allMatchFind=0;
    int sumAllMatchFind=0;
    int sumDeleteallMatchFind=0;
    int numupdatestar=0;
    long IsearchSpace=0,DsearchSpace=0,IdeterminCandite=0,DdeterminCandite=0;

public:
    CSMTopk(Graph& query_graph, Graph& data_grasph,Subgraph& global_subgraph, uint max_num_results,
              bool print_prep, bool print_enum, bool homo);
    ~CSMTopk() override ;

    void Preprocessing() override;
    void InitialMatching(const std::string &path) override;

    void AddEdge(uint v1, uint v2, uint label,float weight) override;
    void RemoveEdge(uint v1, uint v2,uint label) override;
    void AddVertex(uint id, uint label) override;
    void RemoveVertex(uint id) override;
    void InitialTopK(const std::string &path) override;//�õ���ʼ��֮���Top k�������
    void updateTopK() override;
    void deleteUpdateTopK() override;
    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;
    void PrintAverageTime(int len);
private:
    void GenerateMatchingOrder();
    void FindMatches(uint flag,uint order_index, uint depth,
                     std::vector<uint> m, size_t &num_results,float density_s); //flag==0 initial flag==1 update
    int addMatchRecords(MatchRecord* r);//1 ��ʾ�ɹ����� 2��ʾ�ڵ��ظ� 3��ʾ�ڵ��ѱȵ�kthС
    void CreateStarIndex();
    float GetBackWeight(uint order_index,uint depth);
    void updateStarIndex(uint match_index,uint caddidate_v,uint candidate_u,int candidate_v_index);
    vector<int> EdgeisInMatchOrder(uint v1,uint v2,uint v1label,uint v2label,uint velabel);
    void searchMatches(int depth,uint matchorderindex,searchType flag);
    bool LabelFilter(uint data_v,uint query_v);
    void matchVertex(bool isFirstEdge,uint depth,uint data_v,float w);
    void popVertex(uint depth,uint data_v);
    void popVertex(uint data_v,uint matchorderindex,uint depth, const std::vector<uint>&uk_neighbor);
    void densityFilter(uint matchorder_index,uint depth, std::vector<SingleCandidate>&singleVertexCandidate);
    void createLabelToQueryVertex();

   bool updaterightNeighborCandidate(int matchorderindex,uint uk,uint uk_neigh,bool isFirstEdge,uint vk,const std::vector<uint>&uk_neighbor);
   void InitialLocalIndex(int matchorderindex);
   void getIntersetSingleCandidate( std::vector<SingleCandidate>&candidates,int matchorderindex,int depth);
   void createGlobalSubgraph();//����ȫ����ͼ
   bool updateGlobalSubgraph(uint v1,uint v2,uint label,float weight, std::vector<int>&match);
   bool updateGlobalGraphHelp(int m,uint u1,uint u2,uint u1label,uint u2label, uint v1,uint v2,uint v1label,uint v2label,uint elabel,const std::vector<std::vector<uint>>&mcandidate,bool &flag);
   void updateglobalVertexStarIndex(uint u1,uint v1,uint u1label,uint elabel,uint n, const std::vector<std::vector<uint>>&mcandidate);//����v1��ѡ�ڵ㣬������ȫ������
   void deleteGlobalSubgraph(uint v1, uint v2,uint elabel,float weight, std::vector<int> &match);
   void deleteUpdateglobalVertexStarIndex(uint u1,uint v1,uint v2,uint n);
   bool deleteGlobalSubgraphHelp(int m,uint u1,uint u2,uint u1label,uint u2label, uint v1,uint v2,uint v1label,uint v2label,
                                 uint elabel,float weight, std::vector<std::vector<uint>>&mcandidate);
   void deleteGlobalGraphCandidateEdges(uint m,uint u1,uint v1,std::vector<std::vector<uint>>&mcandidate);
   bool deleteMatchRecordWithEdge(uint v1, uint v1label,uint v2, uint v2label,uint label,std::vector<int> &match);
   bool SearchMatchesWithEdge(uint m,uint v1,uint v2,float weight,uint u1,uint u2,searchType type);

};

#endif //MATCHING_CSMTOPk
