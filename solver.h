#pragma once

#include <queue>
#include <vector>
#include <list>
#include <string>

typedef int32_t TVertex;

class UnknownEdgeException : public std::out_of_range {
 public:
  explicit UnknownEdgeException(int32_t edge)
      : std::out_of_range("Edge " + std::to_string(edge) + " doesn't exist") {}
};

class UnknownVertexException : public std::out_of_range {
 public:
  explicit UnknownVertexException(TVertex vertex)
      : std::out_of_range("Vertex " + std::to_string(vertex) +
                          " doesn't exist") {}
};

class FlowMoreThanCapacityException : public std::logic_error {
 public:
  explicit FlowMoreThanCapacityException(int32_t edge)
      : std::logic_error("Edge " + std::to_string(edge) +
                         ": flow more than capacity overflowed") {}
};

template <typename TFlow = long long>
struct Edge {
  TVertex start, finish;
  TFlow capacity, flow;
  Edge(TVertex v, TVertex to, TFlow c, TFlow f)
      : start(v), finish(to), capacity(c), flow(f) {}
};

template <typename TFlow = long long>
class Network {
 public:
  Network(size_t n = 1, TVertex s = 0, TVertex t = 0)
      : nVertices_(n),
        source_(s),
        sink_(t),
        edges_(),
        begins_(n, -1),
        nexts_() {
    if (source_ < 0 || source_ >= nVertices_)
      throw UnknownVertexException(source_);

    if (sink_ < 0 || sink_ >= nVertices_)
      throw UnknownVertexException(sink_);
  }

  void addEdge(TVertex start, TVertex finish, TFlow capacity) {
    addEdgeLocal_(start, finish, capacity, TFlow(0));
    addEdgeLocal_(finish, start, TFlow(0), TFlow(0));
  }

  TVertex getSource() { return source_; }

  TVertex getSink() { return sink_; }

  size_t getNumVertices() { return nVertices_; }
  size_t getNumEdges() { return edges_.size(); }

  class Iterator {
   public:
    explicit Iterator(Network& myNetwork, int32_t numEdge)
        : myNetwork_(myNetwork), numEdge_(numEdge) {}

    Iterator& operator=(Iterator it) {
      this->myNetwork_ = it.myNetwork_;
      this->numEdge_ = it.numEdge_;
      return *this;
    }

    bool isValid() const { return numEdge_ != -1; }

    void goNext() {
      if (numEdge_ < 0 || numEdge_ >= myNetwork_.getNumEdges())
        throw UnknownEdgeException(numEdge_);
      numEdge_ = myNetwork_.nexts_[numEdge_];
    }

    void clear() {
      if (numEdge_ < 0 || numEdge_ >= myNetwork_.getNumEdges())
        throw UnknownEdgeException(numEdge_);
      myNetwork_.clearEdge_(numEdge_);
    }

    TFlow getCapacity() {
      if (numEdge_ < 0 || numEdge_ >= myNetwork_.getNumEdges())
        throw UnknownEdgeException(numEdge_);

      return myNetwork_.edges_[numEdge_].capacity;
    }

    TFlow getFlow() {
      if (numEdge_ < 0 || numEdge_ >= myNetwork_.getNumEdges())
        throw UnknownEdgeException(numEdge_);

      return myNetwork_.edges_[numEdge_].flow;
    }

    TVertex getStart() {
      if (numEdge_ < 0 || numEdge_ >= myNetwork_.getNumEdges())
        throw UnknownEdgeException(numEdge_);

      return myNetwork_.edges_[numEdge_].start;
    }

    TVertex getFinish() {
      if (numEdge_ < 0 || numEdge_ >= myNetwork_.getNumEdges())
        throw UnknownEdgeException(numEdge_);

      return myNetwork_.edges_[numEdge_].finish;
    }

    TFlow getResidualCapacity() {
      if (numEdge_ < 0 || numEdge_ >= myNetwork_.getNumEdges())
        throw UnknownEdgeException(numEdge_);

      return myNetwork_.edges_[numEdge_].capacity -
             myNetwork_.edges_[numEdge_].flow;
    }

    Iterator getBackEdge() {
      if (numEdge_ < 0 || numEdge_ >= myNetwork_.getNumEdges())
        throw UnknownEdgeException(numEdge_);

      return Iterator(myNetwork_, numEdge_ ^ 1);
    }

    void push(TFlow f) {
      if (numEdge_ < 0 || numEdge_ >= myNetwork_.getNumEdges())
        throw UnknownEdgeException(numEdge_);

      myNetwork_.push_(numEdge_, f);
    }

   private:
    Network& myNetwork_;
    size_t numEdge_;
  };

  Iterator begin(TVertex vertex) {
    if (vertex < 0 || vertex >= getNumVertices())
      throw UnknownVertexException(vertex);

    return Iterator(*this, begins_[vertex]);
  }

  TFlow getFlow() {
    TFlow ans = 0;

    for (Network::Iterator it = begin(source_); it.isValid(); it.goNext())
      ans += it.getFlow();

    return ans;
  }

  void clearFlow() {
    for (TVertex vertex = 0; vertex < nVertices_; ++vertex)
      for (Network::Iterator it = begin(vertex); it.isValid(); it.goNext())
        it.clear();
  }

 private:
  size_t nVertices_;
  TVertex source_, sink_;
  std::vector<Edge<TFlow> > edges_;
  std::vector<TVertex> begins_, nexts_;

  void push_(size_t numEdge, TFlow f) {
    edges_[numEdge].flow += f;

    if (edges_[numEdge].flow > edges_[numEdge].capacity)
      throw FlowMoreThanCapacityException(numEdge);

    edges_[numEdge ^ 1].flow -= f;

    if (edges_[numEdge ^ 1].flow > edges_[numEdge ^ 1].capacity)
      throw FlowMoreThanCapacityException(numEdge ^ 1);
  }

  void clearEdge_(size_t numEdge) {
    edges_[numEdge].flow = 0;
    edges_[numEdge ^ 1].flow = 0;
  }

  void addEdgeLocal_(TVertex start,
                     TVertex finish,
                     TFlow capacity,
                     TFlow flow) {
    if (start < 0 || start >= nVertices_)
      throw UnknownVertexException(start);

    if (finish < 0 || finish >= nVertices_)
      throw UnknownVertexException(finish);

    edges_.emplace_back(Edge<TFlow>(start, finish, capacity, flow));

    if (flow > capacity)
      throw FlowMoreThanCapacityException(edges_.size() - 1);

    nexts_.push_back(begins_[start]);
    begins_[start] = edges_.size() - 1;
  }
};

template <typename TFlow = long long>
class MalhotraKumarMaheshwari {
 public:
  MalhotraKumarMaheshwari(Network<TFlow>& Net) : myNetwork_(Net) {}

  TFlow findMaxFlow() {
    TFlow maxFlow = 0;

    while (initializeDistance_()) {
      initializePotentials_();

      while (true) {
        TFlow addFlow = findBlockingFlow_();
        if (addFlow == -1)
          break;

        maxFlow += addFlow;
      }
    }

    return maxFlow;
  }

 private:
  Network<TFlow>& myNetwork_;
  std::vector<int32_t> dist_;
  std::vector<TFlow> pIn_, pOut_;

  void initializePotentials_() {
    pIn_.assign(myNetwork_.getNumVertices(), 0);
    pOut_.assign(myNetwork_.getNumVertices(), 0);

    for (TVertex vertex = 0; vertex < myNetwork_.getNumVertices(); ++vertex) {
      for (auto it = myNetwork_.begin(vertex); it.isValid(); it.goNext()) {
        if (dist_[it.getStart()] + 1 == dist_[it.getFinish()] &&
            dist_[it.getStart()] != -1) {
          pOut_[it.getStart()] += it.getResidualCapacity();
          pIn_[it.getFinish()] += it.getResidualCapacity();
        }
      }
    }
  }

  TFlow getPotential_(TVertex vertex) {
    if (vertex == myNetwork_.getSource())
      return pOut_[vertex];

    if (vertex == myNetwork_.getSink())
      return pIn_[vertex];

    return std::min(pIn_[vertex], pOut_[vertex]);
  }

  bool initializeDistance_() {
    std::queue<TVertex> queue;
    dist_.assign(myNetwork_.getNumVertices(), -1);
    dist_[myNetwork_.getSource()] = 0;
    queue.push(myNetwork_.getSource());

    while (!queue.empty()) {
      TVertex vertex = queue.front();
      queue.pop();

      for (auto it = myNetwork_.begin(vertex); it.isValid(); it.goNext()) {
        if (it.getResidualCapacity() > 0 && dist_[it.getFinish()] == -1) {
          dist_[it.getFinish()] = dist_[it.getStart()] + 1;
          queue.push(it.getFinish());
        }
      }
    }

    return dist_[myNetwork_.getSink()] != -1;
  }

  TVertex getMinPotentialVertex_() {
    TVertex minPotentialVertex = -1;

    for (TVertex vertex = 0; vertex < myNetwork_.getNumVertices(); ++vertex) {
      if (dist_[vertex] != -1 &&
          (minPotentialVertex == -1 ||
           getPotential_(vertex) < getPotential_(minPotentialVertex)))
        minPotentialVertex = vertex;
    }

    return minPotentialVertex;
  }

  void deleteVertex_(TVertex vertex) {
    for (auto it = myNetwork_.begin(vertex); it.isValid(); it.goNext()) {
      if (dist_[it.getStart()] + 1 == dist_[it.getFinish()] &&
          dist_[it.getStart()] != -1) {
        pOut_[it.getStart()] -= it.getResidualCapacity();
        pIn_[it.getFinish()] -= it.getResidualCapacity();
      }

      if (dist_[it.getBackEdge().getStart()] + 1 ==
              dist_[it.getBackEdge().getFinish()] &&
          dist_[it.getBackEdge().getStart()] != -1) {
        pOut_[it.getBackEdge().getStart()] -=
            it.getBackEdge().getResidualCapacity();
        pIn_[it.getBackEdge().getFinish()] -=
            it.getBackEdge().getResidualCapacity();
      }
    }

    dist_[vertex] = -1;
  }

  void tryPush_(TVertex startVertex, TFlow flow, int32_t type) {
    std::queue<TVertex> queue;
    std::vector<TFlow> toPush;

    toPush.assign(myNetwork_.getNumVertices(), 0);
    toPush[startVertex] += flow;
    queue.push(startVertex);

    while (!queue.empty()) {
      TVertex vertex = queue.front();
      queue.pop();

      for (auto it = myNetwork_.begin(vertex); it.isValid(); it.goNext()) {
        if (toPush[vertex] == 0)
          break;

        if (dist_[it.getStart()] + type != dist_[it.getFinish()] ||
            dist_[it.getStart()] == -1 || dist_[it.getFinish()] == -1)
          continue;

        TFlow pushNext = toPush[vertex];

        if (type == 1)
          pushNext = std::min(pushNext, it.getResidualCapacity());
        else
          pushNext = std::min(pushNext, it.getBackEdge().getResidualCapacity());

        if (pushNext == 0)
          continue;

        if (type == 1) {
          it.push(pushNext);
          pOut_[it.getStart()] -= pushNext;
          pIn_[it.getFinish()] -= pushNext;

        } else {
          it.getBackEdge().push(pushNext);
          pOut_[it.getBackEdge().getStart()] -= pushNext;
          pIn_[it.getBackEdge().getFinish()] -= pushNext;
        }

        if (type == 1) {
          if (toPush[it.getFinish()] == 0 &&
              it.getFinish() != myNetwork_.getSink())
            queue.push(it.getFinish());
        } else {
          if (toPush[it.getFinish()] == 0 &&
              it.getFinish() != myNetwork_.getSource())
            queue.push(it.getFinish());
        }

        toPush[vertex] -= pushNext;
        toPush[it.getFinish()] += pushNext;
      }
    }
  }

  TFlow findBlockingFlow_() {
    TVertex minPotentialVertex = getMinPotentialVertex_();

    if (minPotentialVertex == -1)
      return -1;

    if (getPotential_(minPotentialVertex) == 0) {
      deleteVertex_(minPotentialVertex);
      return 0;
    }

    TFlow toPush = getPotential_(minPotentialVertex);
    tryPush_(minPotentialVertex, toPush, 1);
    tryPush_(minPotentialVertex, toPush, -1);
    return toPush;
  }
};

template <typename TFlow = long long>
class Preflow {
 public:
  Preflow(Network<TFlow>& Net)
      : myNetwork_(Net),
        overFlow_(Net.getNumVertices(), 0),
        height_(Net.getNumVertices(), 0) {
    for (auto it = myNetwork_.begin(myNetwork_.getSource()); it.isValid();
         it.goNext()) {
      TFlow capacity = it.getCapacity();
      it.push(capacity);
      overFlow_[it.getStart()] -= capacity;
      overFlow_[it.getFinish()] += capacity;
    }

    height_[myNetwork_.getSource()] = myNetwork_.getNumVertices();
  }

  TFlow findMaxFlow() {
    std::list<TVertex> vertices;

    for (TVertex vertex = 0; vertex < myNetwork_.getNumVertices(); ++vertex) {
      if (vertex != myNetwork_.getSink() && vertex != myNetwork_.getSource())
        vertices.push_back(vertex);
    }

    auto itVertex = vertices.begin();

    while (itVertex != vertices.end()) {
      TVertex vertex = *itVertex;
      auto it = myNetwork_.begin(vertex);
      int32_t heightBefore = height_[vertex];
      discharge_(vertex, it);

      if (height_[vertex] > heightBefore) {
        vertices.erase(itVertex);
        vertices.push_front(vertex);
        itVertex = vertices.begin();
      }

      itVertex++;
    }

    return overFlow_[myNetwork_.getSink()];
  }

 private:
  Network<TFlow>& myNetwork_;
  std::vector<int32_t> height_;
  std::vector<TFlow> overFlow_;
  std::vector<typename Network<TFlow>::Iterator> edgePtr_;

  void push_(typename Network<TFlow>::Iterator it) {
    TFlow addFlow =
        std::min(it.getResidualCapacity(), overFlow_[it.getStart()]);
    it.push(addFlow);
    overFlow_[it.getFinish()] += addFlow;
    overFlow_[it.getStart()] -= addFlow;
  }

  void relabel_(TVertex vertex) {
    int32_t h = 2 * myNetwork_.getNumVertices();

    for (auto it = myNetwork_.begin(vertex); it.isValid(); it.goNext()) {
      if (it.getResidualCapacity() != 0)
        h = std::min(h, height_[it.getFinish()]);
    }
    height_[vertex] = h + 1;
  }

  void discharge_(TVertex vertex, typename Network<TFlow>::Iterator it) {
    while (overFlow_[vertex] > 0) {
      if (!it.isValid()) {
        relabel_(vertex);
        it = myNetwork_.begin(vertex);
      } else {
        if (it.getResidualCapacity() != 0 &&
            height_[vertex] == height_[it.getFinish()] + 1)
          push_(it);
        else
          it.goNext();
      }
    }
  }
};
