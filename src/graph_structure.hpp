#pragma once
#include <pcut/white_graph.hpp>
#include <algorithm>
#include <bit>
#include <numeric>

namespace pcut::detail {
inline void number(std::string& key, std::uint64_t n) {
    char bytes[8];
    for (unsigned i=0;i<8;++i) bytes[i]=static_cast<char>((n>>(8*i))&255);
    key.append(bytes,8);
}
inline void real(std::string& key,double x) { number(key,std::bit_cast<std::uint64_t>(x==0 ? 0.0 : x)); }
inline void string(std::string& key,const std::string& s) { number(key,s.size()); key+=s; }
template<class T> void sequence(std::string& key,const std::vector<T>& s) {
    number(key,s.size()); for (auto x : s) number(key,static_cast<std::uint64_t>(x));
}
inline std::string space_key(const LocalSpace& space) {
    std::string key;
    string(key,space.name); sequence(key,space.charges); sequence(key,space.particles);
    sequence(key,space.parity); real(key,space.vacuum_energy); return key;
}
inline std::string channel_key(const OperatorChannel& op) {
    std::string key;
    number(key,op.fermionic); number(key,static_cast<std::uint64_t>(op.matrix.rows()));
    number(key,static_cast<std::uint64_t>(op.matrix.cols()));
    for (Eigen::Index i=0;i<op.matrix.size();++i) { real(key,op.matrix.data()[i].real()); real(key,op.matrix.data()[i].imag()); }
    return key;
}
inline std::vector<std::size_t> channel_order(const GraphEdge& edge) {
    std::vector<std::size_t> order(edge.channels.size()); std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](auto a,auto b) { return channel_key(edge.channels[a])<channel_key(edge.channels[b]); });
    return order;
}
inline std::string edge_structure(const GraphEdge& edge) {
    std::string key; number(key,edge.legs.size()); number(key,edge.channels.size());
    for (auto c : channel_order(edge)) string(key,channel_key(edge.channels[c]));
    return key;
}
} // namespace pcut::detail
