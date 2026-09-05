#include <pcut/operator.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace pcut {
namespace {
std::size_t matrix_size(std::size_t basis, std::size_t orders, std::size_t budget) {
    if (orders == 0 || (basis && (basis > budget / basis || orders > budget / basis / basis)))
        throw std::length_error("operator matrix-element budget exceeded");
    if (basis > static_cast<std::size_t>(std::numeric_limits<Eigen::Index>::max()))
        throw std::length_error("operator basis exceeds Eigen indexing");
    return basis * basis * orders;
}
OperatorBlock empty_block(const ClusterModel& model, unsigned order,
                          std::optional<int> particles, OperatorOptions options) {
    OperatorBlock block{model.spaces(),zero_charge_basis(model,particles,options.max_basis,options.max_basis_visits),{}};
    matrix_size(block.basis.size(),static_cast<std::size_t>(order)+1,options.max_matrix_elements);
    const auto size=static_cast<Eigen::Index>(block.basis.size());
    block.coefficients.assign(static_cast<std::size_t>(order)+1,Matrix::Zero(size,size));
    return block;
}
Coordinate translate(const Coordinate& a, const Coordinate& b, int sign=1) {
    if (a.size()!=b.size()) throw std::invalid_argument("operator coordinate dimension mismatch");
    Coordinate out(a.size());
    for (std::size_t d=0; d<a.size(); ++d) {
        const auto x=static_cast<long long>(a[d])+static_cast<long long>(sign)*b[d];
        if (x<std::numeric_limits<int>::min() || x>std::numeric_limits<int>::max())
            throw std::overflow_error("operator coordinate overflow");
        out[d]=static_cast<int>(x);
    }
    return out;
}
void validate_block(const OperatorBlock& block, const ClusterModel& model) {
    std::set<State> unique;
    for (auto state : block.basis)
        if (model.charge(state)!=0 || !unique.insert(state).second)
            throw std::invalid_argument("operator basis must contain unique charge-zero states");
    const auto size=static_cast<Eigen::Index>(block.basis.size());
    if (block.coefficients.empty()) throw std::invalid_argument("empty operator series");
    for (const auto& matrix : block.coefficients)
        if (matrix.rows()!=size || matrix.cols()!=size || !matrix.allFinite())
            throw std::invalid_argument("invalid operator matrix");
}
}
std::vector<State> zero_charge_basis(const ClusterModel& model, std::optional<int> particles,
                                     std::size_t max_basis, std::size_t max_visits) {
    if (!max_basis || (particles && *particles<0)) throw std::invalid_argument("invalid external basis limit or particle number");
    if (particles) for (const auto& space : model.spaces()) if (space.particles.empty())
        throw std::invalid_argument("particle selection requires metadata");
    std::vector<int> lower(model.sites()+1), upper(model.sites()+1);
    if (particles) for (std::size_t s=model.sites(); s>0; --s) {
        const auto& space=model.spaces()[s-1];
        int lo=1'000'000, hi=0;
        for (std::size_t l=0;l<space.charges.size();++l) if (space.charges[l]==0) {
            lo=std::min(lo,space.particles[l]); hi=std::max(hi,space.particles[l]);
        }
        lower[s-1]=lower[s]+lo; upper[s-1]=upper[s]+hi;
    }
    std::size_t visits=0;
    std::vector<State> result;
    std::vector<unsigned> local(model.sites());
    std::function<void(std::size_t,int)> visit=[&](std::size_t site,int number) {
        if (++visits>max_visits) throw std::length_error("external basis traversal budget exceeded");
        if (particles && (number+lower[site]>*particles || number+upper[site]<*particles)) return;
        if (site==model.sites()) {
            if (!particles || number==*particles) {
                if (result.size()>=max_basis) throw std::length_error("charge-zero basis budget exceeded");
                result.push_back(model.encode(local));
            }
            return;
        }
        const auto& space=model.spaces()[site];
        for (unsigned l=0; l<space.charges.size(); ++l) if (space.charges[l]==0) {
            const int next=number+(particles ? space.particles[l] : 0);
            if (particles && next>*particles) continue;
            local[site]=l; visit(site+1,next);
        }
    };
    visit(0,0);
    std::sort(result.begin(),result.end());
    return result;
}
OperatorBlock zero_charge_operator(const ClusterModel& model, const EffectiveOperator& effective,
                                   std::optional<int> particles, OperatorOptions options) {
    model.require_operator_grading();
    if (particles && !model.conserves_particles())
        throw std::invalid_argument("selected particle number is not conserved by the model");
    auto result=empty_block(model,effective.order(),particles,options);
    result.coefficients=effective.block(model,result.basis,options.solver);
    return result;
}
Matrix OperatorBlock::evaluate(double lambda) const {
    if (!std::isfinite(lambda)) throw std::invalid_argument("nonfinite expansion parameter");
    const ClusterModel model(spaces,{});
    validate_block(*this,model);
    Matrix result=Matrix::Zero(static_cast<Eigen::Index>(basis.size()),static_cast<Eigen::Index>(basis.size()));
    for (auto it=coefficients.rbegin(); it!=coefficients.rend(); ++it) result=lambda*result+*it;
    return result;
}
void add_embedded_operator(OperatorBlock& parent, const OperatorBlock& child,
                           const std::vector<std::size_t>& map, Complex scale, std::size_t max_work) {
    if (&parent==&child) throw std::invalid_argument("operator embedding must not alias");
    if (map.size()!=child.spaces.size() || parent.coefficients.size()!=child.coefficients.size() ||
        !std::isfinite(scale.real()) || !std::isfinite(scale.imag()))
        throw std::invalid_argument("operator embedding shape or scale mismatch");
    const ClusterModel pm(parent.spaces,{}), cm(child.spaces,{});
    validate_block(parent,pm); validate_block(child,cm);
    std::set<std::size_t> used;
    bool graded=false;
    for (const auto& space : parent.spaces) graded |= !space.parity.empty();
    for (std::size_t s=0; s<map.size(); ++s) {
        if (map[s]>=parent.spaces.size() || !used.insert(map[s]).second)
            throw std::invalid_argument("invalid operator vertex map");
        const auto& a=child.spaces[s]; const auto& b=parent.spaces[map[s]];
        if (a.charges!=b.charges || a.particles!=b.particles || a.parity!=b.parity)
            throw std::invalid_argument("incompatible embedded local spaces");
    }
    if (graded) for (const auto& space : parent.spaces) if (space.parity.empty())
        throw std::invalid_argument("graded embedding needs parity on every site");
    std::size_t expected_basis=1;
    for (const auto& space : child.spaces) {
        const auto zeros=static_cast<std::size_t>(std::count(space.charges.begin(),space.charges.end(),0));
        if (expected_basis>child.basis.size()/zeros)
            throw std::invalid_argument("child operator must retain the full Q=0 basis");
        expected_basis*=zeros;
    }
    // validate_block already checked uniqueness and Q=0 membership, so matching
    // the tensor dimension proves completeness without re-enumerating states.
    if (expected_basis!=child.basis.size() || !std::is_sorted(child.basis.begin(),child.basis.end()))
        throw std::invalid_argument("child operator must retain the full ordered Q=0 basis");
    const auto nc=child.basis.size(), np=parent.basis.size(), orders=child.coefficients.size();
    if (nc && np && (np>max_work/nc || orders>max_work/nc/np))
        throw std::length_error("operator embedding work budget exceeded");
    std::map<State,Eigen::Index> child_index, parent_index;
    std::vector<std::vector<unsigned>> decoded;
    std::vector<unsigned> child_parity;
    for (std::size_t i=0; i<nc; ++i) {
        child_index.emplace(child.basis[i],static_cast<Eigen::Index>(i));
        decoded.push_back(cm.decode(child.basis[i]));
        unsigned p=0;
        if (graded) for (std::size_t s=0; s<map.size(); ++s) p^=child.spaces[s].parity[decoded.back()[s]];
        child_parity.push_back(p);
    }
    // Only even operators admit this identity-spectator embedding rule.
    if (graded) for (std::size_t i=0; i<nc; ++i) for (std::size_t j=0; j<nc; ++j)
        if (child_parity[i]!=child_parity[j]) for (const auto& h : child.coefficients)
            if (h(static_cast<Eigen::Index>(i),static_cast<Eigen::Index>(j))!=Complex{})
                throw std::invalid_argument("cannot embed a parity-odd Hamiltonian");
    for (std::size_t i=0; i<np; ++i) parent_index.emplace(parent.basis[i],static_cast<Eigen::Index>(i));
    for (std::size_t j=0; j<np; ++j) {
        auto local=pm.decode(parent.basis[j]);
        std::vector<unsigned> selected(map.size());
        for (std::size_t s=0; s<map.size(); ++s) selected[s]=local[map[s]];
        const auto col=child_index.at(cm.encode(selected));
        const int input_sign=graded ? graded_gather_sign(parent.spaces,local,map) : 1;
        for (std::size_t i=0; i<nc; ++i) {
            bool nonzero=false;
            for (const auto& h : child.coefficients) nonzero |= h(static_cast<Eigen::Index>(i),col)!=Complex{};
            if (!nonzero) continue;
            for (std::size_t s=0; s<map.size(); ++s) local[map[s]]=decoded[i][s];
            const auto row=parent_index.find(pm.encode(local));
            if (row==parent_index.end()) continue;
            const int sign=input_sign*(graded ? graded_gather_sign(parent.spaces,local,map) : 1);
            for (std::size_t n=0; n<orders; ++n)
                parent.coefficients[n](row->second,static_cast<Eigen::Index>(j))+=
                    scale*static_cast<double>(sign)*child.coefficients[n](static_cast<Eigen::Index>(i),col);
        }
    }
}
LinkedOperator linked_zero_charge(const ClusterCatalog& catalog, const EffectiveOperator& effective,
                                  OperatorOptions options) {
    if (catalog.max_edges()<effective.order()) throw std::invalid_argument("catalog must cover operator order");
    LinkedOperator result{catalog.lattice(),effective.order(),{}};
    std::size_t stored=0;
    for (const auto& entry : catalog.entries()) {
        const auto model=cluster_model(result.lattice,entry.edges);
        auto remaining=options; remaining.max_matrix_elements-=stored;
        auto block=zero_charge_operator(model,effective,{},remaining);
        stored+=matrix_size(block.basis.size(),block.coefficients.size(),remaining.max_matrix_elements);
        block.coefficients[0].setZero(); // only bare reference constants at Q=0
        for (const auto& sub : entry.subclusters)
            add_embedded_operator(block,result.weights.at(sub.index).block,sub.vertex_map,-1,options.max_embedding_work);
        result.weights.push_back({entry.edges,entry.sites,std::move(block)});
    }
    return result;
}
OperatorBlock assemble_operator(const LinkedOperator& linked, const std::vector<Site>& sites,
                                const Cluster& edges, std::optional<int> particles, OperatorOptions options) {
    linked.lattice.validate();
    std::map<Site,std::size_t> indices;
    std::vector<LocalSpace> spaces;
    for (const auto& site : sites) {
        if (site.cell.size()!=linked.lattice.dimension || site.basis>=linked.lattice.cell.size() ||
            !indices.emplace(site,spaces.size()).second) throw std::invalid_argument("invalid assembly site");
        spaces.push_back(linked.lattice.cell[site.basis]);
    }
    if (particles) for (std::size_t type=0;type<linked.lattice.interactions.size();++type)
        if (!cluster_model(linked.lattice,{{type,Coordinate(linked.lattice.dimension,0)}}).conserves_particles())
            throw std::invalid_argument("selected particle number is not conserved by the lattice");
    const std::set<Edge> available(edges.begin(),edges.end());
    if (available.size()!=edges.size()) throw std::invalid_argument("duplicate assembly edge");
    for (const auto& site : vertices(linked.lattice,edges)) if (!indices.contains(site))
        throw std::invalid_argument("assembly edge has missing vertex");
    const ClusterModel model(spaces,{},linked.lattice.gap);
    auto result=empty_block(model,linked.order,particles,options);
    result.coefficients[0].diagonal().setConstant(model.vacuum_energy());
    for (const auto& weight : linked.weights) {
        if (weight.edges.empty() || weight.sites!=vertices(linked.lattice,weight.edges) ||
            weight.block.coefficients.size()!=static_cast<std::size_t>(linked.order)+1)
            throw std::invalid_argument("invalid linked operator weight");
        // Fix one colored edge. Every matching target edge gives one translation,
        // so embeddings retain multiplicity without rotations or automorphism factors.
        const auto& anchor=weight.edges.front();
        for (const auto& target : edges) if (target.type==anchor.type) {
            const auto shift=translate(target.origin,anchor.origin,-1);
            bool exists=true;
            for (const auto& edge : weight.edges)
                exists &= available.contains({edge.type,translate(edge.origin,shift)});
            if (!exists) continue;
            std::vector<std::size_t> map;
            for (auto site : weight.sites) {
                site.cell=translate(site.cell,shift); map.push_back(indices.at(site));
            }
            add_embedded_operator(result,weight.block,map,1,options.max_embedding_work);
        }
    }
    return result;
}
}
