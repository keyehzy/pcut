#include <pcut/operator.hpp>
#include "detail.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace pcut {
namespace {
OperatorBlock empty_block(const ClusterModel& model, unsigned order,
                          std::optional<int> particles) {
    OperatorBlock block{model.spaces(),zero_charge_basis(model,particles),{}};
    block.coefficients=detail::matrix_series(block.basis.size(),static_cast<std::size_t>(order)+1);
    return block;
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
std::vector<State> zero_charge_basis(const ClusterModel& model, std::optional<int> particles) {
    if (particles && *particles<0) throw std::invalid_argument("invalid particle number");
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
    std::vector<State> result;
    std::vector<unsigned> local(model.sites());
    std::function<void(std::size_t,int)> visit=[&](std::size_t site,int number) {
        if (particles && (number+lower[site]>*particles || number+upper[site]<*particles)) return;
        if (site==model.sites()) {
            if (!particles || number==*particles) {
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
                                   std::optional<int> particles) {
    model.require_operator_grading();
    if (particles && !model.conserves_particles())
        throw std::invalid_argument("selected particle number is not conserved by the model");
    OperatorBlock result{model.spaces(),zero_charge_basis(model,particles),{}};
    result.coefficients=effective.block(model,result.basis);
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
namespace {
// Basis setup is independent of an embedding and of subsequent coefficient updates.
struct EmbeddingBasis {
    ClusterModel model;
    std::map<State,Eigen::Index> indices;
    std::vector<std::vector<unsigned>> decoded;
    bool graded=false;
    explicit EmbeddingBasis(const OperatorBlock& block) : model(block.spaces,{}) {
        validate_block(block,model);
        for (const auto& space : block.spaces) graded |= !space.parity.empty();
        if (graded) for (const auto& space : block.spaces) if (space.parity.empty())
            throw std::invalid_argument("graded embedding needs parity on every site");
        for (std::size_t i=0;i<block.basis.size();++i) {
            indices.emplace(block.basis[i],static_cast<Eigen::Index>(i));
            decoded.push_back(model.decode(block.basis[i]));
        }
    }
    void require_child(const OperatorBlock& child) const {
        std::size_t expected=1;
        for (const auto& space : child.spaces) {
            const auto zeros=static_cast<std::size_t>(std::count(space.charges.begin(),space.charges.end(),0));
            if (expected>child.basis.size()/zeros)
                throw std::invalid_argument("child operator must retain the full Q=0 basis");
            expected*=zeros;
        }
        if (expected!=child.basis.size() || !std::is_sorted(child.basis.begin(),child.basis.end()))
            throw std::invalid_argument("child operator must retain the full ordered Q=0 basis");
        if (!graded) return;
        std::vector<unsigned> parity(decoded.size());
        for (std::size_t i=0;i<decoded.size();++i)
            for (std::size_t s=0;s<child.spaces.size();++s) parity[i]^=child.spaces[s].parity[decoded[i][s]];
        for (std::size_t i=0;i<decoded.size();++i) for (std::size_t j=0;j<decoded.size();++j)
            if (parity[i]!=parity[j]) for (const auto& h : child.coefficients)
                if (h(static_cast<Eigen::Index>(i),static_cast<Eigen::Index>(j))!=Complex{})
                    throw std::invalid_argument("cannot embed a parity-odd Hamiltonian");
    }
};
struct EmbeddingPlan {
    const EmbeddingBasis& parent;
    const EmbeddingBasis& child;
    const std::vector<std::size_t>& map;
    std::vector<std::pair<std::size_t,std::size_t>> inversions;
    EmbeddingPlan(const EmbeddingBasis& p, const EmbeddingBasis& c, const std::vector<std::size_t>& m) : parent(p), child(c), map(m) {
        if (map.size()!=child.model.sites()) throw std::invalid_argument("operator embedding shape mismatch");
        inversions=detail::gather_inversions(parent.model.sites(),map);
        for (std::size_t s=0;s<map.size();++s) {
            const auto& a=child.model.spaces()[s]; const auto& b=parent.model.spaces()[map[s]];
            if (a.charges!=b.charges || a.particles!=b.particles || a.parity!=b.parity)
                throw std::invalid_argument("incompatible embedded local spaces");
        }
    }
    void add(OperatorBlock& target, const OperatorBlock& source, Complex scale) const {
        std::vector<unsigned> selected(map.size());
        for (std::size_t j=0;j<parent.decoded.size();++j) {
            auto local=parent.decoded[j];
            for (std::size_t s=0;s<map.size();++s) selected[s]=local[map[s]];
            const auto col=child.indices.at(child.model.encode(selected));
            const int input_sign=parent.graded ? detail::gather_sign(target.spaces,local,inversions) : 1;
            for (std::size_t i=0;i<child.decoded.size();++i) {
                bool nonzero=false;
                for (const auto& h : source.coefficients) nonzero |= h(static_cast<Eigen::Index>(i),col)!=Complex{};
                if (!nonzero) continue;
                for (std::size_t s=0;s<map.size();++s) local[map[s]]=child.decoded[i][s];
                const auto row=parent.indices.find(parent.model.encode(local));
                if (row==parent.indices.end()) continue;
                const int sign=input_sign*(parent.graded ? detail::gather_sign(target.spaces,local,inversions) : 1);
                for (std::size_t n=0;n<source.coefficients.size();++n)
                    target.coefficients[n](row->second,static_cast<Eigen::Index>(j))+=
                        scale*static_cast<double>(sign)*source.coefficients[n](static_cast<Eigen::Index>(i),col);
            }
        }
    }
};
}
void add_embedded_operator(OperatorBlock& parent, const OperatorBlock& child,
                           const std::vector<std::size_t>& map, Complex scale) {
    if (&parent==&child) throw std::invalid_argument("operator embedding must not alias");
    if (parent.coefficients.size()!=child.coefficients.size() ||
        !std::isfinite(scale.real()) || !std::isfinite(scale.imag()))
        throw std::invalid_argument("operator embedding shape or scale mismatch");
    const EmbeddingBasis p(parent), c(child);
    c.require_child(child);
    EmbeddingPlan(p,c,map).add(parent,child,scale);
}
LinkedOperator linked_zero_charge(const WhiteGraphExpansion& catalog, const EffectiveOperator& effective) {
    if (catalog.max_edges()<effective.order()) throw std::invalid_argument("catalog must cover operator order");
    LinkedOperator result{catalog.lattice(),effective.order(),{}};
    for (const auto& entry : catalog.embeddings()) {
        if (entry.edges.size()>effective.order()) break;
        const auto& model=catalog.structural_model(entry);
        model.require_operator_grading();
        OperatorBlock block{model.spaces(),zero_charge_basis(model),{}};
        block.coefficients=catalog.block(entry,effective,block.basis,true);
        block.coefficients[0].setZero(); // only bare reference constants at Q=0
        result.weights.push_back({entry.edges,entry.sites,std::move(block)});
    }
    return result;
}
OperatorBlock assemble_operator(const LinkedOperator& linked, const std::vector<Site>& sites,
                                const Cluster& edges, std::optional<int> particles) {
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
    auto result=empty_block(model,linked.order,particles);
    result.coefficients[0].diagonal().setConstant(model.vacuum_energy());
    const EmbeddingBasis parent(result);
    for (const auto& weight : linked.weights) {
        if (weight.edges.empty() || weight.sites!=vertices(linked.lattice,weight.edges) ||
            weight.block.coefficients.size()!=static_cast<std::size_t>(linked.order)+1)
            throw std::invalid_argument("invalid linked operator weight");
        const EmbeddingBasis child(weight.block);
        child.require_child(weight.block);
        // Fix one colored edge. Every matching target edge gives one translation,
        // so embeddings retain multiplicity without rotations or automorphism factors.
        const auto& anchor=weight.edges.front();
        for (const auto& target : edges) if (target.type==anchor.type) {
            const auto shift=detail::translate(target.origin,anchor.origin,-1);
            bool exists=true;
            for (const auto& edge : weight.edges)
                exists &= available.contains({edge.type,detail::translate(edge.origin,shift)});
            if (!exists) continue;
            std::vector<std::size_t> map;
            for (auto site : weight.sites) {
                site.cell=detail::translate(site.cell,shift); map.push_back(indices.at(site));
            }
            EmbeddingPlan(parent,child,map).add(result,weight.block,1);
        }
    }
    return result;
}
}
