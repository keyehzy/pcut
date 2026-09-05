#include <pcut/sectors.hpp>
#include <algorithm>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace pcut {
namespace {
std::vector<State> sector_basis(const ClusterModel& model, unsigned max_charge, std::size_t max_basis) {
    if (max_charge > 1'000'000) throw std::invalid_argument("requested sector charge is too large");
    std::vector<State> result;
    std::vector<unsigned> state(model.sites());
    std::function<void(std::size_t,unsigned)> enumerate = [&](std::size_t site, unsigned charge) {
        if (site == model.sites()) {
            if (result.size() >= max_basis) throw std::length_error("sector basis budget exceeded");
            result.push_back(model.encode(state)); return;
        }
        const auto& space = model.spaces()[site];
        for (unsigned l=0; l<space.charges.size(); ++l) {
            const auto q=static_cast<unsigned>(space.charges[l]);
            if (q <= max_charge-charge) { state[site]=l; enumerate(site+1,charge+q); }
        }
    };
    enumerate(0,0);
    std::sort(result.begin(),result.end(),[&](State a, State b) {
        const int qa=model.charge(a), qb=model.charge(b);
        return qa==qb ? a<b : qa<qb;
    });
    return result;
}
std::vector<Excitation> excitations(const ClusterModel& model, State state, const std::vector<Site>& sites) {
    const auto local=model.decode(state);
    std::vector<Excitation> result;
    for (std::size_t i=0;i<local.size();++i) if (local[i]!=0) result.push_back({sites[i],local[i]});
    return result;
}
Kernel normalize_kernel(Kernel kernel) {
    Coordinate anchor;
    bool first=true;
    for (const auto* list : {&kernel.output,&kernel.input}) for (const auto& e : *list) {
        if (first || e.site.cell<anchor) { anchor=e.site.cell; first=false; }
    }
    for (auto* list : {&kernel.output,&kernel.input}) for (auto& e : *list)
        for (std::size_t d=0;d<anchor.size();++d) {
            const auto v=static_cast<long long>(e.site.cell[d])-anchor[d];
            if (v<std::numeric_limits<int>::min() || v>std::numeric_limits<int>::max())
                throw std::overflow_error("kernel coordinate overflow");
            e.site.cell[d]=static_cast<int>(v);
        }
    return kernel;
}
}
IrreducibleSectors irreducible_sectors(const ClusterModel& model, const EffectiveOperator& effective,
                                      unsigned max_charge, SectorOptions options) {
    model.require_product_vacuum();
    if (model.fermionic()) throw std::invalid_argument("tensor sector kernels do not support graded fermionic terms");
    IrreducibleSectors result;
    result.basis=sector_basis(model,max_charge,options.max_basis);
    result.kernels=effective.block(model,result.basis,options.solver);
    std::map<State,Eigen::Index> indices;
    std::vector<std::vector<unsigned>> decoded;
    for (std::size_t i=0;i<result.basis.size();++i) {
        indices.emplace(result.basis[i],static_cast<Eigen::Index>(i));
        decoded.push_back(model.decode(result.basis[i]));
    }
    std::vector<unsigned> local(model.sites(),0);
    std::size_t contractions=0;
    // Lower-charge input columns have already been converted to kernels.
    for (std::size_t j=0;j<result.basis.size();++j) for (std::size_t i=0;i<result.basis.size();++i) {
        if (model.charge(result.basis[i])!=model.charge(result.basis[j])) continue;
        std::vector<State> common;
        for (std::size_t s=0;s<model.sites();++s) if (decoded[i][s]!=0 && decoded[i][s]==decoded[j][s]) {
            local[s]=decoded[i][s]; common.push_back(model.encode(local)); local[s]=0;
        }
        std::function<void(std::size_t,State)> subtract = [&](std::size_t start,State removed) {
            for (std::size_t s=start;s<common.size();++s) {
                if (++contractions>options.max_contractions) throw std::length_error("spectator contraction budget exceeded");
                const State next=removed+common[s];
                const auto row=indices.at(result.basis[i]-next);
                const auto col=indices.at(result.basis[j]-next);
                for (auto& matrix : result.kernels)
                    matrix(static_cast<Eigen::Index>(i),static_cast<Eigen::Index>(j))-=matrix(row,col);
                subtract(s+1,next);
            }
        };
        subtract(0,0);
    }
    return result;
}
LinkedSectors linked_expand_sectors(const ClusterCatalog& catalog, const EffectiveOperator& effective,
                                    unsigned max_charge, SectorOptions options) {
    if (catalog.max_edges()<effective.order()) throw std::invalid_argument("catalog does not cover perturbation order");
    const auto& lattice=catalog.lattice();
    const auto order=effective.order();
    for (const auto& space : lattice.cell) space.require_product_vacuum();
    for (const auto& term : lattice.interactions) if (term.fermionic)
        throw std::invalid_argument("tensor sector kernels do not support graded fermionic terms");
    LinkedSectors result{Series(order+1),{}};
    for (std::size_t b=0;b<lattice.cell.size();++b) {
        const auto& space=lattice.cell[b];
        result.energy_per_cell[0]+=space.vacuum_energy;
        for (unsigned l=1;l<space.charges.size();++l) if (static_cast<unsigned>(space.charges[l])<=max_charge) {
            const Excitation e{{Coordinate(lattice.dimension,0),b},l};
            Series value(order+1); value[0]=lattice.gap*space.charges[l];
            result.kernels.emplace(Kernel{{e},{e}},std::move(value));
        }
    }
    if (result.kernels.size()>options.max_kernels) throw std::length_error("linked kernel budget exceeded");
    struct Weight { IrreducibleSectors sectors; std::vector<std::vector<unsigned>> local; };
    std::vector<Weight> weights;
    for (const auto& entry : catalog.entries()) {
        const auto model=cluster_model(lattice,entry.edges);
        Weight weight{irreducible_sectors(model,effective,max_charge,options),{}};
        weight.sectors.kernels[0].setZero(); // embed bare on-site terms separately
        std::map<State,Eigen::Index> indices;
        for (std::size_t i=0;i<weight.sectors.basis.size();++i) {
            indices.emplace(weight.sectors.basis[i],static_cast<Eigen::Index>(i));
            weight.local.push_back(model.decode(weight.sectors.basis[i]));
        }
        for (const auto& sub : entry.subclusters) {
            const auto& child=weights[sub.index];
            std::vector<Eigen::Index> map;
            for (const auto& local : child.local) {
                std::vector<unsigned> parent(model.sites(),0);
                for (std::size_t s=0;s<local.size();++s) parent[sub.vertex_map[s]]=local[s];
                map.push_back(indices.at(model.encode(parent)));
            }
            for (unsigned n=1;n<=order;++n)
                for (std::size_t i=0;i<map.size();++i) for (std::size_t j=0;j<map.size();++j)
                    weight.sectors.kernels[n](map[i],map[j])-=child.sectors.kernels[n](static_cast<Eigen::Index>(i),static_cast<Eigen::Index>(j));
        }
        for (unsigned n=1;n<=order;++n) result.energy_per_cell[n]+=weight.sectors.kernels[n](0,0);
        for (std::size_t i=0;i<weight.sectors.basis.size();++i) for (std::size_t j=0;j<weight.sectors.basis.size();++j) {
            if (i==0 && j==0) continue;
            Series value(order+1);
            bool nonzero=false;
            for (unsigned n=1;n<=order;++n) {
                value[n]=weight.sectors.kernels[n](static_cast<Eigen::Index>(i),static_cast<Eigen::Index>(j));
                nonzero=nonzero || value[n]!=Complex{};
            }
            if (!nonzero) continue;
            const auto key=normalize_kernel({excitations(model,weight.sectors.basis[i],entry.sites),excitations(model,weight.sectors.basis[j],entry.sites)});
            auto& target=result.kernels[key];
            if (target.empty()) target.resize(order+1);
            for (unsigned n=1;n<=order;++n) target[n]+=value[n];
            if (result.kernels.size()>options.max_kernels) throw std::length_error("linked kernel budget exceeded");
        }
        weights.push_back(std::move(weight));
    }
    return result;
}
}
