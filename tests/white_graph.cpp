#include <pcut/pcut.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <algorithm>
#include <limits>
#include <numeric>
#include <set>

using namespace pcut;
namespace {
Coordinate shift(Coordinate a,const Coordinate& b,int sign=1) {
    for (std::size_t d=0;d<a.size();++d) a[d]+=sign*b[d];
    return a;
}
// Independent colored edge-animal enumeration. It knows nothing about abstract
// graph identity, canonical labels, graph mappings or the symbolic evaluator.
std::vector<Cluster> colored_animals(const PeriodicLattice& lattice,unsigned order) {
    std::set<Cluster> level;
    if (order) for (std::size_t t=0;t<lattice.interactions.size();++t)
        level.insert({{t,Coordinate(lattice.dimension)}});
    std::vector<Cluster> result;
    for (unsigned n=1;n<=order;++n) {
        result.insert(result.end(),level.begin(),level.end());
        std::set<Cluster> next;
        if (n<order) for (const auto& c : level) for (const auto& site : vertices(lattice,c))
            for (std::size_t t=0;t<lattice.interactions.size();++t) for (const auto& leg : lattice.interactions[t].legs)
                if (site.basis==leg.basis) {
                    const Edge e{t,shift(site.cell,leg.cell,-1)};
                    if (std::find(c.begin(),c.end(),e)!=c.end()) continue;
                    auto grown=c; grown.push_back(e); next.insert(normalize(grown).cluster);
                }
        level=std::move(next);
    }
    return result;
}
Series colored_vacuum(const PeriodicLattice& lattice,const EffectiveOperator& effective) {
    std::map<Cluster,Series> weights;
    Series result(effective.order()+1);
    for (const auto& space : lattice.cell) result[0]+=space.vacuum_energy;
    for (const auto& c : colored_animals(lattice,effective.order())) {
        auto weight=effective.vacuum(cluster_model(lattice,c)); weight[0]=0;
        // Independently enumerate all proper subsets, including repeated
        // normalized children at different positions in the parent.
        for (unsigned bits=1;bits+1<(1u<<c.size());++bits) {
            Cluster sub;
            for (std::size_t e=0;e<c.size();++e) if (bits&(1u<<e)) sub.push_back(c[e]);
            if (!connected(lattice,sub)) continue;
            const auto& w=weights.at(normalize(sub).cluster);
            for (unsigned n=1;n<=effective.order();++n) weight[n]-=w[n];
        }
        for (unsigned n=1;n<=effective.order();++n) result[n]+=weight[n];
        weights.emplace(c,std::move(weight));
    }
    return result;
}
PeriodicLattice square() {
    auto lattice=models::ising_chain(); lattice.dimension=2;
    const auto channels=lattice.interactions[0].channels;
    lattice.interactions={{{{{0,0},0},{{1,0},0}},channels},{{{{0,0},0},{{0,1},0}},channels}};
    lattice.interactions[0].channels[0].coupling=0.7;
    lattice.interactions[1].channels[0].coupling=-0.2;
    return lattice;
}
WhiteGraph graph(std::size_t n,std::vector<std::vector<std::size_t>> legs) {
    WhiteGraph g; g.spaces.resize(n,models::ising_chain().cell[0]);
    for (auto& support : legs) {
        Matrix op=Matrix::Identity(1u<<support.size(),1u<<support.size());
        g.edges.push_back({std::move(support),{{op,false}}});
    }
    return g;
}
WhiteGraph permute(const WhiteGraph& g,const std::vector<std::size_t>& permutation) {
    auto out=g;
    for (std::size_t v=0;v<g.spaces.size();++v) out.spaces[permutation[v]]=g.spaces[v];
    std::reverse(out.edges.begin(),out.edges.end());
    for (auto& e : out.edges) {
        for (auto& v : e.legs) v=permutation[v];
        std::reverse(e.channels.begin(),e.channels.end());
    }
    return out;
}
void compare_blocks(const WhiteGraphExpansion& expansion,const EffectiveOperator& effective,bool zero=false) {
    for (const auto& e : expansion.embeddings()) {
        const auto model=cluster_model(expansion.lattice(),e.edges);
        std::vector<State> basis;
        if (zero) basis=zero_charge_basis(model);
        else for (State s=0;s<model.dimension();++s) basis.push_back(s);
        const auto expected=effective.block(model,basis), actual=expansion.block(e,effective,basis);
        for (std::size_t n=0;n<actual.size();++n) {
            INFO("order " << n << ", edges " << e.edges.size());
            REQUIRE((actual[n]-expected[n]).norm()<1e-10);
        }
    }
}
}
TEST_CASE("Canonical labels and maps cover arbitrary permutations and automorphisms", "[white]") {
    const auto star=graph(5,{{0,1},{0,2},{0,3},{0,4}});
    const auto a=canonicalize(star);
    REQUIRE(a.vertex_automorphisms==24); // not just lattice rotations/reflections
    auto permutation=std::vector<std::size_t>{0,1,2,3,4};
    do {
        const auto relabeled=permute(star,permutation);
        const auto b=canonicalize(relabeled);
        REQUIRE(a.key==b.key);
        for (std::size_t e=0;e<b.graph.edges.size();++e)
            for (std::size_t leg=0;leg<b.graph.edges[e].legs.size();++leg)
                REQUIRE(b.map.vertices[b.graph.edges[e].legs[leg]]==relabeled.edges[b.map.edges[e]].legs[leg]);
    } while (std::next_permutation(permutation.begin(),permutation.end()));
    REQUIRE(a.key!=canonicalize(graph(5,{{0,1},{1,2},{2,3},{3,4}})).key);
    REQUIRE(canonicalize(graph(4,{{0,1},{1,2},{2,3},{3,0}})).vertex_automorphisms==4);
    REQUIRE(canonicalize(graph(2,{{0,1},{0,1}})).graph.edges.size()==2);
    auto species=star; species.spaces[1].name="different species";
    REQUIRE(canonicalize(species).key!=a.key);
    auto spectrum=star; spectrum.spaces[1].charges={0,2};
    REQUIRE(canonicalize(spectrum).key!=a.key);
    auto op=star; op.edges[0].channels[0].matrix(0,0)=2;
    REQUIRE(canonicalize(op).key!=a.key);
    auto invalid=star; invalid.edges[0].legs={0,0};
    REQUIRE_THROWS_AS(canonicalize(invalid),std::invalid_argument);
}
TEST_CASE("Straight and bent paths share white evaluations with exact infinite embeddings", "[white]") {
    auto lattice=square();
    const EffectiveOperator effective(Coefficients(charge_changes(lattice),3));
    const auto cache=std::make_shared<GraphCache>();
    (void)linked_expand(WhiteGraphExpansion(models::ising_chain(),3,cache),effective);
    const WhiteGraphExpansion expansion(lattice,3,cache);
    const auto colored=colored_animals(lattice,3);
    REQUIRE(expansion.embeddings().size()==colored.size());
    for (std::size_t i=0;i<colored.size();++i) REQUIRE(expansion.embeddings()[i].edges==colored[i]);
    const Cluster straight{{0,{0,0}},{0,{1,0}}}, bent{{0,{0,0}},{1,{1,0}}};
    auto id=[&](const Cluster& edges) {
        return std::find_if(expansion.embeddings().begin(),expansion.embeddings().end(),[&](const auto& e) { return e.edges==edges; })->graph;
    };
    REQUIRE(id(straight)==id(bent));
    const auto actual=linked_expand(expansion,effective);
    const auto expected=colored_vacuum(lattice,effective);
    for (unsigned n=0;n<=3;++n) REQUIRE(std::abs(actual.energy_per_cell[n]-expected[n])<1e-12);
    REQUIRE(expansion.cache()->evaluations()==expansion.graphs().size());
    REQUIRE(expansion.cache()->evaluations()<colored.size());
    const auto count=expansion.cache()->evaluations();
    lattice.interactions[0].channels[0].coupling=0;
    lattice.interactions[1].channels[0].coupling=1.3;
    const auto sweep=linked_expand(expansion.bind(lattice.couplings()),effective);
    const auto independent=colored_vacuum(lattice,effective);
    for (unsigned n=0;n<=3;++n) REQUIRE(std::abs(sweep.energy_per_cell[n]-independent[n])<1e-12);
    REQUIRE(expansion.cache()->evaluations()==count);
    REQUIRE(expansion.cache()->hits()>0);
}
TEST_CASE("Channel monomials and mapped subcluster subtraction preserve multiplicity", "[white]") {
    auto lattice=square();
    auto second=lattice.interactions[0].channels[0];
    second.op.matrix=Matrix::Identity(4,4); second.coupling=0.31;
    for (auto& interaction : lattice.interactions) interaction.channels.push_back(second);
    // Reorder channels in one color; variable maps must follow structure.
    std::reverse(lattice.interactions[1].channels.begin(),lattice.interactions[1].channels.end());
    const WhiteGraphExpansion expansion(lattice,2);
    const EffectiveOperator effective(Coefficients(charge_changes(lattice),2));
    compare_blocks(expansion,effective);
    const ScalarEvaluator scalar([](const WhiteGraph& g,double,unsigned n) {
        SymbolicSeries result(n+1);
        for (std::size_t v=0;v<g.variables();++v) result[1][Monomial{}.multiplied(v)]=1;
        return result;
    });
    const auto sum=linked_scalar(expansion,2,scalar);
    REQUIRE(sum.per_cell[1].real()==Catch::Approx(0.7-0.2+2*0.31));
    for (std::size_t i=0;i<expansion.embeddings().size();++i) if (expansion.embeddings()[i].edges.size()>1)
        REQUIRE(std::abs(sum.weights[i][1])<1e-14);
    for (const auto& entry : expansion.graphs()) for (const auto& sub : entry.subclusters) {
        const auto& child=expansion.graphs()[sub.graph].canonical.graph;
        REQUIRE(sub.map.channels.size()==child.variables());
        REQUIRE(sub.map.edges.size()==child.edges.size());
        for (std::size_t e=0;e<child.edges.size();++e) for (std::size_t l=0;l<child.edges[e].legs.size();++l)
            REQUIRE(sub.map.vertices[child.edges[e].legs[l]]==entry.canonical.graph.edges[sub.map.edges[e]].legs[l]);
    }
    const ScalarEvaluator vacuum([effective](const WhiteGraph& g,double gap,unsigned order) {
        const auto block=effective.symbolic_block(g.model(gap),{0});
        SymbolicSeries series(order+1);
        for (unsigned n=1;n<=order;++n) {
            const auto it=block.coefficients[n].find({0,0});
            if (it!=block.coefficients[n].end()) series[n]=it->second;
        }
        return series;
    });
    const auto linked=linked_scalar(expansion,2,vacuum);
    const auto independent=colored_vacuum(lattice,effective);
    for (unsigned n=0;n<=2;++n) REQUIRE(std::abs(linked.per_cell[n]-independent[n])<1e-12);
    const auto count=expansion.cache()->evaluations();
    (void)linked_scalar(expansion,2,vacuum);
    REQUIRE(expansion.cache()->evaluations()==count);
    Polynomial p{{Monomial{}.multiplied(0).multiplied(1).multiplied(1),2}};
    Polynomial relabeled; const std::vector<std::size_t> map{2,0};
    add_polynomial(relabeled,p,1,&map);
    REQUIRE(substitute(p,{3,4})==substitute(relabeled,{4,0,3}));
    REQUIRE_THROWS_AS(substitute(p,{1}),std::out_of_range);
    REQUIRE_THROWS_AS(Monomial({{0,65}}).degree(),std::invalid_argument);
}
TEST_CASE("Ordered hyperedges parallel templates and species match colored calculations", "[white]") {
    LocalSpace a{{0,1},-0.1,"a"}, b{{0,1,2},0.2,"b"};
    Matrix v=Matrix::Zero(12,12); v(1,4)=Complex(0,0.3); v(4,1)=std::conj(v(1,4)); v(0,11)=v(11,0)=0.2;
    PeriodicLattice lattice{1,{a,b},{},1.7};
    lattice.interactions.push_back({{{{1},0},{{0},1},{{0},0}},{{{v,false},0.7}}});
    lattice.interactions.push_back(lattice.interactions[0]); // intentionally double a template
    lattice.interactions[1].channels[0].coupling=-0.2;
    const WhiteGraphExpansion expansion(lattice,2);
    const EffectiveOperator effective(Coefficients(charge_changes(lattice),2));
    REQUIRE(expansion.embeddings().size()==colored_animals(lattice,2).size());
    compare_blocks(expansion,effective);
    const auto result=linked_expand(expansion,effective,{false});
    const auto expected=colored_vacuum(lattice,effective);
    for (unsigned n=0;n<=2;++n) REQUIRE(std::abs(result.energy_per_cell[n]-expected[n])<1e-12);
    auto combined=lattice; combined.interactions.resize(1); combined.interactions[0].channels[0].coupling=0.5;
    const auto one=linked_expand(WhiteGraphExpansion(combined,2),effective,{false});
    for (unsigned n=0;n<=2;++n) REQUIRE(std::abs(result.energy_per_cell[n]-one.energy_per_cell[n])<1e-12);
    auto reversed=lattice; std::swap(reversed.interactions[0].legs[0],reversed.interactions[0].legs[2]);
    const WhiteGraphExpansion changed(reversed,1);
    // Exchanging equivalent-species legs of an isolated hyperedge is a graph
    // relabeling; the physical embedding must still retain their ordered roles.
    REQUIRE(changed.graphs().size()==1);
    const EffectiveOperator first(Coefficients(charge_changes(reversed),1));
    compare_blocks(changed,first);
}
TEST_CASE("Graded canonical state permutations reproduce Hubbard graph blocks", "[white][fermion]") {
    auto lattice=models::hubbard_square();
    lattice.interactions[0].channels[0].coupling=0.4;
    lattice.interactions[1].channels[0].coupling=-0.7;
    const EffectiveOperator effective(Coefficients(charge_changes(lattice),2));
    const WhiteGraphExpansion expansion(lattice,2);
    compare_blocks(expansion,effective,true);
    REQUIRE(expansion.cache()->evaluations()==expansion.graphs().size());
    const auto count=expansion.cache()->evaluations();
    (void)linked_zero_charge(expansion,effective);
    REQUIRE(expansion.cache()->evaluations()==2*count);
    lattice.interactions[0].channels[0].coupling=0;
    (void)linked_zero_charge(expansion.bind(lattice.couplings()),effective);
    REQUIRE(expansion.cache()->evaluations()==2*count);
}
TEST_CASE("Cache context includes exact program sectors operators and local metadata", "[white]") {
    auto g=canonicalize(graph(2,{{0,1}}));
    GraphCache cache;
    const EffectiveOperator first(Coefficients({0},1)), second(Coefficients({0},2));
    const auto a=cache.block(g,1,first,{0});
    REQUIRE(cache.block(g,1,first,{0})==a);
    REQUIRE(cache.evaluations()==1);
    (void)cache.block(g,2,first,{0});
    (void)cache.block(g,1,second,{0});
    (void)cache.block(g,1,first,{0,1});
    auto changed=g.graph; changed.spaces[0].vacuum_energy=0.1;
    (void)cache.block(canonicalize(changed),1,first,{0});
    changed=g.graph; changed.edges[0].channels[0].matrix(0,0)=2;
    (void)cache.block(canonicalize(changed),1,first,{0});
    changed=g.graph; for (auto& space : changed.spaces) space.parity={0,1};
    changed.edges[0].channels[0].fermionic=true;
    (void)cache.block(canonicalize(changed),1,first,{0});
    REQUIRE(cache.evaluations()==7);
    auto modified=g; modified.graph.spaces[0].charges={0,3};
    (void)cache.block(modified,1,first,{0}); // actual structure, not a stale supplied key
    const EffectiveOperator alphabet(Coefficients({-1,0,1},1));
    (void)cache.block(g,1,alphabet,{0});
    REQUIRE(cache.evaluations()==9);
    REQUIRE(a->coefficients.size()==2);
    REQUIRE_THROWS_AS(cache.block(g,1,first,{0,0}),std::invalid_argument);
    const ScalarEvaluator invalid([](const WhiteGraph&,double,unsigned n) { SymbolicSeries s(n+1); s[1][Monomial{}]=1; return s; });
    REQUIRE_THROWS_AS(cache.scalar(g,1,1,invalid),std::invalid_argument);
}

TEST_CASE("Symbolic evaluation keeps energy units and all reachable intermediate charges", "[white]") {
    // External Q=0 reaches Q=7. This is not an external-sector truncation.
    LocalSpace space{{0,7},-0.3,"ladder"};
    Matrix v(2,2); v << 0,0.4,0.4,0;
    const EffectiveOperator effective(Coefficients({-7,7},4));
    for (double scale : {1e-200,1.0,1e200}) {
        auto local=space; local.vacuum_energy*=scale;
        const ClusterModel model({local},{{{0},scale*v}},scale);
        const auto symbolic=effective.symbolic_block(model,{0});
        const auto numeric=effective.vacuum(model);
        for (unsigned n=0;n<=4;++n) {
            const auto it=symbolic.coefficients[n].find({0,0});
            const auto value=it==symbolic.coefficients[n].end() ? Complex{} : substitute(it->second,{1});
            REQUIRE(std::abs((value-numeric[n])/scale)<1e-12);
        }
    }
}

TEST_CASE("Compact monomials preserve wide indices large powers and immutable copies", "[white]") {
    const auto small=Monomial{}.multiplied(3).multiplied(3);
    auto large=small;
    const auto wide=std::numeric_limits<std::size_t>::max();
    large=large.multiplied(wide);
    REQUIRE(small.degree()==2);
    REQUIRE(large.degree()==3);
    std::size_t largest=0;
    large.for_each_power([&](std::size_t v,unsigned) { largest=std::max(largest,v); });
    REQUIRE(largest==wide);
    auto power=Monomial{};
    for (unsigned n=0;n<64;++n) power=power.multiplied(0);
    REQUIRE(power.degree()==64);
    REQUIRE_THROWS_AS(power.multiplied(0),std::length_error);
    REQUIRE_THROWS_AS(Monomial({{0,0}}),std::invalid_argument);
    REQUIRE_THROWS_AS(Monomial({{1,1},{1,2}}),std::invalid_argument);
    const auto merged=small.relabel({0,0,0,20});
    REQUIRE(merged==Monomial({{20,2}}));
    Polynomial p{{power,2}};
    REQUIRE(substitute(p,{1})==Complex{2});
}

TEST_CASE("Factored symbolic program retains a wide complex channel alphabet", "[white]") {
    const LocalSpace space{{0,1},-0.2,"wide"};
    std::vector<LocalTerm> symbolic,numerical;
    std::vector<double> values;
    for (unsigned c=0;c<40;++c) {
        Matrix matrix=Matrix::Zero(2,2);
        matrix(0,1)=Complex(0.1+0.002*c,0.03);
        matrix(1,0)=std::conj(matrix(0,1));
        matrix(1,1)=0.001*c;
        const double ratio=c%3==0 ? 0 : (c%2==0 ? -0.07 : 0.13);
        symbolic.push_back({{0},matrix,false});
        numerical.push_back({{0},ratio*matrix,false});
        values.push_back(ratio);
    }
    const EffectiveOperator effective(Coefficients({-1,0,1},2));
    const auto actual=effective.symbolic_block(ClusterModel({space},symbolic,1.7),{0,1});
    const auto expected=effective.block(ClusterModel({space},numerical,1.7),{0,1});
    for (unsigned n=0;n<=2;++n) for (State i=0;i<2;++i) for (State j=0;j<2;++j) {
        const auto it=actual.coefficients[n].find({i,j});
        const auto value=it==actual.coefficients[n].end() ? Complex{} : substitute(it->second,values);
        REQUIRE(std::abs(value-expected[n](static_cast<Eigen::Index>(i),static_cast<Eigen::Index>(j)))<1e-12);
    }
}

TEST_CASE("Coupling bindings share plans and reject malformed ratios", "[white]") {
    const WhiteGraphExpansion original(square(),2);
    const auto bound=original.bind({{0},{-0.8}});
    REQUIRE(&bound.graphs()==&original.graphs());
    REQUIRE(&bound.embeddings()==&original.embeddings());
    REQUIRE(original.lattice().interactions[0].channels[0].coupling==0.7);
    REQUIRE(bound.lattice().interactions[0].channels[0].coupling==0);
    REQUIRE_THROWS_AS(original.bind({{1}}),std::invalid_argument);
    REQUIRE_THROWS_AS(original.bind({{1,2},{3}}),std::invalid_argument);
    REQUIRE_THROWS_AS(original.bind({{1},{std::numeric_limits<double>::infinity()}}),std::invalid_argument);
}

TEST_CASE("Hubbard edge support agrees with full mapped subtraction at unequal and zero couplings", "[white][hubbard]") {
    auto lattice=models::hubbard_square();
    const WhiteGraphExpansion plan(lattice,4);
    const EffectiveOperator effective(Coefficients(charge_changes(lattice),4));
    for (const auto& ratios : {std::vector<std::vector<double>>{{0.7},{-0.2}},{{0},{0.9}}}) {
        const auto bound=plan.bind(ratios);
        const auto actual=linked_zero_charge(bound,effective);
        std::vector<OperatorBlock> weights;
        for (std::size_t i=0;i<bound.embeddings().size();++i) {
            const auto& embedding=bound.embeddings()[i];
            auto weight=zero_charge_operator(cluster_model(bound.lattice(),embedding.edges),effective);
            weight.coefficients[0].setZero();
            for (const auto& sub : embedding.subclusters)
                add_embedded_operator(weight,weights[sub.index],sub.map.vertices,-1);
            REQUIRE(weight.basis==actual.weights[i].block.basis);
            for (unsigned n=0;n<=4;++n) {
                const auto& expected=weight.coefficients[n];
                REQUIRE((actual.weights[i].block.coefficients[n]-expected).norm()<1e-10*std::max(1.0,expected.norm()));
            }
            weights.push_back(std::move(weight));
        }
    }
    REQUIRE(plan.cache()->evaluations()==plan.graphs().size());
}
