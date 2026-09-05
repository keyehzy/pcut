#include <pcut/pcut.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <numeric>

namespace {
using namespace pcut;
using Bond=std::pair<unsigned,unsigned>;
// Independent full-Fock construction: elementary creation matrices in GLOBAL
// orbital order. No LocalTerm, cluster embedding, or pCUT application is used.
Matrix annihilator(unsigned sites, unsigned orbital) {
    const unsigned dim=1u<<(2*sites);
    Matrix c=Matrix::Zero(dim,dim);
    for (unsigned state=0; state<dim; ++state) if ((state>>orbital)&1u) {
        int sign=1;
        for (unsigned k=0; k<orbital; ++k) if ((state>>k)&1u) sign=-sign;
        c(state^(1u<<orbital),state)=sign;
    }
    return c;
}
Matrix kinetic(unsigned sites, const std::vector<Bond>& bonds) {
    const unsigned dim=1u<<(2*sites);
    Matrix v=Matrix::Zero(dim,dim);
    for (auto [a,b] : bonds) for (unsigned spin=0; spin<2; ++spin) {
        const Matrix hop=annihilator(sites,2*a+spin).adjoint()*annihilator(sites,2*b+spin);
        v-=hop+hop.adjoint();
    }
    return v;
}
unsigned doublons(State state, unsigned sites) {
    unsigned q=0;
    for (unsigned i=0; i<sites; ++i) q+=((state>>(2*i))&3u)==3u;
    return q;
}
ClusterModel finite(unsigned sites, const std::vector<Bond>& bonds) {
    std::vector<LocalTerm> terms;
    for (auto [a,b] : bonds) terms.push_back({{a,b},models::hubbard_hopping(),true});
    return {std::vector<LocalSpace>(sites,models::hubbard_site()),terms};
}
Matrix restrict_matrix(const Matrix& h, const std::vector<State>& basis) {
    const auto size=static_cast<Eigen::Index>(basis.size());
    Matrix out(size,size);
    for (std::size_t i=0; i<basis.size(); ++i) for (std::size_t j=0; j<basis.size(); ++j)
        out(static_cast<Eigen::Index>(i),static_cast<Eigen::Index>(j))=h(basis[i],basis[j]);
    return out;
}
void near(const Matrix& a, const Matrix& b, double tolerance=2e-11) {
    REQUIRE(a.rows()==b.rows()); REQUIRE(a.cols()==b.cols());
    REQUIRE((a-b).norm()<tolerance);
}
const EffectiveOperator& fourth() {
    static const EffectiveOperator effective(Coefficients({-1,0,1},4));
    return effective;
}
const LinkedOperator& square_linked() {
    static const auto linked=linked_zero_charge(WhiteGraphExpansion(models::hubbard_square(),4),fourth());
    return linked;
}
std::vector<State> spin_basis(unsigned sites) {
    std::vector<State> basis;
    for (unsigned s=0; s<(1u<<sites); ++s) {
        State state=0;
        for (unsigned i=0; i<sites; ++i) state|=State(1+((s>>i)&1u))<<(2*i);
        basis.push_back(state);
    }
    return basis;
}
Matrix spin_block(const OperatorBlock& block,unsigned n) {
    const auto states=spin_basis(static_cast<unsigned>(block.spaces.size()));
    Matrix h(states.size(),states.size());
    for (std::size_t i=0;i<states.size();++i) for (std::size_t j=0;j<states.size();++j) {
        const auto row=std::lower_bound(block.basis.begin(),block.basis.end(),states[i])-block.basis.begin();
        const auto col=std::lower_bound(block.basis.begin(),block.basis.end(),states[j])-block.basis.begin();
        h(static_cast<Eigen::Index>(i),static_cast<Eigen::Index>(j))=block.coefficients[n](row,col);
    }
    return h;
}
Matrix spin_dot(unsigned sites,unsigned a,unsigned b) {
    const unsigned dim=1u<<sites;
    Matrix h=Matrix::Zero(dim,dim);
    for (unsigned s=0;s<dim;++s) {
        const bool opposite=((s>>a)&1u)!=((s>>b)&1u);
        h(s,s)=opposite ? -0.25 : 0.25;
        if (opposite) h(s^(1u<<a)^(1u<<b),s)=0.5;
    }
    return h;
}
}

TEST_CASE("Degenerate charge-zero spaces and strict vacuum drivers", "[hubbard][contracts]") {
    const auto model=finite(2,{{0,1}});
    REQUIRE(model.conserves_particles());
    REQUIRE(model.spaces()[0].charges==std::vector<int>{0,0,0,1});
    REQUIRE(model.charge(model.encode({1,2}))==0);
    REQUIRE(model.particle_number(model.encode({1,2}))==2);
    REQUIRE(zero_charge_basis(model).size()==9);
    REQUIRE(zero_charge_basis(model,1).size()==4);
    REQUIRE(zero_charge_basis(model,2).size()==4);
    REQUIRE_THROWS_AS(fourth().vacuum(model),std::invalid_argument);
    REQUIRE_THROWS_AS(one_particle_basis(model),std::invalid_argument);
    REQUIRE_THROWS_AS(irreducible_sectors(model,fourth(),0),std::invalid_argument);
    const WhiteGraphExpansion catalog(models::hubbard_chain(),4);
    REQUIRE_THROWS_AS(linked_expand(catalog,fourth()),std::invalid_argument);
    REQUIRE_THROWS_AS(linked_expand_sectors(catalog,fourth(),0),std::invalid_argument);
    auto space=models::hubbard_site(); space.parity[1]=0;
    REQUIRE_THROWS_AS(space.validate(),std::invalid_argument);
    space=models::hubbard_site(); space.particles.pop_back();
    REQUIRE_THROWS_AS(space.validate(),std::invalid_argument);
    Matrix odd=Matrix::Zero(4,4); odd(0,1)=odd(1,0)=1;
    REQUIRE_THROWS_AS(ClusterModel({models::hubbard_site()},{{{0},odd,true}}),std::invalid_argument);
    REQUIRE(zero_charge_basis(model,3).empty());
    REQUIRE(zero_charge_basis(finite(20,{}),1000).empty());
    Matrix pair=Matrix::Zero(4,4); pair(0,3)=pair(3,0)=1;
    const ClusterModel nonconserving({models::hubbard_site()},{{{0},pair,true}});
    REQUIRE_FALSE(nonconserving.conserves_particles());
    REQUIRE_THROWS_AS(zero_charge_operator(nonconserving,fourth(),0),std::invalid_argument);
    const ClusterModel ungraded(std::vector<LocalSpace>(2,models::hubbard_site()),
                                {{{0,1},models::hubbard_hopping(),false}});
    REQUIRE_THROWS_AS(zero_charge_operator(ungraded,fourth()),std::invalid_argument);
    REQUIRE_THROWS_AS(linked_zero_charge(WhiteGraphExpansion(models::hubbard_chain(),3),fourth()),std::invalid_argument);
}

TEST_CASE("Fermionic hopping matches global Fock matrices for nonadjacent and reversed legs", "[hubbard][fermion]") {
    const std::vector<Bond> bonds{{2,0},{0,3},{3,1},{1,2}};
    const auto model=finite(4,bonds);
    const Matrix v=kinetic(4,bonds);
    Matrix exact=v*0.17;
    for (unsigned s=0;s<256;++s) exact(s,s)+=doublons(s,4);
    near(model.dense_hamiltonian(0.17),exact);
    REQUIRE(model.changes()==std::vector<int>{-1,0,1});
    // A spectator electron reverses a nonadjacent hopping matrix element.
    const auto hop=finite(3,{{0,2}});
    REQUIRE(hop.apply(0,{{hop.encode({1,0,0}),1}}).at(hop.encode({0,0,1}))==Complex(-1));
    REQUIRE(hop.apply(0,{{hop.encode({1,2,0}),1}}).at(hop.encode({0,2,1}))==Complex(1));
    const auto child=zero_charge_operator(finite(2,{{0,1}}),fourth());
    for (const std::vector<std::size_t>& map : {std::vector<std::size_t>{0,2}, {2,0}, {3,1}}) {
        const auto target=finite(4,{{static_cast<unsigned>(map[0]),static_cast<unsigned>(map[1])}});
        auto embedded=zero_charge_operator(finite(4,{}),fourth());
        add_embedded_operator(embedded,child,map);
        const auto direct=zero_charge_operator(target,fourth());
        for (unsigned n=0;n<=4;++n) near(embedded.coefficients[n],direct.coefficients[n]);
    }
    auto parent=zero_charge_operator(finite(3,{}),fourth());
    REQUIRE_THROWS_AS(add_embedded_operator(parent,child,{1,1}),std::invalid_argument);
    const auto selected=zero_charge_operator(finite(2,{{0,1}}),fourth(),1);
    REQUIRE_THROWS(add_embedded_operator(parent,selected,{0,2}));
}

TEST_CASE("Hubbard dimer analytic singlet and projected one-electron hopping", "[hubbard][analytic]") {
    const auto model=finite(2,{{0,1}});
    const auto half=zero_charge_operator(model,fourth(),2);
    const Matrix exchange=4*(spin_dot(2,0,1)-0.25*Matrix::Identity(4,4));
    near(half.coefficients[2],exchange);
    near(half.coefficients[4],-4*exchange);
    near(half.coefficients[1],Matrix::Zero(4,4));
    near(half.coefficients[3],Matrix::Zero(4,4));
    const auto quarter=zero_charge_operator(model,fourth(),1);
    near(quarter.coefficients[1],restrict_matrix(kinetic(2,{{0,1}}),quarter.basis));
    for (unsigned n=2;n<=4;++n) near(quarter.coefficients[n],Matrix::Zero(4,4));
    double previous=0;
    for (double x : {0.06,0.03}) {
        const Eigen::SelfAdjointEigenSolver<Matrix> ed(half.evaluate(x));
        const double exact=(1-std::sqrt(1+16*x*x))/2;
        const double error=std::abs(ed.eigenvalues()[0]-exact);
        REQUIRE(error<7e-6);
        if (previous) REQUIRE(previous/error>60);
        previous=error;
    }
}

TEST_CASE("Second order contains exchange and all three-site correlated hoppings", "[hubbard][second]") {
    const unsigned sites=3;
    const std::vector<Bond> bonds{{0,1},{1,2}};
    const auto model=finite(sites,bonds);
    const auto block=zero_charge_operator(model,fourth());
    const Matrix v=kinetic(sites,bonds);
    Matrix q=Matrix::Zero(64,64);
    for (unsigned s=0;s<64;++s) if (doublons(s,sites)==1) q(s,s)=1;
    near(block.coefficients[1],restrict_matrix(v,block.basis));
    near(block.coefficients[2],restrict_matrix(-v*q*v,block.basis));
    // Independently expand the t-J + three-site formula in global fermions:
    // -sum_(i!=k neighbors of j,s) [c†is n_j,-s c_ks
    //                         - c†is c†j,-s c_js c_k,-s], with P on both ends.
    std::vector<Matrix> c;
    for (unsigned a=0;a<6;++a) c.push_back(annihilator(sites,a));
    std::vector<Matrix> n;
    for (const auto& op : c) n.push_back(op.adjoint()*op);
    Matrix tj=Matrix::Zero(64,64), three=tj;
    for (auto [a,b] : bonds) {
        const Matrix za=(n[2*a]-n[2*a+1])*0.5, zb=(n[2*b]-n[2*b+1])*0.5;
        const Matrix plusa=c[2*a].adjoint()*c[2*a+1], plusb=c[2*b].adjoint()*c[2*b+1];
        tj+=4*(za*zb+0.5*(plusa*plusb.adjoint()+plusa.adjoint()*plusb)-
               0.25*(n[2*a]+n[2*a+1])*(n[2*b]+n[2*b+1]));
    }
    for (unsigned i : {0u,2u}) for (unsigned spin=0;spin<2;++spin) {
        const unsigned k=2-i, j=1, other=1-spin;
        three-=c[2*i+spin].adjoint()*n[2*j+other]*c[2*k+spin]-
               c[2*i+spin].adjoint()*c[2*j+other].adjoint()*c[2*j+spin]*c[2*k+other];
    }
    REQUIRE(restrict_matrix(three,block.basis).norm()>1);
    near(block.coefficients[2],restrict_matrix(tj+three,block.basis));
    const auto child=zero_charge_operator(finite(2,{{0,1}}),fourth());
    auto weight=block;
    add_embedded_operator(weight,child,{0,1},-1);
    add_embedded_operator(weight,child,{1,2},-1);
    near(weight.coefficients[2],restrict_matrix(three,block.basis));
}

TEST_CASE("Fourth order matches primary canonical formulas in the sign-generator convention", "[hubbard][literature]") {
    // Chernyshev et al. cond-mat/0407255, AM and Canonical_T_explicit:
    // CT1 + gamma [T0,[T0,T- T+]], gamma=1/4 for this continuous generator.
    const Coefficients exact({-1,0,1},4);
    const std::map<Word,Rational> fixtures{
        {{0},1},{{-1,1},-1},{{-1,0,1},1},{{-1,1,0},Rational(-1)/2},{{0,-1,1},Rational(-1)/2},
        {{-1,1,-1,1},1},{{-1,-1,1,1},Rational(-1)/2},{{-1,0,0,1},-1},
        {{-1,0,1,0},1},{{0,-1,0,1},1},{{-1,1,0,0},Rational(-1)/4},
        {{0,0,-1,1},Rational(-1)/4},{{0,-1,1,0},Rational(-1)/2}};
    for (const auto& [word,value] : fixtures) REQUIRE(exact.at(word)==value);
    for (const auto& [word,value] : exact.terms()) {
        int q=0; bool allowed=true;
        for (auto it=word.rbegin();it!=word.rend();++it) { q+=*it; if (q<0) allowed=false; }
        if (allowed) REQUIRE(fixtures.contains(word));
    }
    const std::vector<Bond> bonds{{0,1},{1,3},{3,2},{2,0}};
    const auto model=finite(4,bonds);
    const Matrix v=kinetic(4,bonds);
    Matrix minus=Matrix::Zero(256,256), zero=minus, plus=minus;
    for (unsigned i=0;i<256;++i) for (unsigned j=0;j<256;++j) {
        const int delta=static_cast<int>(doublons(i,4))-static_cast<int>(doublons(j,4));
        if (delta==-1) minus(i,j)=v(i,j);
        if (delta==0) zero(i,j)=v(i,j);
        if (delta==1) plus(i,j)=v(i,j);
    }
    const Matrix a=minus*plus;
    REQUIRE(restrict_matrix(minus*minus*plus*plus,zero_charge_basis(model,4)).norm()>1);
    // Two-doublon intermediate states are essential even though external Q=0.
    const Matrix ct1=a*a-0.5*minus*minus*plus*plus-minus*zero*zero*plus+
        minus*zero*plus*zero+zero*minus*zero*plus-0.5*(a*zero*zero+zero*zero*a);
    const Matrix rotation=0.25*(zero*zero*a+a*zero*zero-2*zero*a*zero);
    const auto all=zero_charge_operator(model,fourth());
    near(all.coefficients[4],restrict_matrix(ct1+rotation,all.basis));
    for (const auto& h : all.coefficients) {
        near(h,h.adjoint());
        for (std::size_t i=0;i<all.basis.size();++i) for (std::size_t j=0;j<all.basis.size();++j)
            if (model.particle_number(all.basis[i])!=model.particle_number(all.basis[j]))
                REQUIRE(std::abs(h(static_cast<Eigen::Index>(i),static_cast<Eigen::Index>(j)))<1e-12);
    }
    REQUIRE(restrict_matrix(rotation,zero_charge_basis(model,2)).norm()>0.1);
    near(restrict_matrix(rotation,zero_charge_basis(model,4)),Matrix::Zero(16,16));
    // Delannoy et al. cond-mat/0412033, Hs4, directly at half filling.
    const auto half=zero_charge_operator(model,fourth(),4);
    near(half.coefficients[4],restrict_matrix(a*a-minus*zero*zero*plus-0.5*minus*minus*plus*plus,half.basis));
}

TEST_CASE("Operator linking reconstructs connected and disconnected systems with arbitrary Fock ordering", "[hubbard][linked]") {
    const auto lattice=models::hubbard_chain();
    const WhiteGraphExpansion catalog(lattice,4);
    const auto linked=linked_zero_charge(catalog,fourth());
    REQUIRE(catalog.embedding_subclusters(1).size()==2);
    for (const Cluster& edges : {Cluster{{0,{0}},{0,{1}},{0,{2}}}, Cluster{{0,{0}},{0,{2}}}}) {
        const auto model=cluster_model(lattice,edges);
        const auto sites=vertices(lattice,edges);
        const auto assembled=assemble_operator(linked,sites,edges);
        const auto direct=zero_charge_operator(model,fourth());
        for (unsigned n=0;n<=4;++n) near(assembled.coefficients[n],direct.coefficients[n]);
    }
    // Disconnected union has zero connected residual, with interleaved components.
    const auto disconnected=finite(4,{{0,2},{1,3}});
    auto residual=zero_charge_operator(disconnected,fourth());
    const auto bond=zero_charge_operator(finite(2,{{0,1}}),fourth());
    add_embedded_operator(residual,bond,{0,2},-1);
    add_embedded_operator(residual,bond,{3,1},-1);
    for (const auto& h : residual.coefficients) near(h,Matrix::Zero(81,81));
    const Cluster plaquette{{0,{0,0}},{1,{1,0}},{0,{0,1}},{1,{0,0}}};
    const auto& square=square_linked();
    auto sites=vertices(square.lattice,plaquette);
    std::swap(sites[0],sites[3]); std::swap(sites[1],sites[2]);
    std::vector<Bond> bonds;
    for (const auto& edge : plaquette) {
        const auto vs=vertices(square.lattice,{edge});
        bonds.emplace_back(static_cast<unsigned>(std::find(sites.begin(),sites.end(),vs[0])-sites.begin()),
                           static_cast<unsigned>(std::find(sites.begin(),sites.end(),vs[1])-sites.begin()));
    }
    for (int electrons : {2,4}) {
        const auto assembled=assemble_operator(square,sites,plaquette,electrons);
        const auto direct=zero_charge_operator(finite(4,bonds),fourth(),electrons);
        REQUIRE(assembled.basis.size()==(electrons==2 ? 24 : 16));
        for (unsigned n=0;n<=4;++n) near(assembled.coefficients[n],direct.coefficients[n]);
    }
}

TEST_CASE("Infinite square-lattice fourth-order spin couplings and plaquette ring exchange", "[hubbard][literature][linked]") {
    const auto& linked=square_linked();
    std::map<Coordinate,std::array<double,5>> exchange;
    unsigned rings=0;
    for (const auto& weight : linked.weights) {
        const unsigned sites=static_cast<unsigned>(weight.sites.size());
        for (unsigned n : {2u,4u}) {
            const Matrix h=spin_block(weight.block,n);
            for (unsigned a=0;a<sites;++a) for (unsigned b=a+1;b<sites;++b) {
                // Orthogonal Pauli projection, not a fit. Tr(H Sz_a Sz_b)/Tr((Sz_a Sz_b)^2).
                double j=0;
                for (unsigned s=0;s<(1u<<sites);++s)
                    j+=h(s,s).real()*(((s>>a)&1u)==((s>>b)&1u) ? 4 : -4)/(1u<<sites);
                Coordinate delta(2);
                for (unsigned d=0;d<2;++d) delta[d]=weight.sites[b].cell[d]-weight.sites[a].cell[d];
                exchange[delta][n]+=j;
            }
        }
        if (weight.edges.size()==4 && sites==4) {
            ++rings;
            // Sorted square vertices: 0=(0,0),1=(0,1),2=(1,0),3=(1,1).
            const Matrix ring=spin_dot(4,0,1)*spin_dot(4,3,2)+spin_dot(4,0,2)*spin_dot(4,1,3)-
                              spin_dot(4,0,3)*spin_dot(4,1,2);
            Matrix expected=80*ring+Matrix::Identity(16,16);
            for (unsigned a=0;a<4;++a) for (unsigned b=a+1;b<4;++b) expected-=4*spin_dot(4,a,b);
            near(spin_block(weight.block,4),expected);
        }
    }
    REQUIRE(rings==1);
    for (const auto& [delta,j] : exchange) {
        const int norm=delta[0]*delta[0]+delta[1]*delta[1];
        REQUIRE(j[2]==Catch::Approx(norm==1 ? 4.0 : 0.0).margin(1e-11));
        REQUIRE(j[4]==Catch::Approx(norm==1 ? -24.0 : (norm==2 || norm==4 ? 4.0 : 0.0)).margin(1e-10));
    }
}

TEST_CASE("Half and quarter finite-system spectra converge to independent Hubbard ED", "[hubbard][ed]") {
    for (const std::vector<Bond>& bonds : {std::vector<Bond>{{0,1},{1,2},{2,3}}, {{0,1},{1,3},{3,2},{2,0}}}) {
        const auto model=finite(4,bonds);
        const Matrix v=kinetic(4,bonds);
        for (int particles : {2,4}) {
            const auto effective=zero_charge_operator(model,fourth(),particles);
            std::vector<State> full;
            for (unsigned s=0;s<256;++s) if (std::popcount(s)==particles) full.push_back(s);
            double previous=0;
            for (double x : {0.04,0.02}) {
                Matrix h=x*v;
                for (unsigned s=0;s<256;++s) h(s,s)+=doublons(s,4);
                const Eigen::SelfAdjointEigenSolver<Matrix> ed(restrict_matrix(h,full)), low(effective.evaluate(x));
                REQUIRE(ed.info()==Eigen::Success); REQUIRE(low.info()==Eigen::Success);
                const double error=(ed.eigenvalues().head(low.eigenvalues().size())-low.eigenvalues()).cwiseAbs().maxCoeff();
                REQUIRE(error<5e-5);
                if (previous) REQUIRE(previous/error>25);
                previous=error;
            }
        }
    }
}

TEST_CASE("First-order Hubbard linking ignores larger catalog animals", "[hubbard]") {
    const auto lattice=pcut::models::hubbard_chain();
    const pcut::EffectiveOperator first(pcut::Coefficients(pcut::charge_changes(lattice),1));
    const auto small=pcut::linked_zero_charge(pcut::WhiteGraphExpansion(lattice,1),first);
    const auto large=pcut::linked_zero_charge(pcut::WhiteGraphExpansion(lattice,2),first);
    REQUIRE(small.weights.size()==1);
    REQUIRE(large.weights.size()==1);
    for (unsigned n=0;n<=1;++n) REQUIRE(large.weights[0].block.coefficients[n]==small.weights[0].block.coefficients[n]);
    const pcut::EffectiveOperator zero(pcut::Coefficients(pcut::charge_changes(lattice),0));
    REQUIRE(pcut::linked_zero_charge(pcut::WhiteGraphExpansion(lattice,2),zero).weights.empty());
}
TEST_CASE("Symbolic Hubbard bindings preserve reversed legs and coupling sweeps", "[hubbard][sweep]") {
    auto lattice=pcut::models::hubbard_chain();
    lattice.interactions[0].legs[1].cell={2};
    std::reverse(lattice.interactions[0].legs.begin(),lattice.interactions[0].legs.end());
    const pcut::EffectiveOperator second(pcut::Coefficients(pcut::charge_changes(lattice),2));
    const pcut::WhiteGraphExpansion initial(lattice,2);
    (void)pcut::linked_zero_charge(initial,second);
    const auto count=initial.cache()->evaluations();
    lattice.interactions[0].channels[0].coupling*=0.7;
    const auto rebound=initial.bind(lattice.couplings());
    for (std::size_t index=0;index<rebound.embeddings().size();++index) {
        const auto& entry=rebound.embeddings()[index];
        const auto model=pcut::cluster_model(lattice,entry.edges);
        const auto expected=pcut::zero_charge_operator(model,second);
        const auto actual=rebound.block(index,second,expected.basis);
        for (unsigned n=0;n<=2;++n) near(actual[n],expected.coefficients[n]);
    }
    // Full raw blocks and edge-support projections have distinct cache contexts.
    REQUIRE(initial.cache()->evaluations()==2*count);
}

TEST_CASE("Charge-zero basis enumeration retains the complete degenerate manifold", "[hubbard]") {
    const auto model=finite(8,{});
    const auto basis=pcut::zero_charge_basis(model);
    REQUIRE(basis.size()==6561); // 3^8, beyond the former 4096-state cap
    REQUIRE(std::is_sorted(basis.begin(),basis.end()));
    REQUIRE(std::adjacent_find(basis.begin(),basis.end())==basis.end());
    for (auto state : basis) REQUIRE(model.charge(state)==0);
    REQUIRE_THROWS_AS(pcut::zero_charge_basis(model,-1),std::invalid_argument);
}
