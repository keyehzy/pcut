#include <pcut/pcut.hpp>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc,char** argv) {
    try {
        std::string geometry="chain", filling="half";
        unsigned order=4;
        double ratio=0.05, U=1;
        for (int i=1;i<argc;++i) {
            const std::string arg=argv[i];
            if (arg=="--help") {
                std::cout<<"pcut_hubbard [--lattice chain|square] [--filling half|quarter]\n"
                         <<"             [--order 1..4] [--ratio t/U] [--U positive-energy]\n"
                         <<"Build infinite-lattice operator weights, then evaluate a four-site\n"
                         <<"open chain or square plaquette at Ne=4 (half) or Ne=2 (quarter).\n";
                return 0;
            }
            if (++i>=argc) throw std::invalid_argument("missing option value");
            const std::string value=argv[i];
            if (arg=="--lattice") geometry=value;
            else if (arg=="--filling") filling=value;
            else {
                std::size_t used=0;
                const double number=std::stod(value,&used);
                if (used!=value.size() || !std::isfinite(number)) throw std::invalid_argument("invalid numeric option");
                if (arg=="--ratio") ratio=number;
                else if (arg=="--U") U=number;
                else if (arg=="--order") {
                    if (number<1 || number>4 || number!=std::floor(number)) throw std::invalid_argument("example order must be 1..4");
                    order=static_cast<unsigned>(number);
                } else throw std::invalid_argument("unknown option: "+arg);
            }
        }
        if ((geometry!="chain" && geometry!="square") || (filling!="half" && filling!="quarter") || U<=0)
            throw std::invalid_argument("invalid geometry, filling, or U");
        const auto lattice=geometry=="chain" ? pcut::models::hubbard_chain() : pcut::models::hubbard_square();
        const pcut::Coefficients coefficients({-1,0,1},order);
        const pcut::EffectiveOperator effective(coefficients);
        const pcut::WhiteGraphExpansion catalog(lattice,order);
        const auto linked=pcut::linked_zero_charge(catalog,effective);
        const pcut::Cluster edges=geometry=="chain" ? pcut::Cluster{{0,{0}},{0,{1}},{0,{2}}} :
            pcut::Cluster{{0,{0,0}},{1,{1,0}},{0,{0,1}},{1,{0,0}}};
        const int electrons=filling=="half" ? 4 : 2;
        const auto block=pcut::assemble_operator(linked,pcut::vertices(lattice,edges),edges,electrons);
        const Eigen::SelfAdjointEigenSolver<pcut::Matrix> ed(U*block.evaluate(ratio));
        if (ed.info()!=Eigen::Success) throw std::runtime_error("finite Hamiltonian diagonalization failed");
        std::cout<<std::setprecision(14)
                 <<"# H/U = Q + (t/U) V; V = -sum(c_dagger c + h.c.)\n"
                 <<"# infinite-lattice connected operator weights: "<<linked.weights.size()<<"\n"
                 <<"# finite open "<<geometry<<", L=4, Ne="<<electrons<<", basis="<<block.basis.size()<<"\n"
                 <<"# coefficients below are dimensionless matrices of H_eff/U\n"
                 <<"order,frobenius_norm\n";
        for (unsigned n=0;n<=order;++n) std::cout<<n<<','<<block.coefficients[n].norm()<<'\n';
        std::cout<<"# finite truncated-Hamiltonian ground energy (units of U input): "<<ed.eigenvalues()[0]<<'\n'
                 <<"# This spectrum is not a thermodynamic energy or a phase prediction.\n";
    } catch (const std::exception& e) {
        std::cerr<<"pcut_hubbard: "<<e.what()<<'\n'; return 1;
    }
}
