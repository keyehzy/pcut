#include <pcut/coefficients.hpp>
#include <algorithm>
#include <cstdlib>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace pcut {
namespace {
using Polynomial = std::vector<Rational>;
using Function = std::map<int, Polynomial>; // exp(-rate*l) * polynomial(l)
int sign(int x) { return (x > 0) - (x < 0); }
void add(Polynomial& p, std::size_t i, const Rational& value) {
    if (p.size() <= i) p.resize(i + 1);
    p[i] += value;
}
void clean(Function& f) {
    for (auto it = f.begin(); it != f.end();) {
        auto& p = it->second;
        while (!p.empty() && p.back() == 0) p.pop_back();
        if (p.empty()) it = f.erase(it); else ++it;
    }
}
void product_add(Function& result, const Function& a, const Function& b, int factor) {
    for (const auto& [r, p] : a) for (const auto& [s, q] : b)
        for (std::size_t i = 0; i < p.size(); ++i)
            for (std::size_t j = 0; j < q.size(); ++j)
                add(result[r+s], i+j, factor * p[i] * q[j]);
}
Function integrate(const Function& source, int decay) {
    Function f;
    for (const auto& [rate, p] : source) {
        const int d = decay - rate;
        if (d == 0) {
            for (std::size_t k = 0; k < p.size(); ++k)
                add(f[rate], k+1, p[k] / (k+1));
        } else {
            Polynomial q(p.size());
            for (std::size_t k = p.size(); k-- > 0;) {
                q[k] = p[k];
                if (k+1 < q.size()) q[k] -= (k+1) * q[k+1];
                q[k] /= d;
                add(f[rate], k, q[k]);
            }
            if (!q.empty()) add(f[decay], 0, -q[0]);
        }
    }
    clean(f);
    return f;
}
}
Coefficients::Coefficients(std::vector<int> changes, unsigned order)
    : order_(order), changes_(std::move(changes)) {
    if (order > 64) throw std::invalid_argument("order exceeds supported recursion depth (64)");
    std::sort(changes_.begin(), changes_.end());
    changes_.erase(std::unique(changes_.begin(), changes_.end()), changes_.end());
    int bound = 0;
    for (int m : changes_) {
        if (m < -1'000'000 || m > 1'000'000)
            throw std::invalid_argument("charge change is too large");
        bound = std::max(bound, std::abs(m));
    }
    std::map<Word, Function> functions;
    for (unsigned n = 1; n <= order; ++n) {
        Word word(n);
        std::function<void(unsigned,int)> enumerate = [&](unsigned pos, int total) {
            if (pos < n) {
                for (int m : changes_) { word[pos] = m; enumerate(pos+1, total+m); }
                return;
            }
            if (std::abs(total) > bound) return; // preserved block band width
            Function f;
            if (n == 1) f[std::abs(total)] = {Rational(1)};
            else {
                Function source;
                int left_sum = 0;
                for (unsigned split = 1; split < n; ++split) {
                    left_sum += word[split-1];
                    const int factor = sign(left_sum) - sign(total-left_sum);
                    if (factor == 0) continue;
                    const auto left = functions.find(Word(word.begin(), word.begin()+split));
                    const auto right = functions.find(Word(word.begin()+split, word.end()));
                    if (left != functions.end() && right != functions.end())
                        product_add(source, left->second, right->second, factor);
                }
                clean(source);
                f = integrate(source, std::abs(total));
            }
            if (total == 0) {
                const auto zero = f.find(0);
                if (zero != f.end()) {
                    if (zero->second.size() != 1) throw std::logic_error("nonconvergent pCUT coefficient");
                    if (zero->second[0] != 0) terms_.emplace(word, zero->second[0]);
                }
            }
            if (!f.empty()) functions.emplace(word, std::move(f));
        };
        enumerate(0, 0);
    }
}
Rational Coefficients::at(const Word& word) const {
    if (word.empty() || word.size() > order_) throw std::out_of_range("word order outside table");
    for (int m : word) if (!std::binary_search(changes_.begin(), changes_.end(), m))
        throw std::invalid_argument("word contains charge change outside table alphabet");
    const auto it = terms_.find(word);
    return it == terms_.end() ? Rational(0) : it->second;
}
}
