#ifndef LOGLOGAPPROX_H_
#define LOGLOGAPPROX_H_

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>

namespace detail {
inline constexpr std::array<double, 1024> make_log2_uint_table() noexcept {
    std::array<double, 1024> table{};
    for (std::size_t i = 1; i <= table.size(); ++i) {
        table[i - 1] = std::log2(static_cast<double>(i));
    }
    return table;
}

inline const std::array<double, 1024>& log2_uint_table() noexcept {
    static constexpr auto table = make_log2_uint_table();
    return table;
}

inline double log2_small_uint(std::uint32_t value) noexcept {
    if (value == 0U) {
        return -std::numeric_limits<double>::infinity();
    }
    if (value <= 1024U) {
        return log2_uint_table()[value - 1U];
    }
    return std::log2(static_cast<double>(value));
}

inline double fast_log2_mantissa_1_to_2_order2(double value) noexcept {
    // value is in [1, 2).  We do not call frexp() again here because we already know
    // the mantissa is value / 2.0 and the exponent contribution is 1.
    const double dM = 0.5 * value;
    const double denom = std::fma(std::fma(0.49463685172392841, dM, 1.426594307123505), dM, 0.2533316901691966);
    const double numer = std::fma(std::fma(1.9127166899499954, dM, -0.68851400593499545), dM, -1.22420645509838);
    const double x = 1.0 / denom;
    return 1.0 + x * numer;
}

inline double fast_log2_mantissa_1_to_2_order3(double value) noexcept {
    const double dM = 0.5 * value;
    const double denom = std::fma(std::fma(std::fma(0.22977948696488379, dM, 1.4961611668393175), dM, 1.071708023446889), dM, 0.084444549259932208);
    const double numer = std::fma(std::fma(std::fma(1.1098414161667869, dM, 1.4491119665946153), dM, -2.0697678829202806), dM, -0.48918550780729392);
    const double x = 1.0 / denom;
    return 1.0 + x * numer;
}

inline double fast_log2_mantissa_1_to_2_order4(double value) noexcept {
    const double dM = 0.5 * value;
    const double denom = std::fma(std::fma(std::fma(std::fma(0.1068562844523792, dM, 2.0062979261642901), dM, 1.2392957064266512), dM, 0.63680961689938775), dM, 0.028211791264274255);
    const double numer = std::fma(std::fma(std::fma(std::fma(0.59329970349044314, dM, 2.3979646338966889), dM, -0.96358966800238843), dM, -1.8439274267589987), dM, -0.18374724264449727);
    const double x = 1.0 / denom;
    return 1.0 + x * numer;
}

inline double fast_log2_mantissa_1_to_2_order5(double value) noexcept {
    const double dM = 0.5 * value;
    const double denom = std::fma(std::fma(std::fma(std::fma(std::fma(0.1636476193217775432220, dM, 2.929373781304793311620), dM, 5.900277160781167040682), dM, 8.339836924158730013801), dM, 1.037334812091420754854), dM, 0.02900560762882882170083);
    const double numer = std::fma(std::fma(std::fma(std::fma(std::fma(1.000000000000000000000, dM, 7.732354604084914484474), dM, -8.646069665182267272030), dM, 3.893816323736673190581), dM, -3.773895413775759433150), dM, -0.2062058488638166542373);
    const double x = 1.0 / denom;
    return 1.0 + x * numer;
}
} // namespace detail

// Fast approximation of log2(log2(x)) for x > 4.
//
// Let x = 2^m * e with e in [1,2), so m = floor(log2(x)). Then
//   log2(log2(x)) = log2(m + log2(e))
//                  = log2(m) + log2(1 + log2(e) / m)
//
// The final log2(1 + y) term is evaluated with the same polynomial-ratio idea
// used by the other approximations in this project.
inline double fastLog2Log2(double x) noexcept {
    if (!std::isfinite(x)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (x < 4.0) {
        return std::log2(std::log2(x));
    }

    int exponent = 0;
    const double mantissa = std::frexp(x, &exponent);
    const int integer_m = exponent - 1;

    // x = 2^m * e, with e in [1, 2).  We know y = log2(e) / m is in [0, 1], so
    // 1 + y is in [1, 2) and can be handled by the same fastLog2pN-style formula,
    // without another frexp() and without extra checks.
    const double e = mantissa * 2.0;
    const double y = detail::fast_log2_mantissa_1_to_2_order2(e) / static_cast<double>(integer_m);

    return detail::log2_small_uint(static_cast<std::uint32_t>(integer_m)) +
           detail::fast_log2_mantissa_1_to_2_order2(1.0 + y);
}

extern "C" {
double loglog(double x) noexcept {
    if (x <= 4.0) {
        return std::log(std::log(x));
    }
    constexpr double ln2 = 0.693147180559945309417232121458176568075500134360255254120L;
    constexpr double lnln2 = -0.3665129205816643; // ln(ln(2))
    const double log2log2x = fastLog2Log2(x);
    return lnln2 + ln2 * log2log2x;
}
}
#endif // LOGLOGAPPROX_H_
