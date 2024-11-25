#ifndef BSM_H
#define BSM_H

#include <cmath>
#include <array>

constexpr long double default_S = 5000;
constexpr long double default_K = 5000;
constexpr long double default_r = 0.049;
constexpr long double default_y = 0.011;
constexpr long double default_sigma = 0.155;
constexpr long double default_T = 1;

constexpr size_t TradingDaysPerYear = 252;
#define DaysToYears(Days) (Days / TradingDaysPerYear)
#define numSensitivities 14

#define d1(S, K, r, y, sigma, T) (std::log(S / K) + r - y ) / (sigma * std::sqrt(T)) + sigma * std::sqrt(T) / 2
#define d2(S, K, r, y, sigma, T) d1(S, K, r, y, sigma, T) - sigma * std::sqrt(T)

#define N(x, mu, sigma, T) (0.5 / (sigma * std::sqrt(T)) * (1 + std::erf((x - mu) / (sigma * std::sqrt(2 * T)))))
#define phi(x, mu, sigma, T) (std::exp(-(x - mu) * (x - mu) / (2 * sigma * sigma * T)) / (sigma * std::sqrt(2 * M_PI * T)))

#define BSM_Call(S, K, r, y, sigma, T) (S * std::exp(-y * T) * N(d1(static_cast<long double>(S), static_cast<long double>(K), r, y, sigma, T), 0, 1, T) - K * std::exp(-r * T) * N(d2(static_cast<long double>(S), static_cast<long double>(K), r, y, sigma, T), 0, 1, T))
#define BSM_Put(S, K, r, y, sigma, T) (K * std::exp(-r * T) * (1 - N(d2(static_cast<long double>(S), static_cast<long double>(K), r, y, sigma, T), 0, 1, T)) - S * std::exp(-y * T) * (1 - N(d1(static_cast<long double>(S), static_cast<long double>(K), r, y, sigma, T), 0, 1, T)))

#define BSM_Call_Delta_S(S, K, r, y, sigma, T) std::exp(-y * T) * N(d1(S, K, r, y, sigma, T), 0, 1, T)
#define BSM_Put_Delta_S(S, K, r, y, sigma, T) std::exp(-y * T) * (N(d1(S, K, r, y, sigma, T), 0, 1, T) - 1)
#define BSM_Call_Delta_K(S, K, r, y, sigma, T) -std::exp(-r * T) * N(d2(S, K, r, y, sigma, T), 0, 1, T)
#define BSM_Put_Delta_K(S, K, r, y, sigma, T) std::exp(-r * T) * N(-d2(S, K, r, y, sigma, T), 0, 1, T)
#define BSM_Gamma_SS(S, K, r, y, sigma, T) std::exp(-y * T) * phi(d1(S, K, r, y, sigma, T), 0, 1, T) / (S * sigma * std::sqrt(T))
#define BSM_Gamma_KK(S, K, r, y, sigma, T) std::exp(-r * T) * phi(d1(S, K, r, y, sigma, T), 0, 1, T) / (K * sigma * std::sqrt(T))
#define BSM_Put_Gamma_SK(S, K, r, y, sigma, T) (-S / K) * BSM_Gamma_SS(S, K, r, y, sigma, T)

#define BSM_Call_Theta(S, K, r, y, sigma, T) -S * phi(d1(S, K, r, y, sigma, T), 0, 1, T) * sigma / (2 * std::sqrt(T)) + y * S * std::exp(-y * T) * N(d1(S, K, r, y, sigma, T), 0, 1, T) - r * K * std::exp(-r * T) * N(d2(S, K, r, y, sigma, T), 0, 1, T)
#define BSM_Put_Theta(S, K, r, y, sigma, T) -S * phi(d1(S, K, r, y, sigma, T), 0, 1, T) * sigma / (2 * std::sqrt(T)) - y * S * std::exp(-y * T) * N(-d1(S, K, r, y, sigma, T), 0, 1, T) + r * K * std::exp(-r * T) * N(-d2(S, K, r, y, sigma, T), 0, 1, T)
#define BSM_Vega(S, K, r, y, sigma, T) S * std::exp(-y * T) * phi(d1(S, K, r, y, sigma, T), 0, 1, T) * std::sqrt(T) / 100 // Re-base to price change per bp instead of per pct change
#define BSM_Call_Rho_r(S, K, r, y, sigma, T) K * T * std::exp(-r * T) * N(d2(S, K, r, y, sigma, T), 0, 1, T) / 100 // Re-base to price change per bp instead of per pct change
#define BSM_Put_Rho_r(S, K, r, y, sigma, T) -K * T * std::exp(-r * T) * (1 - N(d2(S, K, r, y, sigma, T), 0, 1, T)) / 100 // Re-base to price change per bp instead of per pct change
#define BSM_Call_Rho_y(S, K, r, y, sigma, T) -S * T * std::exp(-y * T) * N(d1(S, K, r, y, sigma, T), 0, 1, T) / 100 // Re-base to price change per bp instead of per pct change
#define BSM_Put_Rho_y(S, K, r, y, sigma, T) S * T * std::exp(-y * T) * (1 - N(d1(S, K, r, y, sigma, T), 0, 1, T)) / 100 // Re-base to price change per bp instead of per pct change

#define BSM_CoN_Call(S, K, r, y, sigma, T) std::exp(-r*T)*N(d2(S, K, r, y, sigma, T), 0, 1, T)
#define BSM_CoN_Put(S, K, r, y, sigma, T) std::exp(-r*T)*(1-N(d2(S, K, r, y, sigma, T), 0, 1, T))
#define BSM_AoN_Call(S, K, r, y, sigma, T) S*std::exp(-y*T)*N(d1(S, K, r, y, sigma, T), 0, 1, T) 
#define BSM_AoN_Put(S, K, r, y, sigma, T) S*std::exp(-y*T)*(1-N(d1(S, K, r, y, sigma, T), 0, 1, T))

constexpr long double ATM_Call = default_S*BSM_Call((long double)1.0, (long double)1.0, default_r, default_y, default_sigma, default_T);
constexpr long double ATM_Put = default_S*BSM_Put((long double)1.0, (long double)1.0, default_r, default_y, default_sigma, default_T);
constexpr long double OTM10pct_Call = default_S*BSM_Call((long double)1.0, (long double)1.1, default_r, default_y, default_sigma, default_T);
constexpr long double OTM10pct_Put = default_S*BSM_Put((long double)1.0, (long double)0.9, default_r, default_y, default_sigma, default_T);
constexpr long double ITM10pct_Call = default_S*BSM_Call((long double)1.0, (long double)0.9, default_r, default_y, default_sigma, default_T);
constexpr long double ITM10pct_Put = default_S*BSM_Put((long double)1.0, (long double)1.1, default_r, default_y, default_sigma, default_T);


constexpr long double ATM_CoN_Call = BSM_CoN_Call((long double)1.0, (long double)1.0, default_r, default_y, default_sigma, default_T);
constexpr long double ATM_AoN_Call = default_S*BSM_AoN_Call((long double)1.0, (long double)1.0, default_r, default_y, default_sigma, default_T);
constexpr long double ATM_CoN_Put = BSM_CoN_Put((long double)1.0, (long double)1.0, default_r, default_y, default_sigma, default_T);
constexpr long double ATM_AoN_Put = default_S*BSM_AoN_Put((long double)1.0, (long double)1.0, default_r, default_y, default_sigma, default_T);

// {Call Delta S, Put Delta S, Call Delta K, Put Delta K, Gamma SS, Gamma KK, Gamma SK, Call Theta, Put Theta, Call Vega, Put Vega, Call Rho r, Put Rho r, Call Rho y, Put Rho y}
#define sensitivities(S, K, r, y, sigma, T) std::array<long double, numSensitivities>{BSM_Call_Delta_S(S, K, r, y, sigma, T), BSM_Put_Delta_S(S, K, r, y, sigma, T), BSM_Call_Delta_K(S, K, r, y, sigma, T), BSM_Put_Delta_K(S, K, r, y, sigma, T), BSM_Gamma_SS(S, K, r, y, sigma, T), BSM_Gamma_KK(S, K, r, y, sigma, T), BSM_Put_Gamma_SK(S, K, r, y, sigma, T), BSM_Call_Theta(S, K, r, y, sigma, T), BSM_Put_Theta(S, K, r, y, sigma, T), BSM_Vega(S, K, r, y, sigma, T), BSM_Call_Rho_r(S, K, r, y, sigma, T), BSM_Put_Rho_r(S, K, r, y, sigma, T), BSM_Call_Rho_y(S, K, r, y, sigma, T), BSM_Put_Rho_y(S, K, r, y, sigma, T)}

constexpr std::array<long double, numSensitivities> ATM_Sensitivities = sensitivities(default_S, static_cast<long double>(default_S), default_r, default_y, default_sigma, default_T);
constexpr std::array<long double, numSensitivities> OTM10pct_Sensitivities = sensitivities(default_S, static_cast<long double>(default_S*0.9), default_r, default_y, default_sigma, default_T);
constexpr std::array<long double, numSensitivities> ITM10pct_Sensitivities = sensitivities(default_S, static_cast<long double>(default_S*1.1), default_r, default_y, default_sigma, default_T);



constexpr long double strikeGranularity = (long double)0.001; //0.1% granularity

// Breeden Litzenberger Replicating Portfolios for Binary Options with and without ATM Vol Skew, values becoming largely mismatched with skew
template <typename T, typename M>
constexpr inline long double skewedVol(const T& S, const M& K) {return (long double)0.155 + 0.3*logl(static_cast<long double>(S)/static_cast<long double>(K));} //Volatility skew
constexpr double volATM = skewedVol(default_S, default_S);
constexpr long double skewedVolATMPlus = skewedVol(default_S, default_S*(1+strikeGranularity));
constexpr long double skewedVolATMMinus = skewedVol(default_S, default_S*(1-strikeGranularity));
constexpr long double default_S_minus = static_cast<long double>(default_S*(1-strikeGranularity));
constexpr long double default_S_plus = default_S*(1+strikeGranularity);

constexpr long double BSM_Call_ATM_Minus_skewed = BSM_Call(default_S, default_S_minus, default_r, default_y, skewedVolATMMinus, 1);
constexpr long double BSM_Call_ATM_Plus_skewed = BSM_Call(default_S, default_S_plus, default_r, default_y, skewedVolATMPlus, 1);
constexpr long double BSM_Put_ATM_Minus_skewed = BSM_Put(default_S, default_S_minus, default_r, default_y, skewedVolATMMinus, 1);
constexpr long double BSM_Put_ATM_Plus_skewed = BSM_Put(default_S, default_S_plus, default_r, default_y, skewedVolATMPlus, 1);

constexpr long double BSM_Call_ATM_Minus = BSM_Call(default_S, default_S_minus, default_r, default_y, volATM, 1);
constexpr long double BSM_Call_ATM_Plus = BSM_Call(default_S, default_S_plus, default_r, default_y, volATM, 1);
constexpr long double BSM_Put_ATM_Minus = BSM_Put(default_S, default_S_minus, default_r, default_y, volATM, 1);
constexpr long double BSM_Put_ATM_Plus = BSM_Put(default_S, default_S_plus, default_r, default_y, volATM, 1);


#endif // BSM_H