#ifndef BREEDEN_LITZENBERGER_H
#define BREEDEN_LITZENBERGER_H

#include "bsmUtils.hpp"
#include <iostream>

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


// Breeden Litzenberger Replicating Portfolios for Binary Options with and without ATM Vol Skew
#define normalizer(S)  1/(strikeGranularity*S) //To ensure the BL_CoN contract pays off 1 dollar in payoff region

// Breeden Litzenberger Replicating Portfolios for Binary Options without ATM Vol Skew
constexpr long double ATM_BL_CoN_Call_super = normalizer(default_S)*(BSM_Call_ATM_Minus - ATM_Call);
constexpr long double ATM_BL_CoN_Call_sub = normalizer(default_S)*(ATM_Call - BSM_Call_ATM_Plus);
constexpr long double ATM_BL_CoN_Put_sub = normalizer(default_S)*(-BSM_Put_ATM_Minus + ATM_Put);
constexpr long double ATM_BL_CoN_Put_super = normalizer(default_S)*(-ATM_Put + BSM_Put_ATM_Plus);

constexpr long double ATM_BL_AoN_Call_super = default_S*ATM_BL_CoN_Call_super + ATM_Call;
constexpr long double ATM_BL_AoN_Call_sub = default_S*ATM_BL_CoN_Call_sub + ATM_Call;
constexpr long double ATM_BL_AoN_Put_super = default_S*ATM_BL_CoN_Put_super - ATM_Put;
constexpr long double ATM_BL_AoN_Put_sub = default_S*ATM_BL_CoN_Put_sub - ATM_Put;

// Breeden Litzenberger Replicating Portfolios for Binary Options with ATM Vol Skew
#define BL_CoN_Call_sub(S, K, r, y, T) BSM_Call(S, K*(static_cast<long double>(1-strikeGranularity)), r, y, skewedVol(S, K*(1-strikeGranularity)), 1) - BSM_Call(S, K, r, y, skewedVol(S, K), 1)
#define BL_CoN_Call_super(S, K, r, y, T) normalizer(S)*(BSM_Call(S, K, r, y, skewedVol(static_cast<long double>(S), static_cast<long double>(K)), 1) - BSM_Call(S, K*(1+strikeGranularity), r, y, skewedVol(S, K*(1+strikeGranularity)), 1))
#define BL_CoN_Put_super(S, K, r, y, T) normalizer(S)*(-BSM_Put(S, K*(1-strikeGranularity), r, y, skewedVol(static_cast<long double>(S), static_cast<long double>(K*(1-strikeGranularity))), 1) + BSM_Put(S, K, r, y, skewedVol(static_cast<long double>(S), static_cast<long double>(K)), 1))
#define BL_CoN_Put_sub(S, K, r, y, T) normalizer(S)*(-BSM_Put(S, K, r, y, skewedVol(static_cast<long double>(S), static_cast<long double>(K)), 1) + BSM_Put(S, K*(1+strikeGranularity), r, y, skewedVol(static_cast<long double>(S), static_cast<long double>(K*(1+strikeGranularity))), 1))

#define BL_AoN_Call_super(S, K, r, y, T) K*BL_CoN_Call_super(S, K, r, y, T) + BSM_Call(S, K, r, y, skewedVol(S, K), 1)
#define BL_AoN_Call_sub(S, K, r, y, T) K*BL_CoN_Call_sub(S, K, r, y, T) + BSM_Call(S, K, r, y, skewedVol(S, K), 1)
#define BL_AoN_Put_super(S, K, r, y, T) K*BL_CoN_Put_super(S, K, r, y, T) - BSM_Put(S, K, r, y, skewedVol(S, K), 1)
#define BL_AoN_Put_sub(S, K, r, y, T) K*BL_CoN_Put_sub(S, K, r, y, T) - BSM_Put(S, K, r, y, skewedVol(S, K), 1)

// Constexpr versions for ATM +/- chosen Granularity of Traded Strikes with Volatility Skew
constexpr long double ATM_BL_CoN_Call_super_with_vol_skew = normalizer(default_S)*(BSM_Call(default_S, default_S*(1-strikeGranularity), default_r, default_y, skewedVolATMMinus, 1) - BSM_Call(default_S, default_S, default_r, default_y, volATM, 1));
constexpr long double ATM_BL_CoN_Call_sub_with_vol_skew = normalizer(default_S)*(BSM_Call(default_S, default_S, default_r, default_y, volATM, 1) - BSM_Call(default_S, default_S*(1+strikeGranularity), default_r, default_y, skewedVolATMPlus, 1));
constexpr long double ATM_BL_CoN_Put_sub_with_vol_skew = normalizer(default_S)*(-BSM_Put(default_S, default_S*(1-strikeGranularity), default_r, default_y, skewedVolATMMinus, 1) + BSM_Put(default_S, default_S, default_r, default_y, volATM, 1));
constexpr long double ATM_BL_CoN_Put_super_with_vol_skew = normalizer(default_S)*(-BSM_Put(default_S, default_S, default_r, default_y, vol(default_S, default_S), 1) + BSM_Put(default_S, default_S*(1+strikeGranularity), default_r, default_y, skewedVolATMPlus, 1));

constexpr long double ATM_BL_AoN_Call_super_with_vol_skew = default_S*ATM_BL_CoN_Call_super_with_vol_skew + BSM_Call(default_S, default_S, default_r, default_y, volATM, 1);
constexpr long double ATM_BL_AoN_Call_sub_with_vol_skew = default_S*ATM_BL_CoN_Call_sub_with_vol_skew + BSM_Call(default_S, default_S, default_r, default_y, volATM, 1);
constexpr long double ATM_BL_AoN_Put_super_with_vol_skew = default_S*ATM_BL_CoN_Put_super_with_vol_skew - BSM_Put(default_S, default_S, default_r, default_y, volATM, 1);
constexpr long double ATM_BL_AoN_Put_sub_with_vol_skew = default_S*ATM_BL_CoN_Put_sub_with_vol_skew - BSM_Put(default_S, default_S, default_r, default_y, volATM, 1);





#endif
