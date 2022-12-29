/**
 * @file       zeval.hpp
 * 
 * @author     Alex Skrynnyk (github.com/skrynnyk)
 *
 * @brief      This is a header-only implementation of a utility for evaluation
 *             of Z-domain Transfer Functions. It makes it trivial to set up and 
 *             evaluate the difference-equation form of discrete TFs, simply
 *             requiring the user to specify the numerator and denominator
 *             coefficients.
 *             
 * @detail     At its core, all that Zeval does is take any discrete transfer
 *             function specified in the follwing form
 * 
 * 
 *                          b[0] + b[1] * z^-1 + ... + b[i] * z^-i
 *                 G(z) = ------------------------------------------
 *                             1 + a[0] * z^-1 + ... + a[j] * z^-j
 * 
 *                        
 *             and stood it up and evaluate it in its difference equation form, 
 *             where `Y` is output, and `U` is input
 * 
 *               
 *                 Y(n) =   b[0] * U(n)   + ... + b[i] * U(n-i) 
 *                        - a[0] * Y(n-1) - ... - a[j] * Y(n-j-1)
 * 
 * 
 * @copyright  Copyright 2021 Alex Skrynnyk. All rights reserved.
 * 
 * @license    This project is released under the MIT License.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in 
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/* Includes ------------------------------------------------------------------*/
#include <array>
#include <cstdint>
#include <type_traits>
#include <tuple>

/* Private typedefs ----------------------------------------------------------*/
namespace zeval {
namespace detail {

template <size_t ...I>
struct index_sequence {};

template <size_t N, size_t ...I>
struct make_index_sequence : public make_index_sequence<N - 1, N - 1, I...> {};

template <size_t ...I>
struct make_index_sequence<0, I...> : public index_sequence<I...> {};

/** A wrapper type for coefficient values, because currently C++ does not allow
 *  floating point Non-Type Template Parameters */
struct C {
    float val;
    
    constexpr explicit C(const float& v) 
        : val{ v }
    {}
};

/** Tag type for evaluating the coefficient sign WITHOUT negation */
struct pos_t {};

/** Tag type for evaluating the coefficient sign WITH negation */
struct neg_t {};

/** Tage type for verifying DTF template parameters as `*CSet` classes */
struct cset_t {};

/* Public typedefs -----------------------------------------------------------*/
/**
 * @brief      This class represents a set (i.e. numerator or denominator) of 
 *             distinct DTF coefficients that, given corresponding inputs, can 
 *             compute its own value.
 *
 * @tparam     Cs    Set of individual coefficients, each wrapped by a `C` type
 */
template<typename... Cs>
class HeterogeneousCSet : detail::cset_t {
public:
    /** Number of coefficients in this set */
    constexpr static size_t size = sizeof...(Cs);

public:
    /** Container of individual coefficients */
    std::tuple<decltype(Cs::val)...> _cs;
    // TODO: see if tuple is really necessary as storage (AS 7/12/2021)

public:
    constexpr HeterogeneousCSet(Cs... cs)
        : _cs{ cs.val... }
    {}

    constexpr float operator()(const std::array<float, size>& ins,  pos_t sign)
    {
        return compute(_cs, ins, sign);
    } 

    constexpr float operator()(const std::array<float, size>& ins, neg_t sign)
    {
        return compute(_cs, ins, sign);
    } 

private:

    /* Compile-time recursion for computing coefficient set output */

    template<typename Tuple, typename Sign>
    constexpr float compute(const Tuple& cs,
                            const std::array<float, size> ins,
                            Sign sign)
    {
        using tup_size_t = typename std::tuple_size<Tuple>;
        return compute(cs, detail::make_index_sequence<tup_size_t::value>{}, 
                       ins, sign);
    }

    template<typename Tuple, size_t ... I>
    constexpr float compute(const Tuple& cs, detail::index_sequence<I ...>, 
                            const std::array<float, size> ins, pos_t)
    {
        return ((std::get<I>(cs) * ins[I]) + ...);
    }

    template<typename Tuple, size_t ... I>
    constexpr float compute(const Tuple& cs, detail::index_sequence<I ...>, 
                            const std::array<float, size> ins, neg_t)
    {
        return ((-1.0 * std::get<I>(cs) * ins[I]) + ...);
    }
};

/**
 * @brief      This class represents a set (i.e. numerator or denominator) of 
 *             distinct DTF coefficients that, given corresponding inputs, can 
 *             compute its own value.
 *
 * @tparam     Cs    Set of individual coefficients, each wrapped by a `C` type
 */
template<std::size_t N>
class RuntimeHeterogeneousCSet : detail::cset_t {
public:
    /** Number of coefficients in this set */
    constexpr static size_t size = N;

private:
    /** Container of individual coefficients */
    std::array<float, N> _cs;
    // TODO: see if tuple is really necessary as storage (AS 7/12/2021)

public:
    constexpr RuntimeHeterogeneousCSet(const std::array<float, N>& cs)
        : _cs{ cs }
    {}

    constexpr float operator()(const std::array<float, size>& ins, 
                               detail::pos_t sign)
    {
        float val = 0;
        for (int i = 0; i < ins.size(); i++) {
            val += _cs[i] * ins[i];
        }
        return val;
    }

    constexpr float operator()(const std::array<float, size>& ins, 
                               detail::neg_t sign)
    {
        float val = 0;
        for (int i = 0; i < ins.size(); i++) {
            val += -1.0f * _cs[i] * ins[i];
        }
        return val;
    }
};

/**
 * @brief      This class represents a set (i.e. numerator or denominator) of a 
 *             singular, repeated DTF  coefficient that, given corresponding
 *             inputs, can compute its own value.
 *             
 * @detail     This class serves as a space-optimization -- if the coefficient 
 *             is the same for every input, then why store multiple copies of 
 *             the same coefficient value? :)
 *
 * @tparam     Cs    Set of individual coefficients, each wrapped by a `C` type
 */
template<typename... Cs>
class HomogeneousCSet : detail::cset_t {
public:
    /** Number of coefficients in this set */
    constexpr static size_t size = sizeof...(Cs);

// private:
    /** The singluar coefficient */
    float _c;

public:
    constexpr HomogeneousCSet(Cs... cs)
        : _c( std::get<0>(std::tuple<Cs...>{ cs... }).val )
    {}

    constexpr float operator()(const std::array<float, size>& ins, pos_t sign)
    {
        return compute<sizeof...(Cs)>(_c, ins, sign);
    }

    constexpr float operator()(const std::array<float, size>& ins, neg_t sign)
    {
        return compute<sizeof...(Cs)>(_c, ins, sign);
    }

private:

    /* Compile-time recursion for computing coefficient set output */

    template<size_t Num, typename Sign>
    constexpr float compute(const float& c,
                            const std::array<float, size> ins,
                            Sign sign)
    {
        return compute(c, detail::make_index_sequence<Num>{}, 
                       ins, sign);
    }

    template<size_t... I>
    constexpr float compute(const float& c, detail::index_sequence<I...>, 
                            const std::array<float, size> ins, pos_t)
    {
        return ((c * ins[I]) + ...);
    }

    template<size_t... I>
    constexpr float compute(const float& c, detail::index_sequence<I...>, 
                            const std::array<float, size> ins, neg_t)
    {
        return ((-1.0 * c * ins[I]) + ...);
    }
};

/**
 * @brief      This class represents a set (i.e. numerator or denominator) of a 
 *             singular, repeated DTF  coefficient that, given corresponding
 *             inputs, can compute its own value.
 *             
 * @detail     This class serves as a space-optimization -- if the coefficient 
 *             is the same for every input, then why store multiple copies of 
 *             the same coefficient value? :)
 *
 * @tparam     Cs    Set of individual coefficients, each wrapped by a `C` type
 */
// template<typename... Cs>
template<std::size_t N>
class RuntimeHomogeneousCSet : detail::cset_t {
public:
    /** Number of coefficients in this set */
    constexpr static size_t size = N;

private:
    /** The singluar coefficient */
    float _c;

public:
    constexpr RuntimeHomogeneousCSet(const std::array<float, N>& cs)
        : _c( cs[0] )
    {}

    constexpr float operator()(const std::array<float, size>& ins, 
                               detail::pos_t sign)
    {
        float val = 0;
        for (const auto& in : ins) {
            val += _c * in;
        }
        return val;
    }

    constexpr float operator()(const std::array<float, size>& ins, 
                               detail::neg_t sign)
    {
        float val = 0;
        for (const auto& in : ins) {
            val += -1.0f * _c * in;
        }
        return val;
    }
};
} /* namespace detail */

/**
 * @brief      This class implements the actual DTF in difference equation form
 *
 * @tparam     TAsPolicy    Policy for the As coefficient set
 * @tparam     TBsPolicy    Policy for the Bs coefficient set
 * @tparam     TTickPolicy  Policy for getting system tick time
 */
template<class TAsPolicy, class TBsPolicy, class TTickPolicy>
class DTF {
    static_assert(std::is_base_of<detail::cset_t, TAsPolicy>::value &&
                  std::is_base_of<detail::cset_t, TBsPolicy>::value, 
                  "Not a `HeterogeneousCSet` nor a `HomogeneousCSet` policy");
public:
    // TODO: make type specified via class template
    using Ys_t = std::array<float, TAsPolicy::size>;
    using Us_t = std::array<float, TBsPolicy::size>;
    
    /** DTF sampling interval in ticks (tick units must be managed by user) */
    const uint16_t Ts = 1;

public:
    /** A coefficient set (applies to TF outputs, i.e. Ys) */
    TAsPolicy _As;

    /** B coefficient set (applies to TF inputs, i.e. Us) */
    TBsPolicy _Bs;
    
    /** Last tick at which the TF was evaluated */
    uint32_t _last_tick = 0;

    /** Storage for DTF outputs */
    Ys_t _Ys{ 0 };

    /** Storage for DTF inputs */
    Us_t _Us{ 0 };

public:
    constexpr DTF(const TBsPolicy& bs,const TAsPolicy& as,  uint16_t ts)
        : _As{ as }
        , _Bs{ bs }
        , Ts{ ts }
    {
        // TODO: make C/Y/U type derived from template parameter (AS 7/11/2021)
    }

    /**
     * @brief      Function call operator evaluates the TF given a new input `U`
     * 
     * @details    Uses the `TTickPolicy` class to assess whether or not the next
     *             evaluation interval is met, returns the previous output `Y` 
     *             if not (i.e. Zero-order Hold)
     *
     * @param[in]  U     New latest TF input
     *
     * @return     The result of DTF evaluation
     */
    constexpr float operator()(float U) {
        // Zero-order Hold
        // TODO: make policy-configurable? (AS 7/11/21)
        uint32_t tick = TTickPolicy::getTick();
        if (tick < _last_tick + Ts) { return _Ys[0]; }

        addPoint(U, _Us);
        float Y = _Bs(_Us, detail::pos_t{}) + _As(_Ys, detail::neg_t{});

        addPoint(Y, _Ys);

        _last_tick = tick;
        return Y;
    }

    /**
     * @brief      Returns the latest evaluation value of the TF
     *
     * @return     latest computed output of this TF object
     */
    constexpr float value() const {
        return _Ys[0];
    }

    /**
     * @brief      Resets the state (history of inputs/outputs) of the TF
     */
    constexpr void reset() {
        for (auto& y : _Ys) { y = 0; }
        for (auto& u : _Us) { u = 0; }
    }

    /**
     * @brief      Sets the initial conditions.
     * 
     * @note       Take care to set only when not actively evaluating new inputs
     *
     * @param[in]  Us    The new value of input records
     * @param[in]  Ys    The new value of output records
     */
    constexpr void setInitialConditions(const Us_t& Us, const Ys_t& Ys) {
        for (int i = 0; i < Us.size(); i++) { _Us[i] = Us[i]; }
        for (int i = 0; i < Ys.size(); i++) { _Ys[i] = Ys[i]; }
    }

    /**
     * @brief      Updates the coefficient values
     * 
     * @note       Currently also resets inteernal state to guanratee stability
     *
     * @param[in]  as    The new As coefficients
     * @param[in]  bs    The new Bs coefficients
     */
    constexpr void setCs(const TAsPolicy& as, const TBsPolicy& bs) {
        _As = as;
        _Bs = bs;
        reset();
    }

private:
    /**
     * @brief      Adds a new value to Ys or Us, shifting older values by 1 to
     *             the right in order to maintain same ordered association of
     *             values with their corresponding coefficients
     *
     * @param[in]  val      The value
     * @param      storage  The storage
     *
     * @tparam     N        Number of values in storage
     */
    template<size_t N>
    static constexpr void addPoint(float val, std::array<float, N>& storage) {
        float temp = storage[N - 1], temp1 = 0;
        for (size_t i = 0; i < N; i++) { 
            temp1 = storage[i];
            storage[i] = temp;
            temp = temp1;
        }
        storage[0] = val;
    }
};

/* Public functions ----------------------------------------------------------*/
/**
 * @brief      Helper for creating a discrete TF object
 *
 * @param[in]  As           The A coefficients
 * @param[in]  Bs           The B coefficients
 * @param[in]  <unnamed>    System tick getter policy object
 * @param[in]  Ts           Discretized TF sampling time
 *
 * @tparam     TAsPolicy    Policy for the As coefficient set
 * @tparam     TBsPolicy    Policy for the Bs coefficient set
 * @tparam     TTickPolicy  Policy for getting system tick time
 *
 * @return     The usable TF object
 */
template<class TAsPolicy, class TBsPolicy, class TTickPolicy>
constexpr auto makeDTF(const TBsPolicy& Bs, const TAsPolicy& As, TTickPolicy,
                       uint16_t Ts)
    -> DTF<TAsPolicy, TBsPolicy, TTickPolicy>
{
    // TODO: consume the first A coef as 1 to follow convention (AS 7/12/2021) 
    return { Bs, As, Ts };
}

/**
 * @brief      Makes a DTF coefficient set with unique values
 *
 * @param[in]  cs    The set of coefficient values
 *
 * @tparam     TCs   Type of the individual coefficients
 *
 * @return     A `HeterogeneousCSet` containing the specified coefs
 */
template<class... TCs>
constexpr auto makeCs(const TCs&... cs)
    -> detail::HeterogeneousCSet<decltype(detail::C(TCs{}))...>
{
    using cs_t = typename std::tuple_element<0, std::tuple<TCs...>>::type;
    static_assert(std::is_same<float, cs_t>::value, "Invalid C type");
    return { detail::C{ cs }... };
}

/**
 * @brief      Makes a DTF coefficient set composed of the same value
 *
 * @param[in]  cs    The set of coefficient values
 *
 * @tparam     TCs   Type of the individual coefficients
 *
 * @return     A `HomogeneousCSet` containing the specified coefs
 */
template<class... TCs>
constexpr auto makeHomogeneousCs(const TCs&... cs)
    -> detail::HomogeneousCSet<decltype(detail::C(TCs{}))...>
{
    using cs_t = typename std::tuple_element<0, std::tuple<TCs...>>::type;
    static_assert(std::is_same<float, cs_t>::value, "Invalid C type");
    return { detail::C{ cs }... };
}

namespace util {
template<uint32_t Tsize, size_t ... I>
constexpr auto makeMovingAvg(::zeval::detail::index_sequence<I ...>)
    -> decltype(makeHomogeonousCs((1.0f / Tsize * (I / I))...))
{
    return makeHomogeonousCs(((void) I, (1.0f / Tsize))...);
}

template<std::uint32_t Tsize, size_t ... I>
constexpr auto makeDelay(::zeval::detail::index_sequence<I ...>)
    -> decltype(makeCs(((Tsize - I > 1) ? 0.0f : 1.0f)...))
{
    return makeCs(((Tsize - I > 1) ? 0.0f : 1.0f)...);
}

/**
 * @brief      Makes a moving average (FIR) low pass filter\
 *
 * @param[in]  p            System tick getter policy object
 * @param[in]  Ts           Discretized TF sampling time
 *
 * @tparam     Tsize        Size of the moving average window
 * @tparam     TTickPolicy  Policy for getting system tick time
 *
 * @return     The usable TF object
 */
template<uint32_t Tsize, class TTickPolicy>
constexpr auto makeMovingAvg(TTickPolicy p, uint16_t Ts)
    //-> decltype(makeDTF(makeCs(0.0f), BCs, p, Ts))
{
    auto BCs = makeMovingAvg<Tsize>(
        ::zeval::detail::make_index_sequence<Tsize>{}
    );
    return makeDTF(makeCs(0.0f), BCs, p, Ts);
}

/**
 * @brief      Makes a simple single pole IIR low pass filter
 *
 * @param[in]  alpha        The alpha value ("smoothing factor")
 * @param[in]  p            System tick getter policy object
 * @param[in]  Ts           Discretized TF sampling time
 *
 * @tparam     TTickPolicy  Policy for getting system tick time
 *
 * @return     The usable TF object
 */
template<class TTickPolicy>
constexpr auto makeSinglePoleIIR(float alpha, TTickPolicy p, uint16_t Ts)
    -> decltype(makeDTF(makeCs(-1 * alpha), makeCs(1 - alpha), p, Ts))
{
    return makeDTF(makeCs(-1 * alpha), makeCs(1 - alpha), p, Ts);
}

/**
 * @brief      Makes a discrete differentiator
 *
 * @param[in]  p            System tick getter policy object
 * @param[in]  Ts           Discretized TF sampling time
 *
 * @tparam     TTickPolicy  Policy for getting system tick time
 *
 * @return     The usable TF object
 */
template<class TTickPolicy>
constexpr auto makeDifferentiator(TTickPolicy p, uint16_t Ts)
    -> decltype(makeDTF(makeCs(0.0f), makeCs(1.0f, -1.0f), p, Ts))
{
    return makeDTF(makeCs(0.0f), makeCs(1.0f, -1.0f), p, Ts); // TODO
}

/**
 * @brief      Makes a discrete rectangular integrator.
 *
 * @param[in]  p            System tick getter policy object
 * @param[in]  Ts           Discretized TF sampling time
 *
 * @tparam     TTickPolicy  Policy for getting system tick time
 *
 * @return     The usable TF object
 */
template<class TTickPolicy>
constexpr auto makeIntegrator(TTickPolicy p, uint16_t Ts)
    -> decltype(makeDTF(makeCs(-1.0f), makeCs(1.0f), p, Ts))
{
    // Rectangular integrator
    return makeDTF(makeCs(-1.0f), makeCs(1.0f), p, Ts); // TODO
}

/**
 * @brief      Makes a delay line filter
 *
 * @param[in]  p            System tick getter policy object
 * @param[in]  Ts           Discretized TF sampling time
 *
 * @tparam     Tsize        Length of the delay
 * @tparam     TTickPolicy  Policy for getting system tick time
 *
 * @return     The usable TF object
 */
template<std::uint32_t Tsize, class TTickPolicy>
constexpr auto makeDelay(TTickPolicy p, uint16_t Ts)
    // -> decltype(makeDTF(makeCs(0.0f), makeCs(0.0f, 0.0f, 0.0f), p, Ts))
{
    // Constant delay
    auto BCs = makeDelay<Tsize>(
        ::zeval::detail::make_index_sequence<Tsize>{}
    );
    return makeDTF(makeCs(0.0f), BCs, p, Ts); // TODO
}
} /* namespace util */
} /* namespace zeval */
