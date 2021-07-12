/**
 * @defgroup   DTFEVAL dtfeval
 *
 * @brief      This file implements the Discrete Transfer Function (DTF)
 *             Evaluating utility class `dtfe::DTF` that sets up the
 *             difference-equation form of a DTF given its numerator and
 *             denominator coefficients.
 *
 * @author     Alex Skrynnyk
 * @date       2021
 */

/* Includes ------------------------------------------------------------------*/
#include <array>
#include <cstdint>
#include <type_traits>
#include <tuple>

/* Private typedefs ----------------------------------------------------------*/
namespace dtfe {
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

/** Tage type for verifying DTF template parameters as `CSet` classes */
struct cset_t {};

} /* namespace detail */

/* Public typedefs -----------------------------------------------------------*/
/**
 * @brief      This class represents a set of DTF coefficients (i.e. numerator 
 *             or denominator) that, given corresponding inputs, can compute its
 *             own value.
 *
 * @tparam     Cs    Set of individual coefficients, each wrapped by a `C` type
 */
template<typename... Cs>
class CSet : detail::cset_t {
public:
    /** Number of coefficients in this set */
    constexpr static std::size_t size = sizeof...(Cs);

private:
    /** Container of individual coefficients */
    std::tuple<decltype(Cs::val)...> _cs;

public:
    constexpr CSet(Cs... cs)
        : _cs{ cs.val... }
    {}

    constexpr float operator()(const std::array<float, size>& ins, 
                               detail::pos_t sign)
    {
        return compute(_cs, ins, sign);
    } 

    constexpr float operator()(const std::array<float, size>& ins,
                               detail::neg_t sign)
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
                            const std::array<float, size> ins,
                            detail::pos_t)
    {
        return ((std::get<I>(cs) * ins[I]) + ...);
    }

    template<typename Tuple, size_t ... I>
    constexpr float compute(const Tuple& cs, detail::index_sequence<I ...>, 
                            const std::array<float, size> ins,
                            detail::neg_t)
    {
        return ((-1.0 * std::get<I>(cs) * ins[I]) + ...);
    }
};

/**
 * @brief      This class implements the actual DTF in difference equation form
 *
 * @tparam     TAs          { description }
 * @tparam     TBs          { description }
 * @tparam     TickPolicy   { description }
 */
template<class TAs, class TBs, class TickPolicy>
class DTF {
    static_assert(std::is_base_of<detail::cset_t, TAs>::value, "Not a `CSet`");
    static_assert(std::is_base_of<detail::cset_t, TBs>::value, "Not a `CSet`");

private:
    TAs _a;
    TBs _b;
    
    /** DTF sampling interval in ticks (tick units must be managed by user) */
    const std::uint16_t _ts = 1;
    /** Last tick at which the TF was evaluated */
    std::uint32_t _last_tick = 0;

    /** Storage for DTF outputs */
    std::array<float, TAs::size> _Ys{ 0 };

    /** Storage for DTF inputs */
    std::array<float, TBs::size> _Us{ 0 };

public:
    constexpr DTF(const TAs& as, const TBs& bs, std::uint16_t ts)
        : _a{ as }
        , _b{ bs }
        , _ts{ ts }
    {
        // TODO: make C/Y/U type derived from template parameter (AS 7/11/2021)
    }

    /**
     * @brief      Function call operator evaluates the TF given a new input `U`
     * 
     * @details    Uses the `TickPolicy` class to assess whether or not the next
     *             evaluation interval is met, returns the previous output `Y` 
     *             if not (i.e. Zero-order Hold)
     *
     * @param[in]  U     New latest TF input
     *
     * @return     The result of DTF evaluation
     */
    float operator()(float U) {
        // Zero-order Hold
        // TODO: make policy-configurable? (AS 7/11/21)
        std::uint32_t tick = TickPolicy::getTick();
        if (tick < _last_tick + _ts) { return _Ys[0]; }

        add(U, _Us);
        float Y = _b(_Us, detail::pos_t{}) + _a(_Ys, detail::neg_t{});
        add(Y, _Ys);

        _last_tick = tick;
        return Y;
    }

private:
    /**
     * @brief      { function_description }
     *
     * @param[in]  val      The value
     * @param      storage  The storage
     *
     * @tparam     N        { description }
     */
    template<std::size_t N>
    static constexpr void add(float val, std::array<float, N>& storage) {
        float temp = storage[N - 1], temp1 = 0;
        for (std::size_t i = 0; i < N; i++) { 
            temp1 = storage[i];
            storage[i] = temp;
            temp = temp1;
        }
        storage[0] = val;
    }
};

/* Public functions ----------------------------------------------------------*/
/**
 * @brief      Makes a DTF
 *
 * @param[in]  As           { parameter_description }
 * @param[in]  Bs           { parameter_description }
 * @param[in]  <unnamed>    { parameter_description }
 * @param[in]  Ts           { parameter_description }
 *
 * @tparam     TAs          { description }
 * @tparam     TBs          { description }
 * @tparam     TickPolicy  { description }
 *
 * @return     { description_of_the_return_value }
 */
template<class TAs, class TBs, class TickPolicy>
constexpr auto makeDTF(const TAs& As, const TBs& Bs, TickPolicy,
                       std::uint16_t Ts) -> DTF<TAs, TBs, TickPolicy> {
    return { As, Bs, Ts };
}

/**
 * @brief      Makes coefficients for a DTF
 *
 * @param[in]  cs    The create struct
 *
 * @tparam     TCs   { description }
 *
 * @return     { description_of_the_return_value }
 */
template<class... TCs>
constexpr auto makeCs(const TCs&... cs) -> CSet<decltype(detail::C(TCs{}))...> {
    return { detail::C{ cs }... };
}

/**
 * @brief      Makes a moving average (FIR) filter DTF
 *
 * @param[in]  p            { parameter_description }
 * @param[in]  Ts           { parameter_description }
 *
 * @tparam     Tsize        { description }
 * @tparam     TickPolicy  { description }
 *
 * @return     { description_of_the_return_value }
 */
template<std::uint32_t Tsize, class TickPolicy>
constexpr auto makeMovingAvg(TickPolicy p, std::uint16_t Ts) {
    return makeDTF(makeCs(0.0f), makeCs(1 /* TODO */), p, Ts);
}

/**
 * @brief      Makes a single pole IIR filter DTF
 *
 * @param[in]  alpha        The alpha
 * @param[in]  p            { parameter_description }
 * @param[in]  Ts           { parameter_description }
 *
 * @tparam     Tsize        { description }
 * @tparam     TickPolicy  { description }
 *
 * @return     { description_of_the_return_value }
 */
template<std::uint32_t Tsize, class TickPolicy>
constexpr auto makeSinglePoleIIR(float alpha, TickPolicy p, std::uint16_t Ts) {
    return makeDTF(makeCs(-1 * alpha), makeCs(1 - alpha), p, Ts);
}

} /* namespace dtfe */
