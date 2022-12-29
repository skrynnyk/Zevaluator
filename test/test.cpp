#include <zeval/zeval.hpp>
#include <catch2/catch_all.hpp>
#include <cmath>

//==============================================================================
// Coefficient Set Tests
//==============================================================================

template<std::size_t NumCoefs>
struct CsetTestCase {
    std::array<float, NumCoefs> ins;
    float out;
};


TEST_CASE("Heterogeneous CSet")
{
    SECTION("Single C")
    {
        constexpr std::size_t k_num_coefs = 1;
        auto num     = zeval::makeCs(4.0f);
        using test_t = CsetTestCase<k_num_coefs>;
        auto cases   = std::array<test_t, 4>{
            test_t{ {0}, 0.0f },
            test_t{ {1}, 4.0f },
            test_t{ {2}, 8.0f }
        };

        REQUIRE(num.size == k_num_coefs);

        for (auto& c : cases) {
            REQUIRE(num(c.ins, zeval::detail::pos_t{}) == c.out);
            REQUIRE(num(c.ins, zeval::detail::neg_t{}) == -1.0 * c.out);
        }
    }

    SECTION("Multiple Cs")
    {
        constexpr std::size_t k_num_coefs = 2;
        auto num     = zeval::makeCs(4.0f, 2.5f);
        using test_t = CsetTestCase<k_num_coefs>;
        auto cases   = std::array<test_t, 4>{
            test_t{ {1, 0}, 4.0f },
            test_t{ {0, 1}, 2.5f },
            test_t{ {1, 1}, 6.5f },
            test_t{ {2, 3}, 15.5f }
        };

        REQUIRE(num.size == k_num_coefs);

        for (auto& c : cases) {
            REQUIRE(num(c.ins, zeval::detail::pos_t{}) == c.out);
            REQUIRE(num(c.ins, zeval::detail::neg_t{}) == -1.0 * c.out);
        }
    }
}


TEST_CASE("Homogeneous CSet")
{
    SECTION("Single C")
    {
        constexpr std::size_t k_num_coefs = 1;
        auto num     = zeval::makeHomogeneousCs(4.0f);
        using test_t = CsetTestCase<k_num_coefs>;
        auto cases   = std::array<test_t, 4>{
            test_t{ {0}, 0.0f },
            test_t{ {1}, 4.0f },
            test_t{ {2}, 8.0f }
        };

        REQUIRE(num._c == 4.0f);
        REQUIRE(num.size == k_num_coefs);

        for (auto& c : cases) {
            REQUIRE(num(c.ins, zeval::detail::pos_t{}) == c.out);
            REQUIRE(num(c.ins, zeval::detail::neg_t{}) == -1.0 * c.out);
        }
    }

    SECTION("Multiple Cs")
    {
        constexpr std::size_t k_num_coefs = 2;
        auto num     = zeval::makeHomogeneousCs(4.0f, 4.0f);
        using test_t = CsetTestCase<k_num_coefs>;
        auto cases   = std::array<test_t, 4>{
            test_t{ {1, 0}, 4.0f },
            test_t{ {0, 1}, 4.0f },
            test_t{ {1, 1}, 8.0f },
            test_t{ {2, 3}, 20.0f }
        };

        REQUIRE(num._c == 4.0f);
        REQUIRE(num.size == k_num_coefs);

        for (auto& c : cases) {
            REQUIRE(num(c.ins, zeval::detail::pos_t{}) == c.out);
            REQUIRE(num(c.ins, zeval::detail::neg_t{}) == -1.0 * c.out);
        }
    }
}

//==============================================================================
// Transfer Function Tests
//==============================================================================

struct MockTickPolicy {
    static uint32_t getTick() { static int i = 0; return i++; }
};

struct TFTestCase {
    float in;
    float out;
};

/* Reference -- must be correct by inspection */
class ReferenceImpl {
public:
    float a0 = 0, b0 = 0, b1 = 0;

    float yk_1 = 0;
    float uk_1 = 0;
    bool is_n1 = true;

    ReferenceImpl(float b0, float b1, float a0)
        : a0{ a0 }
        , b0{ b0 }
        , b1{ b1 }
    {}

    float compute(float reading) {
        float yk1 = 0;
        float uk1 = 0;

        if (is_n1) {
            yk1 = 0;
            uk1 = 0;
            is_n1 = false;
        } else {
            yk1 = yk_1;
            uk1 = uk_1;
        }

        float u = reading;
        uk_1 = u;
        float y = b0 * u + b1 * uk1 - a0 * yk1;
        yk_1 = y;
        return y;
    }

};


bool equal_within_perc(float ref, float test, float percent = 0.1) {
    return (std::fabs(ref - test) <= (ref * (percent / 100.0)));
}


TEST_CASE("DTF class")
{
    SECTION("Basic class operation")
    {
        constexpr std::size_t k_ts = 1;
        auto num = zeval::makeCs(4.0f);
        auto den = zeval::makeCs(1.0f);
        auto Gz  = zeval::makeDTF(num, den, MockTickPolicy{}, k_ts);

        REQUIRE(Gz.value() == 0.0f);
        
        (void) Gz(5);
        (void) Gz(2);
        REQUIRE(Gz.value() != 0.0f);
        
        Gz.reset();
        REQUIRE(Gz.value() == 0.0f);
    }

    SECTION("Unity TF output")
    {
        constexpr std::size_t k_num_coefs = 1;
        constexpr std::size_t k_ts = 1;
        auto num     = zeval::makeCs(1.0f);
        auto den     = zeval::makeCs(0.0f);
        auto Gz      = zeval::makeDTF(num, den, MockTickPolicy{}, k_ts);
        using test_t = TFTestCase;
        auto cases   = std::array<test_t, 5>{
            test_t{ 5.0f, 5.0f },
            test_t{ 5.0f, 5.0f },
            test_t{ 5.0f, 5.0f },
            test_t{ 5.0f, 5.0f },
            test_t{ 5.0f, 5.0f }
        };

        float out = 0.0f;
        for (auto& c : cases) {
            out = Gz(c.in);
            CAPTURE(out, c.out);
            REQUIRE(out == c.out);
        }
    }

    constexpr auto Bs = std::array<float, 2>{ 0.038607f, -0.035720f };
    constexpr auto As = std::array<float, 1>{ -0.996750f };

    auto Gz_ref   = ReferenceImpl{ Bs[0], Bs[1],
                                   As[0] };
    auto Gz_zeval = zeval::makeDTF(zeval::makeCs(Bs[0], Bs[1]),
                                   zeval::makeCs(As[0]),
                                   MockTickPolicy{}, 
                                   1);

    SECTION("Reference vs Zeval: Internal State")
    {
        constexpr std::size_t k_num_iters = 20;
        constexpr float       input       = 15;

        float output_ref = 0;
        float output_zeval = 0;
        for (int i = 0; i < k_num_iters; i++) {
            output_ref   = Gz_ref.compute(input);
            output_zeval = Gz_zeval(input);
            // CAPTURE(output_ref, output_zeval);
            REQUIRE(equal_within_perc(Gz_ref.yk_1, Gz_zeval._Ys[0]));
            REQUIRE(equal_within_perc(Gz_ref.uk_1, Gz_zeval._Us[0]));
            // CHECK(equal_within_perc(output_ref, output_zeval));
        }
    }

    SECTION("Reference vs Zeval: Step Response")
    {
        constexpr std::size_t k_num_iters = 2000;
        constexpr float       input       = 15;

        float output_ref = 0;
        float output_zeval = 0;
        for (int i = 0; i < k_num_iters; i++) {
            output_ref   = Gz_ref.compute(input);
            output_zeval = Gz_zeval(input);
            CAPTURE(output_ref, output_zeval);
            REQUIRE(equal_within_perc(output_ref, output_zeval));
        }
    }
}


//==============================================================================

int main( int argc, char* argv[] )
{
  Catch::Session session; // There must be exactly one instance

  // writing to session.configData() here sets defaults
  // this is the preferred way to set them

  int returnCode = session.applyCommandLine( argc, argv );
  if( returnCode != 0 ) // Indicates a command line error
        return returnCode;

  // writing to session.configData() or session.Config() here
  // overrides command line args
  // only do this if you know you need to

  int numFailed = session.run();

  // numFailed is clamped to 255 as some unices only use the lower 8 bits.
  // This clamping has already been applied, so just return it here
  // You can also do any post run clean-up here
  return numFailed;
}
