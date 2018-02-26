#define BOOST_TEST_MODULE "test_index_sequence"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/util/index_sequence.hpp>
#include <type_traits>

BOOST_AUTO_TEST_CASE(test_make_index_sequence)
{
    {
        typedef mjolnir::index_sequence<0> answer;
        const bool make_index_seq_1_result =
            std::is_same<mjolnir::make_index_sequence<1>, answer>::value;
        BOOST_CHECK(make_index_seq_1_result);
    }
    {
        typedef mjolnir::index_sequence<0, 1> answer;
        const bool make_index_seq_2_result =
            std::is_same<mjolnir::make_index_sequence<2>, answer>::value;
        BOOST_CHECK(make_index_seq_2_result);
    }
    {
        typedef mjolnir::index_sequence<0, 1, 2> answer;
        const bool make_index_seq_3_result =
            std::is_same<mjolnir::make_index_sequence<3>, answer>::value;
        BOOST_CHECK(make_index_seq_3_result);
    }
    {
        typedef mjolnir::index_sequence<0, 1, 2, 3> answer;
        const bool make_index_seq_4_result =
            std::is_same<mjolnir::make_index_sequence<4>, answer>::value;
        BOOST_CHECK(make_index_seq_4_result);
    }
    {
        typedef mjolnir::index_sequence<0, 1, 2, 3, 4> answer;
        const bool make_index_seq_5_result =
            std::is_same<mjolnir::make_index_sequence<5>, answer>::value;
        BOOST_CHECK(make_index_seq_5_result);
    }
    {
        typedef mjolnir::index_sequence<0, 1, 2, 3, 4, 5> answer;
        const bool make_index_seq_6_result =
            std::is_same<mjolnir::make_index_sequence<6>, answer>::value;
        BOOST_CHECK(make_index_seq_6_result);
    }
    {
        typedef mjolnir::index_sequence<0, 1, 2, 3, 4, 5, 6> answer;
        const bool make_index_seq_7_result =
            std::is_same<mjolnir::make_index_sequence<7>, answer>::value;
        BOOST_CHECK(make_index_seq_7_result);
    }
    {
        typedef mjolnir::index_sequence<0, 1, 2, 3, 4, 5, 6, 7> answer;
        const bool make_index_seq_8_result =
            std::is_same<mjolnir::make_index_sequence<8>, answer>::value;
        BOOST_CHECK(make_index_seq_8_result);
    }
    {
        typedef mjolnir::index_sequence<0, 1, 2, 3, 4, 5, 6, 7, 8> answer;
        const bool make_index_seq_9_result =
            std::is_same<mjolnir::make_index_sequence<9>, answer>::value;
        BOOST_CHECK(make_index_seq_9_result);
    }
    {
        typedef mjolnir::index_sequence<0, 1, 2, 3, 4, 5, 6, 7, 8, 9> answer;
        const bool make_index_seq_10_result =
            std::is_same<mjolnir::make_index_sequence<10>, answer>::value;
        BOOST_CHECK(make_index_seq_10_result);
    }
}
