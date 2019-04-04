#define BOOST_TEST_MODULE "test_static_string"
#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <mjolnir/util/static_string.hpp>

BOOST_AUTO_TEST_CASE(test_construction)
{
    {
        mjolnir::static_string<8> str;
        BOOST_TEST(str.empty());
        BOOST_TEST(str.size()     == 0);
        BOOST_TEST(str.length()   == 0);
        BOOST_TEST(str.capacity() == 7);
        BOOST_TEST(str.max_size() == 7);
    }
    {
        mjolnir::static_string<8> str("hello");
        BOOST_TEST(!str.empty());
        BOOST_TEST(str.size()     == 5);
        BOOST_TEST(str.length()   == 5);
        BOOST_TEST(str.capacity() == 7);
        BOOST_TEST(str.max_size() == 7);
    }

    {
        mjolnir::static_string<8> str(std::string("hello"));
        BOOST_TEST(!str.empty());
        BOOST_TEST(str.size()     == 5);
        BOOST_TEST(str.length()   == 5);
        BOOST_TEST(str.capacity() == 7);
        BOOST_TEST(str.max_size() == 7);
    }

    {
        mjolnir::static_string<8> str(mjolnir::static_string<6>("hello"));
        BOOST_TEST(!str.empty());
        BOOST_TEST(str.size()     == 5);
        BOOST_TEST(str.length()   == 5);
        BOOST_TEST(str.capacity() == 7);
        BOOST_TEST(str.max_size() == 7);
    }

    {
        bool thrown = false;
        try
        {
            mjolnir::static_string<8> str("hello, world!");
        }
        catch(std::length_error)
        {
            thrown = true;
        }
        BOOST_TEST(thrown);
    }
}

BOOST_AUTO_TEST_CASE(test_substitution)
{
    mjolnir::static_string<8> str;
    BOOST_TEST(str.empty());
    BOOST_TEST(str.size()     == 0);
    BOOST_TEST(str.length()   == 0);
    BOOST_TEST(str.capacity() == 7);
    BOOST_TEST(str.max_size() == 7);

    str = "hello";

    BOOST_TEST(!str.empty());
    BOOST_TEST(str.size()     == 5);
    BOOST_TEST(str.length()   == 5);
    BOOST_TEST(str.capacity() == 7);
    BOOST_TEST(str.max_size() == 7);

    str = std::string("hello");
    BOOST_TEST(!str.empty());
    BOOST_TEST(str.size()     == 5);
    BOOST_TEST(str.length()   == 5);
    BOOST_TEST(str.capacity() == 7);
    BOOST_TEST(str.max_size() == 7);

    str = mjolnir::static_string<6>("hello");
    BOOST_TEST(!str.empty());
    BOOST_TEST(str.size()     == 5);
    BOOST_TEST(str.length()   == 5);
    BOOST_TEST(str.capacity() == 7);
    BOOST_TEST(str.max_size() == 7);
}

BOOST_AUTO_TEST_CASE(test_access)
{
    {
        mjolnir::static_string<8> str("hello");
        BOOST_TEST(str[0] == 'h');
        BOOST_TEST(str[1] == 'e');
        BOOST_TEST(str[2] == 'l');
        BOOST_TEST(str[3] == 'l');
        BOOST_TEST(str[4] == 'o');

        BOOST_TEST(str.at(0) == 'h');
        BOOST_TEST(str.at(1) == 'e');
        BOOST_TEST(str.at(2) == 'l');
        BOOST_TEST(str.at(3) == 'l');
        BOOST_TEST(str.at(4) == 'o');
    }
    {
        mjolnir::static_string<8> str(std::string("hello"));
        BOOST_TEST(str[0] == 'h');
        BOOST_TEST(str[1] == 'e');
        BOOST_TEST(str[2] == 'l');
        BOOST_TEST(str[3] == 'l');
        BOOST_TEST(str[4] == 'o');

        BOOST_TEST(str.at(0) == 'h');
        BOOST_TEST(str.at(1) == 'e');
        BOOST_TEST(str.at(2) == 'l');
        BOOST_TEST(str.at(3) == 'l');
        BOOST_TEST(str.at(4) == 'o');
    }
    {
        mjolnir::static_string<8> str(mjolnir::static_string<6>("hello"));
        BOOST_TEST(str[0] == 'h');
        BOOST_TEST(str[1] == 'e');
        BOOST_TEST(str[2] == 'l');
        BOOST_TEST(str[3] == 'l');
        BOOST_TEST(str[4] == 'o');

        BOOST_TEST(str.at(0) == 'h');
        BOOST_TEST(str.at(1) == 'e');
        BOOST_TEST(str.at(2) == 'l');
        BOOST_TEST(str.at(3) == 'l');
        BOOST_TEST(str.at(4) == 'o');
    }
    {
        mjolnir::static_string<8> str("hello");
        str[0] = 'w';
        str[1] = 'o';
        str[2] = 'r';
        str[3] = 'l';
        str[4] = 'd';

        BOOST_TEST(str[0] == 'w');
        BOOST_TEST(str[1] == 'o');
        BOOST_TEST(str[2] == 'r');
        BOOST_TEST(str[3] == 'l');
        BOOST_TEST(str[4] == 'd');
    }
    {
        mjolnir::static_string<8> str("hello");
        str.at(0) = 'w';
        str.at(1) = 'o';
        str.at(2) = 'r';
        str.at(3) = 'l';
        str.at(4) = 'd';

        BOOST_TEST(str.at(0) == 'w');
        BOOST_TEST(str.at(1) == 'o');
        BOOST_TEST(str.at(2) == 'r');
        BOOST_TEST(str.at(3) == 'l');
        BOOST_TEST(str.at(4) == 'd');
    }
}

BOOST_AUTO_TEST_CASE(test_comparison)
{
    {
        mjolnir::static_string<8> str("hello");
        BOOST_TEST(str == str);
        BOOST_TEST(str == "hello");
        BOOST_TEST(str == std::string("hello"));
    }
    {
        mjolnir::static_string<8> str1("hello");
        mjolnir::static_string<8> str2("world");
        BOOST_TEST(str1 != str2);
        BOOST_TEST(str1 != "world");
        BOOST_TEST(str1 != std::string("world"));

        BOOST_TEST(str1 < str2);
        BOOST_TEST(str1 < "world");
        BOOST_TEST(str1 < std::string("world"));
        BOOST_TEST(str1 <= str2);
        BOOST_TEST(str1 <= "world");
        BOOST_TEST(str1 <= std::string("world"));

        BOOST_TEST(!(str1 > str2));
        BOOST_TEST(!(str1 > "world"));
        BOOST_TEST(!(str1 > std::string("world")));
        BOOST_TEST(!(str1 >= str2));
        BOOST_TEST(!(str1 >= "world"));
        BOOST_TEST(!(str1 >= std::string("world")));
    }
}

BOOST_AUTO_TEST_CASE(test_concatenation)
{
    {
        mjolnir::static_string<16> str("h");
        str += 'e';
        BOOST_TEST(str == "he");
        str += "ll";
        BOOST_TEST(str == "hell");
        str += std::string("o ");
        BOOST_TEST(str == "hello ");
        str += mjolnir::static_string<6>("world");
        BOOST_TEST(str == "hello world");
    }
    {
        mjolnir::static_string<8> str1("hello "), str2("world");
        BOOST_TEST(str1 + str2 == "hello world");
        BOOST_TEST((str1 + str2).capacity() == 15);
    }
    {
        std::string str1("hello");
        mjolnir::static_string<8> str2("world");
        BOOST_TEST(str1 + str2 == "helloworld");
        BOOST_TEST(str2 + str1 == "worldhello");
        str1 += str2;
        BOOST_TEST(str1 == "helloworld");
    }
    {
        mjolnir::static_string<8> str1("hello");
        BOOST_TEST(str1 + "world" == "helloworld");
        BOOST_TEST("world" + str1 == "worldhello");
    }
}
