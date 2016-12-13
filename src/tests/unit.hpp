#ifndef SIMPLE_UNIT_TEST_HPP_
#define SIMPLE_UNIT_TEST_HPP_
#include <functional>
#include <vector>
#include <string>
#include <iostream>

struct UnitTest {
    std::function<bool()> function;
    std::string description;
};

class TestRunner
{
    public:
        void runTests(void) {
            for (auto &test : tests_) {
                std::cout << test.description << ": " <<
                    ((test.function())?("pass"):("fail")) << '\n';
            }
        }
        void addTest(const std::function<bool()> &f, const std::string &desc) {
            tests_.emplace_back(UnitTest{f, desc});
        }
    private:
        std::vector<UnitTest> tests_;
};

extern TestRunner runner;

class Adder
{
    public:
        Adder(const std::function<bool()> &f, const std::string &desc) { runner.addTest(f, desc); }
};

#define UTEST_NAME2(f,l) f##l
#define UTEST_NAME(f,l) UTEST_NAME2(f,l)
#define UTEST(section_name, desc, fun) Adder UTEST_NAME(section_name, __LINE__)(fun, desc);

#endif
