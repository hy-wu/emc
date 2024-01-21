#include <iostream>

int main() {
    double temp;
    long num_test;
    std::cout << "Enter the number of test cases: ";
    std::cin >> temp;
    num_test = static_cast<long>(temp);
    std::cout << "The number of test cases is " << num_test << std::endl;
    std::cout << "The type of num_test is " << typeid(num_test).name() << std::endl;
    return 0;
}
