#include <random>
int test(){
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // std::uniform_real_distribution<double> dis(0, 1);
    for (int i = 0; i < 100; i++) {
    //     std::cout << dis(gen) << std::endl;
    // }
    std::cout << rand() << std::endl;
    }
    return 0;
}