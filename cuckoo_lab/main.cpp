//
// Created by Torin Tung Kwok on 2022/11/20.
//

#include <libcuckoo/cuckoohash_map.hh>

using namespace libcuckoo;

int main()
{
    cuckoohash_map<size_t, double> hash_table;
    {
        auto lt = hash_table.lock_table();
        std::cout << lt.insert(526).first->second << std::endl;
        std::cout << lt.insert(526, 99).first->second << std::endl;
    }
    // hash_table.update_fn()
    hash_table.upsert(526, [](auto &n) { n += 1; });
    std::cout << hash_table.find(526) << std::endl;

    hash_table.upsert(1024, [](auto &n) { n *= 10; });
    hash_table.upsert(1024, [](auto &n) { n *= 10; });
    std::cout << hash_table.find(1024) << std::endl;

    double found;
    // hash_table.find(2048, found);
    hash_table.insert(2048, 985);

    // NOTE(kwok): Deadlock!
    // hash_table.find_fn(2048, [&](const auto &found)
    //                    {
    //                        hash_table.template update_fn(2048, [](auto &found)
    //                                                      {
    //                                                          found = 4096;
    //                                                      }
    //                        );
    //                        std::cout << "🥥" << std::endl;
    //                    }
    // );

    hash_table.template update_fn(2048, [](auto &found)
                                  {
                                      found = 4096;
                                      std::cout << "🥥" << std::endl;
                                  }
    );

    std::cout << hash_table.find(2048) << std::endl;

    auto result = hash_table.find_fn(2048, [](const auto &n)
                                     {
                                         std::cout << "found! " << std::to_string(n) << std::endl;
                                     }
    );
    std::cout << result << std::endl;

    cuckoohash_map<size_t, std::vector<int>> h;
    std::vector<int> vec{14, 4};
    vec.reserve(4);
    h.insert(415, std::vector{5});

    for (auto n: h.find(415)) {
        std::cout << n << std::endl;
    }

    h.update(1, std::vector{15, 5});
    for (auto n: h.find(1)) {
        std::cout << n << std::endl;
    }

    vec[3] = 4;
    for (auto n: vec) {
        std::cout << n << std::endl;
    }

    return 0;
}