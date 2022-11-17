#ifndef BULLDOG_MODULES_ENGINE_SRC_BUCKET_POOL_HPP_
#define BULLDOG_MODULES_ENGINE_SRC_BUCKET_POOL_HPP_

#include "bucket.h"
#include <map>
#include <array>
#include <filesystem>
#include <bulldog/core_util.hpp>

struct BucketMeta
{
    Bucket bucket_;
    int bucket_count_ = 0;
    int post_flop_bucket_count_ = 0;
};

class BucketPool
{
public:
    virtual ~BucketPool()
    {
        for (auto &it: bucket_pool_) {
            for (auto &lossless_it: it) {
                for (auto &[key, val]: lossless_it) {
                    delete val;
                }
            }
        }
    }

    BucketMeta *Get(const std::string &name, int round, bool lossless = false)
    {
        if (!Has(name, round, lossless)) {
            // Loading bucket lazily.
            bucket_pool_[round][lossless][name] = LoadBucket(name, lossless);
            logger::info("🚚cache lazily loaded for name=%s, round=%d, lossless=%d", name, round, lossless);
        } else {
            logger::info("✅cache hit for name=%s, round=%d, lossless=%d", name, round, lossless);
        }
        return bucket_pool_[round][lossless][name];
    }

    bool Has(const std::string &name, int round, bool lossless)
    {
        auto it = bucket_pool_[round][lossless].find(name);
        return it != bucket_pool_[round][lossless].end();
    }

    static BucketMeta *LoadBucket(const std::string &name, bool lossless = false)
    {
        auto meta = new BucketMeta();
        std::filesystem::path dir(BULLDOG_DIR_DATA_ABS);
        if (name.substr(0, 12) == "hierarchical") {
            // Get the round number.
            std::vector<std::string> parsed_str;
            split_string(name, "_", parsed_str); // TODO(kwok): Boost it!
            meta->bucket_.LoadHierarchical(name);
            meta->bucket_count_ = std::atoi(parsed_str[1].c_str()) * std::atoi(parsed_str[2].c_str());
            meta->post_flop_bucket_count_ = std::atoi(parsed_str[2].c_str());
        } else if (name.substr(0, 5) == "waugh") {
            std::filesystem::path bucket_dir(name);
            std::vector<std::string> parsed_str;
            split_string(name, "_", parsed_str);
            // int num_buckets = std::stoi(parsed_str[5]);
            uint8_t r = std::stoul(parsed_str[2]);
            size_t num_buckets = meta->bucket_.LoadClassicFromFlexbuffers(dir / bucket_dir, r, lossless);
            logger::info("number of buckets %lu; round %u", num_buckets, r);
            meta->bucket_count_ = num_buckets;
        } else {
            std::filesystem::path bucket_file(name + ".bin");
            meta->bucket_.LoadClassicFromFile(dir / bucket_file);
            logger::debug("loaded bucket from file %s", dir / bucket_file);
            //get the bucket number // may need to extract this.
            std::vector<std::string> parsed_str;
            split_string(name, "_", parsed_str);
            meta->bucket_count_ = std::stoi(parsed_str[1]);
        }
        return meta;
    }

    void InsertBucket(BucketMeta *meta, std::string name, int r, bool lossless)
    {
        bucket_pool_[r][lossless][name] = meta;
    }

private:
    std::array<
            std::array<
                    std::map<std::string /* name */, BucketMeta *>,
                    2
            >,
            4 /* round */
    > bucket_pool_;
};

#endif //BULLDOG_MODULES_ENGINE_SRC_BUCKET_POOL_HPP_
