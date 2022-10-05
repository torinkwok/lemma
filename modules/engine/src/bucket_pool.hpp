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
            for (auto &[key, val]: it) {
                delete val;
            }
        }
    }

    BucketMeta *Get(const std::string &name, int round)
    {
        if (!Has(name, round)) {
            // Loading bucket lazily.
            bucket_pool_[round][name] = LoadBucket(name);
        }
        return bucket_pool_[round][name];
    }

    bool Has(const std::string &name, int round)
    {
        auto it = bucket_pool_[round].find(name);
        return it != bucket_pool_[round].end();
    }

    static BucketMeta *LoadBucket(const std::string &name)
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
            int num_buckets = std::stoi(parsed_str[5]);
            uint8_t r = std::stoul(parsed_str[2]);
            logger::info("number of buckets %d; round %u", num_buckets, r);
            meta->bucket_.LoadClassicFromFlexbuffers(dir / bucket_dir, r);
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

    void InsertBucket(BucketMeta *meta, std::string name, int r)
    {
        bucket_pool_[r][name] = meta;
    }

private:
    std::array<
            std::map<std::string /* name */, BucketMeta *>,
            4 /* round */
    > bucket_pool_;
};

#endif //BULLDOG_MODULES_ENGINE_SRC_BUCKET_POOL_HPP_
