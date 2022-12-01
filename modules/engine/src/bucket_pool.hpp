#ifndef AUTODIDACT_MODULES_ENGINE_SRC_BUCKET_POOL_HPP_
#define AUTODIDACT_MODULES_ENGINE_SRC_BUCKET_POOL_HPP_

#include "bucket.h"
#include <map>
#include <array>
#include <filesystem>
#include <autodidact/core_util.hpp>

struct sBucketMeta
{
    Bucket bucket;
    size_t bucket_count = 0;
    size_t post_flop_bucket_count = 0;
};

class BucketPool
{
public:
    virtual ~BucketPool()
    {
        for (auto &buckets_it: _bucket_meta_pool) {
            for (auto &lossy_and_lossless_it: buckets_it) {
                for (auto &[key, val]: lossy_and_lossless_it) {
                    delete val;
                }
            }
        }
    }

    sBucketMeta *LazyLoadBucketMeta(const std::string &name, int round, bool lossless = false)
    {
        if (!Has(name, round, lossless)) {
            // Loading bucket lazily.
            _bucket_meta_pool[round][lossless][name] = LoadBucket(name, lossless);
            logger::info("🚛cache lazily loaded for name=%s, round=%d, lossless=%d", name, round, lossless);
        } else {
            logger::info("🎯cache hit for name=%s, round=%d, lossless=%d", name, round, lossless);
        }
        return _bucket_meta_pool[round][lossless][name];
    }

    bool Has(const std::string &name, int round, bool lossless)
    {
        auto it = _bucket_meta_pool[round][lossless].find(name);
        return it != _bucket_meta_pool[round][lossless].end();
    }

    static sBucketMeta *LoadBucket(const std::string &name, bool lossless = false)
    {
        auto *meta = new sBucketMeta();
        std::filesystem::path dir(AUTODIDACT_DIR_DATA_ABS);
        if (name.substr(0, 12) == "hierarchical") {
            // LazyLoadBucketMeta the round number.
            std::vector<std::string> parsed_str;
            split_string(name, "_", parsed_str); // TODO(kwok): Boost it!
            meta->bucket.LoadHierarchical(name);
            meta->bucket_count = std::atoll(parsed_str[1].c_str()) * std::atoll(parsed_str[2].c_str());
            meta->post_flop_bucket_count = std::atoll(parsed_str[2].c_str());
        } else if (name.substr(0, 5) == "waugh") {
            std::filesystem::path bucket_dir(name);
            std::vector<std::string> parsed_str;
            split_string(name, "_", parsed_str);
            // int n_buckets = std::stoi(parsed_str[5]);
            uint8_t r = std::stoul(parsed_str[2]);
            size_t n_buckets = meta->bucket.LoadClassicFromFlexbuffers(dir / bucket_dir, r, lossless);
            logger::info("number of buckets %lu; round %u", n_buckets, r);
            meta->bucket_count = n_buckets;
        } else {
            std::filesystem::path bucket_file(name + ".bin");
            meta->bucket.LoadClassicFromFile(dir / bucket_file);
            logger::debug("loaded bucket from file %s", dir / bucket_file);
            //get the bucket number // may need to extract this.
            std::vector<std::string> parsed_str;
            split_string(name, "_", parsed_str);
            meta->bucket_count = std::stoi(parsed_str[1]);
        }
        return meta;
    }

    void InsertBucket(sBucketMeta *meta, const std::string &name, int r, bool lossless)
    {
        _bucket_meta_pool[r][lossless][name] = meta;
    }

private:
    std::array<
            std::array<
                    std::map<std::string /* name */, sBucketMeta *>,
                    2
            >,
            4 /* round */
    > _bucket_meta_pool;
};

#endif //AUTODIDACT_MODULES_ENGINE_SRC_BUCKET_POOL_HPP_
