// FIXME(kwok): If this directive is put below the import of "node.h", "expected body of lambda expression" error pops up.
#include <flatbuffers/flexbuffers.h>

#include "bucket.h"
#include "node.h"

#include <bulldog/logger.hpp>

#include <cereal/archives/binary.hpp>
#include <filesystem>
#include <bulldog/card_util.h>
#include <utility>
#include <bulldog/core_util.hpp>
#include <pokerstove/peval/CardSetGenerators.h>

extern "C" {
#include <bulldog/game.h>
}

void Bucket::LoadClassicFromFlexbuffers(const std::string &dir, uint8_t r) {
    assert(r >= 0 && r < 4);
    type_ = WAUGH_BUCKET;

    // TODO(kwok): Justify it.
    uint_fast32_t rounds = 0;
    uint8_t *cards_per_round = nullptr;
    hand_indexer_t *out_indexer = nullptr;
    switch (r) {
        case 0: {
            rounds = 1;
            cards_per_round = new uint8_t[]{2};
            out_indexer = &_preflop_indexer;
            break;
        }
        case 1:
            rounds = 2;
            cards_per_round = new uint8_t[]{2, 3};
            out_indexer = &_flop_indexer;
            break;
        case 2:
            rounds = 2;
            cards_per_round = new uint8_t[]{2, 4};
            out_indexer = &_turn_indexer;
            break;
        case 3:
            rounds = 2;
            cards_per_round = new uint8_t[]{2, 5};
            out_indexer = &_river_indexer;
            break;
        default:
            throw std::runtime_error("Invalid round number " + std::to_string(r));
    }

    assert(hand_indexer_init(rounds, cards_per_round, out_indexer));

    auto num_loaded = _LoadClassicFromFlexbuffers(std::filesystem::path(dir), r);
    logger::info("🪣%lu Waugh-indexed buckets loaded for round %u", num_loaded, r);
    assert((r == 0 && num_loaded == 169) ||
           (r == 1 && num_loaded == 1286792) ||
           (r == 2 && num_loaded == 13960050) ||
           (r == 3 && num_loaded == 123156254));

    delete[] cards_per_round;
}

size_t Bucket::_LoadClassicFromFlexbuffers(const std::string &dir, uint8_t r) {
    assert(type_ == WAUGH_BUCKET);

    std::filesystem::path fxb_chunks_dir_path(dir);
    assert(exists(fxb_chunks_dir_path));
    assert(is_directory(fxb_chunks_dir_path));

    using std::filesystem::directory_iterator;
    size_t fxb_num_chunks = std::distance(directory_iterator(fxb_chunks_dir_path), directory_iterator{});
    size_t processed_count = 0;
    for (size_t chunk_id = 0; chunk_id < fxb_num_chunks; chunk_id++) {
        auto fxb_chunk_path = fxb_chunks_dir_path / ("chunk_" + std::to_string(chunk_id) + ".fxb");
        std::ifstream fxb_chunk_file(fxb_chunk_path, std::ios::binary | std::ios::ate);
        if (fxb_chunk_file.is_open()) {
            std::streamsize size = fxb_chunk_file.tellg();
            fxb_chunk_file.seekg(0, std::ios::beg);
            auto buffer = std::vector<uint8_t>(size);
            if (fxb_chunk_file.read(reinterpret_cast<char *>(buffer.data()), size).bad()) {
                // TODO:(kwok) More informative exception.
                throw std::runtime_error("Failed to read from " + fxb_chunk_path.string());
            }
            auto fxb_tmp = flexbuffers::GetRoot(buffer).AsMap();
            assert(!fxb_tmp.IsTheEmptyMap());
            auto global_offset = std::stoul(fxb_tmp.Keys()[0].AsString().str());
            auto vec = fxb_tmp.Values()[0].AsTypedVector();
            for (size_t local_index = 0; local_index < vec.size(); local_index++) {
                auto isomorphic_index = local_index + global_offset;
                master_map_[r][isomorphic_index] = vec[local_index].AsUInt32();
            }
            processed_count += vec.size();
        } else {
            logger::error("unable to open " + dir);
            exit(EXIT_FAILURE);
        }
        fxb_chunk_file.close();
    }
    return processed_count;
}

void Bucket::LoadClassicFromFile(const std::string &ofile) {
    type_ = CLASSIC_BUCKET;

    std::filesystem::path full_path(ofile);
    std::vector<std::string> str;
    split_string(full_path.filename(), "_", str);
    if (str[0] != "buckets") {
        logger::critical("Bucket unable to load file: %s", ofile);
    }

    std::ifstream is(ofile, std::ios::binary);
    if (is.is_open()) {
        cereal::BinaryInputArchive load(is);
        load(cluster_map_);
        //kmeans_map has been stored as std::map
        for (const auto &[key, val]: cluster_map_) {
            master_map_[0][key] = val;
        }
    } else {
        logger::error("unable to open " + ofile);
        exit(EXIT_FAILURE);
    }
    is.close();
}

void Bucket::LoadRangeColex(Board_t *board, int round) {
    type_ = COLEX_BUCKET;

    std::set<Colex> colex_set;
    int bucket_index = 0;
    //for each hand, compute the colex value.
    for (Card_t low = 0; low < HOLDEM_MAX_DECK - 1; low++) {
        for (Card_t high = low + 1; high < HOLDEM_MAX_DECK; high++) {
            auto hand = PrivHand_t{high, low};
            if (board->PrivHandCrash(hand)) {
                continue;
            }
            auto colex = ComputeColexFromAllCards(high, low, *board, round);
            if (colex_set.find(colex) == colex_set.end()) {
                //insert, and increment the value.
                master_map_[0][colex] = bucket_index;
                colex_set.insert(colex);
                bucket_index++;
            }
        }
    }
    logger::trace("colex bucket at round = %d || num of unique colexes = %d || bucket count = %d",
                  round, colex_set.size(), bucket_index);
}

void Bucket::LoadHierarchicalPublic() {
    //find public bucket for flop
    std::filesystem::path dir(BULLDOG_DIR_DATA_ABS);
    std::filesystem::path pub_file("hierarchical_pubcolex_60_2_3.txt");
    std::ifstream pub_is(dir / pub_file, std::ios::binary);
    if (pub_is.is_open()) {
        std::string line;
        std::getline(pub_is, line);
        while (std::getline(pub_is, line)) {
            int board_idx, bucket;
            std::sscanf(line.c_str(), "%d,%d", &board_idx, &bucket);
            pub_colex_bucket_[board_idx] = bucket;
            assert(bucket < 60);
        }
    } else {
        logger::critical("unable to open %s", dir / pub_file);
    };
    pub_is.close();
}

void Bucket::LoadHierarchical(std::string name) {
    //name: hierarchical_60_500_1
    std::vector<std::string> parsed_str;
    split_string(std::move(name), "_", parsed_str);
    unsigned int round = std::atoi(parsed_str[3].c_str());
    unsigned int pub_bucket_count = std::atoi(parsed_str[1].c_str());
    unsigned int priv_bucket_count = std::atoi(parsed_str[2].c_str());

    type_ = HIERARCHICAL_BUCKET;
    //load hierarchical only works after flop
    assert(round > 0);
    std::filesystem::path dir(BULLDOG_DIR_DATA_ABS);

    LoadHierarchicalPublic();

#if DEV > 1
    Bucket_t total_bucket_count = priv_bucket_count * pub_bucket_count;
    std::vector<bool> bucket_tally;
    bucket_tally.reserve(total_bucket_count);
    for (Bucket_t i = 0; i < total_bucket_count; i++) {
        bucket_tally.emplace_back(false);
    }
#endif

    for (unsigned int pub_bucket = 0; pub_bucket < pub_bucket_count; pub_bucket++) {
        //regular pub_bucket load
        // FIXME(kwok): Does the hard-coded number 2 denote the number of players?
        std::filesystem::path
                priv_file("buckets_" + parsed_str[1] + "_" + parsed_str[2] + "_2_" + std::to_string(round + 2) + "_" +
                          std::to_string(pub_bucket) + ".bin");
        std::ifstream priv_is(dir / priv_file, std::ios::binary);
        if (priv_is.is_open()) {
            cereal::BinaryInputArchive load(priv_is);
            load(cluster_map_);
            logger::debug("loaded pub_bucket file %s of size %d", priv_file, cluster_map_.size());
        } else {
            logger::critical("unable to open %s", dir / priv_file);
        }
        priv_is.close();
        // post process to assign pub_bucket into range of 0-29999
        for (auto const &[key, val]: cluster_map_) {
            assert(val < priv_bucket_count);
            // NOTE(kwok): Kind of like indexing a tensor.
            master_map_[pub_bucket][key] = (pub_bucket * priv_bucket_count) + val; //it > 2^16
#if DEV > 1
            bucket_tally[master_map_[pub_bucket][key]] = true;
#endif
        }
        cluster_map_.clear();
    }
#if DEV > 1
    unsigned int empty_bucket_count = 0;
    for (unsigned int i = 0; i < total_bucket_count; i++) {
        if (!bucket_tally[i]) {
            empty_bucket_count++;
        }
    }
    logger::debug("%d/%d hierarchical buckets in round %d are empty!", empty_bucket_count, total_bucket_count, round);
#endif
}

unsigned int Bucket::GetPublicBucket(unsigned int pub_colex) {
    return pub_colex_bucket_[pub_colex];
}

/*
 * it should be board + 2, e.g. 3+2 for flop
 */
void Bucket::LoadHierarchicalColex(Board_t *board, uint8_t r) {
    if (r > HOLDEM_ROUND_RIVER) {
        logger::critical("round %d does not exist in holdem", r);
    }
    if (r == HOLDEM_ROUND_PREFLOP) {
        logger::critical("does not support board hand colex for preflop", r);
    }
    type_ = HIERARCHICAL_COLEX;
    int board_count = r == 1 ? 3 :
                      r == 2 ? 4 : 5;
    auto iso_board_combo = pokerstove::createCardSet(board_count, pokerstove::Card::SUIT_CANONICAL);
    auto raw_board_combo = pokerstove::createCardSet(board_count);
    auto raw_hand_combo = pokerstove::createCardSet(2);

    int iso_board_cursor = 0;
    int bucket_idx_cursor = 0;

    for (const auto &iso_board: iso_board_combo) {
        auto iso_board_colex = iso_board.colex();
        for (const auto &raw_board: raw_board_combo) {
            if (raw_board.canonize().colex() != iso_board_colex) {
                continue;
            }
            //the raw board belongs to the colex board
            for (const auto &raw_hand: raw_hand_combo) {
                // NOTE(kwok): A hand that is not possible given the board.
                if (raw_hand.intersects(raw_board)) {
                    continue;
                }
                // now the hands are valid.
                auto hand_board = pokerstove::CardSet(raw_hand.str() + raw_board.str());
                Colex all_colex = hand_board.canonize().colex();
                // skip duplicated entires
                if (master_map_[iso_board_colex].find(all_colex) == master_map_[iso_board_colex].end()) {
                    master_map_[iso_board_colex][all_colex] = bucket_idx_cursor;
                    bucket_idx_cursor++;
                }
            }
        }
        iso_board_cursor++;
    }
    //  logger::debug("total [%d board colex] [%d hier colex] buckets", iso_board_cursor, bucket_idx_cursor);
}

void Bucket::LoadSubgameColex(Board_t *board, int round) {
    if (round != HOLDEM_ROUND_RIVER) {
        logger::critical("subgame colex only support river. "
                         "other rounds would be too big. "
                         "use other abstraction algorithms instead.");
    }
    type_ = HIERARCHICAL_COLEX;
    //  auto cmd_begin = std::chrono::steady_clock::now();
    int bucket_idx_cursor = 0;
    for (Card_t c = 0; c < HOLDEM_MAX_DECK; c++) {
        if (board->CardCrash(c)) continue;
        //use a local board
        Board_t local_board = *board;
        local_board.cards[4] = c;

        auto board_set = emptyCardset();
        for (unsigned char card: local_board.cards) {
            AddCardTToCardset(&board_set, card);
        }
        auto board_colex = ComputeColex(Canonize(board_set.cards));
        //for each hand, compute the colex value.
        for (Card_t low = 0; low < HOLDEM_MAX_DECK - 1; low++) {
            for (Card_t high = low + 1; high < HOLDEM_MAX_DECK; high++) {
                auto hand = PrivHand_t{high, low};
                if (local_board.PrivHandCrash(hand)) continue;
                auto full_colex = ComputeColexFromAllCards(high, low, local_board, round);
                if (master_map_[board_colex].find(full_colex) == master_map_[board_colex].end()) {
                    //insert, and increment the value.
                    master_map_[board_colex][full_colex] = bucket_idx_cursor;
                    bucket_idx_cursor++;
                }
            }
        }
    }
    //  auto cmd_time =
    //      std::chrono::duration_cast<std::chrono::milliseconds>(
    //          std::chrono::steady_clock::now() - cmd_begin).count();
    //  logger::debug("generate subgame colex takes %d ms", cmd_time);
}

void Bucket::Save(std::map<unsigned int, unsigned short> &entries, const std::string &ofile) {
    std::ofstream os(ofile, std::ios::binary | std::ios::trunc);
    if (os.is_open()) {
        cereal::BinaryOutputArchive archive(os);
        archive(entries);
    } else {
        logger::error("unable to open " + ofile);
    };
    os.close();
}

//todo: all these can be handled on the outside, at the bucket reader level
// TODO(kwok): `board_colex` will be ignored for classic and colex bucket. Make it explicit.
uint32_t Bucket::Get(unsigned long all_colex, unsigned long board_colex) {
    if (type_ == WAUGH_BUCKET) {
        throw std::runtime_error("For Waugh bucket, invoke `uint32_t Bucket::Get(Cardset*, Cardset*)`");
    }

    if (board_colex > 4294967295 && type_ == HIERARCHICAL_COLEX) {
        logger::critical("now the master map is indexed with unsigned int. so the public colex is not safe.");
    }

    if (all_colex > 4294967295) {
        logger::critical("all colex %d too large out of the range of unsigned int", all_colex);
    }

    switch (type_) {
        case HIERARCHICAL_COLEX:
        case HIERARCHICAL_BUCKET: {
            if (master_map_.empty()) {
                logger::critical("must call Bucket::LoadHierarchical first");
                return INVALID_BUCKET;
            }
            unsigned int pub_bucket;
            if (type_ == HIERARCHICAL_BUCKET) {
                pub_bucket = pub_colex_bucket_[board_colex];
            } else {
                pub_bucket = board_colex;
            }

            //check and query
            auto priv_buckets = master_map_.find(pub_bucket);
            if (priv_buckets == master_map_.end()) {
                logger::error("could not find hierarchical private bucket for [pub board colex %d]", pub_bucket);
                return INVALID_BUCKET;
            }
            auto it = priv_buckets->second.find(all_colex);
            if (it == priv_buckets->second.end()) {
                logger::error("could not find private bucket for [board colex %d] [all colex %d]", board_colex,
                              all_colex);
                return INVALID_BUCKET;
            }
            return it->second;
        }
        case CLASSIC_BUCKET:
        case COLEX_BUCKET: {
            return master_map_[0][all_colex];
        }
        default: {
            logger::critical("bucket does not have type_d %d", type_);
            return INVALID_BUCKET;
        }
    }
}


uint32_t Bucket::Get(Cardset *all_cards, Cardset *board_cards) {
    if (type_ == WAUGH_BUCKET) {
        std::set<WaughCard_t> waugh_cards_set = CardsToWaughCards(all_cards->cards);
        // TODO(kwok): Assert the two are identical.
        // logger::debug("☄️%s - %s", WaughCardsToString(waugh_cards_set), CardsToString(all_cards->cards));
        auto waugh_cards_vec = std::vector(waugh_cards_set.begin(), waugh_cards_set.end());
        WaughCard_t *c_arr = &waugh_cards_vec[0];
        uint8_t r = 0;
        hand_index_t index = 0;
        // TODO(kwok): Refactor these deliberately redundant code.
        switch (waugh_cards_vec.size()) {
            case 2:
                r = 0;
                index = hand_index_last(&_preflop_indexer, c_arr);
                break;
            case 5:
                r = 1;
                index = hand_index_last(&_flop_indexer, c_arr);
                break;
            case 6:
                r = 2;
                index = hand_index_last(&_turn_indexer, c_arr);
                break;
            case 7:
                r = 3;
                index = hand_index_last(&_river_indexer, c_arr);
                break;
            default:
                throw std::runtime_error("Invalid size of Waugh cards");
        }
        // TODO(kwok): Validate the bucket to return.
        auto bucket = master_map_[r][index];
        return bucket;
    }
    return Get(ComputeColex(Canonize(all_cards->cards)),
               ComputeColex(Canonize(board_cards->cards)));
}

uint32_t Bucket::Size() {
    if (type_ == HIERARCHICAL_COLEX) {
        //it is alright, because we never use this for subgame solving for flop
        int sum = 0;
        for (auto [key, val]: master_map_) {
            //      logger::debug("ket  = %d", key);
            sum += val.size();
        }
        logger::trace("hierarchical colex with iso flop = %d total = %d", master_map_.size(), sum);
        return sum;
    }
    return master_map_[0].size();
}

std::unordered_map<unsigned int, uint32_t> Bucket::ExtractMap() {
    return master_map_[0];
}