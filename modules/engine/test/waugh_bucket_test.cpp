//
// Created by Torin Tung Kwok on 2022/9/17.
//

#include "../src/bucket.h"

int main() {
    Bucket bucket;
    bucket.LoadClassicFromFlexbuffers();
    auto all_cards_set = CardsetFromString("4s5s3hJhJc");
    bucket.Get(&all_cards_set, NULL);
    return 0;
}
