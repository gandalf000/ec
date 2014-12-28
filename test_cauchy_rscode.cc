// Copyright (c) 2013, The Authors. All rights reserved.

extern "C" {
#include "common/jerasure.h"
#include "common/jerasure_cauchy.h"
}

#include "common/cauchy_rscode.h"

#include "gtest/gtest.h"

namespace {

TEST(TestCauchyRSCoder, GenerateEncodeMatrix)
{
    CauchyRSCoder *coder =  new CauchyRSCoder(8, 4);
    int *matrix = coder->_GenerateEncodeMatrix();
    int *jerasure_matrix = cauchy_good_general_coding_matrix(8, 4, 8);
    for (int i = 0; i < 8 * 4; i++) {
        ASSERT_EQ(matrix[i], jerasure_matrix[i]);
    }

    char *bit_matrix = coder->m_encoding_bit_matrix;
    int *jerasure_bit_matrix = jerasure_matrix_to_bitmatrix(8, 4, 8, jerasure_matrix);
    for (int i = 0; i < 8 * 4 * 8 * 8; i++) {
        ASSERT_EQ(jerasure_bit_matrix[i], bit_matrix[i]);
    }

    int **schedule = coder->m_encoding_schedule;
    int **jerasure_schedule = jerasure_smart_bitmatrix_to_schedule(8, 4, 8, jerasure_bit_matrix);
    int count1 = 0;
    int count2 = 0;
    for (; schedule[count1][0] >= 0; count1++) {}
    for (; jerasure_schedule[count2][0] >= 0; count2++) {}
    ASSERT_EQ(count1, count2);
    for (int i = 0; i < count1; i++)
        for (int j = 0; j < 5; j++) {
            ASSERT_EQ(schedule[i][j], jerasure_schedule[i][j]);
        }

    delete[] matrix;
    jerasure_free_schedule(jerasure_schedule);
    free(jerasure_matrix);
    free(jerasure_bit_matrix);
    delete coder;
}

TEST(TestCauchyRSCoder, TestEncode)
{
    CauchyRSCoder *coder = new CauchyRSCoder(8, 4);
    char *data_ptrs[8];
    char *code_ptrs[4];
    char *jerasure_code_ptrs[4];
    for (int i = 0; i < 8; i++) {
        data_ptrs[i] = new char[1 << 20];
        int *ptr = reinterpret_cast<int *>(data_ptrs[i]);
        for (unsigned i = 0; i < (1 << 20) / sizeof(int); i += sizeof(int), ptr++) { // NOLINT
            *ptr = random();
        }
    }
    for (int i = 0; i < 4; i++) {
        code_ptrs[i] = new char[1 << 20];
        jerasure_code_ptrs[i] = new char[1 << 20];
    }
    int *jerasure_matrix = cauchy_good_general_coding_matrix(8, 4, 8);
    int *jerasure_bit_matrix = jerasure_matrix_to_bitmatrix(8, 4, 8, jerasure_matrix);
    int **jerasure_schedule = jerasure_smart_bitmatrix_to_schedule(8, 4, 8, jerasure_bit_matrix);

    coder->Encode(data_ptrs, code_ptrs, 1 << 20);
    jerasure_schedule_encode(8, 4, 8, jerasure_schedule, data_ptrs,
            jerasure_code_ptrs, 1 << 20, kPacketSize);
    for (int i = 0; i < 4; i++) {
        ASSERT_EQ(memcmp(code_ptrs[i], jerasure_code_ptrs[i], 1 << 20), 0);
    }

    for (int i = 0; i < 8; i++) {
        delete[] data_ptrs[i];
    }
    for (int i = 0; i < 4; i++) {
        delete[] code_ptrs[i];
        delete[] jerasure_code_ptrs[i];
    }
    jerasure_free_schedule(jerasure_schedule);
    free(jerasure_matrix);
    free(jerasure_bit_matrix);
    delete coder;
}

TEST(TestCauchyRSCoder, TestDecode)
{
    CauchyRSCoder *coder = new CauchyRSCoder(8, 4);
    char *data_ptrs[8];
    char *code_ptrs[4];
    char *erased_data_ptrs[8];
    char *erased_code_ptrs[4];

    for (int i = 0; i < 8; i++) {
        data_ptrs[i] = new char[1 << 20];
        erased_data_ptrs[i] = new char[1 << 20];
        int *ptr = reinterpret_cast<int *>(data_ptrs[i]);
        for (unsigned j = 0; j < (1 << 20) / sizeof(int); j += sizeof(int), ptr++) { // NOLINT
            *ptr = random();
        }

        memcpy(erased_data_ptrs[i], data_ptrs[i], 1 << 20);
    }

    for (int i = 0; i < 4; i++) {
        code_ptrs[i] = new char[1 << 20];
        erased_code_ptrs[i] = new char[1 << 20];
    }

    coder->Encode(data_ptrs, code_ptrs, 1 << 20);

    for (int i = 0; i < 4; i++) {
        memcpy(erased_code_ptrs[i], code_ptrs[i], 1 << 20);
    }

    bool erased[12];
    // test no fail
    memset(erased, 0, sizeof(erased) / sizeof(bool)); // NOLINT
    coder->Decode(erased, erased_data_ptrs, erased_code_ptrs, 1 << 20);
    for (int i = 0; i < 4; i++) {
        ASSERT_EQ(memcmp(code_ptrs[i], erased_code_ptrs[i], 1 << 20), 0);
    }
    for (int i = 0; i < 8; i++) {
        ASSERT_EQ(memcmp(data_ptrs[i], erased_data_ptrs[i], 1 << 20), 0);
    }
    // test 1 to 4 device fail
    int fail_index = 0;
    for (int fail = 1; fail <= 4; fail++) {
        for (int fail_data = 0; fail_data <= fail; fail_data++) {
            for (int fail_code = 0; fail_data + fail_code <= fail; fail_code++) {
                if (fail_data == 0 && fail_code == 0) {
                    continue;
                }
                printf("test fail data parts %d fail code parts %d\n", fail_data, fail_code);
                memset(erased, 0, sizeof(erased) / sizeof(bool)); // NOLINT

                for (int i = 0; i < fail_data; i++) {
                    fail_index = rand() % 8; // NOLINT
                    // ake sure there is no duplicate element in erasures
                    while (erased[fail_index]) {
                        fail_index = rand() % 8;    // NOLINT
                    }

                    erased[fail_index] = true;
                    bzero(erased_data_ptrs[fail_index], 1 << 20);
                }
                for (int i = 0; i < fail_code; i++) {
                    fail_index = (rand() % 4) + 8; // NOLINT
                    // ake sure there is no duplicate element in erasures
                    while (erased[fail_index]) {
                        fail_index = (rand() % 4) + 8;    // NOLINT
                    }

                    erased[fail_index] = true;
                    bzero(erased_code_ptrs[fail_index - 8], 1 << 20);
                }
                coder->Decode(erased, erased_data_ptrs, erased_code_ptrs, 1 << 20);
                for (int i = 0; i < 4; i++) {
                    ASSERT_EQ(memcmp(code_ptrs[i], erased_code_ptrs[i], 1 << 20), 0);
                }
                for (int i = 0; i < 8; i++) {
                    ASSERT_EQ(memcmp(data_ptrs[i], erased_data_ptrs[i], 1 << 20), 0);
                }
            }
        }
    }

    // est optimized
    for (int i = 0; i < 3; i++) {
        memset(erased, 0, sizeof(erased) / sizeof(bool)); // NOLINT
        for (int j = 0; j <= i; j++) {
            fail_index = rand() % 8; // NOLINT
            // ake sure there is no duplicate element in erasures
            while (erased[fail_index]) {
                fail_index = rand() % 8;    // NOLINT
            }
            erased[fail_index] = true;
            bzero(erased_data_ptrs[fail_index], 1 << 20);
        }

        for (int j = 3; j > i; j--) {
            bzero(erased_code_ptrs[j], 1 << 20);
        }
        for (int j = 0; j <= i; j++) {
            memcpy(erased_code_ptrs[j], code_ptrs[j], 1 << 20);
        }

        coder->Decode(erased, erased_data_ptrs, erased_code_ptrs, 1 << 20);
        for (int j = 0; j < 8; j++) {
            ASSERT_EQ(memcmp(data_ptrs[j], erased_data_ptrs[j], 1 << 20), 0);
        }

        for (int j = 0; j <= i; j++) {
            ASSERT_EQ(memcmp(code_ptrs[j], erased_code_ptrs[j], 1 << 20), 0);
        }
    }

    for (int i = 0; i < 8; i++) {
        delete[] data_ptrs[i];
        delete[] erased_data_ptrs[i];
    }
    for (int i = 0; i < 4; i++) {
        delete[] code_ptrs[i];
        delete[] erased_code_ptrs[i];
    }
    delete coder;
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    srand(time(NULL));
    return RUN_ALL_TESTS();
}

}
