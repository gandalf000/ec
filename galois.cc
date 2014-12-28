/**
 * Copyright (c) 2014, The Authors. All rights reserved
 * @file common/galois.cpp
 * @brief implement arithmetic operation on Galois Fiedld(GF(2^8))
 */

#include "common/galois.h"
#include <stddef.h>
#include <sys/syscall.h>

static const int word_size = 8;         // the value of w in Galois Field
static const int gf_size = (1 << word_size);  // the size of GF(2^8)
static const int max_in_gf = (1 << word_size) - 1;  // the max number of GF(2^8)
static const int prim_poly = 0435;

/**
 * @brief static function, create log table and anti log table
 */
static void _CreateLogTable(int **log_table_ptr, int **anti_log_table_ptr) {
    // NOLINT
    size_t table_size = sizeof(int) * gf_size;
    *log_table_ptr = new int[table_size];
    *anti_log_table_ptr = new int[table_size];

    int *log_table = *log_table_ptr;
    int *anti_log_table = *anti_log_table_ptr;

    for (int i = 0; i < gf_size ; i++) {
        log_table[i] = max_in_gf;
        anti_log_table[i] = 0;
    }

    int index = 1;
    for (int log = 0; log < max_in_gf ; log++) {
        log_table[index] = log;
        anti_log_table[log] = index;
        index = (index << 1);
        if (index & gf_size)
            index = (index ^ prim_poly) & max_in_gf;
    }
    anti_log_table[max_in_gf] = anti_log_table[0];
}

void GaloisOperator::_CreateTables() {
    // create log table and anti-log table
    int *log_table = NULL;
    int *anti_log_table = NULL;
    _CreateLogTable(&log_table, &anti_log_table);

    // create multi table and divide table
    // NOLINT
    size_t table_size = sizeof(int) * gf_size * gf_size;
    m_mul_table = new int[table_size];
    m_div_table = new int[table_size];

    // Set mult/div tables
    m_mul_table[0] = 0;   // for x = 0 and y = 0
    m_div_table[0] = -1;

    for (int y = 1; y < gf_size; y++) {  // for x = 0 and y > 0
        m_mul_table[y] = 0;
        m_div_table[y] = 0;
    }

    int index = gf_size;
    for (int x = 1; x < gf_size; x++) {  // x > 0
        m_mul_table[index] = 0;  // y = 0
        m_div_table[index] = -1;
        index++;
        int logx = log_table[x];
        for (int y = 1; y < gf_size; y++) {  //  y > 0
            int logy = log_table[y];
            int mul_table_index = logx + logy;
            if (mul_table_index > max_in_gf)
                mul_table_index -= max_in_gf;
            m_mul_table[index] = anti_log_table[mul_table_index];

            int div_table_index = logx - logy;
            if (div_table_index < 0)
                div_table_index += max_in_gf;
            m_div_table[index] = anti_log_table[div_table_index];
            index++;
        }
    }

    delete[] log_table;
    delete[] anti_log_table;
}

/**
 * @brief return x * y based on Galois Field(GF(2^8))
 */
int GaloisOperator::Multiply(int x, int y) {
    if (UNLIKELY(x == 0 || y == 0))
        return 0;

    return m_mul_table[(x << word_size) | y];
}

/**
 * @brief return x / y based on Galois Field(GF(2^8))
 */
int GaloisOperator::Divide(int x, int y) {
    if (UNLIKELY(y == 0))
        return -1;
    if (UNLIKELY(x == 0))
        return 0;

    return m_div_table[(x << word_size) | y];
}

