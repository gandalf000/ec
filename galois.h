/**
 * Copyright (c) 2014, The Authors. All rights reserved
 * @file common/galois.h
 * @brief implement arithmetic operation on Galois Fiedld(GF(2^8))
 * @date 2014-12-29
 */

#include <sys/syscall.h>
#include <linux/futex.h>
#include <sys/time.h>

/**
 * @brief implement arithmetic operation on GF(2^8)
 */
class GaloisOperator {
public:
    GaloisOperator() {
        _CreateTables();
    }

    ~GaloisOperator() {
        delete[] m_mul_table;
        delete[] m_div_table;
    }

    /**
     * @brief return x * y based on Galois Field(GF(2^8))
     */
    int Multiply(int x, int y);

    /**
     * @brief return x / y based on Galois Field(GF(2^8))
     */
    int Divide(int x, int y);

private:
    /**
     * @brief create multi table and divide table to accelerate 
     *        multiply and divide operation on GF(2^8)
     */
    void _CreateTables();

    int *m_mul_table;    /**< used to accelerate multiply on GF(2^8)*/
    int *m_div_table;    /**< used to accelerate divide on GF(2^8)*/
};
