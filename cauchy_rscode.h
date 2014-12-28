/**
 * Copyright (c) 2014, The Authors. All rights reserved
 * @file cauchy_rscode.h
 * @brief Cauchy Reed-Solomon encoding and decoding library
 */

#include <assert.h>
#include <stddef.h>
#include <linux/futex.h>
#include <sys/time.h>
#include "common/galois.h"

static const int kPacketSize = 4096;
static const int kWordBits = 8;
static const int kCodingUnitSize = kPacketSize * kWordBits;

/**
 * @brief Cauchy Reed-Solomon encoding and decoding library
 */
class CauchyRSCoder {
public:
    CauchyRSCoder(int num_data_parts, int num_code_parts) {
        assert(num_data_parts > 0);
        assert(num_code_parts > 0);

        m_num_data_parts = num_data_parts;
        m_num_code_parts = num_code_parts;
        m_galois_operator = new GaloisOperator;
        m_encoding_schedule = NULL;
        m_encoding_bit_matrix = NULL;

        _Init();
    }

    ~CauchyRSCoder() {
        delete m_galois_operator;
        delete[] m_encoding_bit_matrix;
        _FreeSchedule(m_encoding_schedule);
    }

    /**
     * @brief encoding data_parts_n data parts into code_parts_n coding parts
     *
     * @param data_ptrs     Array of num_data_parts pointers to data
     * @param coding_ptrs   Array of num_code_parts pointers to coding data
     * @param size          Size of memory allocated by data_ptrs in bytes.
     */
    void Encode(char **data_ptrs, char **coding_ptrs, int size);

    /**
     * @brief This function recover from any <=m parts failure
     *
     * @param erased Array of indicators which point out whether the device has
     *               been erased. The index of array range frome 0 to k+m-1.
     *               Id's 0 to k-1 are id's of data devices. Id's k to k+m-1 are
     *               id's of coding devices.if erased[i] is true, then the device
     *               i is failed, and should be recovered frome decoding
     * @param data_ptrs     Array of num_data_parts pointers to data
     * @param coding_ptrs   Array of num_code_parts pointers to coding data
     * @param size          Size of memory allocated by data_ptrs in bytes.
     */
    void Decode(bool *erased, char **data_ptrs, char **coding_ptrs, int size);

private:
    void _Init();

    /**
     * @brief Returns the number of ones in the bitmatrix representation of
     *        the number num. The argument num must exist in GF(2^8).
     */
    int _CountCauchyOnes(int num);

    /**
     * @brief generate cauchy encoding matrix and improve it
     */
    int *_GenerateEncodeMatrix();

    /**
     * @brief convert matrix to bitmatrix, to convert multply and divide opertaion on
     *        GF(2^8) to more efficient XOR opertion, thus making encoding and much faster
     */
    char *_MatrixToBitMatrix(int *matrix);

    /**
     *  @brief convert bitmatrix to schedule, to avoid traversing the matrix
     *         during encoding schedule is a list of 5-tuples:
     *                  < op, sd, sb, dd, db >
     *         where op is an operation code: 0 for copy and 1 for XOR, sd is
     *         the id of the source device and sb is the bit of the source device.
     *         The last two elements, dd and db are the destination device and bit.
     */
    int **_BitMatrixToSchedule(int num_data_parts,
                               int num_code_parts,
                               char *bit_matrix);

    /**
     * @brief do operations in the schedule
     */
    void _DoScheduleOperations(int **schedule, char **ptrs, int size);

    void _FreeSchedule(int **schedule);

    int m_num_data_parts;          ///< number of data parts
    int m_num_code_parts;          ///< number of coding parts
    GaloisOperator  *m_galois_operator;  ///< galois filed operator
    char *m_encoding_bit_matrix;            ///< bit matrix used in encoding/decoding
    int **m_encoding_schedule;     ///< coding schedule used for encoding
};
